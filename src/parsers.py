# -*- coding: utf-8 -*-

import re
import sys
from turtle import pos
from rdkit import Chem
from rdkit.Chem import inchi
import gc
from os.path import join
import numpy as np
from reportlab.lib.pdfencrypt import checkU
from tqdm import tqdm
from Bio import SeqIO
import io
from .extras import *
from timeit import default_timer as timer
import xml.etree.ElementTree as ET
from zipfile import ZipFile

import pandas as pd
import tarfile
import os
import sqlite3
import csv
import json
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

link_head = ['InChIKey', 'SMILES', 'Name', 'InChI', 'CAS Number', 'DrugBank ID',
            'KEGG Compound ID', 'KEGG Drug ID', 'PubChem Compound ID',
            'PubChem Substance ID', 'ChEBI ID', 'ChEMBL ID', 'HET ID',
            'ChemSpider ID', 'BindingDB ID']


P_UNIPROT_CODE = re.compile(r"[OPQ][0-9][A-Z0-9][A-Z0-9][A-Z0-9][0-9]")
P_DISEASE_CODE = re.compile(r"MIM:\d+")
P_GO_ANOT_CODE = re.compile(r"GO:\d{7}")
P_PFAM____CODE = re.compile(r"PF:\d{3,6}")
P_MODRES__CODE = re.compile(r"^\w*(\.|\;)")
P_EC_NUM__CODE = re.compile(r"EC=\d+.\d+.\d+.\d+")
PUBMED_ID_CODE = re.compile(r"PubMed=\d+")
SEQ_RANGE_CODE = re.compile(r"\w+\s+\d+\.\.\d+")
SEQ_NOTE__CODE = re.compile(r"\w+\s+\d+\s")


def regular_type(df):
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    df[numeric_cols] = df[numeric_cols].astype("Int64").astype(str).replace("<NA>", np.nan)  # 兼容 Pandas NA
    return df

def export_quads(triplets, file_descriptor):
    """ Export quads to a file

    Parameters
    ----------
    triplets : list
        list of quads
    file_descriptor : file
        file descriptor of an open writable file
    """

    for s, p, o, q in triplets:
        file_descriptor.write("%s\t%s\t%s\t%s\n" % (s, p, o, q))


def export_triplets(triplets, file_descriptor):
    """ Export triplets to a file
    
    Parameters
    ----------
    triplets : list
        list of triplets
    file_descriptor : file
        file descriptor of an open writable file

    """

    for s, p, o in triplets:
        file_descriptor.write("%s\t%s\t%s\n" % (s, p, o))


def sanatize_text(text):
    """ Replace non alphanumeric characters in text with '_'

    Parameters
    ----------
    text : str
        text to sanatize

    Returns
    -------
    text
        the sanatized text
    """
    if text is None:
        return text
    return re.sub('[^a-zA-Z0-9]', '_', text.strip())


def sanatize_se_txt(txt):
    return txt.strip().replace(" ", "_").lower()


class UniProtGPCRParser:
    """
    A UNIPROT GPCR database text file parser class
    """
    def __init__(self):
        self.filepath = ""
        self.output_dp = ""
        self._filename = "gpcr_entries.json"

    @property
    def filenames(self):
        return [self._filename]

    def __parse_gpcr_entry(self, entry):
        """ Process a Uniprot GPCR txt entry

        Parameters
        ----------
        entry : list
            list of str lines representing the entry lines

        Returns
        -------
        dict
            dictionary of extracted data
        """
        entry_data = {
            "uniprot_id": "",
            "history_ids": [],
            "gene_name": [],
            "gene_synonyms": [],
            "protein_synonyms": [],
            "refseq_ids": [],
            "pubmed_ids": [],
            "seq_ranges": [],
            "seq_annotations": [],
            "pdb_ids": []
        }
        
        entry_dictionary = dict()
        for line in entry:
            line_prefix = line[:2]
            if line_prefix == "  ":
                line_prefix = "AA"

            if line_prefix not in entry_dictionary:
                entry_dictionary[line_prefix] = ""

            entry_dictionary[line_prefix] += " " + line[2:].lstrip()

        # Extract UniProt ID and history IDs
        if "AC" in entry_dictionary:
            ac_parts = entry_dictionary["AC"].strip().split(";")
            entry_data["uniprot_id"] = ac_parts[0].strip()
            # History IDs are the remaining ACs
            for i in range(1, len(ac_parts)):
                if ac_parts[i].strip():
                    entry_data["history_ids"].append(ac_parts[i].strip())

        # Extract gene names and synonyms
        if "GN" in entry_dictionary:
            gn_content = entry_dictionary["GN"].strip()
            
            # Extract gene name (Name=)
            gene_name_match = re.search(r'Name=([^;{]+)', gn_content)
            if gene_name_match:
                entry_data["gene_name"].append(gene_name_match.group(1).strip())
            
            # Extract gene synonyms (Synonyms=)
            synonyms_match = re.search(r'Synonyms=([^;{]+)', gn_content)
            if synonyms_match:
                synonyms = [s.strip() for s in synonyms_match.group(1).split(',')]
                entry_data["gene_synonyms"].extend(synonyms)

        # Extract protein names and synonyms from DE section
        if "DE" in entry_dictionary:
            de_content = entry_dictionary["DE"]
            
            # Extract RecName and AltName Full names
            full_names = re.findall(r'(?:RecName|AltName):[^;]*Full=([^;{]+)', de_content)
            for name in full_names:
                entry_data["protein_synonyms"].append(name.strip())
            
            # Extract Short names
            short_names = re.findall(r'Short=([^;{]+)', de_content)
            for name in short_names:
                entry_data["protein_synonyms"].append(name.strip())

        # Extract RefSeq IDs from DR section
        if "DR" in entry_dictionary:
            dr_content = entry_dictionary["DR"]
            refseq_matches = re.findall(r'RefSeq;\s*([^;]+)', dr_content)
            for match in refseq_matches:
                entry_data["refseq_ids"].append(match.strip())
            pdbsum_matches = re.findall(r'PDBsum;\s*([^;]+)', dr_content)
            for match in pdbsum_matches:
                entry_data["pdb_ids"].append(match.strip())

        # Extract PubMed IDs from RX section
        if "RX" in entry_dictionary:
            rx_content = entry_dictionary["RX"]
            pubmed_matches = re.findall(r'PubMed=(\d+)', rx_content)
            for match in pubmed_matches:
                entry_data["pubmed_ids"].append("pubmed:" + match)

        # Extract sequence ranges and annotations from FT section
        # if "FT" in entry_dictionary:
        #     ft_content = entry_dictionary["FT"]
            
        #     # Extract sequence ranges (similar to existing code)
        #     seq_ranges = re.findall(SEQ_RANGE_CODE, ft_content)
        #     for seq_range in seq_ranges:
        #         entry_data["seq_ranges"].append(seq_range.strip())
            
        #     # Extract sequence annotations (similar to existing code)
        #     seq_annotations = re.findall(SEQ_NOTE__CODE, ft_content)
        #     for seq_annotation in seq_annotations:
        #         entry_data["seq_annotations"].append(seq_annotation.strip())

        return entry_data

    def __parse_gpcr_fasta(self, fasta_filepath, entries_dict):
        """ Parse GPCR FASTA file and add sequence information to entries dictionary

        Parameters
        ----------
        fasta_filepath : str
            absolute path to the FASTA file
        entries_dict : dict
            dictionary of GPCR entries to update with sequence information

        Returns
        -------
        dict
            updated dictionary with sequence information
        """
        print_section_header(f"Parsing GPCR FASTA file ({bcolors.OKGREEN + fasta_filepath + bcolors.ENDC})")
        start = timer()
        
        sequences_processed = 0
        matched_entries = 0
        
        try:
            with open(fasta_filepath, 'r', encoding='utf-8') as fasta_file:
                for record in SeqIO.parse(fasta_file, "fasta"):
                    sequences_processed += 1
                    
                    # Parse the header information
                    # Format: >sp|A0A2R9YJI3|GPR22_DANRE G-protein coupled receptor 22 OS=Danio rerio OX=7955 GN=gpr22a PE=2 SV=1
                    header = record.description
                    
                    # Extract UniProt ID (between first and second |)
                    header_parts = header.split('|')
                    if len(header_parts) >= 3:
                        uniprot_id = header_parts[1].strip()
                        entry_name = header_parts[2].split(' ')[0].strip()  # GPR22_DANRE
                        
                        # Extract full description after the entry name
                        remaining_header = ' '.join(header_parts[2].split(' ')[1:])
                        
                        # Extract species (OS=...)
                        species = ""
                        os_match = re.search(r'OS=([^O]+?)(?=\sOX=|\sGN=|\sPE=|\sSV=|$)', remaining_header)
                        if os_match:
                            species = os_match.group(1).strip()
                        
                        # Extract taxon ID (OX=...)
                        taxon_id = ""
                        ox_match = re.search(r'OX=(\d+)', remaining_header)
                        if ox_match:
                            taxon_id = ox_match.group(1).strip()
                        
                        # Extract gene name (GN=...)
                        gene_name = ""
                        gn_match = re.search(r'GN=([^O]+?)(?=\sOX=|\sOS=|\sPE=|\sSV=|$)', remaining_header)
                        if gn_match:
                            gene_name = gn_match.group(1).strip()
                        
                        # Get full name (everything before OS=)
                        full_name = ""
                        full_name_match = re.search(r'^([^O]+?)(?=\sOS=)', remaining_header)
                        if full_name_match:
                            full_name = f"{entry_name} {full_name_match.group(1).strip()}"
                        else:
                            full_name = entry_name
                        
                        # Get sequence and calculate length
                        sequence = str(record.seq)
                        sequence_length = f"{len(sequence)} AA"
                        
                        # Update the entries dictionary
                        if uniprot_id in entries_dict:
                            # Update existing entry
                            entries_dict[uniprot_id].update({
                                "entry_name": entry_name,
                                "full_name": full_name,
                                "species": species,
                                "taxon_id": taxon_id,
                                "sequence": sequence,
                                "sequence_length": sequence_length
                            })
                            
                            # Add gene name if found and not already present
                            if gene_name and gene_name not in entries_dict[uniprot_id]["gene_name"]:
                                entries_dict[uniprot_id]["gene_name"].append(gene_name)
                            
                            matched_entries += 1
                        else:
                            # Create new entry if not found in text parsing
                            entries_dict[uniprot_id] = {
                                "uniprot_id": uniprot_id,
                                "entry_name": entry_name,
                                "full_name": full_name,
                                "species": species,
                                "taxon_id": taxon_id,
                                "gene_name": [gene_name] if gene_name else [],
                                "gene_synonyms": [],
                                "protein_synonyms": [],
                                "history_ids": [],
                                "refseq_ids": [],
                                "pubmed_ids": [],
                                "seq_ranges": [],
                                "seq_annotations": [],
                                "sequence": sequence,
                                "sequence_length": sequence_length
                            }
                    
                    if sequences_processed % 100 == 0:
                        speed = int(sequences_processed / (timer() - start))
                        msg = prc_sym + f"Processing ({speed}) sequences/second => [processed: {sequences_processed}, matched: {matched_entries}]"
                        print("\r" + msg, end="", flush=True)
                        
        except FileNotFoundError:
            print(f"Warning: FASTA file not found: {fasta_filepath}")
            return entries_dict
        except Exception as e:
            print(f"Error parsing FASTA file: {e}")
            return entries_dict
        
        print(done_sym + " Took %1.2f Seconds." % (timer() - start), flush=True)
        print(f"Processed {sequences_processed} sequences, matched {matched_entries} with existing entries")
        
        return entries_dict

    def parse(self, sources_dp, output_dp):
        """ Parse a GPCR textual data file and FASTA file, output findings to JSON file

        Parameters
        ----------
        sources_dp : str
            absolute path to the sources directory
        output_dp : str
            absolute path of the output directory
        """
        txt_filepath = join(sources_dp, "gpcr_entries.txt")
        fasta_filepath = join(sources_dp, "gpcr_sequences.fasta")
        self.filepath = txt_filepath
        self.output_dp = output_dp

        all_entries = {}
        line_index = 0
        
        print_section_header("Parsing GPCR text file (%s)" % (bcolors.OKGREEN + txt_filepath + bcolors.ENDC))
        start = timer()
        
        # Parse text file first
        with open(txt_filepath, 'r', encoding='utf-8') as fd:
            current_entry = []
            eof = False
            
            while not eof:
                raw_line = fd.readline()
                line = raw_line.rstrip()
                
                if line != "//":
                    current_entry.append(line)
                else:
                    if current_entry:
                        entry_data = self.__parse_gpcr_entry(current_entry)
                        if entry_data["uniprot_id"]:
                            all_entries[entry_data["uniprot_id"]] = entry_data
                    current_entry = []

                line_index += 1
                if line_index % 1000 == 0:
                    speed = int(line_index / (timer() - start))
                    msg = prc_sym + "Processing (%d) lines/second => [entries: %d]" % (speed, len(all_entries))
                    print("\r" + msg, end="", flush=True)

                if raw_line == "":
                    eof = True
                    print(done_sym + " Took %1.2f Seconds." % (timer() - start), flush=True)

        # Parse FASTA file and merge with text data
        all_entries = self.__parse_gpcr_fasta(fasta_filepath, all_entries)

        # Save to JSON file
        output_file = join(output_dp, self._filename)
        with open(output_file, 'w', encoding='utf-8') as json_fd:
            json.dump(all_entries, json_fd, indent=2, ensure_ascii=False)
        
        print(f"Saved {len(all_entries)} GPCR entries to {output_file}")

class ChEMBLParser:
    """
    ChEMBL database parser
    """
    def __init__(self):
        """
        Initialise ChEMBL parser class instance
        """
        self._filenames = [
            'chembl_cpi.tsv',
            'chembl_dti.tsv',
            'chembl_cinfo.tsv',
            'chembl_tinfo.tsv',
            'chembl_cpi.txt',
            'chembl_dti.txt'
        ]
        self.init_querys()
        self.db_file_path = "database/processed/chembl/chembl_sqlite.db"

    @property
    def filenames(self):
        """
        Get MedGen filenames

        Returns
        -------
        filename : str
            the names of the MedGen output files
        """
        return self._filenames

    def init_querys(self):
        self.query_cpi = """
            SELECT
                m.chembl_id  AS compound_chembl_id,
                t.chembl_id  AS target_chembl_id,
                act.standard_type,
                act.standard_relation,
                act.standard_value,
                act.standard_units,
                act.activity_comment,
                GROUP_CONCAT(COALESCE(d.pubmed_id, d.doi), ',') AS pubmed_id_or_doi
                
            FROM activities AS act
                JOIN molecule_dictionary m ON act.molregno = m.molregno
                JOIN compound_records r ON m.molregno = r.molregno
                JOIN docs d ON r.doc_id = d.doc_id
                JOIN assays a ON act.assay_id = a.assay_id
                JOIN target_dictionary t ON a.tid = t.tid
            
            WHERE 
                a.assay_type = 'B'
                
            GROUP BY 
                m.chembl_id , t.chembl_id
                
            ORDER BY 
                m.chembl_id , t.chembl_id;
        """

        self.query_dti = """
            SELECT
                m.chembl_id  AS compound_chembl_id,
                dm.action_type,
                t.chembl_id  AS target_chembl_id,
                GROUP_CONCAT(mr.ref_type || '::' || mr.ref_id, ';') AS refs
                
            FROM drug_mechanism AS dm
                JOIN molecule_dictionary m ON dm.molregno = m.molregno
                JOIN target_dictionary t ON dm.tid = t.tid
                JOIN mechanism_refs mr ON dm.mec_id = mr.mec_id

            GROUP BY 
                m.chembl_id , t.chembl_id

            ORDER BY 
                m.chembl_id , t.chembl_id;
            """

        self.query_cinfo = """
            SELECT
                cs.standard_inchi_key AS InChIKey,
                cs.canonical_smiles AS SMILES,
                m.chembl_id  AS "ChEMBL ID",
                cs.standard_inchi AS InChI

            FROM molecule_dictionary AS m
                JOIN compound_structures cs ON m.molregno = cs.molregno

            GROUP BY 
                m.chembl_id

            ORDER BY 
                m.chembl_id;
            """

        self.query_tinfo = """
            SELECT
                t.chembl_id  AS target_chembl_id,
                t.pref_name AS target_name,
                c.accession AS protein_accession,
                c.sequence AS protein_sequence

            FROM target_dictionary AS t
                JOIN target_components tc ON t.tid = tc.tid
                JOIN component_sequences c ON tc.component_id = c.component_id 

            GROUP BY 
                t.chembl_id

            ORDER BY 
                t.chembl_id;
            """

    def unpack_db(self, filepath, output_dp):
        """ Parse ChEMBL database tar.gz file

        Parameters
        ----------
        filepath : str
            Absolute path to the ChEMBL .tar.gz database file
        output_dp : str
            Absolute path of the output directory
        """
        self.filepath = filepath
        self.output_dp = output_dp

        # Ensure the output directory exists
        os.makedirs(output_dp, exist_ok=True)

        # Extract the tar.gz file
        with tarfile.open(filepath, "r:gz") as tar:
            members = tar.getmembers()
            for member in members:
                # Flatten nested directories
                member.name = os.path.basename(member.name)
                if member.name:  # Only process non-empty names
                    tar.extract(member, path=output_dp)

        # Move extracted files directly to output_dp if nested in subfolders
        for root, _, files in os.walk(output_dp):
            for file in files:
                source_path = os.path.join(root, file)
                target_path = os.path.join(output_dp, file)
                if source_path != target_path:
                    os.rename(source_path, target_path)

        # Clean up any remaining empty directories
        for root, dirs, _ in os.walk(output_dp, topdown=False):
            for directory in dirs:
                dir_path = os.path.join(root, directory)
                if not os.listdir(dir_path):
                    os.rmdir(dir_path)
                # Identify and return the .db file path

        db_file_path = ""
        for file in os.listdir(output_dp):
            if file.endswith(".db"):
                db_file_path = os.path.join(output_dp, file)
                break
        new_file_path = os.path.join(output_dp, "chembl_sqlite.db")
        os.rename(db_file_path, new_file_path)
        self.db_file_path = new_file_path

    def db_query(self, query, output_file):
        print(f"\rQuerying {output_file}, please be patient, this may take a while...", end='', flush=True)
        conn = sqlite3.connect(self.db_file_path)
        cursor = conn.cursor()
        cursor.execute(query)
        columns = [description[0] for description in cursor.description]
        rows = cursor.fetchall()
        print(f"\rwriting {output_file}...", end='', flush=True)
        with open(output_file, "w", newline="", encoding="utf-8") as file:
            writer = csv.writer(file, delimiter='\t')  # 指定分隔符为'\t'
            writer.writerow(columns)  # 写入列名
            writer.writerows(rows)  # 写入数据
        conn.close()

    def form_data(self, output_dp):
        # get gpcr uniprot ids
        with open(join(output_dp, '../uniprot/gpcr_entries.json'), 'r', encoding='utf-8') as fd:
            gpcr_entries = json.load(fd)
        gpcr_uniprot_ids = set([sys.intern(x) for x in gpcr_entries.keys()])
        print(f"Loaded {len(gpcr_uniprot_ids)} GPCR UniProt IDs")

        # 读取映射字典
        cinfo = pd.read_csv(join(output_dp, 'chembl_cinfo.tsv'), sep='\t')
        cdict = {v[0]: v[1] for v in cinfo[['ChEMBL ID', 'InChIKey']].values}
        tinfo = pd.read_csv(join(output_dp, 'chembl_tinfo.tsv'), sep='\t')
        tdict = {v[0]: v[1] for v in tinfo[['target_chembl_id', 'protein_accession']].values}
        # release the RAM
        del cinfo; gc.collect()
        del tinfo; gc.collect()

        # 处理CPI数据
        print(f"\rProcessing CPI data...", end='', flush=True)
        start = timer()
        nb_entries = 0
        chunksize = 10000
        cpi_chunks = []
        
        for chunk in pd.read_csv(join(output_dp, 'chembl_cpi.tsv'), sep='\t', chunksize=chunksize, on_bad_lines='skip'):
            # 添加规范化列
            chunk['compound_inchikey'] = chunk['compound_chembl_id'].map(cdict)
            chunk['target_uniprot_id'] = chunk['target_chembl_id'].map(tdict)
            
            # 使用sys.intern过滤GPCR数据
            chunk['target_uniprot_id_intern'] = chunk['target_uniprot_id'].apply(lambda x: sys.intern(x) if pd.notna(x) else x)
            chunk_filtered = chunk[chunk['target_uniprot_id_intern'].isin(gpcr_uniprot_ids)]
            
            # 删除临时列
            chunk_filtered = chunk_filtered.drop('target_uniprot_id_intern', axis=1)
            
            if not chunk_filtered.empty:
                cpi_chunks.append(chunk_filtered)
            
            nb_entries += len(chunk)
            if nb_entries % 10000 == 0:
                speed = nb_entries / (timer() - start)
                msg = prc_sym + "Processed (%d) CPI entries.  Speed: (%1.2f) entries/second" % (
                    nb_entries, speed)
                print("\r" + msg, end="", flush=True)

        # 合并并保存CPI数据
        if cpi_chunks:
            cpi_final = pd.concat(cpi_chunks, ignore_index=True)
            # 重新排列列顺序，将规范化列放在对应ID后面
            cols = list(cpi_final.columns)
            # 找到compound_chembl_id的位置，在其后插入compound_inchikey
            compound_idx = cols.index('compound_chembl_id')
            cols.insert(compound_idx + 1, cols.pop(cols.index('compound_inchikey')))
            # 找到target_chembl_id的位置，在其后插入target_uniprot_id
            target_idx = cols.index('target_chembl_id')
            cols.insert(target_idx + 1, cols.pop(cols.index('target_uniprot_id')))
            cpi_final = cpi_final[cols]
            
            cpi_final.to_csv(join(output_dp, 'chembl_cpi.tsv'), sep='\t', index=False)
            print(f"\rSaved {len(cpi_final)} filtered CPI entries", flush=True)
        else:
            print(f"\rNo GPCR CPI entries found", flush=True)

        # 处理DTI数据
        print(f"\rProcessing DTI data...", end='', flush=True)
        start = timer()
        nb_entries = 0
        dti_chunks = []
        
        for chunk in pd.read_csv(join(output_dp, 'chembl_dti.tsv'), sep='\t', chunksize=chunksize, on_bad_lines='skip'):
            # 添加规范化列
            chunk['compound_inchikey'] = chunk['compound_chembl_id'].map(cdict)
            chunk['target_uniprot_id'] = chunk['target_chembl_id'].map(tdict)
            
            # 使用sys.intern过滤GPCR数据
            chunk['target_uniprot_id_intern'] = chunk['target_uniprot_id'].apply(lambda x: sys.intern(x) if pd.notna(x) else x)
            chunk_filtered = chunk[chunk['target_uniprot_id_intern'].isin(gpcr_uniprot_ids)]
            
            # 删除临时列
            chunk_filtered = chunk_filtered.drop('target_uniprot_id_intern', axis=1)
            
            if not chunk_filtered.empty:
                dti_chunks.append(chunk_filtered)
            
            nb_entries += len(chunk)
            if nb_entries % 10000 == 0:
                speed = nb_entries / (timer() - start)
                msg = prc_sym + "Processed (%d) DTI entries.  Speed: (%1.2f) entries/second" % (
                    nb_entries, speed)
                print("\r" + msg, end="", flush=True)

        # 合并并保存DTI数据
        if dti_chunks:
            dti_final = pd.concat(dti_chunks, ignore_index=True)
            # 重新排列列顺序，将规范化列放在对应ID后面
            cols = list(dti_final.columns)
            # 找到compound_chembl_id的位置，在其后插入compound_inchikey
            compound_idx = cols.index('compound_chembl_id')
            cols.insert(compound_idx + 1, cols.pop(cols.index('compound_inchikey')))
            # 找到target_chembl_id的位置，在其后插入target_uniprot_id
            target_idx = cols.index('target_chembl_id')
            cols.insert(target_idx + 1, cols.pop(cols.index('target_uniprot_id')))
            dti_final = dti_final[cols]
            
            dti_final.to_csv(join(output_dp, 'chembl_dti.tsv'), sep='\t', index=False)
            print(f"\rSaved {len(dti_final)} filtered DTI entries", flush=True)
        else:
            print(f"\rNo GPCR DTI entries found", flush=True)

        print(done_sym + " GPCR filtering and normalization completed.")

    def parse_chem_db(self, source_dp, output_dp):
        """
        Parse ChEMBL files

        Parameters
        ----------
        source_dp : str
            The path to the source directory
        output_dp : str
            The path to the output directory
        """
        print_section_header(
            "Parsing ChEMBL files (%s)" %
            (bcolors.OKGREEN + source_dp + '/chembl_sqlite.tar.gz' + bcolors.ENDC)
        )
        start = timer()
        nb_entries = 0

        input_fp = join(source_dp, 'chembl_sqlite.tar.gz')
        print(f"\runpacking {input_fp}...", end='', flush=True)
        self.db_file_path = join(output_dp, "chembl_sqlite.db")
        if not os.path.exists(self.db_file_path):
            self.unpack_db(input_fp, output_dp)
        nb_entries += 1
        self.db_query(self.query_cpi, os.path.join(output_dp, 'chembl_cpi.tsv')) # 化合物-蛋白质互作
        nb_entries += 1
        self.db_query(self.query_dti, os.path.join(output_dp, 'chembl_dti.tsv')) # 药物-靶标互作
        nb_entries += 1
        self.db_query(self.query_cinfo, os.path.join(output_dp, 'chembl_cinfo.tsv')) # 化合物检索信息
        nb_entries += 1
        self.db_query(self.query_tinfo, os.path.join(output_dp, 'chembl_tinfo.tsv')) # 蛋白质检索信息
        nb_entries += 1
        print(done_sym + "Processed (%d) files. Took %1.2f Seconds." % (nb_entries, timer() - start), flush=True)
        self.form_data(output_dp)


class IUPHARParser:
    """
    IUPHAR database parser for GPCR interactions
    """
    def __init__(self):
        """
        Initialize IUPHAR parser class instance
        """
        self._filenames = [
            'iuphar_cpi.tsv',
            'iuphar_cinfo.tsv'
        ]
        self.filepath = ""
        self.output_dp = ""

    @property
    def filenames(self):
        """
        Get IUPHAR filenames

        Returns
        -------
        filenames : list
            the names of the IUPHAR output files
        """
        return self._filenames

    def _inchi_to_inchikey(self, inchi):
        try:
            # Parse InChI
            mol = Chem.MolFromInchi(inchi)
            if mol is None:
                return None
            
            # Generate InChIKey
            inchikey_str = inchi.MolToInchiKey(mol)
            return inchikey_str if inchikey_str else None
            
        except Exception as e:
            return None

    def _smiles_to_inchi_inchikey(self, smiles):
        try:
            # Parse SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None, None
            
            # Generate InChI
            inchi_str = inchi.MolToInchi(mol)
            if not inchi_str:
                return None, None
            
            # Generate InChIKey
            inchikey_str = inchi.MolToInchiKey(mol)
            if not inchikey_str:
                return inchi_str, None
            
            return inchi_str, inchikey_str
            
        except Exception as e:
            return None, None

    def parse(self, source_dp, output_dp):
        """
        Parse IUPHAR GPCR interactions file

        Parameters
        ----------
        source_dp : str
            The path to the source directory
        output_dp : str
            The path to the output directory
        """
        self.filepath = join(source_dp, 'gpcr_interactions.csv')
        self.output_dp = output_dp
        
        print_section_header(
            "Parsing IUPHAR GPCR interactions file (%s)" %
            (bcolors.OKGREEN + self.filepath + bcolors.ENDC)
        )
        start = timer()
        
        # Read the data file
        try:
            # Read the file, skipping comment lines and handling encoding
            df = pd.read_csv(
                self.filepath,
                quotechar='"',
                encoding='utf-8',
                comment='#',
                engine='python'  # 用 python 引擎能更好处理复杂引号
            )
            
        except Exception as e:
            print(f"Error reading file: {e}")
            return
        
        # Filter for GPCR targets by loading GPCR UniProt IDs
        gpcr_uniprot_ids = set()
        try:
            gpcr_file = join(output_dp, '../uniprot/gpcr_entries.json')
            if os.path.exists(gpcr_file):
                with open(gpcr_file, 'r', encoding='utf-8') as fd:
                    gpcr_entries = json.load(fd)
                gpcr_uniprot_ids = set([sys.intern(x) for x in gpcr_entries.keys()])
            else:
                print("Warning: GPCR entries file not found, processing all targets")
        except Exception as e:
            print(f"Warning: Could not load GPCR filter: {e}")
        
        # Process CPI data
        print("Processing compound-protein interaction data...")
        cpi_data = []
        
        processed_count = 0
        filtered_count = 0
        
        for idx, row in df.iterrows():
            processed_count += 1
            
            # Extract target information
            target_uniprot_ids = str(row.get('Target UniProt ID', None))
            
            # Extract ligand information
            ligand_id = row.get('Ligand ID', None)
            
            # Extract interaction information
            interaction_type = row.get('Type', None)
            action = row.get('Action', None)
            endogenous = row.get('Endogenous', None)
            primary_target = row.get('Primary Target', None)
            approved = row.get('Approved', None)
            
            # Extract affinity information
            affinity_units = row.get('Affinity Units', None)
            affinity_high = row.get('Affinity High', None)
            affinity_median = row.get('Affinity Median', None)
            affinity_low = row.get('Affinity Low', None)
            original_affinity_units = row.get('Original Affinity Units', None)
            original_affinity_relation = row.get('Original Affinity Relation', None)
            original_affinity_low = row.get('Original Affinity Low nm', None)
            original_affinity_high = row.get('Original Affinity High nm', None)
            original_affinity_median = row.get('Original Affinity Median nm', None)

            # Extract additional information
            pubmed_id = row.get('PubMed ID', '')

            for target_uniprot_id in target_uniprot_ids.split('|'):
                # Filter for GPCR targets if filter is available
                if gpcr_uniprot_ids and target_uniprot_id:
                    target_uniprot_intern = sys.intern(target_uniprot_id) if target_uniprot_id else None
                    if target_uniprot_intern not in gpcr_uniprot_ids:
                        continue           
                filtered_count += 1
                # Create CPI entry
                cpi_entry = {
                    'target_uniprot_id': target_uniprot_id,
                    'ligand_id': ligand_id,
                    'interaction_type': interaction_type,
                    'action': action,
                    'endogenous': endogenous,
                    'primary_target': primary_target,
                    'approved': approved,
                    'affinity_units': affinity_units,
                    'affinity_high': affinity_high,
                    'affinity_median': affinity_median,
                    'affinity_low': affinity_low,
                    'original_affinity_units': original_affinity_units,
                    'original_affinity_relation': original_affinity_relation,
                    'original_affinity_low': original_affinity_low,
                    'original_affinity_high': original_affinity_high,
                    'original_affinity_median': original_affinity_median,
                    'pubmed_id': pubmed_id
                }
                cpi_data.append(cpi_entry)
            
            # Progress update
            if processed_count % 1000 == 0:
                speed = processed_count / (timer() - start)
                msg = prc_sym + "Processed (%d) entries, filtered (%d) GPCR entries. Speed: (%1.2f) entries/second" % (
                    processed_count, filtered_count, speed)
                print("\r" + msg, end="", flush=True)

        # Convert to DataFrame and save CPI data
        cpi_df = None
        if cpi_data:
            cpi_df = pd.DataFrame(cpi_data)
            
            # Save CPI data (will be updated later with compound_inchikey)
            cpi_output_file = join(output_dp, 'iuphar_cpi.tsv')
            cpi_df.to_csv(cpi_output_file, sep='\t', index=False, encoding='utf-8')
        else:
            print("No CPI data to save")
        
        # Process ligand lists to create compound info
        print("Processing ligand lists for compound information...")
        self._process_ligand_list(source_dp, output_dp, cpi_df)
        
        print(done_sym + " Took %1.2f Seconds." % (timer() - start), flush=True)

    def _process_ligand_list(self, source_dp, output_dp, cpi_df):
        """
        Process ligand_list.csv to create compound information file

        Parameters
        ----------
        source_dp : str
            The path to the source directory
        output_dp : str
            The path to the output directory
        cpi_df : pandas.DataFrame
            The CPI dataframe to extract unique ligand IDs from
        """
        ligand_list_file = join(source_dp, 'ligand_list.csv')
        
        # Extract unique ligand IDs from CPI data
        unique_ligand_ids = set()
        if cpi_df is not None and not cpi_df.empty:
            unique_ligand_ids = set(cpi_df['ligand_id'].dropna().astype(str).unique())
        else:
            print("No CPI data available, processing all ligands")
        
        try:
            # Read ligand lists file
            ligand_df = pd.read_csv(
                ligand_list_file,
                quotechar='"',
                encoding='utf-8',
                comment='#',
                engine='python'
            )
            
            # Filter for ligands that appear in CPI data (if CPI data exists)
            if unique_ligand_ids:
                ligand_df['Ligand ID'] = ligand_df['Ligand ID'].astype(str)
                filtered_ligand_df = ligand_df[ligand_df['Ligand ID'].isin(unique_ligand_ids)]
            else:
                filtered_ligand_df = ligand_df
            
            if filtered_ligand_df.empty:
                print("No ligands to process for compound information")
                return
            
            # Create compound info data following link_head format
            cinfo_data = []
            valid_ligand_ids = set()  # Track ligands with valid structures
            conversion_stats = {
                'original_inchikey': 0,
                'converted_from_inchi': 0,
                'converted_from_smiles': 0,
                'no_structure': 0
            }
            
            for idx, row in filtered_ligand_df.iterrows():
                ligand_id = row.get('Ligand ID', None)
                original_inchikey = row.get('InChIKey', None)
                original_inchi = row.get('InChI', None)
                original_smiles = row.get('SMILES', None)
                
                # Determine final InChI and InChIKey
                final_inchikey = None
                final_inchi = original_inchi
                
                # Priority 1: Use existing InChIKey if available
                if original_inchikey and isinstance(original_inchikey, str) and original_inchikey.strip():
                    final_inchikey = original_inchikey.strip()
                    conversion_stats['original_inchikey'] += 1
                
                # Priority 2: Convert from InChI if available
                elif original_inchi and isinstance(original_inchi, str) and original_inchi.strip():
                    final_inchikey = self._inchi_to_inchikey(original_inchi.strip())
                    if final_inchikey:
                        conversion_stats['converted_from_inchi'] += 1
                
                # Priority 3: Convert from SMILES if available
                elif original_smiles and isinstance(original_smiles, str) and original_smiles.strip():
                    converted_inchi, converted_inchikey = self._smiles_to_inchi_inchikey(original_smiles.strip())
                    if converted_inchikey:
                        final_inchikey = converted_inchikey
                        final_inchi = converted_inchi or original_inchi  # Use converted InChI if available
                        conversion_stats['converted_from_smiles'] += 1
                
                # Skip ligands without any valid structure representation
                if not final_inchikey:
                    conversion_stats['no_structure'] += 1
                    continue
                
                # Track valid ligand IDs
                if ligand_id:
                    valid_ligand_ids.add(str(ligand_id))
                
                # Map IUPHAR columns to link_head format
                cinfo_entry = {
                    'InChIKey': final_inchikey,
                    'SMILES': original_smiles,
                    'Name': row.get('Name', None),
                    'InChI': final_inchi,
                    'CAS Number': None,  # Not available in IUPHAR data
                    'DrugBank ID': None,  # Not available in IUPHAR data
                    'KEGG Compound ID': None,  # Not available in IUPHAR data
                    'KEGG Drug ID': None,  # Not available in IUPHAR data
                    'PubChem Compound ID': row.get('PubChem CID', None),
                    'PubChem Substance ID': row.get('PubChem SID', None),
                    'ChEBI ID': None,  # Not available in IUPHAR data
                    'ChEMBL ID': row.get('ChEMBL ID', None),
                    'HET ID': None,  # Not available in IUPHAR data
                    'ChemSpider ID': None,  # Not available in IUPHAR data
                    'BindingDB ID': None,  # Not available in IUPHAR data
                    # Additional IUPHAR-specific fields
                    'Ligand ID': ligand_id,
                    'Species': row.get('Species', None),
                    'Type': row.get('Type', None),
                    'Approved': row.get('Approved', None),
                    'IUPAC name': row.get('IUPAC name', None),
                }
                cinfo_data.append(cinfo_entry)
            
            # Print conversion statistics
            print(f"InChIKey conversion statistics:")
            print(f"  - Original InChIKey: {conversion_stats['original_inchikey']}")
            print(f"  - Converted from InChI: {conversion_stats['converted_from_inchi']}")
            print(f"  - Converted from SMILES: {conversion_stats['converted_from_smiles']}")
            print(f"  - No valid structure (dropped): {conversion_stats['no_structure']}")
            
            # Create DataFrame and save compound info
            if cinfo_data:
                cinfo_df = pd.DataFrame(cinfo_data)
                
                # Remove duplicates based on InChIKey
                cinfo_df = cinfo_df.drop_duplicates(subset=['InChIKey'])
                
                # Save compound info data
                cinfo_output_file = join(output_dp, 'iuphar_cinfo.tsv')
                cinfo_df.to_csv(cinfo_output_file, sep='\t', index=False, encoding='utf-8')
                print(f"Saved {len(cinfo_df)} compound entries to {cinfo_output_file}")
                
                # Update CPI data to include compound_inchikey and filter out invalid ligands
                if cpi_df is not None and not cpi_df.empty:
                    # Filter CPI data to only include ligands with valid structures
                    cpi_df_filtered = cpi_df[cpi_df['ligand_id'].astype(str).isin(valid_ligand_ids)].copy()

                    # Create ligand ID to InChIKey mapping
                    ligand_inchikey_map = dict(zip(cinfo_df['Ligand ID'].astype(str), cinfo_df['InChIKey']))
                    
                    # Add compound_inchikey column to filtered CPI data
                    cpi_df_filtered['compound_inchikey'] = cpi_df_filtered['ligand_id'].astype(str).map(ligand_inchikey_map)
                    
                    # Reorder columns to put compound_inchikey after ligand_id
                    cols = list(cpi_df_filtered.columns)
                    ligand_id_idx = cols.index('ligand_id')
                    cols.insert(ligand_id_idx + 1, cols.pop(cols.index('compound_inchikey')))
                    cpi_df_filtered = cpi_df_filtered[cols]
                    
                    # Save updated CPI data
                    cpi_output_file = join(output_dp, 'iuphar_cpi.tsv')
                    cpi_df_filtered.to_csv(cpi_output_file, sep='\t', index=False, encoding='utf-8')
                    
                    dropped_entries = len(cpi_df) - len(cpi_df_filtered)
                    print(f"Updated CPI data: kept {len(cpi_df_filtered)} entries, dropped {dropped_entries} entries without valid structures")
                
            else:
                print("No compound info data to save")
                
        except Exception as e:
            print(f"Error processing ligand_list.csv: {e}")
            return


class BindingDBParser:
    """
    BindingDB database parser
    """
    def __init__(self):
        """
        Initialise BindingDB parser class instance
        """
        self._filenames = [
            'bindingdb_cpi.tsv',
            'bindingdb_cinfo.tsv'
        ]
        self.filepath = ""
        self.output_dp = ""

    @property
    def filenames(self):
        """
        Get BindingDB filenames

        Returns
        -------
        filename : str
            the names of the BindingDB output files
        """
        return self._filenames

    def parse(self, source_dp, output_dp):
        """
        Parse BindingDB files

        Parameters
        ----------
        source_dp : str
            The path to the source directory
        output_dp : str
            The path to the output directory
        """
        self.filepath = join(source_dp, 'BindingDB_All.zip')
        self.output_dp = output_dp
        
        print_section_header(
            "Parsing BindingDB files (%s)" %
            (bcolors.OKGREEN + self.filepath + bcolors.ENDC)
        )
        start = timer()
        
        # Load GPCR UniProt IDs for filtering
        gpcr_uniprot_ids = set()
        try:
            gpcr_file = join(output_dp, '../uniprot/gpcr_entries.json')
            if os.path.exists(gpcr_file):
                with open(gpcr_file, 'r', encoding='utf-8') as fd:
                    gpcr_entries = json.load(fd)
                gpcr_uniprot_ids = set([sys.intern(x) for x in gpcr_entries.keys()])
            else:
                print("Warning: GPCR entries file not found, processing all targets")
        except Exception as e:
            print(f"Warning: Could not load GPCR filter: {e}")
        
        filename = 'BindingDB_All.tsv'
        cpi_data = []
        cinfo_data = []
        processed_count = 0
        filtered_count = 0
        
        with ZipFile(self.filepath, 'r') as dbzip:
            with dbzip.open(filename, 'r', force_zip64=True) as fd:
                chunksize = 1000
                for chunk in pd.read_csv(fd, sep='\t', chunksize=chunksize, on_bad_lines='skip'):
                    for idx, row in chunk.iterrows():
                        processed_count += 1
                        
                        # Extract target UniProt IDs (try multiple chains)
                        target_uniprot_ids = []
                        for chain_num in range(1, 11):  # Up to 10 chains
                            if chain_num == 1:
                                uniprot_col = 'UniProt (SwissProt) Primary ID of Target Chain'
                            else:
                                uniprot_col = f'UniProt (SwissProt) Primary ID of Target Chain.{chain_num}'
                            
                            if uniprot_col in row and not pd.isna(row[uniprot_col]):
                                uniprot_id = str(row[uniprot_col]).strip()
                                if uniprot_id and uniprot_id not in target_uniprot_ids:
                                    target_uniprot_ids.append(uniprot_id)
                        
                        # Extract compound information
                        compound_id = row.get('BindingDB MonomerID', None)
                        compound_inchikey = row.get('Ligand InChI Key', None)
                        
                        # Skip if no valid compound or target
                        if pd.isna(compound_id) or not target_uniprot_ids:
                            continue
                        
                        # Filter for GPCR targets if filter is available
                        valid_targets = []
                        if gpcr_uniprot_ids:
                            for target_id in target_uniprot_ids:
                                target_intern = sys.intern(target_id) if target_id else None
                                if target_intern in gpcr_uniprot_ids:
                                    valid_targets.append(target_id)
                        else:
                            valid_targets = target_uniprot_ids
                        
                        if not valid_targets:
                            continue
                        
                        filtered_count += 1
                        
                        # Extract affinity data
                        ki_value = row.get('Ki (nM)', None)
                        ic50_value = row.get('IC50 (nM)', None)
                        kd_value = row.get('Kd (nM)', None)
                        ec50_value = row.get('EC50 (nM)', None)
                        kon_value = row.get('kon (M-1-s-1)', None)
                        koff_value = row.get('koff (s-1)', None)
                        
                        # Extract experimental conditions
                        ph_value = row.get('pH', None)
                        temp_value = row.get('Temp (C)', None)
                        
                        # Extract publication information
                        pmid = row.get('PMID', None)
                        doi = row.get('Article DOI', None)
                        
                        
                        # Create CPI entries for each valid target
                        for target_id in valid_targets:
                            cpi_entry = {
                                'compound_bindingdb_id': str(compound_id),
                                'compound_inchikey': compound_inchikey,
                                'target_uniprot_id': target_id,
                                'ki_nm': ki_value,
                                'ic50_nm': ic50_value,
                                'kd_nm': kd_value,
                                'ec50_nm': ec50_value,
                                'kon_m1s1': kon_value,
                                'koff_s1': koff_value,
                                'ph': ph_value,
                                'temperature_c': temp_value,
                                'pmid': pmid,
                                'doi': doi,
                                'curation_source': row.get('Curation/DataSource', None),
                                'pdb_id': row.get('PDB ID(s) for Ligand-Target Complex', None)
                            }
                            cpi_data.append(cpi_entry)
                        
                        # Create compound info entry (only once per unique compound)
                        if compound_inchikey and not any(c.get('InChIKey') == compound_inchikey for c in cinfo_data):
                            cinfo_entry = {
                                'InChIKey': compound_inchikey,
                                'SMILES': row.get('Ligand SMILES', None),
                                'Name': row.get('BindingDB Ligand Name', None),
                                'InChI': row.get('Ligand InChI', None),
                                'CAS Number': None,  # Not available in BindingDB
                                'DrugBank ID': row.get('DrugBank ID of Ligand', None),
                                'KEGG Compound ID': row.get('KEGG ID of Ligand', None),
                                'KEGG Drug ID': None,  # Not available in BindingDB
                                'PubChem Compound ID': row.get('PubChem CID', None),
                                'PubChem Substance ID': row.get('PubChem SID', None),
                                'ChEBI ID': row.get('ChEBI ID of Ligand', None),
                                'ChEMBL ID': row.get('ChEMBL ID of Ligand', None),
                                'HET ID': row.get('Ligand HET ID in PDB', None),
                                'ChemSpider ID': None,  # Not available in BindingDB
                                'BindingDB ID': str(compound_id),
                                # Additional BindingDB-specific fields
                                'IUPHAR_GRAC ID': row.get('IUPHAR_GRAC ID of Ligand', None),
                                'ZINC ID': row.get('ZINC ID of Ligand', None)
                            }
                            cinfo_data.append(cinfo_entry)
                        
                        # Progress update
                        if processed_count % 1000 == 0:
                            speed = processed_count / (timer() - start)
                            msg = prc_sym + "Processed (%d) entries, filtered (%d) GPCR entries. Speed: (%1.2f) entries/second" % (
                                processed_count, filtered_count, speed)
                            print("\r" + msg, end="", flush=True)
        
        print(f"\rProcessed {processed_count} total entries, {filtered_count} GPCR entries")
        
        # Convert to DataFrames and save
        if cpi_data:
            cpi_df = pd.DataFrame(cpi_data)
            # Remove duplicates based on compound-target pairs
            cpi_df = cpi_df.drop_duplicates(subset=['compound_bindingdb_id', 'target_uniprot_id'])
            
            # Save CPI data
            cpi_output_file = join(output_dp, 'bindingdb_cpi.tsv')
            cpi_df.to_csv(cpi_output_file, sep='\t', index=False, encoding='utf-8')
            # print(f"Saved {len(cpi_df)} CPI entries to {cpi_output_file}")
        else:
            print("No CPI data to save")
        
        if cinfo_data:
            cinfo_df = pd.DataFrame(cinfo_data)
            # Remove duplicates based on InChIKey
            cinfo_df = cinfo_df.drop_duplicates(subset=['InChIKey'])
            
            # Save compound info data
            cinfo_output_file = join(output_dp, 'bindingdb_cinfo.tsv')
            cinfo_df.to_csv(cinfo_output_file, sep='\t', index=False, encoding='utf-8')
            # print(f"Saved {len(cinfo_df)} compound entries to {cinfo_output_file}")
        else:
            print("No compound info data to save")
        
        print(done_sym + " Took %1.2f Seconds." % (timer() - start), flush=True)

class DrugBankParser:
    """
    DrugBank database parser
    """
    def __init__(self):
        """
        Initialise DrugBank parser class instance
        """
        self._filenames = [
            'drugbank_cpi.tsv',
            'drugbank_cinfo.tsv'
        ]
        self.filepath = ""
        self.output_dp = ""
        self._ns = {'db':'http://www.drugbank.ca'}

    @property
    def filenames(self):
        """
        Get DrugBank filenames

        Returns
        -------
        filename : str
            the names of the DrugBank output files
        """
        return self._filenames

    def __parse_target(self, target_element, drug_id, rel_type):
        """
        Parse a drug target and return target info

        Parameters
        ----------
        target_element : xml.etree.ElementTree.Element
            xml element
        drug_id : string
            id of the drug
        rel_type : string
            type of target relationship

        Returns
        -------
        list
            list of target entries
        """
        target_entries = []
        
        # Link to uniprot
        poly = target_element.find('./db:polypeptide', self._ns)
        if poly is None:
            return target_entries
            
        poly_id = None
        if 'source' in poly.attrib and poly.attrib['source'] == 'Swiss-Prot':
            poly_id = poly.attrib['id']
        else:
            for extern_id in poly.findall('./db:external-identifiers/db:external-identifier', self._ns):
                res = extern_id.find('db:resource', self._ns)
                val = extern_id.find('db:identifier', self._ns)
                if res is not None and val is not None and res.text == 'UniprotKB':
                    poly_id = sanatize_text(val.text)
                    break
                    
        if poly_id is None:
            return target_entries

        # Gather any references to pubmed
        pmids = target_element.findall('./db:references/db:articles/db:article/db:pubmed-id', self._ns)
        pmids = [sanatize_text(pmid.text) for pmid in filter(lambda x: x.text is not None, pmids)]
        ref_string = ','.join(pmids) if pmids else None

        # Gather all actions
        actions = target_element.findall('./db:actions/db:action', self._ns)
        formatted_actions = []
        for action in actions:
            action_text = sanatize_text(action.text)
            # Normalize action text
            if action_text == 'other' or action_text == 'other_unknown':
                action_text = 'unknown'
            if action_text == '':
                continue
            formatted_actions.append(action_text)
       
        # If no action provided set it to unknown
        if len(formatted_actions) == 0:
            formatted_actions = ['unknown']

        # Create entries for each action
        for action in formatted_actions:
            target_entry = {
                'compound_drugbank_id': drug_id,
                'target_uniprot_id': poly_id,
                'relation_type': rel_type,
                'action': action,
                'pmid': ref_string
            }
            target_entries.append(target_entry)
            
        return target_entries

    def __parse_drug(self, drug_element):
        """
        Parse a top level xml drug entry
        
        Parameters
        ----------
        drug_element : xml.etree.ElementTree.Element
            xml element

        Returns
        -------
        list
            list of CPI entries for this drug
        """
        cpi_entries = []
        
        # Get primary drug ID
        drug_id_elem = drug_element.find('./db:drugbank-id[@primary="true"]', self._ns)
        if drug_id_elem is None:
            return cpi_entries
        
        drug_id = drug_id_elem.text

        # Parse drug targets
        for target in drug_element.findall('./db:targets/db:target', self._ns):
            target_entries = self.__parse_target(target, drug_id, 'DRUG_TARGET')
            cpi_entries.extend(target_entries)

        for carrier in drug_element.findall('./db:carriers/db:carrier', self._ns):
            target_entries = self.__parse_target(carrier, drug_id, 'DRUG_CARRIER')
            cpi_entries.extend(target_entries)

        for transporter in drug_element.findall('./db:transporters/db:transporter', self._ns):
            target_entries = self.__parse_target(transporter, drug_id, 'DRUG_TRANSPORTER')
            cpi_entries.extend(target_entries)

        for enzyme in drug_element.findall('./db:enzymes/db:enzyme', self._ns):
            target_entries = self.__parse_target(enzyme, drug_id, 'DRUG_ENZYME')
            cpi_entries.extend(target_entries)

        return cpi_entries

    def __parse_compound_info(self, source_dp, output_dp):
        """
        Parse compound information from structure links.csv

        Parameters
        ----------
        source_dp : str
            The path to the source directory
        output_dp : str
            The path to the output directory

        Returns
        -------
        pandas.DataFrame
            DataFrame containing compound information
        """
        link_path = join(source_dp, 'external_drug.zip')
        
        # Extract the zip file
        with ZipFile(link_path, 'r') as zip_ref:
            zip_ref.extractall(output_dp)

        # Read structure links CSV
        links_df = pd.read_csv(join(output_dp, 'structure links.csv'))
        
        # Map to link_head format, only keeping columns that exist in DrugBank
        cinfo_data = []
        for idx, row in links_df.iterrows():
            cinfo_entry = {
                'InChIKey': row.get('InChIKey', None),
                'SMILES': row.get('SMILES', None),
                'Name': row.get('Name', None),
                'InChI': row.get('InChI', None),
                'CAS Number': row.get('CAS Number', None),
                'DrugBank ID': row.get('DrugBank ID', None),
                'KEGG Compound ID': row.get('KEGG Compound ID', None),
                'KEGG Drug ID': row.get('KEGG Drug ID', None),
                'PubChem Compound ID': row.get('PubChem Compound ID', None),
                'PubChem Substance ID': row.get('PubChem Substance ID', None),
                'ChEBI ID': row.get('ChEBI ID', None),
                'ChEMBL ID': row.get('ChEMBL ID', None),
                'HET ID': row.get('HET ID', None),
                'ChemSpider ID': row.get('ChemSpider ID', None),
                'BindingDB ID': row.get('BindingDB ID', None)
            }
            cinfo_data.append(cinfo_entry)
        
        # Create DataFrame
        cinfo_df = pd.DataFrame(cinfo_data)
        
        # Remove duplicates based on InChIKey
        cinfo_df = cinfo_df.drop_duplicates(subset=['InChIKey'])
        
        # Clean up extracted file
        os.remove(join(output_dp, 'structure links.csv'))
        
        return cinfo_df

    def parse(self, source_dp, output_dp):
        """
        Parse DrugBank files

        Parameters
        ----------
        source_dp : str
            The path to the source directory
        output_dp : str
            The path to the output directory
        """
        self.filepath = join(source_dp, 'drugbank_all_full_database.xml.zip')
        self.output_dp = output_dp
        
        print_section_header(
            "Parsing DrugBank files (%s)" %
            (bcolors.OKGREEN + self.filepath + bcolors.ENDC)
        )
        start = timer()
        
        # Load GPCR UniProt IDs for filtering
        gpcr_uniprot_ids = set()
        try:
            gpcr_file = join(output_dp, '../uniprot/gpcr_entries.json')
            if os.path.exists(gpcr_file):
                with open(gpcr_file, 'r', encoding='utf-8') as fd:
                    gpcr_entries = json.load(fd)
                gpcr_uniprot_ids = set([sys.intern(x) for x in gpcr_entries.keys()])
                print(f"Loaded {len(gpcr_uniprot_ids)} GPCR UniProt IDs for filtering")
            else:
                print("Warning: GPCR entries file not found, processing all targets")
        except Exception as e:
            print(f"Warning: Could not load GPCR filter: {e}")
        
        # Parse compound information first
        print("Processing compound information...")
        cinfo_df = self.__parse_compound_info(source_dp, output_dp)
        
        # Create DrugBank ID to InChIKey mapping
        drugbank_to_inchikey = {}
        for idx, row in cinfo_df.iterrows():
            drugbank_id = row.get('DrugBank ID')
            inchikey = row.get('InChIKey')
            if drugbank_id and inchikey:
                drugbank_to_inchikey[drugbank_id] = inchikey
        
        # Parse XML file for CPI data
        print("Processing compound-protein interactions...")
        cpi_data = []
        processed_count = 0
        filtered_count = 0
        
        with ZipFile(self.filepath, 'r') as dbzip:
            with dbzip.open('full database.xml', force_zip64=True) as xmlfile:
                for event, elem in ET.iterparse(xmlfile):
                    # Check if this is a drug element (not pathway drug elements)
                    if elem.tag == '{http://www.drugbank.ca}drug' and len(elem) > 2:
                        processed_count += 1
                        
                        # Parse drug targets
                        drug_cpi_entries = self.__parse_drug(elem)
                        
                        # Filter for GPCR targets
                        for entry in drug_cpi_entries:
                            target_uniprot_id = entry['target_uniprot_id']
                            
                            # Filter for GPCR targets if filter is available
                            if gpcr_uniprot_ids:
                                target_intern = sys.intern(target_uniprot_id) if target_uniprot_id else None
                                if target_intern not in gpcr_uniprot_ids:
                                    continue
                            
                            filtered_count += 1
                            
                            # Add InChIKey mapping
                            drugbank_id = entry['compound_drugbank_id']
                            entry['compound_inchikey'] = drugbank_to_inchikey.get(drugbank_id, None)
                            
                            cpi_data.append(entry)
                        
                        # Clear element to save memory
                        elem.clear()
                        
                        # Progress update
                        if processed_count % 100 == 0:
                            speed = processed_count / (timer() - start)
                            msg = prc_sym + "Processed (%d) drugs, filtered (%d) GPCR interactions. Speed: (%1.2f) drugs/second" % (
                                processed_count, filtered_count, speed)
                            print("\r" + msg, end="", flush=True)
        
        print(f"\rProcessed {processed_count} drugs, {filtered_count} GPCR interactions")
        
        # Convert to DataFrame and save CPI data
        if cpi_data:
            cpi_df = pd.DataFrame(cpi_data)
            
            # Reorder columns to put compound_inchikey after compound_drugbank_id
            cols = list(cpi_df.columns)
            drugbank_id_idx = cols.index('compound_drugbank_id')
            cols.insert(drugbank_id_idx + 1, cols.pop(cols.index('compound_inchikey')))
            cpi_df = cpi_df[cols]
            
            # Remove duplicates
            cpi_df = cpi_df.drop_duplicates(subset=['compound_drugbank_id', 'target_uniprot_id', 'action'])
            
            # Save CPI data
            cpi_output_file = join(output_dp, 'drugbank_cpi.tsv')
            cpi_df.to_csv(cpi_output_file, sep='\t', index=False, encoding='utf-8')
            print(f"Saved {len(cpi_df)} CPI entries to {cpi_output_file}")
        else:
            print("No CPI data to save")
        
        # Save compound info data
        if not cinfo_df.empty:
            # Only keep compounds that appear in CPI data
            if cpi_data:
                used_drugbank_ids = set([entry['compound_drugbank_id'] for entry in cpi_data])
                cinfo_df = cinfo_df[cinfo_df['DrugBank ID'].isin(used_drugbank_ids)]
            
            cinfo_output_file = join(output_dp, 'drugbank_cinfo.tsv')
            cinfo_df.to_csv(cinfo_output_file, sep='\t', index=False, encoding='utf-8')
            print(f"Saved {len(cinfo_df)} compound entries to {cinfo_output_file}")
        else:
            print("No compound info data to save")
        
        print(done_sym + " Took %1.2f Seconds." % (timer() - start), flush=True)

    def parse_cpi(self, source_dp, output_dp):
        """
        Legacy method for backward compatibility
        """
        self.parse(source_dp, output_dp)


class PDSPParser:
    """
    PDSP database parser
    """
    def __init__(self):
        """
        Initialise PDSP parser class instance
        """
        self._filenames = [
            'pdsp_cpi.tsv',
            'pdsp_cinfo.tsv'
        ]
        self.filepath = ""
        self.output_dp = ""

    @property
    def filenames(self):
        """
        Get PDSP filenames

        Returns
        -------
        filename : str
            the names of the PDSP output files
        """
        return self._filenames

    def _smiles_to_inchi_inchikey(self, smiles):
        """
        Convert SMILES to standardized InChI and InChIKey using RDKit
        
        Parameters
        ----------
        smiles : str
            SMILES string
            
        Returns
        -------
        tuple
            (inchi, inchikey) or (None, None) if conversion fails
        """
        try:
 
            # Parse SMILES
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None, None
            
            # Generate InChI
            inchi_str = inchi.MolToInchi(mol)
            if not inchi_str:
                return None, None
            
            # Generate InChIKey
            inchikey_str = inchi.MolToInchiKey(mol)
            if not inchikey_str:
                return inchi_str, None
            
            return inchi_str, inchikey_str

        except Exception as e:
            # Skip problematic SMILES
            return None, None

    def parse(self, source_dp, output_dp):
        """
        Parse PDSP files

        Parameters
        ----------
        source_dp : str
            The path to the source directory
        output_dp : str
            The path to the output directory
        """
        self.filepath = join(source_dp, 'ki_data.csv')
        self.output_dp = output_dp
        
        print_section_header(
            "Parsing PDSP files (%s)" %
            (bcolors.OKGREEN + self.filepath + bcolors.ENDC)
        )
        start = timer()
        
        # Load GPCR UniProt IDs and build gene name mapping
        gpcr_uniprot_ids = set()
        gene_name_to_uniprot = {}
        try:
            gpcr_file = join(output_dp, '../uniprot/gpcr_entries.json')
            if os.path.exists(gpcr_file):
                with open(gpcr_file, 'r', encoding='utf-8') as fd:
                    gpcr_entries = json.load(fd)
                gpcr_uniprot_ids = set([sys.intern(x) for x in gpcr_entries.keys()])
                
                # Build gene name to UniProt ID mapping
                for uniprot_id, entry_data in gpcr_entries.items():
                    gene_names = entry_data.get('gene_name', [])
                    gene_synonyms = entry_data.get('gene_synonyms', [])
                    
                    # Add all gene names and synonyms to mapping (uppercase for case-insensitive matching)
                    all_names = gene_names + gene_synonyms
                    for gene_name in all_names:
                        if gene_name and isinstance(gene_name, str):
                            gene_name_upper = gene_name.upper()
                            gene_name_to_uniprot[gene_name_upper] = uniprot_id
                
                print(f"Loaded {len(gpcr_uniprot_ids)} GPCR UniProt IDs")
                print(f"Built gene name mapping with {len(gene_name_to_uniprot)} entries")
            else:
                print("Warning: GPCR entries file not found, processing all targets")
        except Exception as e:
            print(f"Warning: Could not load GPCR filter: {e}")
        
        # Read the PDSP data file
        try:
            df = pd.read_csv(
                self.filepath,
                encoding='utf-8',
                on_bad_lines='skip'
            )
        except Exception as e:
            print(f"Error reading file: {e}")
            return
        
        # drop rows with empty 'SMILES'
        df = df.dropna(subset=['SMILES'])
        df = df[df['SMILES'].str.strip() != '']
        
        cpi_data = []
        cinfo_data = []
        processed_count = 0
        filtered_count = 0
        unmatched_genes = set()
        processed_smiles = {}  # Cache for SMILES to InChI/InChIKey conversion
        
        # Process each row
        for idx, row in df.iterrows():
            processed_count += 1
            
            # Extract target information
            unigene = row.get('Unigene', None)
            target_name = row.get('Name', None)
            species = row.get('species', None)
            
            # Extract compound information
            ligand_id = row.get('Ligand ID', None)
            ligand_name = row.get('Ligand Name', None)
            smiles = row.get('SMILES', None)
            cas_number = row.get('CAS', None)
            nsc_number = row.get('NSC', None)
            
            # Extract binding data
            ki_value = row.get('ki Val', None)
            ki_note = row.get('ki Note', None)
            
            # Extract reference information
            reference = row.get('Reference', None)
            link = row.get('Link', None)
            source = row.get('source', None)
            hotligand = row.get('Hotligand', None)
            
            # Skip if no essential information
            if pd.isna(ligand_id) or pd.isna(unigene):
                continue
            
            # Map gene name to UniProt ID
            target_uniprot_id = None
            if unigene and isinstance(unigene, str):
                unigene_upper = unigene.strip().upper()
                target_uniprot_id = gene_name_to_uniprot.get(unigene_upper)
                
                if not target_uniprot_id:
                    unmatched_genes.add(unigene.strip())
                    continue
            
            if not target_uniprot_id:
                continue
            
            # Filter for GPCR targets (already filtered by gene name mapping)
            if gpcr_uniprot_ids:
                target_intern = sys.intern(target_uniprot_id) if target_uniprot_id else None
                if target_intern not in gpcr_uniprot_ids:
                    continue
            
            filtered_count += 1
            
            # Convert SMILES to InChI and InChIKey (with caching)
            inchi_str = None
            inchikey_str = None
            if smiles and isinstance(smiles, str) and smiles.strip():
                smiles_clean = smiles.strip()
                if smiles_clean in processed_smiles:
                    inchi_str, inchikey_str = processed_smiles[smiles_clean]
                else:
                    inchi_str, inchikey_str = self._smiles_to_inchi_inchikey(smiles_clean)
                    processed_smiles[smiles_clean] = (inchi_str, inchikey_str)
            
            # Create CPI entry
            cpi_entry = {
                'compound_pdsp_ligand_id': str(ligand_id),
                'compound_inchikey': inchikey_str,
                'target_uniprot_id': target_uniprot_id,
                'target_gene_name': unigene.strip(),
                'target_name': target_name,
                'species': species,
                'ki_nm': ki_value,
                'ki_note': ki_note,
                'reference': reference,
                'source': source,
                'hotligand': hotligand,
                'link': link
            }
            cpi_data.append(cpi_entry)
            
            # Create compound info entry (only once per unique SMILES)
            if smiles and not any(c.get('SMILES') == smiles for c in cinfo_data if c.get('SMILES')):
                cinfo_entry = {
                    'InChIKey': inchikey_str,
                    'SMILES': smiles,
                    'Name': ligand_name,
                    'InChI': inchi_str,
                    'CAS Number': cas_number,
                    'DrugBank ID': None,  # Not available in PDSP
                    'KEGG Compound ID': None,  # Not available in PDSP
                    'KEGG Drug ID': None,  # Not available in PDSP
                    'PubChem Compound ID': None,  # Not available in PDSP
                    'PubChem Substance ID': None,  # Not available in PDSP
                    'ChEBI ID': None,  # Not available in PDSP
                    'ChEMBL ID': None,  # Not available in PDSP
                    'HET ID': None,  # Not available in PDSP
                    'ChemSpider ID': None,  # Not available in PDSP
                    'BindingDB ID': None,  # Not available in PDSP
                    # PDSP-specific fields
                    'PDSP Ligand ID': str(ligand_id),
                    'NSC Number': nsc_number
                }
                cinfo_data.append(cinfo_entry)
            
            # Progress update
            if processed_count % 1000 == 0:
                speed = processed_count / (timer() - start)
                msg = prc_sym + "Processed (%d) entries, filtered (%d) GPCR entries. Speed: (%1.2f) entries/second" % (
                    processed_count, filtered_count, speed)
                print("\r" + msg, end="", flush=True)
        
        print(f"\rProcessed {processed_count} total entries, {filtered_count} GPCR entries")
        print(f"Converted {len(processed_smiles)} unique SMILES to InChI/InChIKey")
        
        # Log unmatched genes for debugging
        if unmatched_genes:
            print(f"Warning: {len(unmatched_genes)} unique gene names could not be mapped to GPCR UniProt IDs")
            if len(unmatched_genes) <= 20:  # Only print if not too many
                print(f"Unmatched genes: {', '.join(sorted(unmatched_genes))}")
        
        # Convert to DataFrames and save
        if cpi_data:
            cpi_df = pd.DataFrame(cpi_data)
            
            # Reorder columns to put compound_inchikey after compound_pdsp_ligand_id
            cols = list(cpi_df.columns)
            ligand_id_idx = cols.index('compound_pdsp_ligand_id')
            if 'compound_inchikey' in cols:
                cols.insert(ligand_id_idx + 1, cols.pop(cols.index('compound_inchikey')))
                cpi_df = cpi_df[cols]
            
            # Remove duplicates based on compound-target pairs
            cpi_df = cpi_df.drop_duplicates(subset=['compound_pdsp_ligand_id', 'target_uniprot_id'])
            
            # Save CPI data
            cpi_output_file = join(output_dp, 'pdsp_cpi.tsv')
            cpi_df.to_csv(cpi_output_file, sep='\t', index=False, encoding='utf-8')
            print(f"Saved {len(cpi_df)} CPI entries to {cpi_output_file}")
        else:
            print("No CPI data to save")
        
        if cinfo_data:
            cinfo_df = pd.DataFrame(cinfo_data)
            # Remove duplicates based on InChIKey (preferred) or SMILES if InChIKey is not available

            # if 'InChIKey' in cinfo_df.columns:
            #     # First try to deduplicate by InChIKey (for entries that have it)
            #     valid_inchikey_df = cinfo_df[cinfo_df['InChIKey'].notna()]
            #     invalid_inchikey_df = cinfo_df[cinfo_df['InChIKey'].isna()]
                
            #     if not valid_inchikey_df.empty:
            #         valid_inchikey_df = valid_inchikey_df.drop_duplicates(subset=['InChIKey'])
            #     if not invalid_inchikey_df.empty:
            #         invalid_inchikey_df = invalid_inchikey_df.drop_duplicates(subset=['SMILES'])
                
            #     cinfo_df = pd.concat([valid_inchikey_df, invalid_inchikey_df], ignore_index=True)
            # else:
            #     cinfo_df = cinfo_df.drop_duplicates(subset=['SMILES'])
            
            # Save compound info data
            cinfo_output_file = join(output_dp, 'pdsp_cinfo.tsv')
            cinfo_df.to_csv(cinfo_output_file, sep='\t', index=False, encoding='utf-8')
            print(f"Saved {len(cinfo_df)} compound entries to {cinfo_output_file}")
        else:
            print("No compound info data to save")
        
        print(done_sym + " Took %1.2f Seconds." % (timer() - start), flush=True)


class LLMParser:
    """
    LLM database parser
    """
    def __init__(self):
        self._filenames = [
            'llm_cpi.tsv',
            'llm_cinfo.tsv'
        ]
        self.filepath = ""
        self.output_dp = ""

    @property
    def filenames(self):
        return self._filenames

    def parse(self, source_dp, output_dp):
        print_section_header(
            "Parsing LLM files"
        )
        start = timer()

        # 蛋白质映射提取
        with open(join(source_dp, 'merged_protein_llm.json'), 'r', encoding='utf-8') as fd:
            protein_json = json.load(fd)

        protein_dict = {}

        for k, v in protein_json.items():
            is_gpcr = v.get('is_gpcr', 'no')
            if is_gpcr == 'yes' and v['uniprot_id'] is not None:
                protein_dict[sys.intern(k)] = v['uniprot_id']

        # cinfo 构建
        with open(join(source_dp, 'ligand.json'), 'r', encoding='utf-8') as fd:
            ligand_json = json.load(fd)
        
        ligands = []
        ligand_dict = {}
        for name, info in ligand_json.items():
            ligands.append({
                'InChIKey': info.get('inchikey'),
                'SMILES': info.get('smi'),
                'Name': name,
                'PubChem Compound ID': info.get('cid')
            })
            ligand_dict[sys.intern(name)] = info.get('inchikey')



        ligands_df = pd.DataFrame(ligands)
        ligands_df.to_csv(join(output_dp, 'llm_cinfo.tsv'), sep='\t', index=False, encoding='utf-8')

        # cpi 构建
        with open(join(source_dp, 'merged_fiter_triple.json'), 'r', encoding='utf-8') as fd:
            cpi_json = json.load(fd)
        

        def classify_standard_value(value):
            value = value.strip()
            # 1. 纯数字（整数或小数）
            if re.fullmatch(r'\d+(\.\d+)?', value):
                return 'pure_number'
            # 2. 符号 + 数字
            if re.fullmatch(r'(<=|>=|<|>|=|~)\s*\d+(\.\d+)?', value):
                return 'symbol_number'
            # 3. 范围形式（2.2~3、1-4、3 to 8.2）
            if re.fullmatch(r'\d+(\.\d+)?\s*[-~]\s*\d+(\.\d+)?', value):
                return 'range'
            if re.fullmatch(r'\d+(\.\d+)?\s+to\s+\d+(\.\d+)?', value, re.IGNORECASE):
                return 'range'
            return 'unknown'

        def split_symbol_number(value):
            value = value.strip()
            # 匹配符号和数字，比如 <1、>=5.2、~ 0.01
            match = re.fullmatch(r'(<=|>=|<|>|=|~)\s*(\d+(\.\d+)?)', value)
            if match:
                symbol = match.group(1)
                number = float(match.group(2))
                return symbol, number
            else:
                return None, None
        
        def split_range_number(value):
            value = value.strip()
            # 支持形如 "2.2~3", "1-4", "3 to 8.2" 的情况
            match = re.fullmatch(r'(\d+(\.\d+)?)\s*[-~]\s*(\d+(\.\d+)?)', value)
            if match:
                num1 = float(match.group(1))
                num2 = float(match.group(3))
                return num1, num2
            match = re.fullmatch(r'(\d+(\.\d+)?)\s+to\s+(\d+(\.\d+)?)', value, re.IGNORECASE)
            if match:
                num1 = float(match.group(1))
                num2 = float(match.group(3))
                return num1, num2
            return None, None
        
        def extract_nM_multiplier(s):
            s = s.strip()
            # 标准化：替换所有特殊字符
            s = s.replace("⁻", "-").replace("−", "-").replace("–", "-").replace("—", "-")
            s = s.translate(str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹", "0123456789"))
            # 去除多余空格（保留单位前的空格）
            s = re.sub(r'\s+', '', s)
            # 正则匹配各种形式： ×10^-10M, x10(-10)M, ×10−10M 等
            patterns = [
                r'[×x]10\^?(-?\d+)[mM]',         # ×10^-10M 或 x10^-10M 或 ×10−10M
                r'[×x]10\((-?\d+)\)[mM]',        # ×10(-10)M 或 x10(-10)M
            ]
            for pattern in patterns:
                match = re.fullmatch(pattern, s)
                if match:
                    exp = int(match.group(1))
                    return 10 ** (9 + exp)  # 换算为 nM 的倍率
            return None
            
        cpi_data = []

        processed_count = 0
        for paper, content in cpi_json.items():
            pmc_id = paper.replace('.txt', '')
            for triple in content:
                target = sys.intern(triple['target'])
                ligand = sys.intern(triple['ligand'])
                if target in protein_dict and ligand in ligand_dict:
                    # id 处理
                    target_uniprot_id = protein_dict[target]
                    compound_inchikey = ligand_dict[ligand]

                    # 单位处理
                    mult = None
                    standard_units = triple['unit']
                    if re.search(r'[\u4e00-\u9fff]', standard_units) is not None:
                        continue

                    if standard_units in ['μΜ', 'μℳ', 'uM', 'μmol', 'μm', 'um', 'uM']:
                        mult = 1000
                        standard_units = 'nM'
                    elif standard_units in ['nM', 'nmol', 'nm']:
                        mult = 1
                        standard_units = 'nM'
                    else:
                        mult = extract_nM_multiplier(standard_units)
                        if mult is not None:
                            standard_units = 'nM'

                    # 数值符号处理
                    standard_value = triple['affinity_value']
                    value_form = classify_standard_value(standard_value)
                    if value_form == 'pure_number':
                        standard_value = float(standard_value)
                        if mult is not None:
                            standard_value *= mult
                        standard_relation = '='
                    elif value_form == 'symbol_number':
                        standard_relation, standard_value = split_symbol_number(standard_value)
                        if mult is not None:
                            standard_value *= mult

                    elif value_form == 'range':
                        min_v, max_v = split_range_number(standard_value)
                        if mult is not None:
                            min_v *= mult
                            max_v *= mult
                        standard_value = f"({min_v}, {max_v})"
                        standard_relation = 'in'
                    else:
                        standard_relation = None

                    cpi_entry = {
                        'target_uniprot_id': target_uniprot_id,
                        'compound_inchikey': compound_inchikey,
                        'target_name': target,
                        'standard_type': triple['affinity_type'],
                        'standard_relation': standard_relation,
                        'standard_value': standard_value,
                        'standard_units': standard_units,
                        'reference': paper,
                        'pmc': pmc_id,
                        'confidence': triple['confidence']
                    }   
                    cpi_data.append(cpi_entry)
                # Progress update
                processed_count += 1
                if processed_count % 1000 == 0:
                    speed = processed_count / (timer() - start)
                    msg = prc_sym + "Processed (%d) entries, Speed: (%1.2f) entries/second" % (
                        processed_count, speed)
                    print("\r" + msg, end="", flush=True)

        cpi_df = pd.DataFrame(cpi_data)
        cpi_df.to_csv(join(output_dp, 'llm_cpi.tsv'), sep='\t', index=False, encoding='utf-8')
        print(done_sym + " Took %1.2f Seconds." % (timer() - start), flush=True)


class GLASSParser:
    """
    GLASS database parser
    """
    def __init__(self):
        self._filenames = [
            'glass_cpi.tsv',
            'glass_cinfo.tsv'
        ]
        self.filepath = ""
        self.output_dp = ""

    @property
    def filenames(self):
        return self._filenames
    
    def parse(self, source_dp, output_dp):
        print_section_header(
            "Parsing GLASS files" 
        )
        start = timer()

        def split_symbol_number(value):
            if isinstance(value, float):
                return '=', value
            value = value.strip()
            # 匹配符号和数字，比如 <1、>=5.2、~ 0.01
            match = re.fullmatch(r'(<=|>=|<|>|=|~)\s*(\d+(\.\d+)?)', value)
            if match:
                symbol = match.group(1)
                number = float(match.group(2))
                return symbol, number
            # 匹配纯数字（可有小数）
            match = re.fullmatch(r'(\d+(\.\d+)?)', value)
            if match:
                number = float(match.group(1))
                return '=', number
            return None, None

        def form_glass(df):
            df = df.copy()
            # 1. 拆分 'Value' 列
            df[['standard_relation', 'standard_value']] = df['Value'].apply(
                lambda x: pd.Series(split_symbol_number(x))
            )
            
            # 2. 重命名列
            df = df.rename(columns={
                'UniProt ID': 'target_uniprot_id',
                'InChI Key': 'compound_inchikey',
                'Parameter': 'standard_type',
                'Unit': 'standard_units',
                'Database Source': 'source_database',
                'Reference': 'reference'
            })

            
            
            # 3. 筛选并排序所需列
            cols = [
                'compound_inchikey',
                'target_uniprot_id',
                'source_database',
                'standard_type',
                'standard_relation',
                'standard_value',
                'standard_units',
                'reference'
            ]
            df = df[cols]
            df['pmid'] = ''
            def split_reference_and_pmid(ref_list):
                refs = []
                pmids = []
                for item in ref_list:
                    item = item.strip()
                    if item.isdigit():
                        pmids.append(item)
                    elif item:
                        refs.append(item)
                return ', '.join(refs) if refs else None, ', '.join(pmids) if pmids else None

            df['reference'] = df['reference'].fillna('').astype(str).str.split(',')
            df[['reference', 'pmid']] = df['reference'].apply(
                lambda lst: pd.Series(split_reference_and_pmid(lst))
            )
            return df

        # process all
        print("Processing GLASS cpi data...")
        df_all = pd.read_csv(join(source_dp, 'all.tsv'), sep='\t', encoding='utf-8')
        df_all = form_glass(df_all)
        # 区分是否是旧版 glass
        df_all['glass1'] = 'yes'
        df_all.to_csv(join(output_dp, 'glass_cpi.tsv'), sep='\t', index=False, encoding='utf-8')

        print("Processing GLASS cinfo data...")
        df_cinfo = pd.read_csv(join(source_dp, 'ligands.tsv'), sep='\t', encoding='utf-8')

        df_cinfo = df_cinfo.rename(columns={
            'InChI Key': 'InChIKey',
            'Canonical SMILES': 'SMILES',
            'InChI Std. ID': 'InChI',
            'Ligand Name': 'Name',
            'CID': 'PubChem Compound ID',

        })

        df_cinfo = df_cinfo[['InChIKey', 'SMILES', 'Name', 'InChI', 'PubChem Compound ID']]
        df_cinfo.to_csv(join(output_dp, 'glass_cinfo.tsv'), sep='\t', index=False, encoding='utf-8')

        print(done_sym + " Took %1.2f Seconds." % (timer() - start), flush=True)



        