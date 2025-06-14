# -*- coding: utf-8 -*-
from asyncio import protocols
import gzip
from itertools import chain
import re
import sys
from turtle import pos
import zipfile
import requests
import time
import gc
from os.path import join
import numpy as np
from reportlab.lib.pdfencrypt import checkU
from tqdm import tqdm
from biodblinker import KEGGLinker, GeneNameLinker
from Bio import SeqIO
import io
from .extras import *
from timeit import default_timer as timer
import xml.etree.ElementTree as ET
from zipfile import ZipFile
from collections import defaultdict
import pandas as pd
import tarfile
import os
import sqlite3
import csv
import shutil
import json

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

DDI_SIDE_EFFECT_1 = re.compile(r'The risk or severity of (?P<se>.*) can be (?P<mode>\S+)d when .* is combined with .*')
DDI_SIDE_EFFECT_2 = re.compile(r'.* may (?P<mode>\S+) (?P<se>\S+\s?\w*\s?\w*) of .* as a diagnostic agent.')
DDI_SIDE_EFFECT_3 = re.compile(r'The (?P<se>\S+\s?\w*\s?\w*) of .* can be (?P<mode>\S+)d when used in combination with .*')
DDI_SIDE_EFFECT_4 = re.compile(r'The (?P<se>\S+\s?\w*\s?\w*) of .* can be (?P<mode>\S+)d when it is combined with .*')
DDI_SIDE_EFFECT_5 = re.compile(r'.* can cause a decrease in the absorption of .* resulting in a (?P<mode>\S+) (?P<se>\S+\s?\w*\s?\w*) and potentially a decrease in efficacy.')
DDI_SIDE_EFFECT_6 = re.compile(r'.* may decrease the excretion rate of .* which could result in a (?P<mode>\S+) (?P<se>\S+\s?\w*\s?\w*).')
DDI_SIDE_EFFECT_7 = re.compile(r'.* may increase the excretion rate of .* which could result in a (?P<mode>\S+) (?P<se>\S+\s?\w*\s?\w*) and potentially a reduction in efficacy.')
DDI_SIDE_EFFECT_8 = re.compile(r'The (?P<se>\S+\s?\w*\s?\w*) of .* can be (?P<mode>\S+)d when combined with .*')
DDI_SIDE_EFFECT_9 = re.compile(r'.* can cause an increase in the absorption of .* resulting in an (?P<mode>\S+)d (?P<se>\S+\s?\w*\s?\w*) and potentially a worsening of adverse effects.')
DDI_SIDE_EFFECT_10 = re.compile(r'The risk of a (?P<se>\S+\s?\w*\s?\w*) to .* is (?P<mode>\S+)d when it is combined with .*')
DDI_SIDE_EFFECT_11 = re.compile(r'The (?P<se>\S+\s?\w*\s?\w*) of .* can be (?P<mode>\S+)d when combined with .*')
DDI_SIDE_EFFECT_12 = re.compile(r'The (?P<se>\S+\s?\w*\s?\w*) of the active metabolites of .* can be (?P<mode>\S+)d when .* is used in combination with .*')
DDI_SIDE_EFFECT_13 = re.compile(r'The (?P<se>\S+\s?\w*\s?\w*) of .*, an active metabolite of .* can be (?P<mode>\S+)d when used in combination with .*')
DDI_SIDE_EFFECT_14 = re.compile(r'.* may (?P<mode>\S+) the (?P<se>.*) of .*')
DDI_SIDE_EFFECT_15 = re.compile(r'.* may (?P<mode>\S+) the central nervous system depressant (?P<se>\S+\s?\S*\s?\S*) of .*')

DDI_SIDE_EFFECTS = [
    DDI_SIDE_EFFECT_1, DDI_SIDE_EFFECT_2, DDI_SIDE_EFFECT_3, DDI_SIDE_EFFECT_4,
    DDI_SIDE_EFFECT_5, DDI_SIDE_EFFECT_6, DDI_SIDE_EFFECT_7, DDI_SIDE_EFFECT_8,
    DDI_SIDE_EFFECT_9, DDI_SIDE_EFFECT_10, DDI_SIDE_EFFECT_11, DDI_SIDE_EFFECT_12,
    DDI_SIDE_EFFECT_13, DDI_SIDE_EFFECT_14, DDI_SIDE_EFFECT_15
]

DDI_MODE_MAP = {
    'reduced': "decrease",
    'increase': "increase",
    'higher': "increase",
    'decrease': "decrease",
    'reduce': "decrease",
    'lower': "decrease"
}

DDI_SE_NAME_MAP = {
    "central_nervous_system_depressant_(cns_depressant)_activities": 'cns_depression_activities',
    "(cns_depressant)_activities": 'cns_depression_activities',
    "cns_depression": 'cns_depression_activities',
    "cardiotoxic_activities": 'cardiotoxicity',
    "constipating_activities": 'constipation',
    "excretion": 'excretion_rate',
    "hyperkalemic_activities": 'hyperkalemia',
    "hypertensive_activities": 'hypertension',
    "qtc-prolonging_activities": "qtc_prolongation",
    "tachycardic_activities": "tachycardia",
    "hypokalemic_activities": "hypokalemia",
    "hypoglycemic_activities": "hypoglycemia",
    "hypercalcemic_activities": "hypercalcemia",
    "bradycardic_activities": "bradycardia",
    "neutropenic_activities": "neutropenia",
    "orthostatic_hypotensive_activities": "orthostatic_hypotension",
    "neutropenic_activities": "neutropenia",
    "pseudotumor_cerebri_activities": "pseudotumor_cerebri",
    "sedative_activities": "sedation",
    "ototoxic_activities": "ototoxicity",
    "neuromuscular_blocking_activities": "neuromuscular_blockade",
    "nephrotoxic_activities": "nephrotoxicity",
    "myelosuppressive_activities": "myelosuppression",
    "hypotensive_activities": "hypotension",
    "serum_level": "serum_concentration"
}

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



# ...existing code...

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
            # action_comment = row.get('Action comment', None).strip()
            # selectivity = row.get('Selectivity', None).strip()
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

            # concentration_range = row.get('concentration Range', '').strip()
            
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
        
        print(f"\rProcessed {processed_count} total entries, {filtered_count} GPCR entries")
        
        # Convert to DataFrame and save CPI data
        cpi_df = None
        if cpi_data:
            cpi_df = pd.DataFrame(cpi_data)
            # cpi_df = regular_type(cpi_df)  # Apply type regularization
            
            # Save CPI data
            cpi_output_file = join(output_dp, 'iuphar_cpi.tsv')
            cpi_df.to_csv(cpi_output_file, sep='\t', index=False, encoding='utf-8')
            print(f"Saved {len(cpi_df)} CPI entries to {cpi_output_file}")
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
            
            for idx, row in filtered_ligand_df.iterrows():
                # Map IUPHAR columns to link_head format
                cinfo_entry = {
                    'InChIKey': row.get('InChIKey', None),
                    'SMILES': row.get('SMILES', None),
                    'Name': row.get('Name', None),
                    'InChI': row.get('InChI', None),
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
                    'Ligand ID': row.get('Ligand ID', None),
                    'Species': row.get('Species', None),
                    'Type': row.get('Type', None),
                    'Approved': row.get('Approved', None),
                    'IUPAC name': row.get('IUPAC name', None),
                }
                cinfo_data.append(cinfo_entry)
            
            # Create DataFrame and save
            if cinfo_data:
                cinfo_df = pd.DataFrame(cinfo_data)
                
                # Save compound info data
                cinfo_output_file = join(output_dp, 'iuphar_cinfo.tsv')
                cinfo_df.to_csv(cinfo_output_file, sep='\t', index=False, encoding='utf-8')
            else:
                print("No compound info data to save")
        
            
        except Exception as e:
            print(f"Error processing ligand_list.csv: {e}")
            return
        # Add InChIKey mapping to CPI data

        ligand_inchikey_map = dict(zip(cinfo_df['Ligand ID'], cinfo_df['InChIKey']))
    
        # Add compound_inchikey column to CPI data
        cpi_df['compound_inchikey'] = cpi_df['ligand_id'].map(ligand_inchikey_map)
        
        # Reorder columns to put compound_inchikey after ligand_id
        cols = list(cpi_df.columns)
        ligand_id_idx = cols.index('ligand_id')
        cols.insert(ligand_id_idx + 1, cols.pop(cols.index('compound_inchikey')))
        cpi_df = cpi_df[cols]
        
        # Save updated CPI data
        cpi_output_file = join(output_dp, 'iuphar_cpi.tsv')
        cpi_df.to_csv(cpi_output_file, sep='\t', index=False, encoding='utf-8')