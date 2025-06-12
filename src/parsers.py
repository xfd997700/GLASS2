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


class UniProtTxtParser:
    """
    A UNIPROT database text file parser class
    """
    def __init__(self):
        self.filepath = ""
        self.output_dp = ""
        self._interpro_map = {}
        self._filenames = [
            "uniprot_facts.txt", 
            "uniprot_metadata.txt", 
            "uniprot_ppi.txt",
            "kinase_family.txt",
            "specie_map.json"
        ]
        self._specie_map = {}

    @property
    def filenames(self):
        return self._filenames

    def __parse_txt_entry(self, entry):
        """ Process a Uniprot txt entry

        Parameters
        ----------
        entry : list
            list of str lines representing the entry lines

        Returns
        -------
        list
            list of extracted facts
        """
        valid_species = True
        entry_facts = []
        entry_metadata = []
        entry_ppi = []
        # species_list = map(lambda x: x.node, VALID_SPECIES)
        entry_dictionary = dict()
        for line in entry:
            line_prefix = line[:2]
            if line_prefix == "  ":
                line_prefix = "AA"

            if line_prefix not in entry_dictionary:
                entry_dictionary[line_prefix] = ""

            entry_dictionary[line_prefix] += " " + line[2:].lstrip()

        # ------------------------------------------------------------------------
        # Processing IDs and UniProt ACs
        # ------------------------------------------------------------------------
        entry_code = entry_dictionary["AC"].strip().split(";")[0]
        row = entry_dictionary["AC"].split()
        if len(row) > 2:
            for idx in range(2, len(row)):
                entry_metadata.append([entry_code, "OTHER_ID", row[idx][:-1]])

        entry_name = entry_dictionary["ID"].strip().split(" ")[0].split("_")[0]
        entry_species = entry_dictionary["ID"].strip().split(" ")[0].split("_")[-1]

        entry_metadata.append([entry_code, "NAME", entry_name])
        entry_metadata.append([entry_code, "SPECIES", entry_species])
        self._specie_map[entry_code] = entry_species
        # ------------------------------------------------------------------------
        # Processing DE prefix section
        # ------------------------------------------------------------------------
        if "DE" in entry_dictionary:
            full_names = [v[:-1].strip().replace("Full=", "") for v in
                          re.findall("Full=[^;]+;", entry_dictionary["DE"])]
            short_names = [v[:-1].strip().replace("Short=", "") for v in
                           re.findall("Short=[^;]+;", entry_dictionary["DE"])]

            for full_n in full_names:
                entry_metadata.append([entry_code, "FULL_NAME", full_n])
            for short_n in short_names:
                entry_metadata.append([entry_code, "SHORT_NAME", short_n])

        # ------------------------------------------------------------------------
        # Processing OX prefix section
        # ------------------------------------------------------------------------
        if "OX" in entry_dictionary:
            organism_id = entry_dictionary['OX'].strip()
            organism_id = organism_id.split('=')[1]
            if organism_id.endswith(';'):
                organism_id = organism_id[0:-1]

            organism_id = organism_id.split(' ')[0]
            valid_species = True
        # ------------------------------------------------------------------------
        # Processing OC prefix section
        # ------------------------------------------------------------------------
        if "OC" in entry_dictionary:
            organism_classes = [c.strip() for c in entry_dictionary["OC"].strip().split(";")]
            for oc in organism_classes:
                entry_metadata.append([entry_code, "ORGANISM_CLASS", oc])

        # ------------------------------------------------------------------------
        # Processing RX prefix section
        # ------------------------------------------------------------------------
        if "RX" in entry_dictionary:
            pubmed_ids = ["pubmed:"+v.split("=")[-1] for v in re.findall(PUBMED_ID_CODE, entry_dictionary["RX"])]
            for pm in pubmed_ids:
                entry_metadata.append([entry_code, "RELATED_PUBMED_ID", pm])

        # ------------------------------------------------------------------------
        # Processing KW prefix section
        # ------------------------------------------------------------------------
        if "KW" in entry_dictionary:
            keywords = [v.strip() for v in entry_dictionary["KW"].replace(".", "").strip().split(";")]
            for kw in keywords:
                entry_metadata.append([entry_code, "RELATED_KEYWORD", kw])

        # ------------------------------------------------------------------------
        # Processing DR prefix section
        # ------------------------------------------------------------------------
        if "DR" in entry_dictionary:
            links_lines = entry_dictionary["DR"].strip().split(".")
            links_lines_dict = dict()
            for line in links_lines:
                db_name = line.strip().split(";")[0]
                if db_name not in links_lines_dict:
                    links_lines_dict[db_name] = []
                links_lines_dict[db_name].append(line.strip())

            if "GO" in links_lines_dict:
                go_lines = links_lines_dict["GO"]
                for line in go_lines:
                    go_code = line.split(";")[1].strip()
                    if "; F:" in line:
                        go_code_type = "GO_BP"
                    elif "; P:" in line:
                        go_code_type = "GO_MF"
                    else:
                        go_code_type = "GO_CC"
                    entry_facts.append([entry_code, go_code_type, go_code])

            if "HPA" in links_lines_dict:
                hpa_lines = links_lines_dict["HPA"]
                for line in hpa_lines:
                    hpa_code = line.split(";")[1].strip()
                    entry_facts.append([entry_code, "RELATED_ANTIBODY", hpa_code])

            if "Reactome" in links_lines_dict:
                reactome_lines = links_lines_dict["Reactome"]
                for line in reactome_lines:
                    reactome_code = line.split(";")[1].strip()
                    entry_facts.append([entry_code, "RELATED_PATHWAY", reactome_code])
        
            if "DrugBank" in links_lines_dict:
                drugbank_lines = links_lines_dict["DrugBank"]
                for line in drugbank_lines:
                    drugbank_code = line.split(";")[1].strip()
                    entry_facts.append([entry_code, "TARGET_OF_DRUG", drugbank_code])

            if "InterPro" in links_lines_dict:
                interpro_lines = links_lines_dict["InterPro"]
                for line in interpro_lines:
                    interpro_code = line.split(";")[1].strip()
                    if interpro_code in self._interpro_map:
                        entry_facts.append([entry_code, self._interpro_map[interpro_code], interpro_code])

            if 'PROSITE' in links_lines_dict:
                prosite_lines = links_lines_dict["PROSITE"]
                for line in prosite_lines:
                    prosite_code = line.split(";")[1].strip()
                    entry_facts.append([entry_code, "PS_SEQ_ANN", prosite_code])
            
        # # add string cross references   
        #     if "STRING" in links_lines_dict:
        #         string_lines = links_lines_dict["STRING"]
        #         for line in string_lines:
        #             string_code = line.strip()
        #             entry_facts.append([entry_code, "STRING_ID", string_code])
        # ------------------------------------------------------------------------
        # Processing OC prefix section
        # ------------------------------------------------------------------------
        if "CC" in entry_dictionary:
            comments_cats_dict = dict()
            comments_list = [v.strip() for v in entry_dictionary["CC"].strip().split("-!-") if len(v)>3]
            for comment in comments_list:
                comment_cat = comment[:comment.find(":")].strip()
                comment_val = comment[comment.find(":")+1:].strip()
                comments_cats_dict[comment_cat] = comment_val
            if "INTERACTION" in comments_cats_dict:
                interactors_uniprot_acs = re.findall(P_UNIPROT_CODE, comments_cats_dict["INTERACTION"])
                for up_id in interactors_uniprot_acs:
                    entry_ppi.append([entry_code, "INTERACTS_WITH", up_id])
            if "DISEASE" in comments_cats_dict:
                disease_codes = re.findall(P_DISEASE_CODE, comments_cats_dict["DISEASE"])
                for c in disease_codes:
                    entry_facts.append([entry_code, "RELATED_GENETIC_DISORDER", c])

        # ------------------------------------------------------------------------
        # Processing FT prefix section [Sequence annotations]
        # ------------------------------------------------------------------------
        if "FT" in entry_dictionary:
            ft_content = entry_dictionary["FT"]
            seq_ranges = [v.strip().split() for v in re.findall(SEQ_RANGE_CODE, ft_content)]
            seq_annotations = [v.strip().split() for v in re.findall(SEQ_NOTE__CODE, ft_content)]

        # ------------------------------------------------------------------------
        if valid_species:
            return entry_facts, entry_metadata, entry_ppi
        else:
            return [], [], []

    def _parse_pkinfam(self, pkinfam_fp, output_dp):
        pattern = re.compile(r'\((?P<acc>[A-Za-z0-9 ]*)\)')
        kinase_family = {}
        prev_line = ''
        header_started = False
        header_found = False
        header = ''
        with open(pkinfam_fp, 'r') as fd:
            for line in fd:
                if set(line.strip()) == set('='):
                    if header_started:
                        header_started = False
                        header_found = True
                        header = prev_line
                        kinase_family[header] = set()
                    else:
                        header_started = True
                        header_found = False
                elif header_found:
                    prev_line = line
                    accs = [x.strip() for x in pattern.findall(line)]
                    kinase_family[header].update(accs)
                else:
                    prev_line = line
        
        output_fd = open(join(output_dp, 'kinase_family.txt'), 'w')
        for family, kinase_set in kinase_family.items():
            if family.startswith('Atypical'):
                parts = family.strip().split()
                if parts[1] == 'PI3/PI4-kinase':
                    family_str = 'AT_PI3/PI4'
                else:
                    family_str = 'AT_'+parts[1]
            else:
                family_str = family.strip().split()[0]
            for kinase in kinase_set:
                output_fd.write(f'{kinase}\tKINASE_FAMILY\t{family_str}\n')
        output_fd.close()


    def parse(self, sources_dp, output_dp):
        """ Parse a Uniprot textual data file and output findings to a set of files in a specified directory

        Parameters
        ----------
        sources_dp : str
            absolute path to the sources directory

        output_dp : str
            absolute path of the output directory
        """
        filepath = join(sources_dp, "swissprot_entries.txt.gz")
        interpro_fp = join(sources_dp, "interpro_entries.txt")
        pkinfam_fp = join(sources_dp, 'pkinfam.txt')

        self.filepath = filepath
        self.output_dp = output_dp

        facts_fd = open(join(output_dp, "uniprot_facts.txt"), "w")
        metadata_fd = open(join(output_dp, "uniprot_metadata.txt"), "w")
        ppi_fd = open(join(output_dp, "uniprot_ppi.txt"), "w")

        line_index = 0
        nb_facts = 0
        nb_metadata = 0
        nb_ppi = 0
        self._parse_pkinfam(pkinfam_fp, output_dp)

        with open(interpro_fp, 'r') as fd:
            next(fd)
            for line in fd:
                interpro_id, interpro_type = line.strip().split('\t')[:2]
                self._interpro_map[interpro_id] = interpro_type.upper()
        with gzip.open(filepath, 'rt') as fd:
            current_entry = []
            print_section_header("Parsing Uniprot file (%s)" % (bcolors.OKGREEN + filepath + bcolors.ENDC))
            start = timer()
            eof = False
            while not eof:
                raw_line = fd.readline()
                line = raw_line.rstrip()
                if line != "//":
                    current_entry.append(line)
                else:
                    facts, metadata, ppi = self.__parse_txt_entry(current_entry)

                    nb_facts += len(facts)
                    nb_metadata += len(metadata)
                    nb_ppi += len(ppi)

                    export_triplets(facts, facts_fd)
                    export_triplets(metadata, metadata_fd)
                    export_triplets(ppi, ppi_fd)
                    facts_fd.flush()
                    metadata_fd.flush()
                    ppi_fd.flush()
                    current_entry = []

                line_index += 1
                if line_index % 5000 == 0:
                    speed = int(line_index / (timer()-start))
                    msg = prc_sym + "Processing (%d) lines/second => [facts:%d - metadata:%d - ppi:%d - total: %d]" \
                          % (speed, nb_facts, nb_metadata, nb_ppi, nb_facts + nb_metadata + nb_ppi)
                    print("\r" + msg, end="", flush=True)

                if raw_line == "":
                    eof = True
                    print(done_sym + " Took %1.2f Seconds." % (timer()-start), flush=True)

        facts_fd.close()
        metadata_fd.close()
        ppi_fd.close()

        with open(join(output_dp, 'specie_map.json'), 'w') as fd:
            json.dump(self._specie_map, fd)


class HumanProteinAtlasParser:
    """
    A Human Protein Atlas database parser
    """
    def __init__(self):
        """
        initialise new class instance
        """
        self.filepath = ""
        self.output_dp = ""
        # TODO: hpa_cellines_exp.txt 空的

    def __parse_hpa_xml_entry(self, entry):
        """ parse an xml element of an HPA entry

        Parameters
        ----------
        entry : xml.etree.ElementTree.Element
            xml element

        Returns
        -------
        dict
            dictionary of parsed data
        """
        parsed_data = {"rna": [], "tissue": [], "ab": []}

        antibody_list = entry.findall("antibody")
        tissue_expression_list = entry.findall("tissueExpression")
        rna_expression_list = entry.findall("rnaExpression")

        if tissue_expression_list is not None:
            for tissue_expression in tissue_expression_list:
                data_elements = tissue_expression.findall("data")
                for el in data_elements:
                    cell_line = el.find("cellLine")
                    if cell_line is not None:
                        cell_line = cell_line.text
                        cell_line_level = el.find("level").text
                        parsed_data["tissue"].append(["CELLINE", cell_line, cell_line_level])

                    tissue = el.find("tissue").text.replace(" ", "_")
                    tissue_level = el.find("level").text
                    t_cells = el.findall("tissueCell")
                    parsed_data["tissue"].append(["TISSUE", tissue, tissue_level])
                    for tc in t_cells:
                        cell_type = tc.find("cellType").text.replace(" ", "_")
                        cell_type_level = tc.find("level").text
                        parsed_data["tissue"].append(["TISSUE", tissue + "__" + cell_type, cell_type_level])
        if rna_expression_list is not None:
            for rna_expression in rna_expression_list:
                data_elements = rna_expression.findall("data")
                for el in data_elements:
                    cell_line = el.find("cellLine")
                    if cell_line is not None:
                        cell_line_name = cell_line.text.replace(" ", "_")
                        cell_line_level = "NormRNA_EXP:" + el[1].attrib["expRNA"]
                        cellosaurus_id = cell_line.attrib["cellosaurusID"]
                        cellosaurus_id = "NA" if cellosaurus_id == "" else cellosaurus_id
                        organ = cell_line.attrib["organ"]
                        suffix = "%s#%s" % (organ, cellosaurus_id)
                        parsed_data["rna"].append(["CELLINE", cell_line_name, cell_line_level, suffix])

                    tissue = el.find("tissue")
                    if tissue is not None:
                        tissue = tissue.text.replace(" ", "_")
                        tissue_level = el.find("level").text
                        t_cells = el.findall("tissueCell")
                        parsed_data["rna"].append(["RNA", tissue, tissue_level])
                        for tc in t_cells:
                            cell_type = tc.find("cellType").text.replace(" ", "_")
                            cell_type_level = tc.find("level").text
                            parsed_data["rna"].append(["RNA", tissue + "__" + cell_type, cell_type_level])

        if antibody_list is not None:
            for antibody_elem in antibody_list:
                ab_id = antibody_elem.attrib["id"]
                antigen = antibody_elem.find("antigenSequence").text
                parsed_data["ab"].append(["ANTIBODY", ab_id, antigen])
        return parsed_data

    def parse_database_xml(self, filepath, output_dp):
        """ Parse HPA xml file

        Parameters
        ----------
        filepath : str
            absolute file path of the hpa xml file
        output_dp : str
            path of the output directory
        """
        self.filepath = filepath
        self.output_dp = output_dp
        xml_fd = gzip.open(filepath)
        cl_exp_fp = join(output_dp, "hpa_cellines_exp.txt")
        tissue_exp_fp = join(output_dp, "hpa_tissues_exp.txt")
        ab_data_fp = join(output_dp, "hpa_antibodies.txt")

        cl_exp_fd = open(cl_exp_fp, "w")
        tissue_exp_fd = open(tissue_exp_fp, "w")
        ab_data_fd = open(ab_data_fp, "w")
        
        print_section_header("Parsing HPA XML file (%s)" % (bcolors.OKGREEN + filepath + bcolors.ENDC))
        start = timer()
        nb_entries = 0
        for event, entry in ET.iterparse(xml_fd, events=('start', 'end')):
            if entry.tag == "entry" and event == "end" and len(list(entry)) > 2:
                nb_entries += 1
                if nb_entries % 5 == 0:
                    speed = nb_entries / (timer() - start)
                    msg = prc_sym + "Processed (%d) entries.  Speed: (%1.5f) entries/second" % (nb_entries, speed)
                    print("\r" + msg, end="", flush=True)

                id_elem = entry.find("identifier")
                if len(id_elem) >= 1:
                    db_name = id_elem[0].attrib["db"]
                    if db_name == "Uniprot/SWISSPROT":
                        entry_id = id_elem[0].attrib["id"]
                        entry_data = self.__parse_hpa_xml_entry(entry)
                        # export antibody data
                        for _, ab, ag in entry_data["ab"]:
                            ag = ag if ag is not None else "-"
                            ab_data_fd.write("%s\t%s\t%s\n" % (ab, entry_id, ag))
                        # export tissue data
                        for context, ts, level in entry_data["tissue"]:
                            if level is not None:
                                if context == "TISSUE":
                                    tissue_exp_fd.write("%s\t%s\t%s\n" % (entry_id, ts, level))
                                else:
                                    cl_exp_fd.write("%s\t%s\t%s\n" % (entry_id, ts, level))
                        # export rna data
                        for rna_data in entry_data["rna"]:
                            if len(rna_data) == 3:
                                _, cl, level = rna_data
                                tissue_exp_fd.write("%s\t%s\t%s\n" % (entry_id, cl, level))
                            if len(rna_data) == 4:
                                _, cl, level, organ = rna_data
                                cl_exp_fd.write("%s\t%s\t%s\t%s\n" % (entry_id, cl, organ, level))
                entry.clear()

        print(done_sym + " Took %1.2f Seconds." % (timer() - start), flush=True)
        ab_data_fd.close()
        tissue_exp_fd.close()
        cl_exp_fd.close()

