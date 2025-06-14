import os
from configparser import RawConfigParser

from bs4 import BeautifulSoup
import requests
from .file_io import download_file_md5_check
from .extras import print_section_header, print_bold_line
from os.path import join
import sys
from shutil import copyfile
from os.path import join, isfile
# from .parsers import VALID_SPECIES

def download_uniprot_files(sources_dp, srcs_cp):
    """ Download uniprot gpcr files

    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading UniProt GPCR data files")
    gpcr_entries_fp = join(sources_dp, "gpcr_entries.txt")
    download_file_md5_check(srcs_cp["uniprot"]["swissprot_entries"], gpcr_entries_fp)
    gpcr_seq_fp = join(sources_dp, "gpcr_sequences.fasta")
    download_file_md5_check(srcs_cp["uniprot"]["swissprot_fasta"], gpcr_seq_fp)

    print_bold_line()

def download_iuphar_data(sources_dp, srcs_cp):
    """ Download iuphar gpcr files

    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading IUPHAR GPCR data files")
    iuphar_fp = join(sources_dp, "gpcr_interactions.csv")
    download_file_md5_check(srcs_cp["iuphar"]["gpcr_interactions"], iuphar_fp)
    ligand_fp = join(sources_dp, "ligand_list.csv")
    download_file_md5_check(srcs_cp["iuphar"]["ligand_list"], ligand_fp)
    print_bold_line()

def download_chembl_data(sources_dp, srcs_cp):
    """ Download chembl database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading ChEMBL data files")
    chembl_sql = join(sources_dp, "chembl_sqlite.tar.gz")
    download_file_md5_check(srcs_cp["chembl"]["chembl_sql"], chembl_sql)
    print_bold_line()

def download_drugbank_data(sources_dp, srcs_cp, username, password):
    """ Download drugbank database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading Drugbank data files")

    full_database_fp = join(sources_dp, "drugbank_all_full_database.xml.zip")
    download_file_md5_check(
        srcs_cp["drugbank"]["drugbank_all_full_database"],
        full_database_fp,
        username = username,
        password = password
    )

    download_file_md5_check(
        srcs_cp["drugbank"]["external_drug"],
        join(sources_dp, "external_drug.zip"),
        username=username,
        password=password
    )

    print_bold_line()

def download_bindingdb_data(sources_dp, srcs_cp):
    """ Download bindingdb database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Moving BindingDB data files")
    if isfile(join("manually_prepared", "BindingDB_All.zip")):
        copyfile(join("manually_prepared", "BindingDB_All.zip"), join(sources_dp, "BindingDB_All.zip"))
    else:
        raise FileNotFoundError("File BindingDB_All.zip not found in manually_prepared folder,\
                                please check manually_prepared/README.md to prepare the data")
    print_bold_line()

def download_pdsp_data(sources_dp, srcs_cp):
    """ Download PDSP file

    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)

    if isfile(join("manually_prepared", "KiDatabase.csv")):
        copyfile(join("manually_prepared", "KiDatabase.csv"), join(sources_dp, 'ki_data.csv'))
    else:
        raise FileNotFoundError("File KiDatabase.csv not found in manually_prepared folder,\
                                please check manually_prepared/README.md to prepare the data")
    print_bold_line()
    print_bold_line()
