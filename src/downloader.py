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

def download_unichem_data(sources_dp, srcs_cp):
    """ Download unichem database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading UniChem data files")
    reference = join(sources_dp, "reference.tsv.gz")
    download_file_md5_check(srcs_cp["unichem"]["reference"], reference)
    source = join(sources_dp, "source.tsv.gz")
    download_file_md5_check(srcs_cp["unichem"]["source"], source)
    structure = join(sources_dp, "structure.tsv.gz")
    download_file_md5_check(srcs_cp["unichem"]["structure"], structure)
    print_bold_line()


def download_glass_data(sources_dp, srcs_cp):
    """ Download GLASS database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading GLASS data files")
    cpi_fp = join(sources_dp, "cpi.csv")
    download_file_md5_check(srcs_cp["glass"]["cpi"], cpi_fp)
    ligands_fp = join(sources_dp, "ligands.tsv")
    download_file_md5_check(srcs_cp["glass"]["ligands"], ligands_fp)

    print_bold_line()

def download_llm_data(sources_dp, srcs_cp):
    """ Download llm database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading llm data files")
    mesh_fp = join(sources_dp, "mesh_supp.xml")
    download_file_md5_check(srcs_cp["ref"]["supp"], mesh_fp)
    cid2mesh_fp = join(sources_dp, "cid2mesh.txt")
    download_file_md5_check(srcs_cp["ref"]["cid2mesh"], cid2mesh_fp)

    if isfile(join("manually_prepared", "llm.jsonl")):
        copyfile(join("manually_prepared", "llm.jsonl"), join(sources_dp, "llm.jsonl"))
    else:
        raise FileNotFoundError("File llm.jsonl not found in manually_prepared folder,\
                                please check manually_prepared/README.md to prepare the data")
    print_bold_line()

def download_pubchem_data(sources_dp, srcs_cp):
    """ Download pubchem database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading PubChem data files")

    url = srcs_cp["pubchem"]["inchikey"]
    response = requests.get(url)
    response.raise_for_status()  # 确保请求成功
    # 使用 BeautifulSoup 解析 HTML
    soup = BeautifulSoup(response.text, "html.parser")
    # 提取所有 .gaf.gz 文件链接

    for link in soup.find_all("a"):
        href = link.get("href")
        if href and 'pc_inchikey2compound' in href:  # 筛选出 .gaf.gz 文件
            filename = href.split("/")[-1]
            save_path = join(sources_dp, filename)
            cur_url = join(url, href)
            download_file_md5_check(cur_url, save_path)
    
    title = join(sources_dp, "source.tsv.gz")
    download_file_md5_check(srcs_cp["pubchem"]["title"], title)
    iupac = join(sources_dp, "iupac.tsv.gz")
    download_file_md5_check(srcs_cp["pubchem"]["iupac"], iupac)
    print_bold_line()
