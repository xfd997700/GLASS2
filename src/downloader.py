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

def download_iuphar_data(source_dp, srcs_cp):
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
    iuphar_fp = join(sources_dp, "gpcr_interactions.txt")
    download_file_md5_check(srcs_cp["iuphar"]["gpcr_interactions"], iuphar_fp)
    print_bold_line()



def download_cellosaurus_data(sources_dp, srcs_cp):
    """ Download cellosaurus database files

    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading Cellosaurus data files")
    cellosaurus_data_fp = join(sources_dp, "cellosaurus_data.txt")
    download_file_md5_check(srcs_cp["cellosaurus"]["cellosaurus_txt"], cellosaurus_data_fp)
    print_bold_line()


def download_hpa_data(sources_dp, srcs_cp):
    """ Download cellosaurus database files

    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading Human Protein Atlas data files")
    hpa_data_fp = join(sources_dp, "proteinatlas.xml.gz")
    normal_tissues_fp = join(sources_dp, "normal_tissues.tsv.gz")
    rna_cellines_fp = join(sources_dp, "rna_cellines.tsv.gz")
    pathology_fp = join(sources_dp, "pathology.tsv.gz")

    download_file_md5_check(srcs_cp["hpa"]["xml_database"], hpa_data_fp)
    download_file_md5_check(srcs_cp["hpa"]["normal_tissues"], normal_tissues_fp)
    download_file_md5_check(srcs_cp["hpa"]["rna_cellines"], rna_cellines_fp)
    download_file_md5_check(srcs_cp["hpa"]["pathology"], pathology_fp)

    print_bold_line()


def download_reactome_data(sources_dp, srcs_cp):
    """ Download reactome database files

    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading Reactome data files")

    reactome_ppi_fp = join(sources_dp, "reactome_ppi.txt")
    reactome_pathway_rels_fp = join(sources_dp, "reactome_pathway_rels.txt")
    reactome_protein_complex_rels_fp = join(sources_dp, "reactome_protein_complex_rels.txt")
    reactome_complex_pathway_rels_fp = join(sources_dp, "reactome_complex_pathway_rels.txt")
    reactome_go_mapping_fp = join(sources_dp, "reactome_go_mapping.txt.gz")
    reactome_omim_mapping_fp = join(sources_dp, "reactome_omim_mapping.txt")
    reactome_pathway_list_fp = join(sources_dp, "reactome_pathway_list.txt")
    reactome_protein_pathway_rels_fp = join(sources_dp, "reactome_protein_pathway_rels.txt")
    reactome_go_rels_fp = join(sources_dp, "reactome_go_rels.txt")
    download_file_md5_check(srcs_cp["reactome"]["reactome_ppi"], reactome_ppi_fp)
    download_file_md5_check(srcs_cp["reactome"]["reactome_pathway_rels"], reactome_pathway_rels_fp)
    download_file_md5_check(srcs_cp["reactome"]["reactome_protein_complex_rels"], reactome_protein_complex_rels_fp)
    download_file_md5_check(srcs_cp["reactome"]["reactome_complex_pathway_rels"], reactome_complex_pathway_rels_fp)
    download_file_md5_check(srcs_cp["reactome"]["reactome_go_mapping"], reactome_go_mapping_fp)
    download_file_md5_check(srcs_cp["reactome"]["reactome_omim_mapping"], reactome_omim_mapping_fp)
    download_file_md5_check(srcs_cp["reactome"]["reactome_pathway_list"], reactome_pathway_list_fp)
    download_file_md5_check(srcs_cp["reactome"]["reactome_protein_pathway_rels"], reactome_protein_pathway_rels_fp)
    download_file_md5_check(srcs_cp["reactome"]["reactome_go_rels"], reactome_go_rels_fp)
    print_bold_line()

def download_ctd_data(sources_dp, srcs_cp):
    """ Download ctd database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading CTD data files")

    ctd_chemical_id_mapping = join(sources_dp,'CTD_chemicals.tsv.gz')
    ctd_gene_id_mapping =join(sources_dp,'CTD_genes.tsv.gz')
    ctd_chemical_gene_interactions = join(sources_dp,'CTD_chem_gene_ixns.tsv.gz')
    ctd_chemical_disase_association =join(sources_dp,'CTD_chemicals_diseases.tsv.gz')
    ctd_chemical_pathway_association = join(sources_dp,'CTD_chem_pathways_enriched.tsv.gz')
    ctd_gene_disease_association = join(sources_dp,'CTD_genes_diseases.tsv.gz')
    ctd_gene_pathway_association = join(sources_dp,'CTD_genes_pathways.tsv.gz')
    ctd_disease_pathway_assiciation = join(sources_dp,'CTD_diseases_pathways.tsv.gz')
    ctd_disease_molecular_function = join(sources_dp,'CTD_disease_molecular_function.tsv.gz')
    ctd_disease_cellular_component = join(sources_dp,'CTD_disease_cellular_component.tsv.gz')
    ctd_disease_biological_process = join(sources_dp,'CTD_disease_biological_process.tsv.gz')
    ctd_chemical_phenotype = join(sources_dp,'CTD_chemical_phenotype.tsv.gz')

    try:
        download_file_md5_check(srcs_cp["ctd"]["gene_disease_association"], ctd_gene_disease_association, bypass_md5=True)
    except Exception as _:
        print()
        print(f'Unable to download {srcs_cp["ctd"]["gene_disease_association"]}')
        print('Please download manually through a browser and place in the data/sources folder')
        sys.exit(-1)

    download_file_md5_check(srcs_cp["ctd"]["chemical_id_mapping"], ctd_chemical_id_mapping)
    download_file_md5_check(srcs_cp["ctd"]["gene_id_mapping"], ctd_gene_id_mapping)
    download_file_md5_check(srcs_cp["ctd"]["chemical_gene_interactions"], ctd_chemical_gene_interactions)
    download_file_md5_check(srcs_cp["ctd"]["chemical_disase_association"], ctd_chemical_disase_association)
    download_file_md5_check(srcs_cp["ctd"]["chemical_pathway_association"], ctd_chemical_pathway_association)
    download_file_md5_check(srcs_cp["ctd"]["gene_pathway_association"], ctd_gene_pathway_association)
    download_file_md5_check(srcs_cp["ctd"]["disease_pathway_assiciation"], ctd_disease_pathway_assiciation)
    download_file_md5_check(srcs_cp["ctd"]["disease_molecular_function"], ctd_disease_molecular_function)
    download_file_md5_check(srcs_cp["ctd"]["disease_cellular_component"], ctd_disease_cellular_component)
    download_file_md5_check(srcs_cp["ctd"]["disease_biological_process"], ctd_disease_biological_process)
    download_file_md5_check(srcs_cp["ctd"]["chemical_phenotype"], ctd_chemical_phenotype)

    # the following file should be prepared manually on https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi
    manual_files = ["CIDs.txt.gz", "InChIKeys.txt.gz", "InChIs.txt.gz", "SMILES.txt.gz", "ctd_name.txt"]
    for manual_file in manual_files:
        if isfile(join("manually_prepared", manual_file)):
            copyfile(join("manually_prepared", manual_file), join(sources_dp, manual_file))
        else:
            raise FileNotFoundError(f"File {manual_file} not found in manually_prepared folder,\
                                    please check manually_prepared/README.md to prepare the data")
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


def download_phosphosite_data(sources_dp, srcs_cp):
    """ Download phosphositeplus database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading PhosphoSitePlus data files")

    phosphorylation_site_fp = join(sources_dp, "phosphorylation_site.txt.gz")
    kinase_substrate_fp = join(sources_dp, "kinase_substrate.txt.gz")

    download_file_md5_check(srcs_cp["phosphositeplus"]["phosphorylation_site"], phosphorylation_site_fp)
    download_file_md5_check(srcs_cp["phosphositeplus"]["kinase_substrate"], kinase_substrate_fp)

    print_bold_line()


def download_intact_data(sources_dp, srcs_cp):
    """ Download intact database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading Intact data files")

    intact_zip_fp = join(sources_dp, "intact.zip")

    download_file_md5_check(srcs_cp["intact"]["intact_zip"], intact_zip_fp)

    print_bold_line()


def download_sider_data(sources_dp, srcs_cp):
    """ Download sider database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading SIDER data files")

    sider_interactions_fp = join(sources_dp, "sider_interactions.tsv.gz")
    sider_side_effects_fp = join(sources_dp, "sider_side_effects.tsv.gz")

    download_file_md5_check(srcs_cp["sider"]["indications"], sider_interactions_fp)
    download_file_md5_check(srcs_cp["sider"]["side_effects"], sider_side_effects_fp)

    print_bold_line()


def download_kegg_data(sources_dp, srcs_cp):
    """ Download kegg database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading KEGG data files")

    disease_fp = join(sources_dp, "diseases.txt")

    download_file_md5_check(srcs_cp["kegg"]["diseases"], disease_fp)

    print_bold_line()


def download_mesh_data(sources_dp, srcs_cp):
    """ Download mesh database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading MeSH data files")

    disease_fp = join(sources_dp, "mesh_diseases.txt")
    supp_fp = join(sources_dp, "mesh_supp_concepts.xml")
    download_file_md5_check(srcs_cp["mesh"]["mesh_diseases"], disease_fp)
    download_file_md5_check(srcs_cp["mesh"]["mesh_supp_concepts"], supp_fp)
    print_bold_line()


def download_medgen_data(sources_dp, srcs_cp):
    """ Download medgen database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading MedGen data files")

    disease_fp = join(sources_dp, "medgen_omim_mappings.txt.gz")
    download_file_md5_check(srcs_cp["medgen"]["medgen_omim_mappings"], disease_fp)
    disease_fp = join(sources_dp, "medgen_mappings.txt.gz")
    download_file_md5_check(srcs_cp["medgen"]["medgen_mappings"], disease_fp)
    print_bold_line()


def download_hijazi20_data(sources_dp, srcs_cp):
    """ Download Hijazi20 database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading Hijazi20 data files")

    pdts_fp = join(sources_dp, "hijazi20.xlsm")
    download_file_md5_check(srcs_cp["hijazi20"]["hijazi20_pdt"], pdts_fp)
    print_bold_line()


def download_smpdb_data(sources_dp, srcs_cp):
    """ Download Cutillas20 database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading Cutillas data files")

    pathways_fp = join(sources_dp, "smpdb_pathways.csv.zip")
    download_file_md5_check(srcs_cp["smpdb"]["smpdb_pathways"], pathways_fp)

    metabolites_fp = join(sources_dp, "smpdb_metabolites.csv.zip")
    download_file_md5_check(srcs_cp["smpdb"]["smpdb_metabolites"], metabolites_fp)

    protein_fp = join(sources_dp, "smpdb_proteins.csv.zip")
    download_file_md5_check(srcs_cp["smpdb"]["smpdb_proteins"], protein_fp)

    biopax_fp = join(sources_dp, "smpdb_biopax.zip")
    download_file_md5_check(srcs_cp["smpdb"]["smpdb_biopax"], biopax_fp)
    print_bold_line()


def download_kegg_data(sources_dp, srcs_cp):
    """ Download kegg database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading KEGG data files")
    disease_fp = join(sources_dp, "diseases.txt")
    download_file_md5_check(srcs_cp["kegg"]["diseases"], disease_fp)
    disease_fp = join(sources_dp, "map.txt")
    download_file_md5_check(srcs_cp["kegg"]["map"], disease_fp)
    disease_fp = join(sources_dp, "hsa.txt")
    download_file_md5_check(srcs_cp["kegg"]["hsa"], disease_fp)
    print_bold_line()


def download_foodb_data(sources_dp, srcs_cp):
    """ Download foodb database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading FOODB data files")
    foodb_fp = join(sources_dp, "foodb.tar.gz")
    download_file_md5_check(srcs_cp["foodb"]["foodb_csv"], foodb_fp)
    print_bold_line()

def download_tcmbank_data(sources_dp, srcs_cp):
    """ Download tcmbank database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading TCMBank data files")

    target_fp = join(sources_dp, "gene.xlsx")
    disease_fp = join(sources_dp, "disease.xlsx")
    ingredient_fp = join(sources_dp, "ingredient.xlsx")
    herb_fp = join(sources_dp, "herb.xlsx")

    download_file_md5_check(srcs_cp["tcmbank"]["tcm_target"], target_fp)
    download_file_md5_check(srcs_cp["tcmbank"]["tcm_ingredient"], ingredient_fp)
    download_file_md5_check(srcs_cp["tcmbank"]["tcm_herb"], herb_fp)
    download_file_md5_check(srcs_cp["tcmbank"]["tcm_disease"], disease_fp)

    print_bold_line()

def download_hgnc_data(sources_dp, srcs_cp):
    """ Download hgnc database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading HGNC data files")
    hgnc_complete_set = join(sources_dp, "hgnc_complete_set.txt")
    download_file_md5_check(srcs_cp["hgnc"]["hgnc_complete_set"], hgnc_complete_set)
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

def download_hmdb_data(sources_dp, srcs_cp):
    """ Download hmdb database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading HMDB data files")
    metabolites_fp = join(sources_dp, "hmdb_metabolites.zip")
    download_file_md5_check(srcs_cp["hmdb"]["all_metabo"], metabolites_fp)

    protein_fp = join(sources_dp, "hmdb_proteins.zip")
    download_file_md5_check(srcs_cp["hmdb"]["all_protein"], protein_fp)
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

def download_go_data(sources_dp, srcs_cp):
    """ Download GO database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading GO data files")
    go_all = join(sources_dp, "ontology.json")
    download_file_md5_check(srcs_cp["go"]["go_ontology"], go_all)

    go_annotation = join(sources_dp, "goa")
    os.makedirs(go_annotation, exist_ok=True)

    print_section_header("Moving QuickGO data files")
    # the following file should be prepared manually on https://pubchem.ncbi.nlm.nih.gov/idexchange/idexchange.cgi
    manual_files = ["MF.gpad", "CC.gpad", "BP.gpad"]
    for manual_file in manual_files:
        if isfile(join("manually_prepared", manual_file)):
            copyfile(join("manually_prepared", manual_file), join(sources_dp, manual_file))
        else:
            raise FileNotFoundError(f"File {manual_file} not found in manually_prepared folder,\
                                    please check manually_prepared/README.md to prepare the data")
        

    # url = srcs_cp["go"]["go_annotation"]
    # response = requests.get(url)
    # response.raise_for_status()  # 确保请求成功
    # # 使用 BeautifulSoup 解析 HTML
    # soup = BeautifulSoup(response.text, "html.parser")
    # # 提取所有 .gaf.gz 文件链接
    # gaf_links = []
    # for link in soup.find_all("a"):
    #     href = link.get("href")
    #     if href and href.endswith(".gaf.gz"):  # 筛选出 .gaf.gz 文件
    #         filename = href.split("/")[-1]
    #         if not filename.startswith("filtered_goa_uniprot"):
    #             goa_path = join(go_annotation, filename)
    #             download_file_md5_check(href, goa_path)
    # print_bold_line()

def download_string_data(sources_dp, srcs_cp, version="12.0"):
    """ Download STRING database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading STRING data files")
    base_url = srcs_cp["string"]["base_url"]
    for specie in VALID_SPECIES:
        node = specie.node
        specie_links = join(sources_dp, f"{specie.code}.txt.gz")
        url = base_url + f"/protein.links.detailed.v{version}/{node}.protein.links.detailed.v{version}.txt.gz"
        try:
            download_file_md5_check(url, specie_links)
        except:
            print(f"{specie.scientific_name} - {specie.node} NOT FOUND IN STRING, SKIPPING")

    aliases = join(sources_dp, "aliases.txt.gz")
    aliases_link = srcs_cp["string"]["base_url"] + f"/protein.aliases.v{version}.txt.gz"
    download_file_md5_check(aliases_link, aliases)

def download_tcdb_data(sources_dp, srcs_cp):
    """ Download tcdb database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading TCDB data files")
    substrate = join(sources_dp, "substrate.txt")
    download_file_md5_check(srcs_cp["tcdb"]["substrate_map"], substrate)
    protein = join(sources_dp, "protein.txt")
    download_file_md5_check(srcs_cp["tcdb"]["protein_map"], protein)
    go = join(sources_dp, "go.txt")
    download_file_md5_check(srcs_cp["tcdb"]["go_map"], go)
    pfam = join(sources_dp, "pfam.txt")
    download_file_md5_check(srcs_cp["tcdb"]["pfam_map"], pfam)
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
    print_bold_line()

def download_chebi_data(sources_dp, srcs_cp):
    """ Download chebi database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading ChEBI data files")
    sdf = join(sources_dp, "ChEBI_complete.sdf.gz")
    download_file_md5_check(srcs_cp["chebi"]["sdf"], sdf)
    obo = join(sources_dp, "chebi_lite.obo")
    download_file_md5_check(srcs_cp["chebi"]["obo"], obo)

    print_bold_line()


def download_hpo_data(sources_dp, srcs_cp):
    """ Download hpo database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading HPO data files")
    hpoa_fp = join(sources_dp, "phenotype.hpoa")
    download_file_md5_check(srcs_cp["hpo"]["annotation"], hpoa_fp)
    
    ontology_fp = join(sources_dp, "hp.obo")
    download_file_md5_check(srcs_cp["hpo"]["ontology"], ontology_fp)

    g2p_fp = join(sources_dp, "genes_to_phenotype.txt")
    download_file_md5_check(srcs_cp["hpo"]["genes2phenotype"], g2p_fp)

    g2d_fp = join(sources_dp, "genes_to_disease.txt")
    download_file_md5_check(srcs_cp["hpo"]["genes2disease"], g2d_fp)
    print_bold_line()

def download_do_data(sources_dp, srcs_cp):
    """ Download do database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading DO data files")

    ontology_fp = join(sources_dp, "doid.obo")
    download_file_md5_check(srcs_cp["do"]["ontology"], ontology_fp)

    print_bold_line()

# TODO: 将来进行统一的版本管理
def download_stitch_data(sources_dp, srcs_cp, version="5.0"):
    """ Download stitch database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading STITCH data files")

    chemicals_inchikeys_fp = join(sources_dp, "chemicals_inchikeys.tsv.gz")
    download_file_md5_check(srcs_cp["stitch"]["chem_inchikey"].replace('$VER$', 'v' + version), chemicals_inchikeys_fp)

    chemicals_smiles_fp = join(sources_dp, "chemicals_smiles.tsv.gz")
    download_file_md5_check(srcs_cp["stitch"]["chem_smiles"].replace('$VER$', 'v' + version), chemicals_smiles_fp)

    chem_chem_fp = join(sources_dp, "chemical_chemical_links.tsv.gz")
    download_file_md5_check(srcs_cp["stitch"]["chem_chem"].replace('$VER$', 'v' + version), chem_chem_fp)

    base_url = srcs_cp["stitch"]["base_url"]
    for specie in VALID_SPECIES:
        node = specie.node
        specie_links = join(sources_dp, f"{specie.code}_actions.tsv.gz")
        url = base_url + f"/actions.v{version}/{node}.actions.v{version}.tsv.gz"
        try:
            download_file_md5_check(url, specie_links)
        except:
            print(f"{specie.scientific_name} - {specie.node} NOT FOUND IN STITCH, SKIPPING")

    print_bold_line()

def download_umls_data(sources_dp, srcs_cp):
    """ Download umls database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Moving UMLS data files")
    # MRREL.RRF 暂时用不到
    # manual_files = ["MRDEF.RRF", "MRCONSO.RRF", "MRRANK.RRF", "MRREL.RRF"]
    manual_files = ["MRDEF.RRF", "MRCONSO.RRF", "MRRANK.RRF"]
    for manual_file in manual_files:
        if isfile(join("manually_prepared", manual_file)):
            copyfile(join("manually_prepared", manual_file), join(sources_dp, manual_file))
        else:
            raise FileNotFoundError(f"File {manual_file} not found in manually_prepared folder,\
                                    please check manually_prepared/README.md to prepare the data")
    
    print_bold_line()

def download_compath_data(sources_dp, srcs_cp):
    """ Download do database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading ComPath data files")

    compath_fp = join(sources_dp, "compath_mappings.tsv")
    download_file_md5_check(srcs_cp["compath"]["mappings"], compath_fp)

    print_bold_line()

def download_ncbigene_data(sources_dp, srcs_cp):
    """ Download NCBI Gene database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading NCBI Gene data files")

    gene2refseq_fp = join(sources_dp, "gene2refseq.gz")
    download_file_md5_check(srcs_cp["ncbigene"]["gene2refseq"], gene2refseq_fp)
    refseq2uniprot_fp = join(sources_dp, "refseq2uniprot.gz")
    download_file_md5_check(srcs_cp["ncbigene"]["refseq2uniprot"], refseq2uniprot_fp)

    print_bold_line()

def download_rnainter_data(sources_dp, srcs_cp):
    """ Download RNAInter database files
    Parameters
    ----------
    sources_dp : str
        the sources directory path
    srcs_cp : RawConfigParser
        source urls config parser
    """
    os.makedirs(sources_dp, exist_ok=True)
    print_section_header("Downloading RNA Interactome data files")

    # rna_rna = http://www.rnainter.org/raidMedia/download/Download_data_RR.tar.gz
    # rna_protein = http://www.rnainter.org/raidMedia/download/Download_data_RP.tar.gz
    # rna_dna = http://www.rnainter.org/raidMedia/download/Download_data_RD.tar.gz
    # rna_compound = http://www.rnainter.org/raidMedia/download/Download_data_RC.tar.gz

    rna_rna_fp = join(sources_dp, "rna_rna.tar.gz")
    download_file_md5_check(srcs_cp["rnainter"]["rna_rna"], rna_rna_fp)
    rna_protein_fp = join(sources_dp, "rna_protein.tar.gz")
    download_file_md5_check(srcs_cp["rnainter"]["rna_protein"], rna_protein_fp)
    rna_dna_fp = join(sources_dp, "rna_dna.tar.gz")
    download_file_md5_check(srcs_cp["rnainter"]["rna_dna"], rna_dna_fp)
    rna_compound_fp = join(sources_dp, "rna_compound.tar.gz")
    download_file_md5_check(srcs_cp["rnainter"]["rna_compound"], rna_compound_fp)

    print_bold_line()


# TODO: FooDB  FoodOn
# TODO: 中药-TCMSP
# TODO: 次要目标：通过命名学串联食物、中草药 NCBI Taxonomy

