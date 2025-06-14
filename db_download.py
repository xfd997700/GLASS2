#%%
from configparser import RawConfigParser
import os
from src.downloader import *
from src.parsers import *
from src.file_io import export_file_md5, file_has_valid_md5
from src.extras import inf_sym, done_sym
from src.info import ART_LOGO, ART_VERSION, PRINT_DOWNLOAD

config = RawConfigParser()
config.read('src/sources.ini')

username = config['default']['username']
password = config['default']['password']

sources_dp = "database/sources"
preprocessed_dp = "database/processed"
os.makedirs(sources_dp, exist_ok=True)


print("\033[96m")
print(ART_LOGO)
print(ART_VERSION)
print(PRINT_DOWNLOAD)
print("\033[0m")

uniprot_db = os.path.join(sources_dp, 'uniprot/')
download_uniprot_files(uniprot_db, config)
# chembl_db = os.path.join(sources_dp, 'chembl/')
# download_chembl_data(chembl_db, config)
# drugbank_db = os.path.join(sources_dp, 'drugbank/')
# download_drugbank_data(drugbank_db, config, username, password)
# bindingdb_db = os.path.join(sources_dp, 'bindingdb/')
# download_bindingdb_data(bindingdb_db, config)
iuphar_db = os.path.join(sources_dp, 'iuphar/')
download_iuphar_data(iuphar_db, config)
pdsp_db = os.path.join(sources_dp, 'pdsp/')
download_pdsp_data(pdsp_db, config)



reactome_db = os.path.join(sources_dp, 'reactome/')
download_reactome_data(reactome_db, config)
ctd_db = os.path.join(sources_dp, 'ctd/')
download_ctd_data(ctd_db, config)
phosphosite_db = os.path.join(sources_dp, 'phosphosite/')
download_phosphosite_data(phosphosite_db, config)
intact_db = os.path.join(sources_dp, 'intact/')
download_intact_data(intact_db, config)
sider_db = os.path.join(sources_dp, 'sider/')
download_sider_data(sider_db, config)
kegg_db = os.path.join(sources_dp, 'kegg/')
download_kegg_data(kegg_db, config)
mesh_db = os.path.join(sources_dp, 'mesh/')
download_mesh_data(mesh_db, config)
medgen_db = os.path.join(sources_dp, 'medgen/')
download_medgen_data(medgen_db, config)
hijazi20_db = os.path.join(sources_dp, 'hijazi20/')
download_hijazi20_data(hijazi20_db, config)
smpdb_db = os.path.join(sources_dp, 'smpdb/')
download_smpdb_data(smpdb_db, config)
hpa_db = os.path.join(sources_dp, 'hpa/')
download_hpa_data(hpa_db, config)
cellosaurus_db = os.path.join(sources_dp, 'cellosaurus/')
download_cellosaurus_data(cellosaurus_db, config)


foodb_db = os.path.join(sources_dp, 'foodb/')
download_foodb_data(foodb_db, config)
# tcmbank_db = os.path.join(sources_dp, 'tcmbank/')
# download_tcmbank_data(tcmbank_db, config)
hgnc_db = os.path.join(sources_dp, 'hgnc/')
download_hgnc_data(hgnc_db, config)
chembl_db = os.path.join(sources_dp, 'chembl/')
download_chembl_data(chembl_db, config)
hmdb_db = os.path.join(sources_dp, 'hmdb/')
download_hmdb_data(hmdb_db, config)
bindingdb_db = os.path.join(sources_dp, 'bindingdb/')
download_bindingdb_data(bindingdb_db, config)
go_db = os.path.join(sources_dp, 'go/')
download_go_data(go_db, config)
string_db = os.path.join(sources_dp, 'string/')
download_string_data(string_db, config)
tcdb_db = os.path.join(sources_dp, 'tcdb/')
download_tcdb_data(tcdb_db, config)
unichem_db = os.path.join(sources_dp, 'unichem/')
download_unichem_data(unichem_db, config)
pubchem_db = os.path.join(sources_dp, 'pubchem/')
download_pubchem_data(pubchem_db, config)
chebi_db = os.path.join(sources_dp, 'chebi/')
download_chebi_data(chebi_db, config)
hpo_db = os.path.join(sources_dp, 'hpo/')
download_hpo_data(hpo_db, config)
do_db = os.path.join(sources_dp, 'do/')
download_do_data(do_db, config)
stitch_db = os.path.join(sources_dp, 'stitch/')
download_stitch_data(stitch_db, config)
umls_db = os.path.join(sources_dp, 'umls/')
download_umls_data(umls_db, config)
compath_db = os.path.join(sources_dp, 'compath/')
download_compath_data(compath_db, config)
ncbigene_db = os.path.join(sources_dp, 'ncbigene/')
download_ncbigene_data(ncbigene_db, config)
rnainter_db = os.path.join(sources_dp, 'rnainter/')
download_rnainter_data(rnainter_db, config)

#%%
# ----------------------------------------------------------------------
# processing uniprot entries file
uniprot_parser = UniProtGPCRParser()
uniprot_dp = join(preprocessed_dp, 'uniprot')
os.makedirs(uniprot_dp, exist_ok=True)
uniprot_output_fps = [join(uniprot_dp, fn) for fn in uniprot_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in uniprot_output_fps]))

if invalid_md5:
    uniprot_parser.parse(uniprot_db, uniprot_dp)
    for ofp in uniprot_output_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "Uniprot processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing ChEMBL files
chembl_parser = ChEMBLParser()
chembl_dp = join(preprocessed_dp, 'chembl')
os.makedirs(chembl_dp, exist_ok=True)
chembl_fps = [join(chembl_dp, fn) for fn in chembl_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in chembl_fps]))

if invalid_md5:
    chembl_parser.parse_chem_db(chembl_db, chembl_dp)
    for ofp in chembl_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "ChEMBL processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)


# ----------------------------------------------------------------------
# processing ncbi gene entries file
ncbigene_parser = NCBIGeneParser()
ncbigene_dp = join(preprocessed_dp, 'ncbigene')
os.makedirs(ncbigene_dp, exist_ok=True)
ncbigene_output_fps = [join(ncbigene_dp, fn) for fn in ncbigene_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in ncbigene_output_fps]))

if invalid_md5:
    ncbigene_parser.parse(ncbigene_db, ncbigene_dp)
    for ofp in ncbigene_output_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "ncbigene processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing chebi entries file
chebi_parser = CHEBIParser()
chebi_dp = join(preprocessed_dp, 'chebi')
os.makedirs(chebi_dp, exist_ok=True)
chebi_output_fps = [join(chebi_dp, fn) for fn in chebi_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in chebi_output_fps]))

if invalid_md5:
    chebi_parser.parse(chebi_db, chebi_dp)
    for ofp in chebi_output_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "Chebi processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing compath entries file
compath_parser = ComPathParser()
compath_dp = join(preprocessed_dp, 'compath')
os.makedirs(compath_dp, exist_ok=True)
compath_output_fps = [join(compath_dp, fn) for fn in compath_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in compath_output_fps]))

if invalid_md5:
    compath_parser.parse(compath_db, compath_dp)
    for ofp in compath_output_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "ComPath processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing umls entries file
umls_parser = UMLSParser()
umls_dp = join(preprocessed_dp, 'umls')
os.makedirs(umls_dp, exist_ok=True)
umls_output_fps = [join(umls_dp, fn) for fn in umls_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in umls_output_fps]))

if invalid_md5:
    umls_parser.parse(umls_db, umls_dp)
    for ofp in umls_output_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "UMLS processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing pubchem entries file
pubchem_parser = PubChemParser()
pubchem_dp = join(preprocessed_dp, 'pubchem')
os.makedirs(pubchem_dp, exist_ok=True)
pubchem_output_fps = [join(pubchem_dp, fn) for fn in pubchem_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in pubchem_output_fps]))

if invalid_md5:
    pubchem_parser.parse(pubchem_db, pubchem_dp)
    for ofp in pubchem_output_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "PubChem processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)


# ----------------------------------------------------------------------
# processing HPA entries file
hpa_parser = HumanProteinAtlasParser()
hpa_dp = join(preprocessed_dp, 'hpa')
os.makedirs(hpa_dp, exist_ok=True)
hpa_files = ["hpa_antibodies.txt", "hpa_cellines_exp.txt", "hpa_tissues_exp.txt"]
hpa_fps = [join(hpa_dp, fn) for fn in hpa_files]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in hpa_fps]))
if invalid_md5:
    hpa_parser.parse_database_xml(join(hpa_db, "proteinatlas.xml.gz"), hpa_dp)
    for ofp in hpa_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "HPA processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)
# ----------------------------------------------------------------------
# processing Cellosaurus database file
cell_parser = CellosaurusParser()
cell_dp = join(preprocessed_dp, 'cellosaurus')
os.makedirs(cell_dp, exist_ok=True)
cell_parser.parse_db_file(join(cellosaurus_db, "cellosaurus_data.txt"), cell_dp)

cell_fps = [join(cell_dp, fn) for fn in cell_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in cell_fps]))
if invalid_md5:
    cell_parser.parse_db_file(join(cellosaurus_db, "cellosaurus_data.txt"), cell_dp)
    for ofp in cell_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "Cellosaurus processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)
# TODO: export md5 hashes of resulting files as in other databases

# ----------------------------------------------------------------------
# processing DrugBank entries file
drugbank_source_fp = join(drugbank_db, "drugbank_all_full_database.xml.zip")
drugbank_parser = DrugBankParser()
db_dp = join(preprocessed_dp, 'drugbank')
os.makedirs(db_dp, exist_ok=True)
drugbank_fps = [join(db_dp, fn) for fn in drugbank_parser.filelist]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in drugbank_fps]))
if invalid_md5:
    # Skip drugbank processing if the file is not in the source folder
    if os.path.exists(drugbank_source_fp):
        drugbank_parser.parse_drugbank_xml(drugbank_source_fp, db_dp)
        for ofp in drugbank_fps:
            export_file_md5(ofp)
    else:
        print(fail_sym + "Drugbank source not available >>> Skipping Drugbank processing")
else:
    print(inf_sym + "DrugBank processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing MESH
mesh_parser = MESHParser()
mesh_diseases_fp = join(mesh_db, 'mesh_diseases.txt')
mesh_supp_fp = join(mesh_db, "mesh_supp_concepts.xml")
mesh_dp = join(preprocessed_dp, 'mesh')
os.makedirs(mesh_dp, exist_ok=True)
mesh_fps = [join(mesh_dp, fn) for fn in mesh_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in mesh_fps]))
if invalid_md5:
    mesh_parser.parse_mesh(mesh_diseases_fp, mesh_supp_fp, mesh_dp)
    for ofp in mesh_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "MESH processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing Reactome entries file
reactome_parser = ReactomeParser()
reactome_dp = join(preprocessed_dp, 'reactome')
os.makedirs(reactome_dp, exist_ok=True)
reactome_fps = [join(reactome_dp, fn) for fn in reactome_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in reactome_fps]))

if invalid_md5:
    reactome_parser.parse_reactome(reactome_db, reactome_dp)
    for ofp in reactome_fps:
        export_file_md5(ofp)
else:
    print(
        inf_sym + "Reactome processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing CTD entries file
ctd_parser = CTDParser()
ctd_dp = join(preprocessed_dp, 'ctd')
os.makedirs(ctd_dp, exist_ok=True)
ctd_fps = [join(ctd_dp, fn) for fn in ctd_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in ctd_fps]))

if invalid_md5:
    ctd_parser.parse_ctd(ctd_db, join(uniprot_dp, 'uniprot_metadata.txt'), ctd_dp)
    for ofp in ctd_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "CTD processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing Phosphosite entries file
phosphosite_parser = PhosphositeParser()
phosphosite_dp = join(preprocessed_dp, 'phosphosite')
os.makedirs(phosphosite_dp, exist_ok=True)
phosphosite_fps = [join(phosphosite_dp, fn) for fn in phosphosite_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in phosphosite_fps]))

if invalid_md5:
    phosphosite_parser.parse_phosphosite(phosphosite_db, phosphosite_dp)
    for ofp in phosphosite_fps:
        export_file_md5(ofp)
else:
    print(
        inf_sym + "PhosphoSitePlus processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing Intact zip file
intact_parser = IntactParser()
intact_dp = join(preprocessed_dp, 'intact')
os.makedirs(intact_dp, exist_ok=True)
intact_fps = [join(intact_dp, fn) for fn in intact_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in intact_fps]))

if invalid_md5:
    intact_parser.parse_intact(intact_db, intact_dp, join(uniprot_dp, 'uniprot_ppi.txt'))
    for ofp in intact_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "Intact processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing Sider files
sider_parser = SiderParser()
sider_dp = join(preprocessed_dp, 'sider')
os.makedirs(sider_dp, exist_ok=True)
sider_fps = [join(sider_dp, fn) for fn in sider_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in sider_fps]))

if invalid_md5:
    sider_parser.parse_sider(sider_db, sider_dp)
    for ofp in sider_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "Sider processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing MedGen files
medgen_parser = MedgenParser()
medgen_dp = join(preprocessed_dp, 'medgen')
os.makedirs(medgen_dp, exist_ok=True)
medgen_fps = [join(medgen_dp, fn) for fn in medgen_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in medgen_fps]))

if invalid_md5:
    medgen_parser.parse_medgen(medgen_db, medgen_dp)
    for ofp in medgen_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "MedGen processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing Hijazi20 files
hijazi20_parser = Hijazi20Parser()
hijazi20_dp = join(preprocessed_dp, 'hijazi20')
os.makedirs(hijazi20_dp, exist_ok=True)
hijazi20_fps = [join(hijazi20_dp, fn) for fn in hijazi20_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in hijazi20_fps]))

if invalid_md5:
    hijazi20_parser.parse_phosphorylation(hijazi20_db, hijazi20_dp)
    for ofp in hijazi20_fps:
        export_file_md5(ofp)
else:
    print(
        inf_sym + "Hijazi20 processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)


# ----------------------------------------------------------------------
# processing FOODB files
foodb_parser = FOODBParser()
foodb_dp = join(preprocessed_dp, 'foodb')
os.makedirs(foodb_dp, exist_ok=True)
foodb_fps = [join(foodb_dp, fn) for fn in foodb_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in foodb_fps]))

if invalid_md5:
    foodb_parser.parse_foodb(foodb_db, foodb_dp)
    for ofp in foodb_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "foodb processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing HMDB files
hmdb_parser = HMDBParser()
hmdb_dp = join(preprocessed_dp, 'hmdb')
os.makedirs(hmdb_dp, exist_ok=True)
hmdb_fps = [join(hmdb_dp, fn) for fn in hmdb_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in hmdb_fps]))

if invalid_md5:
    hmdb_parser.parse_hmdb(hmdb_db, hmdb_dp)
    for ofp in hmdb_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "hmdb processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing SMPDB files
smpdb_parser = SmpdbParser()
smpdb_dp = join(preprocessed_dp, 'smpdb')
os.makedirs(smpdb_dp, exist_ok=True)
smpdb_fps = [join(smpdb_dp, fn) for fn in smpdb_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in smpdb_fps]))

if invalid_md5:
    smpdb_parser.parse_pathways(smpdb_db, smpdb_dp)
    for ofp in smpdb_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "SMPDB processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)
    
# ----------------------------------------------------------------------
# processing BindingDB files
bindingdb_parser = BindingDBParser()
bindingdb_dp = join(preprocessed_dp, 'bindingdb')
os.makedirs(bindingdb_dp, exist_ok=True)
bindingdb_fps = [join(bindingdb_dp, fn) for fn in bindingdb_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in bindingdb_fps]))

bindingdb_db = os.path.join(sources_dp, 'bindingdb/')
if invalid_md5:
    bindingdb_parser.parse_cpi(bindingdb_db, bindingdb_dp)
    for ofp in bindingdb_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "bindingdb processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing GO files
go_parser = GOParser()
go_dp = join(preprocessed_dp, 'go')
os.makedirs(go_dp, exist_ok=True)
go_fps = [join(go_dp, fn) for fn in go_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in go_fps]))

go_db = os.path.join(sources_dp, 'go/')
if invalid_md5:
    go_parser.parse(go_db, go_dp)
    for ofp in go_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "go processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing STRING files
string_parser = STRINGParser()
string_dp = join(preprocessed_dp, 'string')
os.makedirs(string_dp, exist_ok=True)
string_fps = [join(string_dp, fn) for fn in string_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in string_fps]))

string_db = os.path.join(sources_dp, 'string/')

# string_parser._parse_aliases(string_db, string_dp)

if invalid_md5:
    string_parser.parse(string_db, string_dp)
    for ofp in string_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "string processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing TCDB files
tcdb_parser = TCDBParser()
tcdb_dp = join(preprocessed_dp, 'tcdb')
os.makedirs(tcdb_dp, exist_ok=True)
tcdb_fps = [join(tcdb_dp, fn) for fn in tcdb_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in tcdb_fps]))

tcdb_db = os.path.join(sources_dp, 'tcdb/')
if invalid_md5:
    tcdb_parser.parse(tcdb_db, tcdb_dp)
    for ofp in tcdb_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "tcdb processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing KEGG links
kegg_parser = KeggParser()
kegg_dp = join(preprocessed_dp, 'kegg')
os.makedirs(kegg_dp, exist_ok=True)
kegg_fps = [join(kegg_dp, fn) for fn in kegg_parser.filelist]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in kegg_fps]))
if invalid_md5:
    kegg_parser.parse_kegg(kegg_db, kegg_dp)
    for ofp in kegg_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "KEGG processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing HPO links
hpo_parser = HPOParser()
hpo_dp = join(preprocessed_dp, 'hpo')
os.makedirs(hpo_dp, exist_ok=True)
hpo_fps = [join(hpo_dp, fn) for fn in hpo_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in hpo_fps]))

hpo_db = os.path.join(sources_dp, 'hpo/')
if invalid_md5:
    hpo_parser.parse(hpo_db, hpo_dp)
    for ofp in hpo_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "HPO processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing DO links
do_parser = DOParser()
do_dp = join(preprocessed_dp, 'do')
os.makedirs(do_dp, exist_ok=True)
do_fps = [join(do_dp, fn) for fn in do_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in do_fps]))

do_db = os.path.join(sources_dp, 'do/')
if invalid_md5:
    do_parser.parse(do_db, do_dp)
    for ofp in do_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "DO processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing STITCH links
stitch_parser = STITCHParser()
stitch_dp = join(preprocessed_dp, 'stitch')
os.makedirs(stitch_dp, exist_ok=True)
stitch_fps = [join(stitch_dp, fn) for fn in stitch_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in stitch_fps]))

stitch_db = os.path.join(sources_dp, 'stitch/')
if invalid_md5:
    stitch_parser.parse(stitch_db, stitch_dp)
    for ofp in stitch_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "STITCH processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# ----------------------------------------------------------------------
# UniKG chemical and protein parsing
unikg_parser = UniKGParser()
data_root = 'database/processed'
source_root = 'database/sources'

chem_files = [join(data_root, fn) for fn in unikg_parser.chem_filenames]
prot_files = [join(data_root, fn) for fn in unikg_parser.prot_filenames]

do_chem = False
for files in chem_files:
    if not os.path.exists(files):
        do_chem = True
        break

do_prot = False
for files in prot_files:
    if not os.path.exists(files):
        do_prot = True
        break

if do_chem:
    unikg_parser._parse_compounds(source_root, data_root)
if do_prot:
    unikg_parser._parse_proteins(source_root, data_root)