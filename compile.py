#%%
from configparser import RawConfigParser
import os
from src.downloader import *
from src.parsers import *
from src.file_io import export_file_md5, file_has_valid_md5
from src.extras import inf_sym, done_sym
from src.info import ART_LOGO, ART_VERSION, PRINT_INFO
from src.compile_func import merge_all_compound_info, fetch_xref, merge_all_cpi_data, reg_postprocess
import shutil
config = RawConfigParser()
config.read('src/sources.ini')

username = config['default']['username']
password = config['default']['password']

sources_dp = "database/sources"
processed_dp = "database/processed"
os.makedirs(sources_dp, exist_ok=True)


print("\033[96m")
print(ART_LOGO)
print(ART_VERSION)
print(PRINT_INFO)
print("\033[0m")

uniprot_db = os.path.join(sources_dp, 'uniprot/')
download_uniprot_files(uniprot_db, config)
chembl_db = os.path.join(sources_dp, 'chembl/')
# download_chembl_data(chembl_db, config)
drugbank_db = os.path.join(sources_dp, 'drugbank/')
# download_drugbank_data(drugbank_db, config, username, password)
bindingdb_db = os.path.join(sources_dp, 'bindingdb/')
# download_bindingdb_data(bindingdb_db, config)
iuphar_db = os.path.join(sources_dp, 'iuphar/')
download_iuphar_data(iuphar_db, config)
pdsp_db = os.path.join(sources_dp, 'pdsp/')
download_pdsp_data(pdsp_db, config)
unichem_db = os.path.join(sources_dp, 'unichem/')
# download_unichem_data(unichem_db, config)

pubchem_db = os.path.join(sources_dp, 'pubchem/')
# download_pubchem_data(pubchem_db, config)

llm_db = os.path.join(sources_dp, 'llm/')
download_llm_data(llm_db, config)

glass_db = os.path.join(sources_dp, 'glass/')
# download_glass_data(glass_db, config)

#%%
# ----------------------------------------------------------------------
# processing pubchem entries file
pubchem_parser = PubChemParser()
pubchem_dp = join(processed_dp, 'pubchem')
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
# processing uniprot entries file
uniprot_parser = UniProtGPCRParser()
uniprot_dp = join(processed_dp, 'uniprot')
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
chembl_dp = join(processed_dp, 'chembl')
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
# processing DrugBank entries file
drugbank_source_fp = join(drugbank_db, "drugbank_all_full_database.xml.zip")
drugbank_parser = DrugBankParser()
db_dp = join(processed_dp, 'drugbank')
os.makedirs(db_dp, exist_ok=True)
drugbank_fps = [join(db_dp, fn) for fn in drugbank_parser.filelist]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in drugbank_fps]))
if invalid_md5:
    # Skip drugbank processing if the file is not in the source folder
    if os.path.exists(drugbank_source_fp):
        drugbank_parser.parse(drugbank_source_fp, db_dp)
        for ofp in drugbank_fps:
            export_file_md5(ofp)
    else:
        print(fail_sym + "Drugbank source not available >>> Skipping Drugbank processing")
else:
    print(inf_sym + "DrugBank processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing BindingDB files
bindingdb_parser = BindingDBParser()
bindingdb_dp = join(processed_dp, 'bindingdb')
os.makedirs(bindingdb_dp, exist_ok=True)
bindingdb_fps = [join(bindingdb_dp, fn) for fn in bindingdb_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in bindingdb_fps]))

bindingdb_db = os.path.join(sources_dp, 'bindingdb/')
if invalid_md5:
    bindingdb_parser.parse(bindingdb_db, bindingdb_dp)
    for ofp in bindingdb_fps:
        export_file_md5(ofp)
else:
    print(inf_sym + "bindingdb processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing PDSP entries file
pdsp_parser = PDSPParser()
db_dp = join(processed_dp, 'pdsp')
os.makedirs(db_dp, exist_ok=True)
pdsp_fps = [join(db_dp, fn) for fn in pdsp_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in pdsp_fps]))
if invalid_md5:
    # Skip pdsp processing if the file is not in the source folder
    if os.path.exists(pdsp_db):
        pdsp_parser.parse(pdsp_db, db_dp)
        for ofp in pdsp_fps:
            export_file_md5(ofp)
    else:
        print(fail_sym + "PDSP source not available >>> Skipping PDSP processing")
else:
    print(inf_sym + "PDSP processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing IUPHAR entries file
iuphar_parser = IUPHARParser()
db_dp = join(processed_dp, 'iuphar')
os.makedirs(db_dp, exist_ok=True)
iuphar_fps = [join(db_dp, fn) for fn in iuphar_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in iuphar_fps]))
if invalid_md5:
    # Skip pdsp processing if the file is not in the source folder
    if os.path.exists(iuphar_db):
        iuphar_parser.parse(iuphar_db, db_dp)
        for ofp in iuphar_fps:
            export_file_md5(ofp)
    else:
        print(fail_sym + "IUPHAR source not available >>> Skipping IUPHAR processing")
else:
    print(inf_sym + "IUPHAR processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing LLM entries file
llm_parser = LLMParser()
db_dp = join(processed_dp, 'llm')
os.makedirs(db_dp, exist_ok=True)
llm_fps = [join(db_dp, fn) for fn in llm_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in llm_fps]))
if invalid_md5:
    # Skip pdsp processing if the file is not in the source folder
    if os.path.exists(llm_db):
        llm_parser.parse(llm_db, db_dp)
        for ofp in llm_fps:
            export_file_md5(ofp)
    else:
        print(fail_sym + "LLM source not available >>> Skipping LLM processing")
else:
    print(inf_sym + "LLM processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)

# ----------------------------------------------------------------------
# processing GLASS entries file
glass_parser = GLASSParser()
db_dp = join(processed_dp, 'glass')
os.makedirs(db_dp, exist_ok=True)
glass_fps = [join(db_dp, fn) for fn in glass_parser.filenames]
invalid_md5 = bool(sum([not file_has_valid_md5(ofp) for ofp in glass_fps]))
if invalid_md5:
    # Skip pdsp processing if the file is not in the source folder
    if os.path.exists(glass_db):
        glass_parser.parse(glass_db, db_dp)
        for ofp in glass_fps:
            export_file_md5(ofp)
    else:
        print(fail_sym + "GLASS source not available >>> Skipping GLASS processing")
else:
    print(inf_sym + "GLASS processed files exists with valid md5 hashes %s. >>> Parsing not required." % done_sym)
    
# %%
# ----------------------------------------------------------------------
# Combine all processed data
output_dp = 'database/processed/'
compile_dp = 'database/glass2/'
os.makedirs(compile_dp, exist_ok=True)

result = merge_all_cpi_data(output_dp, compile_dp)

print(inf_sym + "Merged CPI data: %s" % done_sym)
merged_cinfo = merge_all_compound_info(processed_dp)

unique_inchikeys = set(sys.intern(x) for x in result['full']['compound_inchikey'].unique())
filtered_cinfo = merged_cinfo[merged_cinfo['InChIKey'].isin(unique_inchikeys)].reset_index(drop=True)
filtered_cinfo.to_csv(join(compile_dp, 'ligands.tsv'), sep='\t', index=False, encoding='utf-8')

print(inf_sym + f"Filtered compound info saved to {join(compile_dp, 'ligands.tsv')} {done_sym}")

print(inf_sym + "Copying extra data...")
unichem_src = os.path.join(processed_dp, 'unichem', 'inchikey2uci.json')
unichem_dst = os.path.join(compile_dp, 'inchikey2uci.json')
if os.path.exists(unichem_src):
    shutil.copy2(unichem_src, unichem_dst)
else:
    print(fail_sym + f"Unichem mapping file not found at {unichem_src} >>> Skipping copy")

# Copy uniprot gpcr_entries.json to compile_dp as protein.json
uniprot_src = os.path.join(processed_dp, 'uniprot', 'gpcr_entries.json')
protein_dst = os.path.join(compile_dp, 'protein.json')
if os.path.exists(uniprot_src):
    shutil.copy2(uniprot_src, protein_dst)
else:
    print(fail_sym + f"Uniprot GPCR entries file not found at {uniprot_src} >>> Skipping copy")
fetch_xref(filtered_cinfo, compile_dp)


# %%
import json
import pandas as pd
from os.path import join, isfile

compile_dp = 'database/glass2/'
with open('database/glass2/inchikey2uci.json', 'r', encoding='utf-8') as f:
    inchikey2uci = json.load(f)
    # Convert inchikey2uci to filtered_cinfo with ['InChIKey', 'UCI'] columns
filtered_cinfo = pd.DataFrame([
    {'InChIKey': k, 'UCI': v} for k, v in inchikey2uci.items()
])

result = {
    'cls': pd.read_csv(join(compile_dp, 'glass2_cls.csv'), encoding='utf-8'),
    'reg_act': pd.read_csv(join(compile_dp, 'glass2_reg_act.csv'), encoding='utf-8'),
    'reg_inact': pd.read_csv(join(compile_dp, 'glass2_reg_inact.csv'), encoding='utf-8'),
    'full': pd.read_csv(join(compile_dp, 'glass2_full.tsv'), sep='\t', encoding='utf-8')
}


print('adding uci to cpi data...')
for k, df in result.items():
    if isinstance(df, pd.DataFrame) and 'compound_inchikey' in df.columns:
        result[k] = df.merge(filtered_cinfo[['InChIKey', 'UCI']], left_on='compound_inchikey', right_on='InChIKey', how='left').rename(columns={'UCI': 'uci'})
        cols = list(result[k].columns)
        if 'uci' in cols:
            cols.insert(1, cols.pop(cols.index('uci')))
            result[k] = result[k][cols]
        if 'InChIKey' in result[k].columns:
            result[k] = result[k].drop(columns=['InChIKey'])
        
result['cls'].to_csv(join(compile_dp, 'glass2_cls.csv'), index=False, encoding='utf-8')
result['reg_act'].to_csv(join(compile_dp, 'glass2_reg_act.csv'), index=False, encoding='utf-8')
result['reg_inact'].to_csv(join(compile_dp, 'glass2_reg_inact.csv'), index=False, encoding='utf-8')
result['full'].to_csv(join(compile_dp, 'glass2_full.tsv'), sep='\t', index=False, encoding='utf-8')

reg_postprocess(join(compile_dp, 'glass2_reg_act.csv'), join(compile_dp, 'glass2_reg_inact.csv'))
# %%
