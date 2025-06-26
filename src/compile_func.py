from collections import defaultdict
import gzip
import json
import os
from os.path import join
import re
import sys
import numpy as np
import pandas as pd
from tqdm import tqdm
from .extras import *
from timeit import default_timer as timer

def merge_cinfo(output_dp):
    """
    Simplified version of merge compound information - only keep essential fields
    
    Parameters
    ----------
    output_dp : str
        The path to the output directory containing all database folders
    
    Returns
    -------
    pandas.DataFrame
        Simplified merged compound information DataFrame
    """
    print_section_header("Merging compound information (Simplified)")
    start = timer()
    
    # Database file paths
    database_files = {
        'ChEMBL': join(output_dp, 'chembl/chembl_cinfo.tsv'),
        'IUPHAR': join(output_dp, 'iuphar/iuphar_cinfo.tsv'),
        'BindingDB': join(output_dp, 'bindingdb/bindingdb_cinfo.tsv'),
        'DrugBank': join(output_dp, 'drugbank/drugbank_cinfo.tsv'),
        'PDSP': join(output_dp, 'pdsp/pdsp_cinfo.tsv'),
        'GLASS': join(output_dp, 'glass/glass_cinfo.tsv')
    }
    
    # Load all dataframes
    dataframes = []
    total_entries = 0
    
    print("Loading dataframes...")
    for db_name, file_path in database_files.items():
        if os.path.exists(file_path):
            try:
                df = pd.read_csv(file_path, sep='\t', encoding='utf-8')
                print(f"Loaded {len(df)} entries from {db_name}")
                dataframes.append(df)
                total_entries += len(df)
            except Exception as e:
                print(f"Warning: Could not load {db_name} data: {e}")
        else:
            print(f"Warning: {db_name} compound info file not found: {file_path}")
    
    if not dataframes:
        print("No compound information files found!")
        return pd.DataFrame()
    
    print(f"Total entries loaded: {total_entries}")
    
    # Combine all dataframes
    print("Combining dataframes...")
    combined_df = pd.concat(dataframes, ignore_index=True, sort=False)
    
    # Filter out entries without InChIKey
    print("Filtering entries with valid InChIKey...")
    valid_df = combined_df[combined_df['InChIKey'].notna() & 
                          (combined_df['InChIKey'].str.strip() != '')].copy()
    
    invalid_count = len(combined_df) - len(valid_df)
    if invalid_count > 0:
        print(f"Removed {invalid_count} entries without valid InChIKey")
    
    # Define the columns we want to keep
    target_columns = ["Name", "InChIKey", "InChI", "SMILES", "IUPAC name"]
    
    # Only keep columns that exist in the dataframe
    available_columns = [col for col in target_columns if col in valid_df.columns]
    
    # Keep only the target columns
    simplified_df = valid_df[available_columns].copy()
    
    print(f"Keeping columns: {available_columns}")
    
    # Clean the data
    print("Cleaning data...")
    for col in available_columns:
        if col != 'InChIKey':  # Don't modify InChIKey
            # Replace NaN-like values with actual NaN
            simplified_df[col] = simplified_df[col].astype(str)
            simplified_df[col] = simplified_df[col].replace(['nan', 'NaN', '<NA>', 'None'], pd.NA)
            simplified_df[col] = simplified_df[col].str.strip()
            simplified_df[col] = simplified_df[col].replace('', pd.NA)
    
    # 关键优化：排序后去重
    print("Sorting and deduplicating by InChIKey...")
    
    # 创建一个排序键：有Name的条目排前面，Name越长越前面
    print("Creating sort key for Name priority...")
    simplified_df['_name_length'] = simplified_df['Name'].fillna('').astype(str).str.len()
    simplified_df['_has_name'] = simplified_df['Name'].notna() & (simplified_df['Name'].astype(str).str.strip() != '')
    
    # 排序：先按InChIKey分组，然后按是否有Name（有Name的在前），再按Name长度（长的在前）
    print("Sorting by InChIKey and Name")
    sorted_df = simplified_df.sort_values(['InChIKey', '_has_name', '_name_length'], 
                                            ascending=[True, False, False])

    # 去重：保留每个InChIKey的第一条记录（即Name最好的那条）
    print("Dropping duplicates")
    final_df = sorted_df.drop_duplicates(subset=['InChIKey'], keep='first')

    
    # 删除辅助列
    final_df = final_df.drop(columns=['_name_length', '_has_name'])
    
    # 重新按InChIKey排序
    final_df = final_df.sort_values('InChIKey').reset_index(drop=True)
    
    # Save simplified merged compound info
    output_file = join(output_dp, 'merged_cinfo.tsv')
    final_df.to_csv(output_file, sep='\t', index=False, encoding='utf-8')
    
    # Print statistics
    print(f"\nSimplified Merging Statistics:")
    print(f"  - Total input entries: {total_entries}")
    print(f"  - Valid entries with InChIKey: {len(valid_df)}")
    print(f"  - Unique compounds (by InChIKey): {len(final_df)}")
    print(f"  - Columns kept: {len(available_columns)}")
    
    # Print field coverage
    print(f"\nField Coverage in Final Dataset:")
    for col in available_columns:
        if col != 'InChIKey':
            non_empty_count = (final_df[col].notna() & 
                             (final_df[col].astype(str).str.strip() != '') &
                             (final_df[col].astype(str).str.lower() != 'nan')).sum()
            coverage_pct = (non_empty_count / len(final_df)) * 100
            print(f"  - {col}: {non_empty_count}/{len(final_df)} ({coverage_pct:.1f}%)")
    
    print(f"\nSaved simplified merged compound information to: {output_file}")
    print(done_sym + f" Took {timer() - start:.2f} Seconds.")
    
    return final_df


def analyze_cinfo_coverage(df):
    """
    Analyze the coverage of different fields across databases
    """
    print_section_header("Analyzing compound information coverage")
    
    print(f"Analyzing {len(df)} compounds...")
    print(f"\nField Coverage Analysis:")
    print("-" * 50)
    
    # Analyze coverage for each field
    coverage_stats = []
    for column in df.columns:
        if column in ['InChIKey', 'Source_Databases', 'Source_Count']:
            continue
            
        non_null_count = df[column].notna().sum()
        non_empty_count = (df[column].notna() & 
                          (df[column].astype(str).str.strip() != '') &
                          (df[column].astype(str).str.lower() != 'nan')).sum()
        
        coverage_stats.append({
            'Field': column,
            'Non_Null': non_null_count,
            'Non_Empty': non_empty_count,
            'Coverage_Percent': (non_empty_count / len(df)) * 100
        })
    
    # Sort by coverage percentage
    coverage_df = pd.DataFrame(coverage_stats)
    coverage_df = coverage_df.sort_values('Coverage_Percent', ascending=False)
    
    for _, row in coverage_df.iterrows():
        print(f"{row['Field']:<25}: {row['Non_Empty']:>6}/{len(df)} ({row['Coverage_Percent']:>6.1f}%)")
    
    # Save coverage analysis
    # coverage_file = join(output_dp, 'cinfo_coverage_analysis.tsv')
    # coverage_df.to_csv(coverage_file, sep='\t', index=False, encoding='utf-8')
    # print(f"\nCoverage analysis saved to: {coverage_file}")


def merge_all_compound_info(output_dp):
    """
    Complete workflow to merge compound information and analyze coverage
    
    Parameters
    ----------
    output_dp : str
        The path to the output directory
    """
    # Merge compound information
    merged_df = merge_cinfo(output_dp)
    # merged_df = pd.read_csv('database/processed/merged_cinfo.tsv', sep='\t', encoding='utf-8')
    unichem_data = join(output_dp, 'unichem')
    unichem_source = 'database/sources/unichem'
    os.makedirs(unichem_data, exist_ok=True)
    sub_structure = join(unichem_data, 'structure.tsv.gz')
    if not os.path.exists(sub_structure):
        unichem_structure = join(unichem_source, 'structure.tsv.gz')
        unique_inchikeys = set(sys.intern(x) for x in merged_df['InChIKey'].unique())
        print("******************************************************************************")
        print("Start scanning related UniChem UCI, this may take a while...")
        print(f"To track the progress, you can manually check the unzipped file size at:\n\033[95m{unichem_structure}\033[0m\n")

        with gzip.open(unichem_structure, "rb") as f_in, gzip.open(sub_structure, "wb") as f_out:
            # 读取第一行（表头）并写入输出文件
            header = f_in.readline()
            f_out.write(header)
            
            # 直接遍历文件行，并使用 tqdm 显示进度条
            with tqdm(unit="B", unit_scale=True,
                    desc="Scanning...") as pbar:
                for line in f_in:
                    # 解码并去掉换行符
                    decoded_line = line.decode("utf-8").strip()
                    # 按制表符分割
                    columns = decoded_line.split("\t")
                    
                    # 确保行有足够的列，且第 3 列在 target_keys 中
                    if len(columns) > 2 and sys.intern(columns[2]) in unique_inchikeys:
                        f_out.write(line)
                    # 更新 tqdm 进度条
                    pbar.update(len(line))

        print("\033[94mDone scanning unichem structure!\033[0m")
        print(f"Saved {sub_structure}")
        print("******************************************************************************")

    unichem_df = pd.read_csv(sub_structure, sep='\t')
    unichem_df.rename(columns={'STANDARDINCHIKEY': 'InChIKey'}, inplace=True)
    print("merging unichem data...")
    merged_df = merged_df.merge(unichem_df[['UCI', 'STANDARDINCHI', 'InChIKey']], on='InChIKey', how='left', suffixes=('', '_unichem'))
    merged_df['CP'] = np.where(merged_df['InChI'] == merged_df['STANDARDINCHI'], 1, 0)
    merged_df = merged_df.dropna(subset=['UCI'])
    merged_df.rename(columns={'STANDARDINCHI': 'InChI_unichem'}, inplace=True)
    merged_df['UCI'] = merged_df['UCI'].astype(np.int64)
    merged_df['InChIKey'] = merged_df['InChIKey'].astype(str)
    # all_chem_root = join(output_dp, "merged_cinfo.tsv")
    merged_df = merged_df.dropna(subset=['InChIKey'])
    # merged_df.to_csv(all_chem_root, sep='\t', index=False)

    print(f"Done merging, total {len(merged_df)} compounds")
    
    # extract the InChIKey and UCI from merged_df and make a dict{InChIKey: UCI}
    inchikey_uci = dict(zip(merged_df['InChIKey'], 'UCI:' + merged_df['UCI'].astype(str)))
    inchikey_uci_json = join(unichem_data, "inchikey2uci.json")
    with open(inchikey_uci_json, 'w') as f:
        json.dump(inchikey_uci, f)
    
    # analyze_cinfo_coverage(all_chem_root)
    return merged_df


def fetch_xref(merged_df, output_dp):
    """
    Fetch external references for compounds from UniChem reference data
    
    Parameters
    ----------
    merged_df : pandas.DataFrame
        Merged compound information DataFrame containing UCI column
    unichem_source_dp : str
        Path to the UniChem source directory containing source.tsv.gz and reference.tsv.gz
    
    Returns
    -------
    dict
        Dictionary mapping UCI to external database references
        Format: {'UCI:123': {'chembl': ['CHEMBL5546', ...], 'drugbank': []}, ...}
    """
    print_section_header("Fetching external references from UniChem")
    unichem_source_dp = 'database/sources/unichem'
    start = timer()
    
    # File paths
    source_file = join(unichem_source_dp, 'source.tsv.gz')
    reference_file = join(unichem_source_dp, 'reference.tsv.gz')
    
    # Check if files exist
    if not os.path.exists(source_file):
        print(f"Source file not found: {source_file}")
        return {}
    
    if not os.path.exists(reference_file):
        print(f"Reference file not found: {reference_file}")
        return {}
    
    # Step 1: Load source mapping (SRC_ID to database name)
    print("Loading source mapping...")
    source_df = pd.read_csv(source_file, sep='\t', encoding='utf-8')
    
    # Create SRC_ID to database name mapping
    src_id_to_name = {}
    for _, row in source_df.iterrows():
        src_id = row['SRC_ID']
        name = row['NAME'].lower()  # Convert to lowercase for consistency
        src_id_to_name[src_id] = name
    
    print(f"Loaded {len(src_id_to_name)} source databases")
    
    # Step 2: Extract unique UCIs from merged_df
    print("Extracting UCIs from merged data...")
    unique_ucis = set()
    for uci in merged_df['UCI'].dropna().unique():
        unique_ucis.add(int(uci))  # Store as integer for faster lookup
    
    print(f"Found {len(unique_ucis)} unique UCIs to process")
    
    # Step 3: Process reference file to build external reference mapping
    print("Processing reference file...")
    
    # Initialize the result dictionary defaultdict(defaultdict(list))
    uci_xref_dict = defaultdict(lambda: defaultdict(list))
    # Process reference file in chunks to handle large files
    chunk_size = 100000
    processed_rows = 0
    
    print("Reading reference file...")
    ref_df = pd.read_csv(reference_file, sep='\t', encoding='utf-8')
    ref_df = ref_df[ref_df['ASSIGNMENT'] == 1]  # Filter for valid assignments
    for i, row in tqdm(ref_df.iterrows(),
                       desc="Fetching external references",
                       total=len(ref_df),
                       unit="rows"):
        uci = int(row['UCI'])
        if uci in unique_ucis:
            src_id = row['SRC_ID']
            src_compound_id = str(row['SRC_COMPOUND_ID'])
            
            # Get database name from SRC_ID
            db_name = src_id_to_name.get(src_id)
            if db_name:
                uci_key = f"UCI:{uci}"
                uci_xref_dict[uci_key][db_name].append(src_compound_id)



    # with pd.read_csv(reference_file, sep='\t', encoding='utf-8', chunksize=chunk_size) as reader:
    #     for chunk_num, chunk in enumerate(reader):
    #         # Filter for our UCIs and valid assignments only
    #         filtered_chunk = chunk[
    #             (chunk['UCI'].isin(unique_ucis)) & 
    #             (chunk['ASSIGNMENT'] == 1)
    #         ]
            
    #         # Process each row in the filtered chunk
    #         for _, row in filtered_chunk.iterrows():
    #             uci = int(row['UCI'])
    #             src_id = row['SRC_ID']
    #             src_compound_id = str(row['SRC_COMPOUND_ID'])
                
    #             # Get database name from SRC_ID
    #             db_name = src_id_to_name.get(src_id)
    #             if db_name:
    #                 uci_key = f"UCI:{uci}"
    #                 if uci_key in uci_xref_dict:
    #                     uci_xref_dict[uci_key][db_name].append(src_compound_id)
            
    #         processed_rows += len(chunk)
            
    #         # Progress update
    #         if chunk_num % 10 == 0:
    #             speed = processed_rows / (timer() - start)
    #             print(f"\r{prc_sym}Processed {processed_rows} reference rows. Speed: {speed:.0f} rows/second", end="", flush=True)
    
    # print(f"{done_sym} Finished processing reference file. Total rows processed: {processed_rows}", flush=True)
    
    # Step 4: Clean up empty lists and provide statistics
    print("Cleaning up and generating statistics...")
    
    total_xrefs = 0
    db_stats = {}
    
    for uci_key, db_dict in uci_xref_dict.items():
        for db_name, compound_list in db_dict.items():
            if compound_list:
                # Remove duplicates while preserving order
                unique_compounds = list(dict.fromkeys(compound_list))
                uci_xref_dict[uci_key][db_name] = unique_compounds
                total_xrefs += len(unique_compounds)
                
                # Update database statistics
                if db_name not in db_stats:
                    db_stats[db_name] = 0
                db_stats[db_name] += len(unique_compounds)
    
    # Print statistics
    compounds_with_xrefs = sum(1 for uci_dict in uci_xref_dict.values() 
                              if any(db_list for db_list in uci_dict.values()))
    
    print(f"\nExternal Reference Statistics:")
    print(f"  - Total compounds: {len(uci_xref_dict)}")
    print(f"  - Compounds with external references: {compounds_with_xrefs}")
    print(f"  - Total external references: {total_xrefs}")
    
    print(f"\nDatabase-wise Reference Distribution:")
    for db_name, count in sorted(db_stats.items(), key=lambda x: x[1], reverse=True):
        if count > 0:
            print(f"  - {db_name}: {count} references")
    
    print(done_sym + f" Took {timer() - start:.2f} Seconds.")
    save_xref_data(uci_xref_dict, output_dp)
    # return uci_xref_dict


def save_xref_data(uci_xref_dict, unichem_data_dp):
    print("Saving external reference data...")
    
    # unichem_data_dp = join(output_dp, 'unichem')
    # os.makedirs(unichem_data_dp, exist_ok=True)
    
    # Save as JSON
    xref_json_file = join(unichem_data_dp, 'uci_xref.json')
    with open(xref_json_file, 'w', encoding='utf-8') as f:
        json.dump(uci_xref_dict, f, ensure_ascii=False, indent=2)
    
    print(f"Saved external references as JSON: {xref_json_file}")

def split_chembl_pmid_doi(reference: str):
    entries = reference.split(',')
    pmids = []
    dois = []
    
    for entry in entries:
        entry = entry.strip()
        if re.fullmatch(r'\d+', entry):  # 全是数字：PMID
            pmids.append(entry)
        else:  # 含有字母或符号：DOI
            dois.append(entry)
    
    pmid_str = ','.join(pmids)
    doi_str = ','.join(dois)
    
    return pmid_str, doi_str

class CPIMerger:
    """
    Compound-Protein Interaction (CPI) data merger class
    """
    def __init__(self, output_dp):
        self.output_dp = output_dp

        
    
    def _form_chembl_cpi(self):
        """
        Process ChEMBL CPI data with activity filtering
        
        Returns
        -------
        dict
            Dictionary containing:
            - 'formed_cpi': Valid CPI data
            - 'true_negative': Entries with negative activity comments
            - 'uncertain': Entries without activity comments and experimental data
        """
        print("Processing ChEMBL CPI data...")
        
        # Load ChEMBL CPI data
        chembl_file = join(self.output_dp, 'chembl/chembl_cpi.tsv')
        if not os.path.exists(chembl_file):
            print("ChEMBL CPI file not found")
            return {'formed_cpi': pd.DataFrame(), 'true_negative': pd.DataFrame(), 'uncertain': pd.DataFrame()}
        
        chembl_df = pd.read_csv(chembl_file, sep='\t', encoding='utf-8')
        print(f"Loaded {len(chembl_df)} ChEMBL CPI entries")
        
        # Initialize result containers
        formed_cpi = []
        true_negative = []
        uncertain = []
        
        # Define negative activity comment patterns (following the SQLite logic)

        def is_negative_comment(comment: str) -> bool:
            # if comment is None:
            #     return True
            if pd.isna(comment):
                return False
            
            comment_lower = str(comment).lower()
            # 精确匹配不区分大小写
            if comment_lower == "inconclusive":
                return True
            # 含有下列否定或不活性提示词的，视为无效
            keywords = ['non ', 'no ', 'not ', 'weak', 'nactive']
            for kw in keywords:
                if kw in comment_lower:
                    return True
            return False
        
        def has_experimental_data(row):
            """Check if row has experimental data (standard_value and related fields)"""
            standard_value = row.get('standard_value')
            standard_type = row.get('standard_type')
            
            # Check if we have valid experimental data
            if pd.notna(standard_value) and pd.notna(standard_type):
                try:
                    float(standard_value)  # Check if value is numeric
                    return True
                except (ValueError, TypeError):
                    pass
            return False
        
        # Process each entry
        for idx, row in tqdm(chembl_df.iterrows(), total=len(chembl_df), desc="Forming ChEMBL entries"):
            activity_comment = row.get('activity_comment')
            pmid, doi = None, None
            reference = row.get('pubmed_id_or_doi')
            if not pd.isna(reference):
                # Split PubMed ID and DOI
                pmid, doi = split_chembl_pmid_doi(reference)

            # Check for negative activity comments
            if is_negative_comment(activity_comment):
                # True negative entry
                entry = {
                    'compound_inchikey': row['compound_inchikey'],
                    'target_uniprot_id': row['target_uniprot_id'],
                    'source_database': 'chembl',
                    'source_compound_id': row['compound_chembl_id'],
                    'standard_type': row.get('standard_type'),
                    'standard_relation': row.get('standard_relation'),
                    'standard_value': row.get('standard_value'),
                    'standard_units': row.get('standard_units'),
                    'activity_comment': activity_comment,
                    'pmid': pmid,
                    'doi': doi,
                    'classification': 'inactivate',
                    'neg_reason': 'negative chembl comment'
                }
                true_negative.append(entry)
                
            elif has_experimental_data(row):
                # Valid CPI entry with experimental data
                entry = {
                    'compound_inchikey': row['compound_inchikey'],
                    'target_uniprot_id': row['target_uniprot_id'],
                    'source_database': 'chembl',
                    'source_compound_id': row['compound_chembl_id'],
                    'standard_type': row.get('standard_type'),
                    'standard_relation': row.get('standard_relation'),
                    'standard_value': row.get('standard_value'),
                    'standard_units': row.get('standard_units'),
                    'activity_comment': activity_comment,
                    'pmid': pmid,
                    'doi': doi,
                    'classification': 'active',
                }
                formed_cpi.append(entry)
                
            else:
                # Uncertain entry (no activity comment and no experimental data)
                entry = {
                    'compound_inchikey': row['compound_inchikey'],
                    'target_uniprot_id': row['target_uniprot_id'],
                    'source_database': 'chembl',
                    'source_compound_id': row['compound_chembl_id'],
                    'standard_type': row.get('standard_type'),
                    'standard_relation': row.get('standard_relation'),
                    'standard_value': row.get('standard_value'),
                    'standard_units': row.get('standard_units'),
                    'activity_comment': activity_comment,
                    'pmid': pmid,
                    'doi': doi,
                    'classification': 'uncertain',
                }
                uncertain.append(entry)
        
        # Convert to DataFrames
        formed_cpi_df = pd.DataFrame(formed_cpi)
        true_negative_df = pd.DataFrame(true_negative)
        uncertain_df = pd.DataFrame(uncertain)
        # true_negative_df保留has_experimental_data的行，其他放入uncertian
        # Print statistics
        print(f"\nChEMBL CPI Processing Results:")
        print(f"  - Valid CPI entries: {len(formed_cpi_df)}")
        print(f"  - True negative entries: {len(true_negative_df)}")
        print(f"  - Uncertain entries: {len(uncertain_df)}")
        
        # Show standard_type distribution for formed_cpi
        if not formed_cpi_df.empty:
            print(f"\nStandard Type Distribution in Valid CPI:")
            type_counts = formed_cpi_df['standard_type'].value_counts()
            for stype, count in type_counts.head(10).items():
                print(f"  - {stype}: {count}")
        
        return formed_cpi_df, true_negative_df, uncertain_df

    
    def _form_iuphar_cpi(self):
        """
        Process IUPHAR CPI data
        
        Returns
        -------
        dict
            Dictionary containing:
            - 'formed_cpi': Valid CPI data with affinity measurements
            - 'uncertain': Entries without valid affinity data or with '-' units
        """
        print("Processing IUPHAR CPI data...")
        
        # Load IUPHAR CPI data
        iuphar_file = join(self.output_dp, 'iuphar/iuphar_cpi.tsv')
        if not os.path.exists(iuphar_file):
            print("IUPHAR CPI file not found")
            return {'formed_cpi': pd.DataFrame(), 'uncertain': pd.DataFrame()}
        
        iuphar_df = pd.read_csv(iuphar_file, sep='\t', encoding='utf-8')
        print(f"Loaded {len(iuphar_df)} IUPHAR CPI entries")
        
        # Skip entries without valid InChIKey
        valid_df = iuphar_df[iuphar_df['compound_inchikey'].notna() & 
                            (iuphar_df['compound_inchikey'].str.strip() != '')].copy()
        
        skipped_count = len(iuphar_df) - len(valid_df)
        if skipped_count > 0:
            print(f"Skipped {skipped_count} entries without valid InChIKey")
        
        def process_affinity_data(row):
            affinity_units = row.get('affinity_units')

            # Check if affinity_units is '-' (uncertain)
            if str(affinity_units).strip() == '-':
                return None, None, None, None, False
            
            affinity_units_str = str(affinity_units).strip()
            
            # Get affinity values (log values)
            affinity_high = row.get('affinity_high')
            affinity_median = row.get('affinity_median')
            affinity_low = row.get('affinity_low')
            
            # Get original values (nM units)
            original_high = row.get('original_affinity_high')
            original_median = row.get('original_affinity_median')
            original_low = row.get('original_affinity_low')
            
            # Map affinity units to standard types
            affinity_type_map = {
                'pKi': 'Ki',
                'pIC50': 'IC50',
                'pKd': 'Kd',
                'pKB': 'KB',
                'pEC50': 'EC50',
                'pA2': 'pA2'
            }
            
            standard_type = affinity_type_map.get(affinity_units_str, affinity_units_str)
            
            # Determine if we should use log values (prefer log values, fallback to original)
            use_original = False
            
            # Check if we have log values first
            if (pd.notna(affinity_high) or pd.notna(affinity_median) or pd.notna(affinity_low)):
                # Use log values (no units for log values)
                high_val = affinity_high if pd.notna(affinity_high) else None
                median_val = affinity_median if pd.notna(affinity_median) else None
                low_val = affinity_low if pd.notna(affinity_low) else None
                units = None
            elif (pd.notna(original_high) or pd.notna(original_median) or pd.notna(original_low)):
                # Fallback to original values (nM units)
                high_val = original_high if pd.notna(original_high) else None
                median_val = original_median if pd.notna(original_median) else None
                low_val = original_low if pd.notna(original_low) else None
                units = 'nM'
                use_original = True
            else:
                # No values available
                return None, None, None, None, False
            
            # Determine relation and value based on available data
            if high_val is not None and low_val is not None:
                # Both high and low values available
                relation = 'in'
                value = f"({low_val}, {high_val})"
            elif median_val is not None:
                # Only median available
                relation = '='
                value = median_val
            elif low_val is not None:
                # Only low value available
                relation = '>'
                value = low_val
            elif high_val is not None:
                # Only high value available
                relation = '<'
                value = high_val
            else:
                # No valid values
                return None, None, None, None, False
            
            standard_type = standard_type if use_original else affinity_units_str
            return standard_type, relation, value, units, use_original
        
        def process_interaction_type(interaction_type):
            if pd.isna(interaction_type):
                return None
            
            interaction_str = str(interaction_type).strip()
            if not interaction_str or interaction_str.lower() == 'nan':
                return None
            
            return f"drug ({interaction_str.lower()})"
        
        def process_pubmed_ids(pubmed_id_str):
            if pd.isna(pubmed_id_str):
                return None
            
            pubmed_str = str(pubmed_id_str).strip()
            if not pubmed_str or pubmed_str.lower() == 'nan':
                return None
            
            # Split by | and join with ,
            pmids = [pmid.strip() for pmid in pubmed_str.split('|') if pmid.strip()]
            return ','.join(pmids) if pmids else None
        
        # Initialize result containers
        formed_cpi = []
        uncertain = []
        
        # Process each entry
        for idx, row in tqdm(valid_df.iterrows(), total=len(valid_df), desc="Processing IUPHAR entries"):
            
            # Process affinity data
            standard_type, relation, value, units, _ = process_affinity_data(row)
            
            # Process interaction type
            action_type = process_interaction_type(row.get('interaction_type'))
            
            # Process PMID
            pmid = process_pubmed_ids(row.get('pubmed_id'))
            
            if standard_type is not None and relation is not None and value is not None:
                # Valid CPI entry with affinity data
                entry = {
                    'compound_inchikey': row['compound_inchikey'],
                    'target_uniprot_id': row['target_uniprot_id'],
                    'source_database': 'iuphar',
                    'source_compound_id': row.get('ligand_id', ''),
                    'standard_type': standard_type,
                    'standard_relation': relation,
                    'standard_value': value,
                    'standard_units': units,
                    'activity_comment': None,  # IUPHAR doesn't have activity comments
                    'action_type': action_type,
                    'pmid': pmid,
                    'classification': 'active',
                }
                formed_cpi.append(entry)
            else:
                # Uncertain entry (no valid affinity data or '-' units)
                entry = {
                    'compound_inchikey': row['compound_inchikey'],
                    'target_uniprot_id': row['target_uniprot_id'],
                    'source_database': 'iuphar',
                    'source_compound_id': row.get('ligand_id', ''),
                    'standard_type': None,
                    'standard_relation': None,
                    'standard_value': None,
                    'standard_units': None,
                    'activity_comment': None,
                    'action_type': action_type,
                    'pmid': pmid,
                    'classification': 'uncertain',
                }
                uncertain.append(entry)
        
        # Convert to DataFrames
        formed_cpi_df = pd.DataFrame(formed_cpi)
        uncertain_df = pd.DataFrame(uncertain)
        
        # Print statistics
        print(f"\nIUPHAR CPI Processing Results:")
        print(f"  - Valid CPI entries: {len(formed_cpi_df)}")
        print(f"  - Uncertain entries: {len(uncertain_df)}")
        
        # Show standard_type distribution for formed_cpi
        if not formed_cpi_df.empty:
            print(f"\nStandard Type Distribution in Valid CPI:")
            type_counts = formed_cpi_df['standard_type'].value_counts()
            for stype, count in type_counts.head(10).items():
                print(f"  - {stype}: {count}")
            
            # Show relation distribution
            print(f"\nRelation Distribution:")
            relation_counts = formed_cpi_df['standard_relation'].value_counts()
            for relation, count in relation_counts.items():
                print(f"  - {relation}: {count}")
            
        return formed_cpi_df, uncertain_df

    def _form_bindingdb_cpi(self):
        """
        Process BindingDB CPI data
        
        Returns
        -------
        dict
            Dictionary containing:
            - 'formed_cpi': Valid CPI data with experimental measurements
            - 'uncertain': Entries without experimental data
        """
        print("Processing BindingDB CPI data...")
        
        # Load BindingDB CPI data
        bindingdb_file = join(self.output_dp, 'bindingdb/bindingdb_cpi.tsv')
        if not os.path.exists(bindingdb_file):
            print("BindingDB CPI file not found")
            return {'formed_cpi': pd.DataFrame(), 'uncertain': pd.DataFrame()}
        
        bindingdb_df = pd.read_csv(bindingdb_file, sep='\t', encoding='utf-8')
        print(f"Loaded {len(bindingdb_df)} BindingDB CPI entries")
        
        # Skip entries without valid InChIKey
        valid_df = bindingdb_df[bindingdb_df['compound_inchikey'].notna() & 
                            (bindingdb_df['compound_inchikey'].str.strip() != '')].copy()
        
        skipped_count = len(bindingdb_df) - len(valid_df)
        if skipped_count > 0:
            print(f"Skipped {skipped_count} entries without valid InChIKey")
        
        # Define experimental data columns
        exp_data_columns = ['ki_nm', 'ic50_nm', 'kd_nm', 'ec50_nm', 'kon_m1s1', 'koff_s1']
        form_dict = {'ki_nm': ('Ki', 'nM'),
                     'ic50_nm': ('IC50', 'nM'),
                     'kd_nm': ('Kd', 'nM'),
                     'ec50_nm': ('EC50', 'nM'),
                     'kon_m1s1': ('Kon', 'M-1s-1'),
                     'koff_s1': ('Koff', 's-1')}


        def extract_value_relation(value_str):
            """
            Extract numeric value and relation from string like '>10', '<=5.2', '100'
            
            Returns
            -------
            tuple
                (relation, numeric_value) or (None, None) if invalid
            """
            if pd.isna(value_str):
                return None, None
            
            value_str = str(value_str).strip()
            if not value_str:
                return None, None
            
            # Check for relation symbols
            if value_str.startswith('>='):
                relation = '>='
                numeric_part = value_str[2:].strip()
            elif value_str.startswith('<='):
                relation = '<='
                numeric_part = value_str[2:].strip()
            elif value_str.startswith('>'):
                relation = '>'
                numeric_part = value_str[1:].strip()
            elif value_str.startswith('<'):
                relation = '<'
                numeric_part = value_str[1:].strip()
            elif value_str.startswith('='):
                relation = '='
                numeric_part = value_str[1:].strip()
            else:
                # No relation symbol, assume equals
                relation = '='
                numeric_part = value_str
            
            # Try to convert to float
            try:
                numeric_value = float(numeric_part)
                return relation, numeric_value
            except (ValueError, TypeError):
                return None, None

        
        def has_experimental_data(row):
            """Check if row has any experimental data"""
            for col in exp_data_columns:
                value = row.get(col)
                if pd.notna(value) and str(value).strip():
                    relation, numeric_value = extract_value_relation(value)
                    if relation is not None and numeric_value is not None:
                        return True
            return False
        
        # Initialize result containers
        formed_cpi = []
        uncertain = []
        
        # Process each entry
        for idx, row in tqdm(valid_df.iterrows(), total=len(valid_df), desc="Processing BindingDB entries"):
            
            if has_experimental_data(row):
                # Process each experimental data column
                for col in exp_data_columns:
                    value = row.get(col)
                    relation, numeric_value = extract_value_relation(value)
                    
                    if relation is not None and numeric_value is not None:
                        standard_type, units = form_dict[col]
                        
                        # Process PMID - convert to int then str
                        pmid = row.get('pmid')
                        pmid_str = None
                        if pd.notna(pmid):
                            try:
                                pmid_int = int(float(pmid))  # Handle potential float representation
                                pmid_str = str(pmid_int)
                            except (ValueError, TypeError):
                                pmid_str = str(pmid) if pmid else None
                        
                        # Process temperature
                        temperature = row.get('temperature_c')
                        if pd.notna(temperature):
                            try:
                                temperature = float(temperature)
                            except (ValueError, TypeError):
                                temperature = None
                        
                        # Process pH
                        ph = row.get('ph')
                        if pd.notna(ph):
                            try:
                                ph = float(ph)
                            except (ValueError, TypeError):
                                ph = None
                        
                        # Create standardized entry
                        entry = {
                            'compound_inchikey': row['compound_inchikey'],
                            'target_uniprot_id': row['target_uniprot_id'],
                            'source_database': 'bindingdb',
                            'source_compound_id': row['compound_bindingdb_id'],
                            'standard_type': standard_type,
                            'standard_relation': relation,
                            'standard_value': numeric_value,
                            'standard_units': units,
                            'activity_comment': None,  # BindingDB doesn't have activity comments
                            'pmid': pmid_str,
                            'doi': row.get('doi'),
                            'pdb_id': row.get('pdb_id'),
                            'ph': ph,
                            'temperature': temperature,
                            'classification': 'active',
                        }
                        formed_cpi.append(entry)
            else:
                # Uncertain entry (no experimental data)
                # Process PMID
                pmid = row.get('pmid')
                pmid_str = None
                if pd.notna(pmid):
                    try:
                        pmid_int = int(float(pmid))
                        pmid_str = str(pmid_int)
                    except (ValueError, TypeError):
                        pmid_str = str(pmid) if pmid else None
                
                entry = {
                    'compound_inchikey': row['compound_inchikey'],
                    'target_uniprot_id': row['target_uniprot_id'],
                    'source_database': 'bindingdb',
                    'source_compound_id': row['compound_bindingdb_id'],
                    'standard_type': None,
                    'standard_relation': None,
                    'standard_value': None,
                    'standard_units': None,
                    'activity_comment': None,
                    'pmid': pmid_str,
                    'doi': row.get('doi'),
                    'pdb_id': row.get('pdb_id'),
                    'ph': row.get('ph'),
                    'temperature': row.get('temperature_c'),
                    'classification': 'uncertain',
                }
                uncertain.append(entry)
        
        # Convert to DataFrames
        formed_cpi_df = pd.DataFrame(formed_cpi)
        uncertain_df = pd.DataFrame(uncertain)
        
        # Print statistics
        print(f"\nBindingDB CPI Processing Results:")
        print(f"  - Valid CPI entries: {len(formed_cpi_df)}")
        print(f"  - Uncertain entries: {len(uncertain_df)}")
        
        # Show standard_type distribution for formed_cpi
        if not formed_cpi_df.empty:
            print(f"\nStandard Type Distribution in Valid CPI:")
            type_counts = formed_cpi_df['standard_type'].value_counts()
            for stype, count in type_counts.head(10).items():
                print(f"  - {stype}: {count}")
            
            # Show relation distribution
            print(f"\nRelation Distribution:")
            relation_counts = formed_cpi_df['standard_relation'].value_counts()
            for relation, count in relation_counts.items():
                print(f"  - {relation}: {count}")
        
        return formed_cpi_df, uncertain_df


    def _form_drugbank_cpi(self):
        """
        Process DrugBank CPI data
        
        Returns
        -------
        pandas.DataFrame
            Standardized DrugBank CPI data
        """
        print("Processing DrugBank CPI data...")
        
        # Load DrugBank CPI data
        drugbank_file = join(self.output_dp, 'drugbank/drugbank_cpi.tsv')
        if not os.path.exists(drugbank_file):
            print("DrugBank CPI file not found")
            return pd.DataFrame()
        
        drugbank_df = pd.read_csv(drugbank_file, sep='\t', encoding='utf-8')
        print(f"Loaded {len(drugbank_df)} DrugBank CPI entries")
        
        # Skip entries without InChIKey
        valid_df = drugbank_df[drugbank_df['compound_inchikey'].notna() & 
                            (drugbank_df['compound_inchikey'].str.strip() != '')].copy()
        
        skipped_count = len(drugbank_df) - len(valid_df)
        if skipped_count > 0:
            print(f"Skipped {skipped_count} entries without valid InChIKey")
        
        # Initialize result container
        formed_cpi = []
        
        # Process each entry
        for idx, row in tqdm(valid_df.iterrows(), total=len(valid_df), desc="Forming DrugBank entries"):
            # Get action and format as drug action
            action = row.get('action')
            standard_type = None
            
            if pd.notna(action) and str(action).strip():
                action_str = str(action).strip()
                standard_type = f"drug ({action_str})"
            
            # Create standardized entry
            entry = {
                'compound_inchikey': row['compound_inchikey'],
                'target_uniprot_id': row['target_uniprot_id'],
                'source_database': 'drugbank',
                'source_compound_id': row['compound_drugbank_id'],
                'action_type': standard_type,
                # 'standard_relation': None,  # DrugBank doesn't have relation
                # 'standard_value': None,     # DrugBank doesn't have experimental values
                # 'standard_units': None,     # DrugBank doesn't have units
                # 'activity_comment': None,   # DrugBank doesn't have activity comments
                'pmid': row.get('pmid'),
                'classification': 'active',
            }
            formed_cpi.append(entry)
        
        # Convert to DataFrame
        formed_cpi_df = pd.DataFrame(formed_cpi)
        
        # Print statistics
        print(f"\nDrugBank CPI Processing Results:")
        print(f"  - Valid CPI entries: {len(formed_cpi_df)}")
        
        # Show standard_type distribution
        if not formed_cpi_df.empty:
            print(f"\nStandard Type Distribution in DrugBank CPI:")
            type_counts = formed_cpi_df['action_type'].value_counts()
            for stype, count in type_counts.head(10).items():
                print(f"  - {stype}: {count}")
        
        return formed_cpi_df

    def _form_pdsp_cpi(self):
        """
        Process PDSP CPI data

        Returns
        -------
        dict
            Dictionary containing:
            - 'formed_cpi': Valid CPI data with Ki measurements
            - 'uncertain': Entries without Ki data
        """
        print("Processing PDSP CPI data...")

        # Load PDSP CPI data
        pdsp_file = join(self.output_dp, 'pdsp/pdsp_cpi.tsv')
        if not os.path.exists(pdsp_file):
            print("PDSP CPI file not found")
            return {'formed_cpi': pd.DataFrame(), 'uncertain': pd.DataFrame()}

        pdsp_df = pd.read_csv(pdsp_file, sep='\t', encoding='utf-8')
        print(f"Loaded {len(pdsp_df)} PDSP CPI entries")

        # Skip entries without valid InChIKey
        valid_df = pdsp_df[pdsp_df['compound_inchikey'].notna() & 
                            (pdsp_df['compound_inchikey'].str.strip() != '')].copy()

        skipped_count = len(pdsp_df) - len(valid_df)
        if skipped_count > 0:
            print(f"Skipped {skipped_count} entries without valid InChIKey")

        def extract_ki_data(ki_nm, ki_note):
            """
            Extract Ki value and relation from ki_nm and ki_note columns
            
            Returns
            -------
            tuple
                (relation, numeric_value) or (None, None) if invalid
            """
            # Check if ki_nm has valid data
            if pd.isna(ki_nm):
                return None, None
                
            try:
                ki_value = float(ki_nm)
            except (ValueError, TypeError):
                return None, None
            
            # Determine relation from ki_note (default to '=' if empty/null)
            if pd.isna(ki_note) or str(ki_note).strip() == '':
                relation = '='
            else:
                ki_note_str = str(ki_note).strip()
                # Map common symbols
                if ki_note_str in ['>', 'gt']:
                    relation = '>'
                elif ki_note_str in ['<', 'lt']:
                    relation = '<'
                elif ki_note_str in ['>=', 'gte']:
                    relation = '>='
                elif ki_note_str in ['<=', 'lte']:
                    relation = '<='
                elif ki_note_str in ['=', 'eq', '']:
                    relation = '='
                else:
                    # For any other note, assume equals but keep the note
                    relation = '='
            
            return relation, ki_value

        # Initialize result containers
        formed_cpi = []
        # Process each entry
        for idx, row in tqdm(valid_df.iterrows(), total=len(valid_df), desc="Processing PDSP entries"):
            
            # Extract Ki data
            ki_nm = row.get('ki_nm')
            ki_note = row.get('ki_note')
            relation, ki_value = extract_ki_data(ki_nm, ki_note)
            
            # Valid CPI entry with Ki data
            entry = {
                'compound_inchikey': row['compound_inchikey'],
                'target_uniprot_id': row['target_uniprot_id'],
                'source_database': 'pdsp',
                'source_compound_id': row.get('compound_pdsp_ligand_id', ''),
                'standard_type': 'Ki',
                'standard_relation': relation,
                'standard_value': ki_value,
                'standard_units': 'nM',
                'activity_comment': None,  # PDSP doesn't have activity comments
                'reference': row.get('reference'),
                'ki_note': ki_note,  # Keep original note for reference
                'classification': 'active',
            }
            formed_cpi.append(entry)


        # Convert to DataFrames
        formed_cpi_df = pd.DataFrame(formed_cpi)

        # Print statistics
        print(f"\nPDSP CPI Processing Results:")
        print(f"  - Valid CPI entries: {len(formed_cpi_df)}")

        # Show relation distribution for formed_cpi
        if not formed_cpi_df.empty:
            print(f"\nRelation Distribution in Valid CPI:")
            relation_counts = formed_cpi_df['standard_relation'].value_counts()
            for relation, count in relation_counts.items():
                print(f"  - {relation}: {count}")
            
            # Show Ki value statistics
            ki_values = formed_cpi_df['standard_value']
            print(f"\nKi Value Statistics:")
            print(f"  - Mean: {ki_values.mean():.2f} nM")
            print(f"  - Median: {ki_values.median():.2f} nM")
            print(f"  - Min: {ki_values.min():.2f} nM")
            print(f"  - Max: {ki_values.max():.2f} nM")

        return formed_cpi_df

    def _form_glass_cpi(self):
        """
        Process GLASS CPI data
        
        Returns
        -------
        pandas.DataFrame
            Standardized GLASS CPI data
        """
        print("Processing GLASS CPI data...")
        
        # Load GLASS CPI data
        glass_file = join(self.output_dp, 'glass/glass_cpi.tsv')
        if not os.path.exists(glass_file):
            print("GLASS CPI file not found")
            raise FileNotFoundError(f"GLASS CPI file not found at {glass_file}")
        
        formed_cpi_df = pd.read_csv(glass_file, sep='\t', encoding='utf-8')
        return formed_cpi_df



    def _merge_duplicate_interactions(self, cpi_df):
        """
        Merge duplicate compound-target interactions from the same or different databases
        
        Parameters
        ----------
        cpi_df : pandas.DataFrame
            Combined CPI data
            
        Returns
        -------
        pandas.DataFrame
            Merged CPI data with duplicate handling
        """
        print("Merging duplicate interactions...")
        
        # Group by compound_uci and target_uniprot_id
        merged_data = []
        
        grouped = cpi_df.groupby(['compound_uci', 'target_uniprot_id'])
        
        for (compound_uci, target_uniprot_id), group in tqdm(grouped, desc="Merging interactions"):
            # Get all source databases for this interaction
            source_databases = group['source_database'].unique().tolist()
            source_ids = {}
            
            # Collect source IDs from each database
            for db in source_databases:
                db_rows = group[group['source_database'] == db]
                source_ids[f'{db}_id'] = '; '.join(db_rows['source_compound_id'].unique().astype(str))
            
            # Merge interaction data (prioritize by database preference)
            merged_entry = {
                'compound_uci': compound_uci,
                'target_uniprot_id': target_uniprot_id,
                'source_databases': '; '.join(sorted(source_databases)),
                'source_count': len(source_databases),
            }
            
            # Add source IDs
            merged_entry.update(source_ids)
            
            # Merge measurement values (take first non-null value, prioritize by database)
            measurement_columns = ['ki_nm', 'ic50_nm', 'kd_nm', 'ec50_nm', 'standard_value']
            for col in measurement_columns:
                if col in group.columns:
                    non_null_values = group[col].dropna()
                    if len(non_null_values) > 0:
                        merged_entry[col] = non_null_values.iloc[0]
            
            # Merge other important fields
            other_fields = ['interaction_type', 'action', 'standard_type', 'standard_units']
            for col in other_fields:
                if col in group.columns:
                    unique_values = group[col].dropna().unique()
                    if len(unique_values) > 0:
                        merged_entry[col] = '; '.join(unique_values.astype(str))
            
            # Merge references
            references = group['reference'].dropna().unique()
            if len(references) > 0:
                merged_entry['references'] = '; '.join(references.astype(str))
            
            merged_data.append(merged_entry)
        
        return pd.DataFrame(merged_data)
    
    def merge_cpi(self):
        """
        Main function to merge all CPI data

        Returns
        -------
        dict
            Dictionary containing merged DataFrames:
            - 'formed_cpi': Valid CPI data with experimental measurements
            - 'inactive_cpi': Negative activity data  
            - 'uncertain_cpi': Uncertain activity data
        """
        print_section_header("Merging Compound-Protein Interaction (CPI) data")
        start = timer()
        
        # Initialize result containers
        all_formed_cpi = []
        all_inactive_cpi = []
        all_uncertain_cpi = []
        
        # Process ChEMBL (returns 3 DataFrames)
        print("=" * 60)
        chembl_formed, chembl_inactive, chembl_uncertain = self._form_chembl_cpi()
        
        if not chembl_formed.empty:
            all_formed_cpi.append(chembl_formed)
        
        if not chembl_inactive.empty:
            all_inactive_cpi.append(chembl_inactive)
        
        if not chembl_uncertain.empty:
            all_uncertain_cpi.append(chembl_uncertain)
        
        # Process IUPHAR (returns 2 DataFrames)
        print("=" * 60)
        iuphar_formed, iuphar_uncertain = self._form_iuphar_cpi()
        
        if not iuphar_formed.empty:
            all_formed_cpi.append(iuphar_formed)
        
        if not iuphar_uncertain.empty:
            all_uncertain_cpi.append(iuphar_uncertain)
        
        # Process BindingDB (returns 2 DataFrames)
        print("=" * 60)
        bindingdb_formed, bindingdb_uncertain = self._form_bindingdb_cpi()
        
        if not bindingdb_formed.empty:
            all_formed_cpi.append(bindingdb_formed)
        
        if not bindingdb_uncertain.empty:
            all_uncertain_cpi.append(bindingdb_uncertain)
        
        # Process DrugBank (returns 1 DataFrame - all are formed_cpi)
        print("=" * 60)
        drugbank_cpi = self._form_drugbank_cpi()
        
        if not drugbank_cpi.empty:
            all_formed_cpi.append(drugbank_cpi)
        
        # Process PDSP (returns 1 DataFrame - all are formed_cpi)
        print("=" * 60)
        pdsp_cpi = self._form_pdsp_cpi()
        
        if not pdsp_cpi.empty:
            all_formed_cpi.append(pdsp_cpi)

        # Process GLASS (returns 1 DataFrame - all are formed_cpi)
        print("=" * 60)
        glass_cpi = self._form_glass_cpi()
        
        if not glass_cpi.empty:
            all_formed_cpi.append(glass_cpi)
        
        # Combine results
        print("=" * 60)
        print("Combining all CPI data...")
        
        results = {}
        
        # Combine formed_cpi
        if all_formed_cpi:
            results['formed_cpi'] = pd.concat(all_formed_cpi, ignore_index=True, sort=False)
            
            # Sort by source database and compound_inchikey
            results['formed_cpi'] = results['formed_cpi'].sort_values(
                ['source_database', 'compound_inchikey'], 
                ascending=[True, True]
            ).reset_index(drop=True)
            
            # Save merged formed CPI
            results['formed_cpi'].to_csv(join(self.output_dp, 'merged_cpi_formed.tsv'), 
                                    sep='\t', index=False, encoding='utf-8')
        else:
            results['formed_cpi'] = pd.DataFrame()
        
        # Combine inactive_cpi  
        if all_inactive_cpi:
            results['inactive_cpi'] = pd.concat(all_inactive_cpi, ignore_index=True, sort=False)
            
            # Sort by source database and compound_inchikey
            results['inactive_cpi'] = results['inactive_cpi'].sort_values(
                ['source_database', 'compound_inchikey'], 
                ascending=[True, True]
            ).reset_index(drop=True)
            
            # Save merged inactive CPI
            results['inactive_cpi'].to_csv(join(self.output_dp, 'merged_cpi_inactive.tsv'), 
                                        sep='\t', index=False, encoding='utf-8')
        else:
            results['inactive_cpi'] = pd.DataFrame()
        
        # Combine uncertain_cpi
        if all_uncertain_cpi:
            results['uncertain_cpi'] = pd.concat(all_uncertain_cpi, ignore_index=True, sort=False)
            
            # Sort by source database and compound_inchikey
            results['uncertain_cpi'] = results['uncertain_cpi'].sort_values(
                ['source_database', 'compound_inchikey'], 
                ascending=[True, True]
            ).reset_index(drop=True)
            
            # Save merged uncertain CPI
            results['uncertain_cpi'].to_csv(join(self.output_dp, 'merged_cpi_uncertain.tsv'), 
                                        sep='\t', index=False, encoding='utf-8')
        else:
            results['uncertain_cpi'] = pd.DataFrame()
        
        # Print final statistics
        print(f"\n" + "=" * 60)
        print(f"Final CPI Merging Statistics:")
        print(f"  - Total formed CPI entries: {len(results['formed_cpi'])}")
        print(f"  - Total inactive CPI entries: {len(results['inactive_cpi'])}")
        print(f"  - Total uncertain CPI entries: {len(results['uncertain_cpi'])}")
        
        # Database coverage statistics for formed_cpi
        if not results['formed_cpi'].empty:
            print(f"\nDatabase Coverage in Formed CPI:")
            db_counts = results['formed_cpi']['source_database'].value_counts()
            for db, count in db_counts.items():
                print(f"  - {db.upper()}: {count} interactions")
        
        # Show standard_type distribution in formed_cpi
        if not results['formed_cpi'].empty:
            print(f"\nStandard Type Distribution in Formed CPI:")
            type_counts = results['formed_cpi']['standard_type'].value_counts()
            for stype, count in type_counts.head(15).items():
                print(f"  - {stype}: {count}")
        
        print(f"\nSaved merged CPI data files:")
        print(f"  - Formed CPI: {join(self.output_dp, 'merged_cpi_formed.tsv')}")
        print(f"  - Inactive CPI: {join(self.output_dp, 'merged_cpi_inactive.tsv')}")
        print(f"  - Uncertain CPI: {join(self.output_dp, 'merged_cpi_uncertain.tsv')}")
        
        print(done_sym + f" Took {timer() - start:.2f} Seconds.")
        
        return results


def id_cpi_pairs(cpi_df):
    start = timer()
    if cpi_df.empty:
        print(f"No CPI data to process")
        return pd.DataFrame()
    
    # Create a copy to avoid modifying original data
    df = cpi_df.copy()
    
    # Create compound-protein pair identifier and assign unique IDs
    df['compound_protein_pair'] = df['compound_inchikey'] + '|' + df['target_uniprot_id']
    unique_pairs = df['compound_protein_pair'].unique()
    pair_to_id = {pair: f"PAIR_{i+1:08d}" for i, pair in enumerate(unique_pairs)}
    df['pair_id'] = df['compound_protein_pair'].map(pair_to_id)
    
    # Sort by pair_id
    df = df.sort_values('pair_id').reset_index(drop=True)
    # move pair_id to the front
    cols = ['pair_id'] + [col for col in df.columns if col != 'pair_id']
    df = df[cols]
    # Remove temporary column
    df = df.drop(columns=['compound_protein_pair'])
    
    print(f"Found {len(unique_pairs)} unique compound-protein pairs")
    print(f"Took {timer() - start:.2f} Seconds.")
    return df

def deduplicate_cpi_pairs(cpi_df):
    print("Deduplicating CPI pairs...")
    start = timer()
    
    if cpi_df.empty:
        print("No CPI data to process")
        return pd.DataFrame()
    
    original_count = len(cpi_df)
    
    # Create grouping key for identical experimental records
    # Fill NaN with empty string for consistent grouping
    df = cpi_df.copy()
    key_cols = ['pair_id', 'standard_type', 'standard_relation', 'standard_value', 'standard_units']
    for col in key_cols:
        if col in df.columns:
            df[col] = df[col].fillna('')
    
    # Group by experimental record identity
    grouped = df.groupby(key_cols, dropna=False)
    
    def merge_string_fields(series, sep=','):
        """Merge comma-separated string fields into unique set"""
        all_values = set()
        for val in series.dropna():
            if pd.notna(val) and str(val).strip():
                # Split by comma and add to set
                all_values.update([v.strip() for v in str(val).split(',') if v.strip()])
        return sep.join(sorted(all_values)) if all_values else None
    
    # Aggregate function for each group
    agg_dict = {}
    
    # Keep first value for non-mergeable fields
    single_value_cols = [col for col in df.columns if col not in key_cols + ['source_database', 'pmid', 'doi']]
    for col in single_value_cols:
        agg_dict[col] = 'first'
    
    # Merge source databases, pmids, and dois
    merge_cols = ['source_database', 'pmid', 'doi']
    for col in merge_cols:
        if col in df.columns:
            agg_dict[col] = lambda x, col=col: merge_string_fields(x)
    
    # Apply aggregation
    deduplicated_df = grouped.agg(agg_dict).reset_index()
    
    # Restore NaN for empty strings in key columns
    for col in key_cols:
        if col in deduplicated_df.columns:
            deduplicated_df[col] = deduplicated_df[col].replace('', pd.NA)
    
    # Reorder columns to match original
    original_cols = [col for col in cpi_df.columns if col in deduplicated_df.columns]
    deduplicated_df = deduplicated_df[original_cols]
    
    # Sort by pair_id
    deduplicated_df = deduplicated_df.sort_values('pair_id').reset_index(drop=True)
    
    duplicates_removed = original_count - len(deduplicated_df)
    
    print(f"Deduplication Results:")
    print(f"  - Original entries: {original_count}")
    print(f"  - Final entries: {len(deduplicated_df)}")
    print(f"  - Duplicates removed: {duplicates_removed}")
    print(f"  - Deduplication rate: {duplicates_removed/original_count*100:.1f}%")
    
    print(done_sym + f" Took {timer() - start:.2f} Seconds.")
    
    return deduplicated_df

def deduplicate_cpi_pairs(cpi_df):
    print("Deduplicating CPI pairs...")
    start = timer()
    
    if cpi_df.empty:
        print("No CPI data to process")
        return pd.DataFrame()
    
    original_count = len(cpi_df)

    # Create grouping key for identical experimental records
    # Fill NaN with empty string for consistent grouping
    print("Preparing data for grouping...")
    df = cpi_df.copy()

    key_cols = ['pair_id', 'standard_type', 'standard_value', 'standard_units']  # 移除standard_relation
    for col in key_cols:
        if col in df.columns:
            df[col] = df[col].fillna('')
    
    # Group by experimental record identity (without relation)
    print("Grouping identical records...")
    grouped = df.groupby(key_cols, dropna=False)
    num_groups = len(grouped)
    print(f"Found {num_groups} unique experimental records")
    
    def merge_string_fields(series, sep=','):
        """Merge comma-separated string fields into unique set"""
        all_values = set()
        for val in series.dropna():
            if pd.notna(val) and str(val).strip():
                # Split by comma and add to set
                all_values.update([v.strip() for v in str(val).split(',') if v.strip()])
        return sep.join(sorted(all_values)) if all_values else None
    
    def merge_relations(series):
        """Merge relations: if different relations exist, use '~', otherwise keep the original"""
        relations = series.fillna('').astype(str)
        unique_relations = set(rel.strip() for rel in relations if rel.strip())
        
        if len(unique_relations) > 1:
            return '~'  # Different relations -> approximately equal
        elif len(unique_relations) == 1:
            return list(unique_relations)[0]  # Same relation
        else:
            return pd.NA  # No relations
        
    def merge_glass1(series):
        """如果有yes就返回yes，否则返回no或原值"""
        vals = series.dropna().astype(str).str.lower()
        if 'yes' in vals.values:
            return 'yes'
        elif 'no' in vals.values:
            return 'no'
        elif len(vals) > 0:
            return vals.iloc[0]
        else:
            return pd.NA
        
    # Aggregate function for each group
    print("Setting up aggregation functions...")
    agg_dict = {}
    
    # Keep first value for non-mergeable fields
    single_value_cols = [col for col in df.columns if col not in key_cols + ['standard_relation', 'source_database', 'pmid', 'doi']]
    for col in single_value_cols:
        agg_dict[col] = 'first'
    
    # Special handling for standard_relation
    if 'standard_relation' in df.columns:
        agg_dict['standard_relation'] = merge_relations
    
    # Special handling for glass1
    if 'glass1' in df.columns:
        agg_dict['glass1'] = merge_glass1
    
    # Merge source databases, pmids, and dois
    merge_cols = ['source_database', 'pmid', 'doi']
    for col in merge_cols:
        if col in df.columns:
            agg_dict[col] = lambda x, col=col: merge_string_fields(x)
    
    # Apply aggregation with progress tracking
    print("Applying deduplication...")
    deduplicated_df = grouped.agg(agg_dict).reset_index()

    
    # Restore NaN for empty strings in key columns
    print("Cleaning up data...")
    for col in key_cols:
        if col in deduplicated_df.columns:
            deduplicated_df[col] = deduplicated_df[col].replace('', pd.NA)
    
    # Reorder columns to match original
    original_cols = [col for col in cpi_df.columns if col in deduplicated_df.columns]
    deduplicated_df = deduplicated_df[original_cols]
    
    # Sort by pair_id
    print("Final sorting...")

    deduplicated_df = deduplicated_df.sort_values('pair_id').reset_index(drop=True)
    
    duplicates_removed = original_count - len(deduplicated_df)
    
    print(f"\nDeduplication Results:")
    print(f"  - Original entries: {original_count}")
    print(f"  - Final entries: {len(deduplicated_df)}")
    print(f"  - Duplicates removed: {duplicates_removed}")
    print(f"  - Deduplication rate: {duplicates_removed/original_count*100:.1f}%")
    
    # Show relation statistics
    if not deduplicated_df.empty and 'standard_relation' in deduplicated_df.columns:
        print(f"\nRelation Distribution after Deduplication:")
        relation_counts = deduplicated_df['standard_relation'].value_counts()
        for relation, count in relation_counts.items():
            print(f"  - {relation}: {count}")
    
    print(done_sym + f" Took {timer() - start:.2f} Seconds.")
    
    return deduplicated_df

def id_cpi_pairs(cpi_df, prefix='PAIR_'):
    """
    Process CPI data to assign unique pair IDs
    
    Parameters
    ----------
    cpi_df : pandas.DataFrame
        CPI DataFrame
    
    Returns
    -------
    pandas.DataFrame
        CPI DataFrame with pair_id added and sorted by pairs
    """
    print("Processing CPI pairs...")
    start = timer()
    
    if cpi_df.empty:
        print(f"No CPI data to process")
        return pd.DataFrame()
    
    # Create a copy to avoid modifying original data
    df = cpi_df.copy()
    
    # Create compound-protein pair identifier and assign unique IDs
    print("Creating pair identifiers...")
    df['compound_protein_pair'] = df['compound_inchikey'] + '|' + df['target_uniprot_id']
    unique_pairs = df['compound_protein_pair'].unique()
    pair_to_id = {pair: f"{prefix}{i+1:08d}" for i, pair in enumerate(unique_pairs)}
    df['pair_id'] = df['compound_protein_pair'].map(pair_to_id)

    # Sort by pair_id
    print("Sorting by pair_id...")
    df = df.sort_values('pair_id').reset_index(drop=True)

    # Move pair_id to the front
    cols = ['pair_id'] + [col for col in df.columns if col != 'pair_id']
    df = df[cols]
    
    # Remove temporary column
    df = df.drop(columns=['compound_protein_pair'])
    
    print(f"Found {len(unique_pairs)} unique compound-protein pairs")
    print(done_sym + f" Took {timer() - start:.2f} Seconds.")
    
    return df


def postprocess_cpi(df):
    df = df.copy()

    # 1. 处理 'Ki high' 和 'Ki low'
    mask_high = df['standard_type'].fillna('').str.lower() == 'ki high'
    mask_low = df['standard_type'].fillna('').str.lower() == 'ki low'
    df.loc[mask_high, 'standard_type'] = 'Ki'
    df.loc[mask_high, 'standard_relation'] = '<'
    df.loc[mask_low, 'standard_type'] = 'Ki'
    df.loc[mask_low, 'standard_relation'] = '>'

    # 替换 standard_units 列中为 '/nM' 的为 'nM'
    df['standard_units'] = df['standard_units'].replace('/nM', 'nM')

    # 2. 处理 standard_relation 为 'in' 的区间
    mask_in = df['standard_relation'].fillna('').str.lower() == 'in'
    interval_rows = df[mask_in]
    new_rows = []
    for _, row in interval_rows.iterrows():
        val = str(row['standard_value']).replace('(', '').replace(')', '').replace('[', '').replace(']', '')
        x, y = [v.strip() for v in val.split(',')]
        # 新行1: >
        row1 = row.copy()
        row1['standard_relation'] = '>'
        row1['standard_value'] = x
        # 新行2: <
        row2 = row.copy()
        row2['standard_relation'] = '<'
        row2['standard_value'] = y
        new_rows.extend([row1, row2])
    # 删除原有的'in'行，添加新行
    df = df[~mask_in]
    if new_rows:
        df = pd.concat([df, pd.DataFrame(new_rows)], ignore_index=True)

    # 3. 处理 pKi/pKd 还原
    def pki_to_nM(pval):
        try:
            return 10 ** (9 - float(pval))
        except Exception:
            return np.nan

    mask_pki = df['standard_type'].fillna('').str.lower().str.contains('pki')
    mask_pkd = df['standard_type'].fillna('').str.lower().str.contains('pkd')

    df.loc[mask_pki, 'standard_value'] = df.loc[mask_pki, 'standard_value'].apply(pki_to_nM)
    df.loc[mask_pki, 'standard_type'] = 'Ki'
    df.loc[mask_pki, 'standard_units'] = 'nM'

    df.loc[mask_pkd, 'standard_value'] = df.loc[mask_pkd, 'standard_value'].apply(pki_to_nM)
    df.loc[mask_pkd, 'standard_type'] = 'Kd'
    df.loc[mask_pkd, 'standard_units'] = 'nM'

    return df



def mark_inactive_by_threshold(df, threshold=10000):
    df = df.copy()
    df['standard_relation'] = df['standard_relation'].fillna('=')
    mask = (
        df['standard_type'].isin(['IC50', 'EC50', 'Ki', 'Kd']) &
        (df['standard_units'] == 'nM')
    )

    # 只将其他类型中 classification 为空的填为 'uncertain'
    df.loc[~mask & df['classification'].isna(), 'classification'] = 'uncertain'

    sub_df = df[mask].copy()
    sub_df['standard_value_num'] = pd.to_numeric(sub_df['standard_value'], errors='coerce')

    def classify(row):
        rel = str(row['standard_relation']).strip()
        val = row['standard_value_num']
        if pd.isna(val):
            return row.get('classification', None)
        if rel in ['=', '==', '', '~']:
            if val > threshold:
                return 'inactive'
            else:
                return 'active'
        elif rel in ['>', '>=', '>>']:
            if val >= threshold:
                return 'inactive'
            else:
                return 'active'
        elif rel in ['<', '<=', '<<']:
            if val > threshold:
                return 'uncertain'
            else:
                return 'active'
        return 'uncertain'

    # 应用分类
    idx = sub_df.index
    df.loc[idx, 'classification'] = sub_df.apply(classify, axis=1)
    return df

def normalize_uniprot_ids(df,
                          current_ids, history2current,
                          id_col='target_uniprot_id'):
    # 构建一个映射函数
    def map_to_current(uid):
        if uid in current_ids:
            return uid
        return history2current.get(uid, None)

    # 应用映射
    df = df.copy()
    df[id_col] = df[id_col].map(map_to_current)

    # 过滤掉无法映射到现用ID的行
    df = df[df[id_col].notna() & df[id_col].isin(current_ids)].reset_index(drop=True)

    return df

# Usage function to be called from main workflow
def merge_all_cpi_data(output_dp, comile_dp='database/glass2'):
    """
    Merge all CPI data using the CPIMerger class
    """
    # ============== uniprot ID normalization ==============
    # 读取json
    gpcr_json_path='database/processed/uniprot/gpcr_entries.json'
    with open(gpcr_json_path, 'r') as f:
        gpcr_dict = json.load(f)

    # 现用ID集合
    current_ids = set(gpcr_dict.keys())

    # 历史ID到现用ID的映射
    history2current = {}
    for cur_id, info in gpcr_dict.items():
        for hist_id in info.get('history_ids', []):
            history2current[hist_id] = cur_id
    # ======================================================

    cpi_merger = CPIMerger(output_dp)
    result = cpi_merger.merge_cpi()
    general_df, inactive_df, uncertain_df = result['formed_cpi'], result['inactive_cpi'], result['uncertain_cpi']

    # 添加 旧版数据区分
    general_df['glass1'] = general_df['glass1'].fillna('no')

    # 处理 general_df
    general_df = id_cpi_pairs(general_df) # 为每个CPI对分配唯一ID
    general_df = postprocess_cpi(general_df) # 后处理CPI数据
    general_df = deduplicate_cpi_pairs(general_df) # 去重CPI对
    general_df = normalize_uniprot_ids(general_df, current_ids, history2current) # 归一化Uniprot ID
    general_df = mark_inactive_by_threshold(general_df, threshold=10000) # 标记不活跃的CPI对

    # 按照 classification 划分
    active_df = general_df[general_df['classification'] == 'active'].copy()
    general_inactive_df = general_df[general_df['classification'] == 'inactive'].copy()
    general_uncertain_df = general_df[general_df['classification'] == 'uncertain'].copy()

    # 处理 inactive_df
    inactive_df = pd.concat([inactive_df, general_inactive_df], ignore_index=True, sort=False) # 合并不活跃的CPI对
    inactive_df = id_cpi_pairs(inactive_df) # 为每个CPI对分配唯一ID
    inactive_df = postprocess_cpi(inactive_df) # 后处理CPI数据
    inactive_df = deduplicate_cpi_pairs(inactive_df) # 去重 CPI对
    inactive_df = normalize_uniprot_ids(inactive_df, current_ids, history2current) # 归一化Uniprot ID

    # 处理 uncertain_df
    uncertain_df = pd.concat([uncertain_df, general_uncertain_df], ignore_index=True, sort=False) # 合并不确定的CPI对
    uncertain_df = id_cpi_pairs(uncertain_df) # 为每个CPI
    uncertain_df = postprocess_cpi(uncertain_df) # 后处理CPI数据
    uncertain_df = deduplicate_cpi_pairs(uncertain_df) # 去重 CPI对
    uncertain_df = normalize_uniprot_ids(uncertain_df, current_ids, history2current) # 归一化Uniprot ID

    # ============== 分类与回归数据集处理 ==============

    # 1. 分类数据集（cls_df）
    # 合并active和inactive
    cls_df = pd.concat([active_df, inactive_df], ignore_index=True, sort=False)
    # classification转为二元值
    cls_df = cls_df.copy()
    cls_df['classification'] = cls_df['classification'].map({'active': 1, 'inactive': 0})
    # 重新打pair_id
    cls_df = id_cpi_pairs(cls_df)
    # 只保留需要的列
    # TODO: 加上 neg_reason
    keep_cols = ['pair_id', 'compound_inchikey', 'target_uniprot_id', 'classification', 'pmid']
    cls_df = cls_df[keep_cols]
    # 合并pmid并去重
    def merge_pmids(series):
        pmids = set()
        for val in series.dropna():
            pmids.update([v.strip() for v in str(val).split(',') if v.strip()])
        return ','.join(sorted(pmids)) if pmids else None
    cls_df = cls_df.groupby(['pair_id', 'compound_inchikey', 'target_uniprot_id', 'classification'], dropna=False).agg({'pmid': merge_pmids}).reset_index()
    # 检查是否有同一个pair_id对应多个标签（1和0），如果有则去掉
    pair_label_counts = cls_df.groupby('pair_id')['classification'].nunique()
    conflict_pairs = pair_label_counts[pair_label_counts > 1].index
    cls_df = cls_df[~cls_df['pair_id'].isin(conflict_pairs)].reset_index(drop=True)
    # 'compound_inchikey' dropna
    cls_df = cls_df.dropna(subset=['compound_inchikey'])

    # 2. 回归数据集（reg_df）
    reg_mask = lambda df: df['standard_type'].isin(['IC50', 'EC50', 'Ki', 'Kd']) & (df['standard_units'] == 'nM')
    reg_act_df = active_df[reg_mask(active_df)].copy()
    reg_act_df = reg_act_df[['pair_id', 'compound_inchikey',
                            'target_uniprot_id', 'standard_type', 'standard_relation', 'standard_value', 
                            'standard_units', 'pmid', 'doi', 'reference']]
    
    reg_inact_df = inactive_df[reg_mask(inactive_df)].copy()
    reg_inact_df = reg_inact_df[['pair_id', 'compound_inchikey',
                                 'target_uniprot_id', 'standard_type', 'standard_relation', 'standard_value', 
                                 'standard_units', 'pmid', 'doi', 'reference']]
    # 'compound_inchikey' dropna
    reg_act_df = reg_act_df.dropna(subset=['compound_inchikey'])
    reg_inact_df = reg_inact_df.dropna(subset=['compound_inchikey'])

    #  ============== 合并完整数据，重编号并同步 ==============
    full_cpi_df = pd.concat([active_df, inactive_df, uncertain_df], ignore_index=True, sort=False)
    full_cpi_df = id_cpi_pairs(full_cpi_df, prefix='GLASS')  # 为每个CPI对分配唯一ID
    # 'compound_inchikey' dropna
    full_cpi_df = full_cpi_df.dropna(subset=['compound_inchikey'])

    # 构建 (compound_inchikey, target_uniprot_id) 到新pair_id的映射
    pair_map = dict(zip(
        zip(full_cpi_df['compound_inchikey'], full_cpi_df['target_uniprot_id']),
        full_cpi_df['pair_id']
    ))

    # 定义一个同步pair_id的函数
    def sync_pair_id(df):
        df = df.copy()
        df['pair_id'] = df.apply(
            lambda row: pair_map.get((row['compound_inchikey'], row['target_uniprot_id']), pd.NA),
            axis=1
        )
        return df

    # 同步更新cls_df, reg_act_df, reg_inact_df的pair_id
    cls_df = sync_pair_id(cls_df)
    reg_act_df = sync_pair_id(reg_act_df)
    reg_inact_df = sync_pair_id(reg_inact_df)

    cls_df.to_csv(join(comile_dp, 'glass2_cls.csv'), index=False, encoding='utf-8')
    reg_act_df.to_csv(join(comile_dp, 'glass2_reg_act.csv'), index=False, encoding='utf-8')
    reg_inact_df.to_csv(join(comile_dp, 'glass2_reg_inact.csv'), index=False, encoding='utf-8')
    full_cpi_df.to_csv(join(comile_dp, 'glass2_full_cpi.tsv'), sep='\t', index=False, encoding='utf-8')

    # 返回所有结果
    return {
        'cls': cls_df,
        'reg_act': reg_act_df,
        'reg_inact': reg_inact_df,
        'full': full_cpi_df,
    }

