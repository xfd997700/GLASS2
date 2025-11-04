# %%

import gzip
import pandas as pd
import os
import json
import multiprocessing as mp
from functools import partial
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED, rdMolDescriptors, AllChem, SDWriter
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import numpy as np
from tqdm import tqdm


def complete_ligands_data(ligands_path='database/glass2/ligands.tsv', 
                         xref_path='database/glass2/uci_xref.json',
                         output_path=None):
    """
    补全ligands.tsv文件中缺失的InChI、Name和IUPAC name字段
    
    Args:
        ligands_path: ligands.tsv文件路径
        xref_path: uci_xref.json文件路径  
        output_path: 输出文件路径，如果为None则覆盖原文件
    """
    
    print("📖 读取ligands.tsv文件...")
    df = pd.read_csv(ligands_path, sep='\t')
    original_count = len(df)
    print(f"原始数据: {original_count} 条记录")
    
    # 统计缺失情况
    missing_inchi = df['InChI'].isna().sum()
    missing_name = df['Name'].isna().sum()
    missing_iupac = df['IUPAC name'].isna().sum()
    
    print(f"缺失统计: InChI={missing_inchi}, Name={missing_name}, IUPAC={missing_iupac}")
    
    # 1. 补全InChI字段
    print("\n🔧 补全InChI字段...")
    inchi_filled = 0
    
    # 用InChI_unichem填补空缺的InChI
    mask_empty_inchi = df['InChI'].isna()
    mask_has_unichem = ~df['InChI_unichem'].isna()
    fill_mask = mask_empty_inchi & mask_has_unichem
    
    df.loc[fill_mask, 'InChI'] = df.loc[fill_mask, 'InChI_unichem']
    inchi_filled += fill_mask.sum()
    
    # 从SMILES计算剩余的InChI
    still_missing_inchi = df['InChI'].isna()
    if still_missing_inchi.sum() > 0:
        print(f"从SMILES计算InChI: {still_missing_inchi.sum()} 个化合物")
        for idx in tqdm(df[still_missing_inchi].index, desc="计算InChI"):
            smiles = df.loc[idx, 'SMILES']
            if pd.notna(smiles):
                try:
                    mol = Chem.MolFromSmiles(str(smiles))
                    if mol:
                        inchi = Chem.MolToInchi(mol)
                        df.loc[idx, 'InChI'] = inchi
                        inchi_filled += 1
                except Exception as e:
                    print(f"计算InChI失败 (索引{idx}): {e}")
    
    print(f"✅ InChI补全完成: {inchi_filled} 个")
    
    # 2. 准备PubChem数据补全Name和IUPAC
    if missing_name > 0 or missing_iupac > 0:
        print("\n📚 读取UCI交叉引用数据...")
        with open(xref_path, 'r') as f:
            uci_xref = json.load(f)
        
        # 收集需要查找的PubChem IDs
        pubchem_ids_for_names = []
        pubchem_ids_for_iupac = []
        uci_to_idx = {}  # UCI到DataFrame索引的映射
        
        for idx, row in df.iterrows():
            uci = 'UCI:' + str(row['UCI'])
            need_name = pd.isna(row['Name'])
            need_iupac = pd.isna(row['IUPAC name'])
            
            if (need_name or need_iupac) and uci in uci_xref:
                xref_data = uci_xref[uci]
                if 'pubchem' in xref_data and xref_data['pubchem']:
                    pubchem_id = int(xref_data['pubchem'][0])
                    uci_to_idx[uci] = idx
                    
                    if need_name:
                        pubchem_ids_for_names.append((uci, pubchem_id))
                    if need_iupac:
                        pubchem_ids_for_iupac.append((uci, pubchem_id))
        
        # 按PubChem ID排序以优化查找
        pubchem_ids_for_names.sort(key=lambda x: x[1])
        pubchem_ids_for_iupac.sort(key=lambda x: x[1])
        
        print(f"需要查找Name的PubChem IDs: {len(pubchem_ids_for_names)}")
        print(f"需要查找IUPAC的PubChem IDs: {len(pubchem_ids_for_iupac)}")
        
        # 3. 查找Names
        if pubchem_ids_for_names:
            found_names = read_pubchem_titles_optimized(pubchem_ids_for_names)
            name_filled = 0
            for uci, name in found_names.items():
                idx = uci_to_idx[uci]
                df.loc[idx, 'Name'] = name
                name_filled += 1
            print(f"✅ Name补全完成: {name_filled} 个")
        
        # 4. 查找IUPAC names
        if pubchem_ids_for_iupac:
            found_iupac = read_pubchem_iupac_optimized(pubchem_ids_for_iupac)
            print('found_iupac:', len(found_iupac))
            iupac_filled = 0
            for uci, iupac in found_iupac.items():
                idx = uci_to_idx[uci]
                df.loc[idx, 'IUPAC name'] = iupac
                iupac_filled += 1
                if iupac_filled % 100 == 0:
                    print(f"已补全 {iupac_filled} 个 IUPAC name...")
            print(f"✅ IUPAC补全完成: {iupac_filled} 个")
        
        # 5. 用IUPAC填补空缺的Name
        mask_empty_name = df['Name'].isna()
        mask_has_iupac = ~df['IUPAC name'].isna()
        iupac_to_name_mask = mask_empty_name & mask_has_iupac
        
        df.loc[iupac_to_name_mask, 'Name'] = df.loc[iupac_to_name_mask, 'IUPAC name']
        iupac_to_name_count = iupac_to_name_mask.sum()
        print(f"✅ 用IUPAC填补Name: {iupac_to_name_count} 个")
    
    # 6. 保存结果
    if output_path is None:
        output_path = ligands_path
    
    df.to_csv(output_path, sep='\t', index=False)
    
    # 统计最终结果
    final_missing_inchi = df['InChI'].isna().sum()
    final_missing_name = df['Name'].isna().sum()
    final_missing_iupac = df['IUPAC name'].isna().sum()
    
    print(f"\n📊 补全结果统计:")
    print(f"InChI: {missing_inchi} -> {final_missing_inchi} (补全 {missing_inchi - final_missing_inchi})")
    print(f"Name: {missing_name} -> {final_missing_name} (补全 {missing_name - final_missing_name})")
    print(f"IUPAC: {missing_iupac} -> {final_missing_iupac} (补全 {missing_iupac - final_missing_iupac})")
    print(f"✅ 结果已保存至: {output_path}")

def read_pubchem_titles_optimized(uci_pubchem_ids):
    """优化的PubChem titles读取函数"""
    if not uci_pubchem_ids:
        return {}
    
    pubchem_to_uci = {pid: uci for uci, pid in uci_pubchem_ids}
    target_pubchem_ids = [pid for _, pid in uci_pubchem_ids]
    
    found_names = {}
    target_idx = 0
    
    file_path = '../UniBioMap/database/sources/pubchem/CID-Title.gz'
    if not os.path.exists(file_path):
        print(f"⚠️  文件不存在: {file_path}")
        return {}
    
    file_size = os.path.getsize(file_path)
    
    with gzip.open(file_path, 'rt') as f:
        with tqdm(total=file_size, unit='B', unit_scale=True, desc="读取PubChem Titles") as pbar:
            while target_idx < len(target_pubchem_ids):
                line = f.readline()
                if not line:
                    break
                
                pbar.update(len(line.encode('utf-8')))
                
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split('\t', 1)
                if len(parts) != 2:
                    continue
                
                try:
                    cid = int(parts[0])
                    title = parts[1]
                    
                    current_target = target_pubchem_ids[target_idx]
                    
                    if cid == current_target:
                        uci_id = pubchem_to_uci[cid]
                        found_names[uci_id] = title
                        target_idx += 1
                    elif cid > current_target:
                        # 跳过不存在的目标
                        while target_idx < len(target_pubchem_ids) and target_pubchem_ids[target_idx] < cid:
                            target_idx += 1
                        # 检查当前CID是否匹配新目标
                        if target_idx < len(target_pubchem_ids) and cid == target_pubchem_ids[target_idx]:
                            uci_id = pubchem_to_uci[cid]
                            found_names[uci_id] = title
                            target_idx += 1
                except ValueError:
                    continue
    
    return found_names

def read_pubchem_iupac_optimized(uci_pubchem_ids):
    """优化的PubChem IUPAC读取函数"""
    if not uci_pubchem_ids:
        return {}
    
    pubchem_to_uci = {pid: uci for uci, pid in uci_pubchem_ids}
    target_pubchem_ids = [pid for _, pid in uci_pubchem_ids]
    
    found_iupac = {}
    target_idx = 0
    
    file_path = '../UniBioMap/database/sources/pubchem/CID-IUPAC.gz'
    if not os.path.exists(file_path):
        print(f"⚠️  文件不存在: {file_path}")
        return {}
    
    file_size = os.path.getsize(file_path)
    
    with gzip.open(file_path, 'rt') as f:
        with tqdm(total=file_size, unit='B', unit_scale=True, desc="读取PubChem IUPAC") as pbar:
            while target_idx < len(target_pubchem_ids):
                line = f.readline()
                if not line:
                    break
                
                pbar.update(len(line.encode('utf-8')))
                
                line = line.strip()
                if not line:
                    continue
                
                parts = line.split('\t', 1)
                if len(parts) != 2:
                    continue
                
                try:
                    cid = int(parts[0])
                    iupac = parts[1]
                    
                    current_target = target_pubchem_ids[target_idx]
                    
                    if cid == current_target:
                        uci_id = pubchem_to_uci[cid]
                        found_iupac[uci_id] = iupac
                        target_idx += 1
                    elif cid > current_target:
                        # 跳过不存在的目标
                        while target_idx < len(target_pubchem_ids) and target_pubchem_ids[target_idx] < cid:
                            target_idx += 1
                        # 检查当前CID是否匹配新目标
                        if target_idx < len(target_pubchem_ids) and cid == target_pubchem_ids[target_idx]:
                            uci_id = pubchem_to_uci[cid]
                            found_iupac[uci_id] = iupac
                            target_idx += 1
                except ValueError:
                    continue
    
    return found_iupac

def calc_lipinski(mol):
    props = {
        'MW': Descriptors.MolWt(mol),
        'LogP': Descriptors.MolLogP(mol),
        'HBA': rdMolDescriptors.CalcNumHBA(mol),
        'HBD': rdMolDescriptors.CalcNumHBD(mol),
        'RotB': Lipinski.NumRotatableBonds(mol),
    }
    failed = []
    if props['MW'] > 500: failed.append('MW')
    if props['LogP'] > 5: failed.append('LogP')
    if props['HBA'] > 10: failed.append('HBA')
    if props['HBD'] > 5: failed.append('HBD')
    return failed, props

def calc_qed_details(mol):
    props = {}
    qed_props = QED.properties(mol)
    for k, v in qed_props._asdict().items():
        props[f'QED_{k}'] = v
    props['QED'] = QED.qed(mol)
    return props

def smiles_to_sdf(mol, filename, title=None):
    mol_withH = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol_withH, randomSeed=0xf00d)
    if title:
        mol_withH.SetProp('_Name', title)
    writer = SDWriter(filename)
    writer.write(mol_withH)
    writer.close()

def mol_to_svg(mol, filename, width=300, height=300):
    """将分子转换为SVG格式"""
    try:
        drawer = rdMolDraw2D.MolDraw2DSVG(width, height)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()
        svg = drawer.GetDrawingText()
        with open(filename, 'w') as f:
            f.write(svg)
        return True
    except Exception as e:
        print(f"SVG generation failed: {e}")
        return False

def process_single_row(args):
    """处理单行数据的函数"""
    row_data, output_sdf_dir, output_svg_dir = args
    idx, row = row_data
    
    smiles = str(row['SMILES'])
    uci = str(row['UCI'])
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'UCI': uci, 'SMILES': smiles, 'MolFormula': None, 'MolWt': None, 
                   'LipinskiFail': 'InvalidSmiles', 'QED': None, 'SVG_Generated': False}
        
        # 分子式和分子量
        mol_formula = rdMolDescriptors.CalcMolFormula(mol)
        mol_wt = Descriptors.MolWt(mol)
        
        # Lipinski
        lip_fails, lip_props = calc_lipinski(mol)
        
        # QED
        qed_props = calc_qed_details(mol)
        
        # UCI ID处理
        uci_id = uci.replace('UCI:', '')
        
        # SDF生成
        sdf_filename = os.path.join(output_sdf_dir, f'{uci_id}.sdf')
        smiles_to_sdf(mol, sdf_filename, title=str(uci))
        
        # SVG生成
        svg_filename = os.path.join(output_svg_dir, f'{uci_id}.svg')
        svg_success = mol_to_svg(mol, svg_filename)
        
        # 汇总结果
        result = {
            'UCI': uci,
            'SMILES': smiles,
            'MolFormula': mol_formula,
            'MolWt': mol_wt,
            'LipinskiFail': ','.join(lip_fails) if lip_fails else 'Pass',
            'MW': lip_props['MW'],
            'LogP': lip_props['LogP'],
            'HBA': lip_props['HBA'],
            'HBD': lip_props['HBD'],
            'RotB': lip_props['RotB'],
            'SVG_Generated': svg_success
        }
        result.update(qed_props)
        return result
        
    except Exception as e:
        print(f"Error processing {uci}: {e}")
        return {'UCI': uci, 'SMILES': smiles, 'MolFormula': None, 'MolWt': None, 
               'LipinskiFail': 'ProcessingError', 'QED': None, 'SVG_Generated': False}

def process_df_parallel(df, output_sdf_dir='./sdf', output_svg_dir='./svg', n_processes=None):
    """并行处理DataFrame（不分块，直接全量处理）"""
    if n_processes is None:
        n_processes = min(mp.cpu_count(), 8)
    print(f"Using {n_processes} processes for parallel processing...")

    os.makedirs(output_sdf_dir, exist_ok=True)
    os.makedirs(output_svg_dir, exist_ok=True)

    args_list = [(row_data, output_sdf_dir, output_svg_dir) for row_data in df.iterrows()]

    with mp.Pool(processes=n_processes) as pool:
        results = list(tqdm(
            pool.imap(process_single_row, args_list),
            total=len(args_list),
            desc="Parallel Processing"
        ))

    return pd.DataFrame(results)

def process_df(df, output_sdf_dir='./sdf', output_svg_dir='./svg'):
    """保持原有的串行处理函数，用于小数据集"""
    os.makedirs(output_sdf_dir, exist_ok=True)
    os.makedirs(output_svg_dir, exist_ok=True)
    
    results = []
    for idx, row in tqdm(df.iterrows(), total=len(df), desc="Processing molecules"):
        smiles = str(row['SMILES'])
        uci = str(row['UCI'])
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            results.append({'UCI': uci, 'SMILES': smiles, 'MolFormula': None, 'MolWt': None, 
                            'LipinskiFail': 'InvalidSmiles', 'QED': None, 'SVG_Generated': False})
            continue
        
        # 分子式和分子量
        mol_formula = rdMolDescriptors.CalcMolFormula(mol)
        mol_wt = Descriptors.MolWt(mol)
        
        # Lipinski
        lip_fails, lip_props = calc_lipinski(mol)
        
        # QED
        qed_props = calc_qed_details(mol)
        
        # UCI ID处理
        uci_id = uci.replace('UCI:', '')
        
        # SDF生成
        sdf_filename = os.path.join(output_sdf_dir, f'{uci_id}.sdf')
        smiles_to_sdf(mol, sdf_filename, title=str(uci))
        
        # SVG生成
        svg_filename = os.path.join(output_svg_dir, f'{uci_id}.svg')
        svg_success = mol_to_svg(mol, svg_filename)
        
        # 汇总
        row_data = {
            'UCI': uci,
            'SMILES': smiles,
            'MolFormula': mol_formula,
            'MolWt': mol_wt,
            'LipinskiFail': ','.join(lip_fails) if lip_fails else 'Pass',
            'MW': lip_props['MW'],
            'LogP': lip_props['LogP'],
            'HBA': lip_props['HBA'],
            'HBD': lip_props['HBD'],
            'RotB': lip_props['RotB'],
            'SVG_Generated': svg_success
        }
        row_data.update(qed_props)
        results.append(row_data)
    
    return pd.DataFrame(results)

if __name__ == "__main__":
    cinfo_path = 'database/glass2/ligands.tsv'
    output_sdf_dir = 'database/glass2/sdf'
    output_svg_dir = 'database/glass2/svg'
    output_csv_path = 'database/glass2/ligands_prop.tsv'

    complete_ligands_data(
        ligands_path=cinfo_path,
        xref_path='database/glass2/uci_xref.json',
        output_path=None  # 或者None来覆盖原文件
    )
    # 创建输出目录
    os.makedirs(output_sdf_dir, exist_ok=True)
    os.makedirs(output_svg_dir, exist_ok=True)

    # 读取数据
    print("Loading data...")
    df = pd.read_csv(cinfo_path, sep='\t', usecols=['UCI', 'SMILES'])
    print(f"Loaded {len(df)} molecules")

    # 根据数据量选择处理方式
    if len(df) > 10000:  # 大于1万条记录使用并行处理
        print("Using parallel processing for large dataset...")
        chem_prop_df = process_df_parallel(
            df, 
            output_sdf_dir=output_sdf_dir, 
            output_svg_dir=output_svg_dir,
            n_processes=24,
        )
    else:
        print("Using serial processing for small dataset...")
        chem_prop_df = process_df(df, output_sdf_dir, output_svg_dir)

    # 保存结果
    print("Saving results...")
    chem_prop_df.to_csv(output_csv_path, index=False, sep='\t', encoding='utf-8')
    
    # 统计结果
    svg_success_count = chem_prop_df['SVG_Generated'].sum()
    total_count = len(chem_prop_df)
    print(f"Done! Results saved to {output_csv_path}")
    print(f"SDFs saved in {output_sdf_dir}/")
    print(f"SVGs saved in {output_svg_dir}/ ({svg_success_count}/{total_count} successful)")