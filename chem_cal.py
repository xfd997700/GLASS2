# %%

import pandas as pd
import os
import multiprocessing as mp
from functools import partial
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED, rdMolDescriptors, AllChem, SDWriter
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
import numpy as np
from tqdm import tqdm

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
            n_processes=16,
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