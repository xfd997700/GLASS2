# %%
import pandas as pd
from openbabel import pybel
from multiprocessing import Pool, cpu_count
from tqdm import tqdm

def smiles_to_iupac(smiles):
    try:
        mol = pybel.readstring("smi", smiles)
        name = mol.write("iupac").strip()
        return name if name else None
    except Exception:
        return None

if __name__ == "__main__":

    cinfo_path = 'database/glass2/ligands.tsv'
    df = pd.read_csv(cinfo_path, sep='\t', encoding='utf-8')

    # 只处理IUPAC name为空的行
    mask = df['IUPAC name'].isnull() | (df['IUPAC name'].astype(str).str.strip() == '')
    smiles_to_calc = df.loc[mask, 'SMILES']

    num_workers = min(8, cpu_count())
    with Pool(processes=num_workers) as pool:
        iupac_list = list(tqdm(pool.imap(smiles_to_iupac, smiles_to_calc), total=len(smiles_to_calc), desc="Converting SMILES"))

    # 只更新为空的行
    df.loc[mask, 'IUPAC name'] = iupac_list

    print(df.head())
# %%
import pandas as pd
import subprocess

def smiles_to_iupac(smiles):
    try:
        result = subprocess.run(
            ["obabel", "-ismi", "-", "-oname"],  # "-oname"输出IUPAC名称
            input=smiles.encode(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            timeout=5
        )
        # 输出一般就是IUPAC名
        iupac = result.stdout.decode().strip()
        if iupac:
            return iupac
        else:
            return None
    except Exception as e:
        return None

df = pd.DataFrame({'smiles': ['CC(=O)OC1=CC=CC=C1C(=O)O', 'C1CCCCC1', 'invalid_smiles']})
df['IUPAC'] = df['smiles'].apply(smiles_to_iupac)
print(df)
