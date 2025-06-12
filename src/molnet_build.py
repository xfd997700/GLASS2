import os

from rdkit import Chem
from rdkit.Chem import AllChem, RDConfig, BRICS, FragmentCatalog, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold

fg_without_ca_smart = ['[N;D2]-[C;D3](=O)-[C;D1;H3]', 'C(=O)[O;D1]', 'C(=O)[O;D2]-[C;D1;H3]',
                       'C(=O)-[H]', 'C(=O)-[N;D1]', 'C(=O)-[C;D1;H3]', '[N;D2]=[C;D2]=[O;D1]',
                       '[N;D2]=[C;D2]=[S;D1]', '[N;D3](=[O;D1])[O;D1]', '[N;R0]=[O;D1]', '[N;R0]-[O;D1]',
                       '[N;R0]-[C;D1;H3]', '[N;R0]=[C;D1;H2]', '[N;D2]=[N;D2]-[C;D1;H3]', '[N;D2]=[N;D1]',
                       '[N;D2]#[N;D1]', '[C;D2]#[N;D1]', '[S;D4](=[O;D1])(=[O;D1])-[N;D1]',
                       '[N;D2]-[S;D4](=[O;D1])(=[O;D1])-[C;D1;H3]', '[S;D4](=O)(=O)-[O;D1]',
                       '[S;D4](=O)(=O)-[O;D2]-[C;D1;H3]', '[S;D4](=O)(=O)-[C;D1;H3]', '[S;D4](=O)(=O)-[Cl]',
                       '[S;D3](=O)-[C;D1]', '[S;D2]-[C;D1;H3]', '[S;D1]', '[S;D1]', '[#9,#17,#35,#53]',
                       '[C;D4]([C;D1])([C;D1])-[C;D1]',
                       '[C;D4](F)(F)F', '[C;D2]#[C;D1;H]', '[C;D3]1-[C;D2]-[C;D2]1', '[O;D2]-[C;D2]-[C;D1;H3]',
                       '[O;D2]-[C;D1;H3]', '[O;D1]', '[O;D1]', '[N;D1]', '[N;D1]', '[N;D1]']

fName = os.path.join(RDConfig.RDDataDir, 'FunctionalGroups.txt')
fparams = FragmentCatalog.FragCatParams(1, 6, fName)
fg_without_ca_list = [Chem.MolFromSmarts(smarts) for smarts in fg_without_ca_smart]
fg_with_ca_list = [fparams.GetFuncGroup(i) for i in range(39)]



def extract_substructures_all(smiles: str):
    """
    输入一个smiles，输出一个子结构集合（官能团 + BRICS碎片 + Murcko骨架碎片）
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return set()
    
    substructures = set()

    ### 1. 官能团
    for i, fg_mol in enumerate(fg_with_ca_list):  # 你的代码里的 fg_with_ca_list
        if mol.HasSubstructMatch(fg_mol):
            substructures.add(f"FG_{i}")  # 记录命中的官能团编号，也可以换成更具体名字

    ### 2. BRICS碎片
    brics_fragments = set()
    for bond_atoms, brics_labels in BRICS.FindBRICSBonds(mol):
        # 提取连接的原子索引
        frag = Chem.FragmentOnBonds(mol, [mol.GetBondBetweenAtoms(*bond_atoms).GetIdx()])
        frags = Chem.GetMolFrags(frag, asMols=True)
        for f in frags:
            smiles_frag = Chem.MolToSmiles(f, isomericSmiles=True)
            brics_fragments.add(f"BRICS_{smiles_frag}")
    substructures |= brics_fragments

    # ### 3. Murcko骨架碎片
    # murcko_substructures = return_murcko_leaf_structure(smiles)  # 你的函数
    # for frag_name, atom_indices in murcko_substructures.items():
    #     if frag_name != 0:  # 不是整个分子
    #         substructures.add(f"MURCKO_{frag_name}")

    return substructures