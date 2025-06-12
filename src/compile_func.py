from collections import defaultdict
import json
import sys
import time
from biodblinker import KEGGLinker, SiderLinker, MESHLinker, MimLinker
import os
from os import makedirs, listdir, remove, walk
from os.path import join, isdir
from shutil import copy
import gzip
import pandas as pd
import requests
from tqdm import tqdm
import numpy as np
from rdkit import Chem

kegg_linker = KEGGLinker()
sider_linker = SiderLinker()
mesh_linker = MESHLinker()
mim_linker = MimLinker()

data_root = 'database/processed'
source_root = 'database/sources'
output_root = 'database/output'
core_root = 'database/core_kg'
links_root = join(output_root, 'links')
sublinks_root = join(output_root, 'sublinks')
meta_root = join(output_root, 'metadata')
properties_root = join(output_root, 'properties')
drug_properties_root = join(properties_root, 'drug')
protein_properties_root = join(properties_root, 'protein')
pathway_properties_root = join(properties_root, 'pathway')
disease_properties_root = join(properties_root, 'disease')
cell_properties_root = join(properties_root, 'cell')
mim_properties_root = join(properties_root, 'genetic_disorders')
other_root = join(output_root, 'other')

dirs = [
    output_root, links_root, sublinks_root, meta_root, properties_root, drug_properties_root,
    protein_properties_root, pathway_properties_root, disease_properties_root,
    cell_properties_root, mim_properties_root, other_root, core_root
]
for d in dirs:
    makedirs(d, exist_ok=True)

def get_unique_inter_src(inters,
                         label=None):
    """
    Get the set of unique interactions with sources

    Parameters
    ----------
    inters: list
        the list of interactions to filter
    label: str
        the relation label to use for the interaction
    """
    grouped = defaultdict(set)
    for s, o, src in inters:
        grouped[(s, o)].add(src)
    if label:
        inters = [(s, label, o, ';'.join(src_list)) for (s, o), src_list in grouped.items()]
    else:
        inters = [(s, o, ';'.join(src_list)) for (s, o), src_list in grouped.items()]
    return inters

def get_all_proteins():
    """
    Get the set of uniprot proteins to include in the kg

    Parameters
    ----------

    Returns
    -------
    set
        the set of proteins to include in the kg
    """
    set_path = join(data_root, 'uniprot', 'protein_set.json')
    assert os.path.exists(set_path), f"Protein set not found at {set_path},\
        please check if UniKGParser._parse_proteins() has been run correctly in previous steps."
    with open(set_path, 'r') as fd:
        protein_set = json.load(fd)
    return protein_set


def get_proteins_by_metadata(filter_field, accpeted_values):
    """
    Get the set of uniprot proteins for which the metadata value of field
    is ont of the accpeted_values

    Parameters
    ----------
    filter_field: str
        The metadata field to filter by
    accepted_values: list
        The set of accpetable values for the metadata field

    Returns
    -------
    set
        the set of proteins to include in the kg
    """
    uniprot_meta = join(data_root, 'uniprot', 'uniprot_metadata.txt')
    protein_set = set()
    with open(uniprot_meta, 'r') as fd:
        for line in fd:
            protein, field, value = line.strip().split('\t')
            if field == filter_field and value in accpeted_values:
                protein_set.add(protein)

    return protein_set

def get_all_drugs():
    """
    Get the set of drugs to include in the kg

    Parameters
    ----------

    Returns
    -------
    set
        the set of drugs to include in the kg
    """
    with open(join(data_root, 'drugbank', 'db2inchikey.json'), 'r') as fd:
        db2inchikey = json.load(fd)

    db_meta = join(data_root, 'drugbank', 'db_meta.txt')
    db_drugs = set()

    with open(db_meta, 'r') as fd:
        for line in fd:
            drug = line.split('\t')[0]
            if drug in db2inchikey:  
                db_drugs.add(db2inchikey[drug])

    return db_drugs

def get_all_compounds():
    """
    Get the set of compounds to include in the kg
    Parameters
    ----------
    Returns
    -------
    set
        the set of compounds to include in the kg
    """
    set_path = join(data_root, 'unichem', 'inchikey2uci.json')
    assert os.path.exists(set_path), f"Compound set not found at {set_path},\
        please check if UniKGParser._parse_compounds() has been run correctly in previous steps."
    with open(set_path, 'r') as fd:
        chem_dict = json.load(fd)
    chem_set = set(chem_dict.keys())
    return chem_dict, chem_set


def regular_type(df):
    numeric_cols = df.select_dtypes(include=[np.number]).columns
    df[numeric_cols] = df[numeric_cols].astype("Int64").astype(str).replace("<NA>", np.nan)  # 兼容 Pandas NA
    return df


def merge_compound_df(dfs, key_column='InChIKey'):
    df = pd.concat(dfs)

    exclude_cols = {key_column, 'SMILES', 'Name'}

    agg_funcs = {col: lambda x: ';'.join(sorted(set(x.dropna()))) for col in df.columns if col not in exclude_cols}
    agg_funcs.update({col: lambda x: x.iloc[0] for col in df.columns if col == 'SMILES'})

    tqdm.pandas()
    merged_df = df.groupby(key_column, as_index=False).progress_apply(lambda x: x.agg(agg_funcs))

    return merged_df


def get_all_mesh_diseases():
    """
    Get the set of mesh diseases to include in the kg

    Parameters
    ----------

    Returns
    -------
    set
        the set of mesh diseases to include in the kg
    """
    files = [
        join(data_root, 'mesh', 'mesh_disease_meta.txt'),
        join(data_root, 'mesh', 'mesh_scr_disease_meta.txt')
    ]

    mesh_diseases = set()
    for fp in files:
        with open(fp, 'r') as fd:
            for line in fd:
                disease, meta, value = line.strip().split('\t')
                if meta == 'TYPE':
                    mesh_diseases.add(disease)

    return mesh_diseases

def get_umls_disease_mapping():
    """
    Get the mapping from MeSH OMIM ORPHANET to UMLS CUI
    """
    with open(join(data_root, 'umls', 'mesh2cui.json'), 'r') as f:
        mesh2cui = json.load(f)
    with open(join(data_root, 'umls', 'omim2cui.json'), 'r') as f:
        omim2cui = json.load(f)
    with open(join(data_root, 'umls', 'orpha2cui.json'), 'r') as f:
        orpha2cui = json.load(f)

    def intern_dict(d):
        return {sys.intern(k): 'CUI:'+v for k, v in d.items()}

    mesh2cui = intern_dict(mesh2cui)
    omim2cui = intern_dict(omim2cui)
    orpha2cui = intern_dict(orpha2cui)

    umls_maps = dict({
        'mesh': mesh2cui,
        'omim': omim2cui,
        'orpha': orpha2cui
    })
    return umls_maps



def get_all_genetic_disorders():
    """
    Get the set of uniprot proteins for which the metadata value of field
    is ont of the accpeted_values

    Parameters
    ----------
    filter_field: str
        The metadata field to filter by
    accepted_values: list
        The set of accpetable values for the metadata field

    Returns
    -------
    set
        the set of proteins to include in the kg
    """
    # uniprot_meta = join(data_root, 'uniprot', 'uniprot_facts.txt')
    # genetic_disorders = set()
    # with open(uniprot_meta, 'r') as fd:
    #     for line in fd:
    #         protein, field, value = line.strip().split('\t')
    #         if field == 'RELATED_GENETIC_DISORDER':
    #             genetic_disorders.add(value)

    # # TODO: HPO疾病暂且按biokg的做法记录为genetic disorder，后续考虑使用 UMLS_CUI统一表示疾病
    # hpo_disease = join(data_root, 'hpo', 'phenotype_disease.txt')
    # with open(hpo_disease, 'r') as fd:
    #     for line in fd:
    #         _, _, disease = line.strip().split('\t')
    #         genetic_disorders.add(disease)

    genetic_disorders = set()
    med_mim = join(data_root, 'medgen', 'mim_categories.txt')
    with open(med_mim, 'r') as fd:
        for line in fd:
            disease, _, _ = line.strip().split('\t')
            genetic_disorders.add(disease)

    return genetic_disorders


def get_all_unique_ppi(protein_set):
    """
    Get the set of protein-protein interactions to include in the kg with optimized processing.
    """
    # Intern strings in protein_set for faster lookups
    protein_set = {sys.intern(p) for p in protein_set}
    
    files = {
        'uniprot': join(data_root, 'uniprot', 'uniprot_ppi.txt'),
        'reactome': join(data_root, 'reactome', 'reactome_ppi.txt'),
        'intact': join(data_root, 'intact', 'intact_ppi.txt'),
        'string': join(data_root, 'string', 'string_ppi.txt'),
    }
    ppis = set()

    for src, fp in files.items():
        with open(fp, 'r') as fd:
            for line in fd:
                # Efficiently split the line into first and third columns
                parts = line.strip().split('\t')
                if len(parts) < 3:
                    continue  # Skip malformed lines
                s_raw = parts[0].strip()
                o_raw = parts[2].strip()
                # Intern strings to speed up comparisons and set lookups
                s = sys.intern(s_raw)
                o = sys.intern(o_raw)
                # Early filtering: skip if either protein is not in the set
                if s not in protein_set or o not in protein_set:
                    continue
                # Ensure consistent order to avoid duplicates
                if s > o:
                    ppis.add((o, s, src))
                else:
                    ppis.add((s, o, src))
    
    # Directly convert filtered PPI entries to the final format
    unique_triples = get_unique_inter_src(ppis, 'PPI')
    # unique_triples = [(s, 'PPI', o, src) for s, o, src in ppis]
    
    return unique_triples


def get_species_map():
    species_map = {}
    with open(join(data_root, 'uniprot', 'uniprot_metadata.txt'), 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            if p == 'SPECIES':
                species_map[s] = o

    return species_map


# def write_ppi_by_species(ppis):
#     """
#     Get the set of protein proteins interactions to include in the kg

#     Parameters
#     ----------
#     protein_set: set
#         the set of proteins used to filter the ppis
#     Returns
#     -------
#     list
#         the list of unique protein protein interactions
#     """
#     species_map = get_species_map()
#     unique_species = list(set(species_map.values()))
#     species_root = join(links_root, 'PPI_SPECIES')
#     makedirs(species_root) if not isdir(species_root) else None

#     species_file_names = {}
#     for species in set(species_map.values()):
#         species_file_names[species] = join(species_root, f'{species}_ppi.txt')

#     #species_files['OTHER'] = open(join(species_root, 'INTERSPECIES_ppi.txt'), 'w')
    
#     # files = [
#     #     join(data_root, 'uniprot', 'uniprot_ppi.txt'),
#     #     join(data_root, 'reactome', 'reactome_ppi.txt'),
#     #     join(data_root, 'intact', 'intact_ppi.txt')
#     # ]
#     # ppis = set()

#     # for fp in files:
#     #     with open(fp, 'r') as fd:
#     #         for line in fd:
#     #             s, _, o = line.strip().split('\t')[:3]
#     #             if s > o:
#     #                 ppis.add((o, s))
#     #             else:
#     #                 ppis.add((s, o))

#     # unique_triples = []
#     species_files = {}
#     for s, _, o, _ in ppis:
#         if s not in species_map or o not in species_map:
#             continue
#         s_species = species_map[s]
#         o_species = species_map[o]
#         if s_species == o_species:
#             if s_species not in species_files:
#                 species_files[s_species] = open(species_file_names[s_species], 'w')
#             species_files[s_species].write(f'{s}\tPPI\t{o}\n')
#         else:
#             if 'OTHER' not in species_files:
#                 species_files['OTHER'] = open(join(species_root, 'INTERSPECIES_ppi.txt'), 'w')
#             species_files['OTHER'].write(f'{s}\tPPI\t{o}\n')

#     for f in species_files.values():
#         f.close()


def get_all_protein_sequence_annotations(protein_set):
    """
    Get the set of protein sequence annotations to include in the kg

    Parameters
    ----------
    protein_set: set
        the set of proteins used to filter the ppis
    Returns
    -------
    list
        the list of interpro protein sequence annotations
    list
        the list of prosite protein sequence annotations
    """
    protein_set = {sys.intern(p) for p in protein_set}
    seq_ann_root = join(protein_properties_root, 'sequence_annotations')
    makedirs(seq_ann_root) if not isdir(seq_ann_root) else None
    annotation_output_files = {
        'ACTIVE_SITE': open(join(seq_ann_root, 'protein_active_site.txt'), 'w'),
        'BINDING_SITE': open(join(seq_ann_root, 'protein_binding_site.txt'), 'w'),
        'CONSERVED_SITE': open(join(seq_ann_root, 'protein_conserved_site.txt'), 'w'),
        'DOMAIN': open(join(seq_ann_root, 'protein_domain.txt'), 'w'),
        'FAMILY': open(join(seq_ann_root, 'protein_family.txt'), 'w'),
        'HOMOLOGOUS_SUPERFAMILY': open(join(seq_ann_root, 'protein_homologous_superfamily.txt'), 'w'),
        'PTM': open(join(seq_ann_root, 'protein_ptm.txt'), 'w'),
        'REPEAT': open(join(seq_ann_root, 'protein_repeat.txt'), 'w'),
        'GO_BP': open(join(protein_properties_root, 'protein_go_biological_process.txt'), 'w'),
        'GO_CC': open(join(protein_properties_root, 'protein_go_cellular_component.txt'), 'w'),
        'GO_MF': open(join(protein_properties_root, 'protein_go_molecular_function.txt'), 'w'),
        'RELATED_GENETIC_DISORDER': open(join(links_root, 'protein_disease.txt'), 'w')
    }

    with open(join(data_root, 'uniprot', 'uniprot_facts.txt'), 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            s = sys.intern(s)
            if p in annotation_output_files and s in protein_set:
                annotation_output_files[p].write(line)


    for fd in annotation_output_files.values():
        fd.close()

# TODO: add reactome
def get_all_protein_go(protein_set):
    """
    Get the set of protein go annotations to include in the kg

    Parameters
    ----------
    protein_set: set
        the set of proteins used to filter the ppis
    Returns
    -------
    list
        the list of unique protein go annotations
    """
    protein_set = {sys.intern(p) for p in protein_set}
    go_root = join(protein_properties_root, 'go')
    makedirs(go_root) if not isdir(go_root) else None
    go_types = ['GO_BP', 'GO_CC', 'GO_MF']
    pgi = set()
    with open(join(data_root, 'uniprot', 'uniprot_facts.txt'), 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            s = sys.intern(s)
            if p in go_types and s in protein_set:
                pgi.add((s, p, o, 'go'))

    with open(join(data_root, 'go', 'protein_go.txt'), 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            s = sys.intern(s)
            p = p.split('||')[-1]
            if p in go_types and s in protein_set:
                pgi.add((s, p, o, 'go'))
    triples = list(pgi)
    return triples

# TODO: 加入其他数据库，然后考虑是否要从关系角度细分GO类型
def get_all_pathway_go():
    go_types = ['GO_BP', 'GO_CC', 'GO_MF']

    pgi = set()
    # with open(join(data_root, 'reactome', 'reactome_gene_go_mapping.txt'), 'r') as fd:
    #     for line in fd:
    #         protein, map_type, goid, pathway_id, species = line.strip().split('\t')
    #         pgi.add((pathway_id, goid, 'reactome'))
    
    with open(join(data_root, 'reactome', 'reactome_go_rels.txt'), 'r') as fd:
        for line in fd:
            pathway_id, map_type, goid = line.strip().split('\t')
            pgi.add((pathway_id, goid, 'reactome'))

    with open(join(data_root, 'kegg', 'pathway_go.txt'), 'r') as fd:
        for line in fd:
            pathway_id, map_type, goid = line.strip().split('\t')
            pgi.add((pathway_id, goid, 'kegg'))
    
    triples = get_unique_inter_src(pgi, 'PATHWAY_GO')

    return triples


# TODO: 临时函数，将MIM作为 genetic_disorders
def get_all_protein_genetic_disorders(protein_set, genetic_disorder_set):
    """
    Get the set of protein genetic disorders to include in the kg

    Parameters
    ----------
    protein_set: set
        the set of proteins used to filter the ppis
    Returns
    -------
    list
        the list of unique protein genetic disorders
    """
    protein_set = {sys.intern(p) for p in protein_set}
    genetic_disorder_set = {sys.intern(d) for d in genetic_disorder_set}
    protein_genetic_disorders = set()
    
    with open(join(data_root, 'uniprot', 'uniprot_facts.txt'), 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            s = sys.intern(s)
            o = sys.intern(o)
            if p == 'RELATED_GENETIC_DISORDER' and s in protein_set and o in genetic_disorder_set:
                protein_genetic_disorders.add((s, o, 'uniprot'))

    triples = get_unique_inter_src(protein_genetic_disorders, 'PROTEIN_GENETIC_DISORDER')
    return triples

# TODO： 弃用
def get_all_protein_drug_interactions(protein_set, chem_dict):
    """
    Get the set of protein drug interactions to include in the kg

    Parameters
    ----------
    protein_set: set
        the set of proteins used to filter the protein drug interactions
    chem_set: set
        the set of drug used to filter the protein drug interactions

    Returns
    -------
    list
        the list of unique protein drug interactions
    """
    kegg_links = join(data_root, 'kegg', 'gene_drug.txt')
    uniprot_facts = join(data_root, 'uniprot', 'uniprot_facts.txt')
    db_targets = join(data_root, 'drugbank', 'db_targets.txt')
    ctd_drug_protein = join(data_root, 'ctd', 'ctd_drug_protein_interactions')
    # ctd_drug_protein = join(data_root, 'ctd', 'ctd_drug_protein_interactions.txt')
    pdis = set()
    carriers = set()
    enzymes = set()
    transporters = set()
    targets = set()
    with open(kegg_links, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            db_drugs = kegg_linker.convert_drugid_to_drugbank([o])
            protein_acc = kegg_linker.convert_geneid_to_uniprot([s])
            if len(db_drugs[0]) == 0:
                continue
            if len(protein_acc[0]) == 0:
                continue
            for d in db_drugs[0]:
                try:
                    d = chem_dict[d]
                except:
                    continue
                for p in protein_acc[0]:
                    pdis.add((p, d, 'kegg'))

    with open(uniprot_facts, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            if p == 'TARGET_OF_DRUG':
                try:
                    o = chem_dict[o]
                except:
                    continue
                pdis.add((s, o, 'uniprot'))

    with open(db_targets, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')[:3]
            try:
                o = chem_dict[o]
            except:
                continue
            pdis.add((o, s, 'drugbank'))

            if p == 'DRUG_CARRIER':
                carriers.add((o, s))
            elif p == 'DRUG_ENZYME':
                enzymes.add((o, s))
            elif p == 'DRUG_TRANSPORTER':
                transporters.add((o, s))
            else:
                targets.add((o, s))

    unique_triples = []
    function_triples = []
    for protein, drug, src in pdis:
        if protein in protein_set: # and drug in drug_set
            unique_triples.append((drug, 'DPI', protein, src))

    for protein, drug in targets:
        if protein in protein_set:
            function_triples.append((drug, 'DRUG_TARGET', protein))

    for protein, drug in carriers:
        if protein in protein_set:
            function_triples.append((drug, 'DRUG_CARRIER', protein))

    for protein, drug in enzymes:
        if protein in protein_set:
            function_triples.append((drug, 'DRUG_ENZYME', protein))

    for protein, drug in transporters:
        if protein in protein_set:
            function_triples.append((drug, 'DRUG_TRANSPORTER', protein))

    return unique_triples, function_triples

def get_all_chem_protein_interactions(protein_set, chem_set, chem_dict):
    """
    Get the set of compound protein interactions to include in the kg

    Returns
    -------
    list
        the list of unique protein drug interactions
    """
    chem_set = {sys.intern(c) for c in chem_set}
    protein_set = {sys.intern(p) for p in protein_set}

    with open(join(data_root, 'drugbank', 'db2inchikey.json'), 'r') as fd:
        db2inchikey = json.load(fd)

    cpis = set()
    carriers = set()
    enzymes = set()
    transporters = set()
    targets = set()
    chembl_cpi = join(data_root, 'chembl', 'chembl_cpi.txt')
    bindingdb_cpi = join(data_root, 'bindingdb', 'bindingdb_cpi.txt')
    kegg_links = join(data_root, 'kegg', 'gene_drug.txt')
    uniprot_facts = join(data_root, 'uniprot', 'uniprot_facts.txt')
    db_targets = join(data_root, 'drugbank', 'db_dti.txt')
    stitch_cpi = join(data_root, 'stitch', 'stitch_cpi.txt')
    ctd_drug_protein = join(data_root, 'ctd', 'ctd_drug_protein_interactions.txt')

    # TODO: 这里可以细分关系类型
    with open(ctd_drug_protein, 'r') as fd:
        for line in fd:
            s, p, o, prov = line.strip().split('\t')[:4]
            if prov == 'CURATED':
                o = sys.intern(o)
                s = sys.intern(s)
                if s in chem_set and o in protein_set:
                    s = chem_dict[s]
                    cpis.add((s, o, 'ctd'))

    with open(kegg_links, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            db_drugs = kegg_linker.convert_drugid_to_drugbank([o])
            protein_acc = kegg_linker.convert_geneid_to_uniprot([s])
            if len(db_drugs[0]) == 0:
                continue
            if len(protein_acc[0]) == 0:
                continue
            for d in db_drugs[0]:
                d = sys.intern(d)
                if d in chem_set:
                    d = chem_dict[d]
                    for p in protein_acc[0]:
                        p = sys.intern(p)
                        if p in protein_set:
                            cpis.add((d, p, 'kegg'))

    with open(uniprot_facts, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            if p == 'TARGET_OF_DRUG':
                o, s = sys.intern(o), sys.intern(s)
                try:
                    o = db2inchikey[o]
                except:
                    continue
                o = sys.intern(o)
                if o in chem_set:
                    o = chem_dict[o]
                    cpis.add((o, s, 'uniprot'))

    with open(db_targets, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')[:3]
            o, s = sys.intern(o), sys.intern(s)
            if o in protein_set and s in chem_set:
                s = chem_dict[s]
                cpis.add((s, o, 'drugbank'))

                if p == 'DRUG_CARRIER':
                    carriers.add((s, 'DRUG_CARRIER', o))
                elif p == 'DRUG_ENZYME':
                    enzymes.add((s, 'DRUG_ENZYME', o))
                elif p == 'DRUG_TRANSPORTER':
                    transporters.add((s, 'DRUG_TRANSPORTER', o))
                else:
                    targets.add((s, 'DRUG_TARGET', o))

    with open(chembl_cpi, 'r') as fd:
        for line in fd:
            c, _, p = line.strip().split('\t')
            c, p = sys.intern(c), sys.intern(p)
            if c in chem_set and p in protein_set:
                c = chem_dict[c]
                cpis.add((c, p, 'chembl'))

    with open(bindingdb_cpi, 'r') as fd:
        for line in fd:
            row = line.strip().split('\t')
            if len(row) != 3:
                continue
            c, _, p = row
            c, p = sys.intern(c), sys.intern(p)
            if c in chem_set and p in protein_set:
                c = chem_dict[c]
                cpis.add((c, p, 'bindingdb'))

    # TODO: 后面细分CPI关系类型
    with open(stitch_cpi, 'r') as fd:  
        for line in fd:
            c, _, p = line.strip().split('\t')
            c, p = sys.intern(c), sys.intern(p)
            if c in chem_set and p in protein_set:
                c = chem_dict[c]
                cpis.add((c, p, 'stitch'))

    unique_triples = get_unique_inter_src(cpis, 'CPI')

    function_triples = list(targets) + list(carriers) + list(enzymes) + list(transporters)
    # grouped = defaultdict(set)
    # for s, o, src in cpis:
    #     grouped[(s, o)].add(src)
    # cpis = [(s, 'CPI', o, ';'.join(src_list)) for (s, o), src_list in grouped.items()]
    return unique_triples, function_triples


def get_all_chem_pathway_associations(chem_set, chem_dict):
    """
    Get the set of compound pathway associations to include in the kg
    Returns
    -------
    list
        the list of unique protein pathway associations
    """
    chem_set = {sys.intern(c) for c in chem_set}

    kegg_links = join(data_root, 'kegg', 'drug_pathway.txt')
    ctd_pathways = [
        join(data_root, 'ctd', 'ctd_drug_reactome_pathway_association.txt'),
        join(data_root, 'ctd', 'ctd_drug_kegg_pathway_association.txt')
    ]
    db_pathways = join(data_root, 'drugbank', 'db_pathways.txt')
    hmdb_pathways = [join(data_root, 'hmdb', 'ckegg.txt'),
                     join(data_root, 'hmdb', 'csmpdb.txt')]
    
    foodb_pathways = join(data_root, 'foodb', 'compound_pathway.txt')
    smpdb_pathways = join(data_root, 'smpdb', 'smpdb_compound_pathway.txt')

    cpas = set()

    with open(kegg_links, 'r') as fd:
        for line in fd:
            c, _, o = line.strip().split('\t')
            c, o = sys.intern(c), sys.intern(o)
            if c in chem_set:
                c = chem_dict[c]
                cpas.add((c, o, 'kegg'))

    with open(db_pathways, 'r') as fd:
        for line in fd:
            c, _, o = line.strip().split('\t')
            c, o = sys.intern(c), sys.intern(o)
            if c in chem_set:
                c = chem_dict[c]
                cpas.add((c, o, 'drugbank'))

    for f in ctd_pathways:
        with open(f, 'r') as fd:
            for line in fd:
                c, _, o, prov = line.strip().split('\t')
                if prov == 'CURATED':
                    c, o = sys.intern(c), sys.intern(o)
                    if c in chem_set:
                        c = chem_dict[c]
                        cpas.add((c, o, 'ctd'))

    for f in hmdb_pathways:
        with open(f, 'r') as fd:
            for line in fd:
                c, _, o = line.strip().split('\t')
                c, o = sys.intern(c), sys.intern(o)
                if c in chem_set:
                    c = chem_dict[c]
                    cpas.add((c, o, 'hmdb'))

    with open(foodb_pathways, 'r') as fd:
        for line in fd:
            c, _, o = line.strip().split('\t')
            c, o = sys.intern(c), sys.intern(o)
            if c in chem_set:
                c = chem_dict[c]
                cpas.add((c, o, 'foodb'))
    
    with open(smpdb_pathways, 'r') as fd:
        for line in fd:
            c, _, o = line.strip().split('\t')
            c, o = sys.intern(c), sys.intern(o)
            if c in chem_set:
                c = chem_dict[c]
                cpas.add((c, o, 'smpdb'))
    cpas = get_unique_inter_src(cpas, 'CHEM_PATHWAY_ASSOCIATION')

    return cpas

def get_all_chem_go_associations(chem_set, chem_dict):
    """
    Get the set of compound GO to include in the kg
    """
    ctd_drug_go = join(data_root, 'ctd', 'ctd_drug_phenotype.txt')
    chem_set = set(sys.intern(x) for x in chem_set)
    cgis = set()
    # TODO: 这里可以细分关系类型
    with open(ctd_drug_go, 'r') as fd:
        for line in fd:
            s, _, o, _, prov, _ = line.strip().split('\t')
            if prov == 'CURATED':
                s = sys.intern(s)
                if s in chem_set:
                    s = chem_dict[s]
                    cgis.add((s, o, 'ctd'))

    triples = get_unique_inter_src(cgis, 'CHEM_GO_ASSOCIATION')

    return triples 

def get_all_metabo_interactions(protein_set, chem_set, chem_dict):
    chem_set = {sys.intern(c) for c in chem_set}
    protein_set = {sys.intern(p) for p in protein_set}

    pmis = set()
    metabos = join(data_root, 'hmdb', 'metabo.txt')

    with open(metabos, 'r') as fd:
        for line in fd:
            p, _, c = line.strip().split('\t')
            c, p = sys.intern(c), sys.intern(p)
            if c in chem_set and p in protein_set:
                c = chem_dict[c]
                pmis.add((p, 'HAS_METABOLITE', c, 'hmdb'))

    return list(pmis)


def get_all_protein_disease_associations(protein_set, umls_disease_map):
    """
    Get the set of protein disease associations to include in the kg

    Parameters
    ----------
    protein_set: set
        the set of proteins used to filter the protein disease associations
    disease_set: set
        the set of mesh diseases used to filter the protein disease associations

    Returns
    -------
    list
        the list of unique protein disease associations
    """
    # disease_set = {sys.intern(d) for d in disease_set}
    protein_set = {sys.intern(p) for p in protein_set}

    kegg_links = join(data_root, 'kegg', 'gene_disease.txt')
    uniprot_facts = join(data_root, 'uniprot', 'uniprot_facts.txt')
    ctd_protein_disease = join(data_root, 'ctd', 'ctd_protein_disease_association.txt')
    hpo_links = join(data_root, 'hpo', 'protein_disease.txt')
    pdis = set()

    with open(uniprot_facts, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            if p == 'RELATED_GENETIC_DISORDER':
                o = o.split(':')[1]
                o = umls_disease_map['omim'].get(sys.intern(o), None)
                if o is None: continue
                pdis.add((s, o, 'uniprot'))

    with open(kegg_links, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            mesh_diseases = kegg_linker.convert_disease_to_mesh([o])
            omim_diseases = kegg_linker.convert_disease_to_omim([o])
            protein_acc = kegg_linker.convert_geneid_to_uniprot([s])
            if len(protein_acc[0]) == 0:
                continue
            if len(mesh_diseases[0]) == 0 and len(omim_diseases[0]) == 0:
                continue

            for p in protein_acc[0]:
                for d in mesh_diseases[0]:
                    d = umls_disease_map['mesh'].get(sys.intern(d), None)
                    if d is not None:
                        pdis.add((p, d, 'kegg'))

                for d in mesh_diseases[0]:
                    d = umls_disease_map['omim'].get(sys.intern(d), None)
                    if d is not None:
                        pdis.add((p, d, 'kegg'))

    with open(ctd_protein_disease, 'r') as fd:
        for line in fd:
            s, p, o, prov = line.strip().split('\t')[:4]
            if prov == 'CURATED':
                o = umls_disease_map['mesh'].get(sys.intern(o), None)
                if o is not None:
                    pdis.add((s, o, 'ctd'))

    
    with open(hpo_links, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            d_src, d = o.split(':')
            if d_src in ['OMIM', 'MIM']:
                o = umls_disease_map['omim'].get(sys.intern(o), None)
                if o is not None:
                    pdis.add((s, o, 'hpo'))
            elif d_src in ['ORPHA', 'ORPHANET']:
                o = umls_disease_map['orpha'].get(sys.intern(o), None)
                if o is not None:
                    pdis.add((s, o, 'hpo'))


    pdis = get_unique_inter_src(pdis)

    unique_triples = []
    for protein, disease, src in pdis:
        protein = sys.intern(protein)
        if protein in protein_set:
            unique_triples.append((protein, 'PROTEIN_DISEASE_ASSOCIATION', disease, src))

    return unique_triples

def get_all_protein_phenotype_associations(protein_set):
    """
    Get the set of protein phenotype associations to include in the kg
    Returns
    -------
    list
        the list of unique protein phenotype associations
    """
    protein_set = {sys.intern(p) for p in protein_set}

    hpo_links = join(data_root, 'hpo', 'protein_phenotype.txt')
    pdis = set()

    with open(hpo_links, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            pdis.add((s, o, 'hpo'))

    pdis = get_unique_inter_src(pdis)

    unique_triples = []
    for protein, phenotype, src in pdis:
        protein = sys.intern(protein)
        if protein in protein_set:
            unique_triples.append((protein, 'PROTEIN_PHENOTYPE_ASSOCIATION', phenotype, src))

    return unique_triples

def get_all_phenotype_disease_associations(umls_disease_map):
    """
    Get the set of phenotype genetic_disorder associations to include in the kg
    Returns
    -------
    list
        the list of unique phenotype genetic_disorder associations
    """

    hpo_links = join(data_root, 'hpo', 'phenotype_disease.txt')
    pdis = set()

    with open(hpo_links, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            d_src, o = o.split(':')
            if d_src in ['OMIM', 'MIM']:
                o = umls_disease_map['omim'].get(sys.intern(o), None)
                if o is not None:
                    pdis.add((s, o, 'hpo'))
            elif d_src in ['ORPHA', 'ORPHANET']:
                o = umls_disease_map['orpha'].get(sys.intern(o), None)
                if o is not None:
                    pdis.add((s, o, 'hpo'))
            

    unique_triples = get_unique_inter_src(pdis, label='PHENOTYPE_GENETIC_DISORDER')

    return unique_triples

def get_all_phenotype_phenotype_rels():
    """
    Get the set of phenotype phenotype relations to include in the kg
    Returns
    -------
    list
        the list of unique phenotype phenotype relations
    """
    hpo_links = join(data_root, 'hpo', 'phenotype_rels.txt')
    pprs = set()

    with open(hpo_links, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            pprs.add((s, o, 'hpo'))

    pprs = get_unique_inter_src(pprs, 'PHENOTYPE_IS_A_PHENOTYPE')
    return pprs

def get_all_protein_pathway_associations(protein_set):
    """
    Get the set of protein pathway associations to include in the kg

    Parameters
    ----------
    protein_set: set
        the set of proteins used to filter the protein pathway associations

    Returns
    -------
    list
        the list of unique protein pathway associations
    """
    protein_set = {sys.intern(p) for p in protein_set}
    kegg_links = join(data_root, 'kegg', 'gene_pathway.txt')
    uniprot_facts = join(data_root, 'uniprot', 'uniprot_facts.txt')
    ctd_pathways = [
        join(data_root, 'ctd', 'ctd_protein_reactome_pathway_association.txt'),
        join(data_root, 'ctd', 'ctd_protein_kegg_pathway_association.txt')
    ]
    hmdb_pathways = [join(data_root, 'hmdb', 'pkegg.txt'),
                     join(data_root, 'hmdb', 'psmpdb.txt')]
    db_pathways = join(data_root, 'drugbank', 'db_pathways.txt')
    reactome_pathways = join(data_root, 'reactome', 'reactome_protein_pathway_rels.txt')

    ppis = set()
    with open(uniprot_facts, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            if p == 'RELATED_PATHWAY':
                ppis.add((s, o, 'uniprot'))

    with open(kegg_links, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')

            protein_acc = kegg_linker.convert_geneid_to_uniprot([s])
            if len(protein_acc[0]) == 0:
                protein_acc = [[s]]
            for acc in protein_acc[0]:
                ppis.add((acc, o, 'kegg'))

    with open(db_pathways, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            if p == 'PATHWAY_ENZYME':
                ppis.add((o, s, 'drugbank'))
    
    with open(reactome_pathways, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            ppis.add((s, o, 'reactome'))

    for f in ctd_pathways:
        with open(f, 'r') as fd:
            for line in fd:
                s, p, o, prov = line.strip().split('\t')
                if prov == 'CURATED':
                    ppis.add((s, o, 'ctd'))

    for f in hmdb_pathways:
        with open(f, 'r') as fd:
            for line in fd:
                s, p, o = line.strip().split('\t')
                ppis.add((s, o, 'hmdb'))

    ppis = get_unique_inter_src(ppis)
    # grouped = defaultdict(set)
    # for s, o, src in ppis:
    #     grouped[(s, o)].add(src)
    # ppis = [(s, o, ';'.join(src_list)) for (s, o), src_list in grouped.items()]

    unique_triples = []
    for protein, pathway, src in ppis:
        protein = sys.intern(protein)
        if protein in protein_set:
            unique_triples.append((protein, 'PROTEIN_PATHWAY_ASSOCIATION', pathway, src))

    return unique_triples


def get_all_drug_drug_interactions(chem_set, chem_dict):
    """
    Get the set of drug drug interactions to include in the kg

    Parameters
    ----------
    drug_set: set
        the set of drugs used to filter the drug drug interactions

    Returns
    -------
    list
        the list of unique drug drug interactions
    """
    chem_set = {sys.intern(c) for c in chem_set}
    db_interactions = join(data_root, 'drugbank', 'db_ddi.txt')
    stitch_cci = join(data_root, 'stitch', 'stitch_cci.txt')
    ddis = set()

    # TODO: New added stitch cci，后续区分 CCI 和 DDI
    with open(db_interactions, 'r') as fd:
        for line in fd:
            s, _, o = line.strip().split('\t')[:3]
            if s in chem_set and o in chem_set:
                s = chem_dict[s]
                o = chem_dict[o]
                if s > o:
                    ddis.add((o, s, 'drugbank'))
                else:
                    ddis.add((s, o, 'drugbank'))

   
    with open(stitch_cci, 'r') as fd:
        for line in fd:
            s, _, o = line.strip().split('\t')
            s, o = sys.intern(s), sys.intern(o)
            if s in chem_set and o in chem_set:
                s = chem_dict[s]
                o = chem_dict[o]
                if s > o:
                    ddis.add((o, s, 'stitch'))
                else:
                    ddis.add((s, o, 'stitch'))

    # unique_triples = []
    # for d1, d2 in ddis:
    #     d1, d2 = sys.intern(d1), sys.intern(d2)
    #     if d1 in chem_set and d2 in chem_set:
    #         d1, d2 = chem_dict[d1], chem_dict[d2]
    #         unique_triples.append((d1, 'DDI', d2))


    # for c1, c2 in ddis: 
    #     c1, c2 = sys.intern(c1), sys.intern(c2)
    #     if c1 in chem_set and c2 in chem_set:
    #         c1, c2 = chem_dict[c1], chem_dict[c2]
    #         unique_triples.append((c1, 'CCI', c2))

    unique_triples = get_unique_inter_src(ddis, 'CCI')
    return unique_triples


def get_all_drug_pathway_associations(drug_set):
    """
    Get the set of drug pathway associations to include in the kg

    Parameters
    ----------
    drug_set: set
        the set of drugs used to filter the drug pathway associations

    Returns
    -------
    list
        the list of unique drug pathway associations
    """
    kegg_links = join(data_root, 'kegg', 'drug_pathway.txt')

    ctd_pathways = [
        join(data_root, 'ctd', 'ctd_drug_reactome_pathway_association.txt'),
        join(data_root, 'ctd', 'ctd_drug_kegg_pathway_association.txt')
    ]
    db_pathways = join(data_root, 'drugbank', 'db_pathways.txt')
    dpas = set()

    with open(kegg_links, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            db_drugs = kegg_linker.convert_drugid_to_drugbank([s])
            if len(db_drugs[0]) == 0:
                continue

            for d in db_drugs[0]:
                dpas.add((d, o))

    with open(db_pathways, 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            if p == 'DRUG_PATHWAY':
                dpas.add((s, o))

    for f in ctd_pathways:
        with open(f, 'r') as fd:
            for line in fd:
                s, p, o, prov = line.strip().split('\t')
                if prov == 'CURATED':
                    dpas.add((s, o))

    unique_triples = []
    for drug, pathway in dpas:
        if drug in drug_set:
            unique_triples.append((drug, 'DRUG_PATHWAY_ASSOCIATION', pathway))

    return unique_triples

def get_all_chem_disease_associations(chem_set, chem_dict, umls_disease_map):
    """
    Get the set of compound pathway associations to include in the kg

    Parameters
    ----------
    drug_set: set
        the set of compounds used to filter the drug pathway associations

    Returns
    -------
    list
        the list of unique compound pathway associations
    """
    chem_set = {sys.intern(c) for c in chem_set}
    kegg_links = join(data_root, 'kegg', 'disease_drug.txt')
    ctd_links = join(data_root, 'ctd', 'ctd_drug_diseases.txt')
    hmdb_links = join(data_root, 'hmdb', 'cdisease.txt')
    ddis = set()

    with open(kegg_links, 'r') as fd:
        for line in fd:
            s, _, o = line.strip().split('\t')
            db_drugs = kegg_linker.convert_drugid_to_drugbank([o])
            if len(db_drugs[0]) == 0:
                continue

            mesh_diseases = kegg_linker.convert_disease_to_mesh([s])
            omim_diseases = kegg_linker.convert_disease_to_omim([s])
            if len(mesh_diseases[0]) == 0 and len(omim_diseases[0]) == 0:
                continue
            
            # if sys.intern(o) in chem_set:
            #     o = chem_dict[o]
            #     ddis.add((o, s, 'kegg'))
            for dr in db_drugs[0]:
                for ds in mesh_diseases[0]:
                    ds = umls_disease_map['mesh'].get(sys.intern(ds), None)
                    if ds is not None:
                        ddis.add((dr, ds, 'kegg'))
                for ds in omim_diseases[0]:
                    ds = umls_disease_map['omim'].get(sys.intern(ds), None)
                    if ds is not None:
                        ddis.add((dr, ds, 'kegg'))
                    

    with open(ctd_links, 'r') as fd:
        for line in fd:
            s, _, o, prov = line.strip().split('\t')[:4]
            o = umls_disease_map['mesh'].get(sys.intern(o), None)
            if o is None: continue
            if sys.intern(s) in chem_set:
                s = chem_dict[s]
                if prov == 'CURATED':
                    ddis.add((s, o, 'ctd'))
    
    with open(hmdb_links, 'r') as fd:
        for line in fd:
            s, _, o = line.strip().split('\t')
            o = o.split(":")[1]
            o = umls_disease_map['omim'].get(sys.intern(o), None)
            if o is None: continue
            if sys.intern(s) in chem_set:
                s = chem_dict[s]
                ddis.add((s, o, 'hmdb'))

    unique_triples = get_unique_inter_src(ddis, 'CHEM_DISEASE_ASSOCIATION')
    # for drug, disease in ddis:
    #     if drug in drug_set and disease in disease_set:
    #         unique_triples.append((drug, 'DRUG_DISEASE_ASSOCIATION', disease))
    return unique_triples

# TODO: fix this func
def get_all_drug_side_effects(chem_set, chem_dict):
    """
    Get the set of drug side effect associations to include in the kg

    Parameters
    ----------
    drug_set: set
        the set of drugs used to filter the drug side effect associations

    Returns
    -------
    list
        the list of unique drug side effect associations
    """
    chem_set = {sys.intern(c) for c in chem_set}
    
    with open('database/processed/drugbank/db2inchikey.json', 'r') as fd:
        db2inchikey = json.load(fd)

    sider_effects = join(data_root, 'sider', 'sider_effects.txt')
    sider_indications = join(data_root, 'sider', 'sider_indications.txt')

    side_effects = set()
    indications = set()
    with open(sider_effects, 'r') as fd:
        for line in fd:
            d, _, o = line.strip().split('\t')
            db_ids = sider_linker.convert_drugs_to_drugbank([d])
            if len(db_ids[0]) == 0:
                continue

            for db_id in db_ids[0]:
                inchikey = db2inchikey.get(db_id, None)
                if inchikey is None:
                    continue
                inchikey = sys.intern(inchikey)
                if inchikey in chem_set:
                    s = chem_dict[inchikey]
                    side_effects.add((s, o, 'sider'))

    with open(sider_indications, 'r') as fd:
        for line in fd:
            d, _, o = line.strip().split('\t')
            db_ids = sider_linker.convert_drugs_to_drugbank([d])
            if len(db_ids[0]) == 0:
                continue

            for db_id in db_ids[0]:
                inchikey = db2inchikey.get(db_id, None)
                if inchikey is None:
                    continue
                inchikey = sys.intern(inchikey)
                if inchikey in chem_set:
                    s = chem_dict[inchikey]
                    indications.add((s, o, 'sider'))

    

    unique_triples = get_unique_inter_src(side_effects, 'DRUG_SIDEEFFECT_ASSOCIATION')
    unique_indications = get_unique_inter_src(indications, 'DRUG_INDICATION_ASSOCIATION')

    return unique_triples, unique_indications


def get_all_drug_atc_codes(chem_set, chem_dict):
    """
    Get the set of drug atc codes to include in the kg

    Parameters
    ----------
    drug_set: set
        the set of drugs used to filter the drug drug interactions

    Returns
    -------
    list
        the list of unique drug atc codes
    """
    chem_set = {sys.intern(c) for c in chem_set}
    db_interactions = join(data_root, 'drugbank', 'db_atc.txt')
    datc = set()

    with open(db_interactions, 'r') as fd:
        for line in fd:
            s, _, o = line.strip().split('\t')[:3]
            s = sys.intern(s)
            if s in chem_set:
                s = chem_dict[s]
                datc.add((s, o, 'drugbank'))

    unique_triples = get_unique_inter_src(datc, 'DRUG_ATC_CODE')

    return unique_triples


def get_all_disease_pathway_associations(umls_disease_map):
    """
    Get the set of disease pathway associations to include in the kg

    Parameters
    ----------

    Returns
    -------
    list
        the list of unique disease pathway associations
    """
    kegg_links = join(data_root, 'kegg', 'disease_pathway.txt')

    ctd_links = [
        join(data_root, 'ctd', 'ctd_disease_reactome_pathway_association.txt'),
        join(data_root, 'ctd', 'ctd_disease_kegg_pathway_association.txt'),
    ]

    dpis = set()

    with open(kegg_links, 'r') as fd:
        for line in fd:
            s, _, o = line.strip().split('\t')

            mesh_diseases = kegg_linker.convert_disease_to_mesh([s])
            omim_diseases = kegg_linker.convert_disease_to_omim([o])
            if len(mesh_diseases[0]) == 0 and len(omim_diseases[0]) == 0:
                continue
            # for ds in mesh_diseases[0]:
            #     dpis.add((ds, o, 'kegg'))
            for ds in mesh_diseases[0]:
                    ds = umls_disease_map['mesh'].get(sys.intern(ds), None)
                    if ds is not None:
                        dpis.add((ds, o, 'kegg'))
            for ds in mesh_diseases[0]:
                    ds = umls_disease_map['mesh'].get(sys.intern(ds), None)
                    if ds is not None:
                        dpis.add((ds, o, 'kegg'))

    # TODO: 这里似乎没有CURATED，全是INFERRED，考虑到底要不要加入
    for f in ctd_links:
        with open(f, 'r') as fd:
            for line in fd:
                s, _, o, prov = line.strip().split('\t')
                # if prov == 'CURATED':
                s = umls_disease_map['mesh'].get(sys.intern(s), None)
                if s is not None:
                    dpis.add((s, o, 'ctd'))
    
    unique_triples = get_unique_inter_src(dpis, label = 'DISEASE_PATHWAY_ASSOCIATION')

    # unique_triples = []
    # for disease, pathway, src in dpis:
    #     # if disease in disease_set:
    #     unique_triples.append((disease, 'DISEASE_PATHWAY_ASSOCIATION', pathway, src))
    return unique_triples

# TODO: 这里似乎没有CURATED，全是INFERRED，考虑到底要不要加入
def get_all_disease_go_associations(umls_disease_map):
    """
    Get the set of compound GO to include in the kg
    """
    ctd_disease_gos = [
        join(data_root, 'ctd', 'ctd_disease_biological_process.txt'),
        join(data_root, 'ctd', 'ctd_disease_cellular_component.txt'),
        join(data_root, 'ctd', 'ctd_disease_molecular_function.txt')
    ]

    dgis = set()
    for file in ctd_disease_gos:
        with open(file, 'r') as fd:
            for line in fd:
                s, _, o, prov = line.strip().split('\t')[:4]
                s = umls_disease_map['mesh'].get(sys.intern(s), None)
                if s is not None:
                    dgis.add((s, o, 'ctd'))

    triples = get_unique_inter_src(dgis, 'DISEASE_GO_ASSOCIATION')

    return triples 


def get_all_pathway_pathway_associations():
    """
    Get the set of pathway pathway associations to include in the kg

    Parameters
    ----------

    Returns
    -------
    list
        the list of unique pathway pathway associations
    """

    files = {
        'smpdb': join(data_root, 'smpdb', 'smpdb_pathway_pathway.txt'),
        'reactome': join(data_root, 'reactome', 'reactome_pathway_rels.txt'),
        'compath': join(data_root, 'compath', 'pathway_mappings.txt'),
        'kegg': join(data_root, 'kegg', 'pathway_pathway.txt'),
    }
    # TODO: 关系类型可以细分分层但是感觉没必要
    ppas = set()
    for src, fp in files.items():
        with open(fp, 'r') as fd:
            for line in fd:
                s, r, o = line.strip().split('\t')
                if s > o:
                    ppas.add((o, s, src))
                else:
                    ppas.add((s, o, src))

    
    unique_triples = get_unique_inter_src(ppas, 'PATHWAY_PATHWAY_ASSOCIATION')

    return unique_triples


def get_disease_tree(umls_disease_map):
    """
    Get the set of disease disease associations to include in the kg

    Parameters
    ----------

    Returns
    -------
    list
        the list of unique disease disease associations
    """
    disease_tree = join(data_root, 'mesh', 'mesh_disease_tree.txt')
    distree = set()

    with open(disease_tree, 'r') as fd:
        for line in fd:
            s, _, o = line.strip().split('\t')      
            s = umls_disease_map['mesh'].get(sys.intern(s), None)
            if s is not None:
                distree.add((s, o))     
            # distree.add((s, o))

    unique_triples = []
    for disease, category in distree:
        unique_triples.append((disease, 'DISEASE_SUPERGRP', category))
    return unique_triples


# def get_pathway_rels():
#     pathway_rels = set()
#     with open(join(data_root, 'reactome', 'reactome_pathway_rels.txt'), 'r') as fd:
#         for line in fd:
#             s, _, o = line.strip().split('\t')
#             pathway_rels.add((s, o))

#     unique_triples = []
#     for pathway, parent in pathway_rels:
#         unique_triples.append((pathway, 'HAS_PARENT_PATHWAY', parent))

#     return unique_triples

def get_all_go_go_rels():
    triples = set()
    with open(join(data_root, 'go', 'go_go.txt'), 'r') as fd:
        for line in fd:
            s, r, o = line.strip().split('\t')
            # capitalize r
            r = 'GO_' + r.upper() + '_GO'
            if s.startswith('GO:') and o.startswith('GO:'):
                triples.add((s, r, o, 'go'))
    return list(triples)

def get_complex_pathway_rels():
    top_level_pathways = set()
    complex_pathways = set()
    with open(join(data_root, 'reactome', 'reactome_complex_pathway_rels.txt'), 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')
            if p == 'COMPLEX_PATHWAY':
                complex_pathways.add((s, o))
            else:
                top_level_pathways.add((s, o))

    unique_triples = []
    for _complex, pathway in complex_pathways:
        unique_triples.append((_complex, 'COMPLEX_IN_PATHWAY', pathway, 'reactome'))

    tl_triples = []
    for _complex, pathway in top_level_pathways:
        tl_triples.append((_complex, 'COMPLEX_TOP_LEVEL_PATHWAY', pathway, 'reactome'))
    return unique_triples, tl_triples


def get_protein_complex_rels(protein_set):
    protein_set = {sys.intern(p) for p in protein_set}
    protein_complex = set()
    with open(join(data_root, 'reactome', 'reactome_protein_complex_rels.txt'), 'r') as fd:
        for line in fd:
            s, _, o = line.strip().split('\t')[:3]
            protein_complex.add((s, o))

    unique_triples = []
    for protein, _complex in protein_complex:
        if sys.intern(protein) in protein_set:
            unique_triples.append((protein, 'MEMBER_OF_COMPLEX', _complex, 'reactome'))

    return unique_triples


def get_all_protein_expressions(protein_set):
    """
    Get the set of protein tissue expressions to include in the kg

    Parameters
    ----------
    protein_set: set
        the set of proteins used to filter the protein tissue expressions
    Returns
    -------
    list
        the list of unique protein tissue expressions
    """
    protein_set = {sys.intern(p) for p in protein_set}
    hpa_tissue_expression = join(data_root, 'hpa', 'hpa_tissues_exp.txt')
    high_expression = set()
    medium_expression = set()
    low_expression = set()
    expression = set()
    tissue_structure = set()
    with open(hpa_tissue_expression, 'r') as fd:
        for line in fd:
            p, t, level = line.strip().split('\t')
            parts = t.split('__')
            if len(parts) == 2:
                tissue = parts[0]
                cell = parts[1]
                tissue_structure.add((cell, tissue))
            if level in ['low', 'medium', 'high']:
                expression.add((p, t))
            if level == 'low':
                low_expression.add((p, t))
            elif level == 'medium':
                medium_expression.add((p, t))
            elif level == 'high':
                high_expression.add((p, t))

    unique_triples = []
    # for protein, tissue in expression:
    #     if protein in protein_set:
    #         unique_triples.append((protein, 'Protein_Expression', tissue))
    hpa_other = join(other_root, 'hpa')
    makedirs(hpa_other) if not isdir(hpa_other) else None
    level_output = open(join(hpa_other, 'protein_expression_level.txt'), 'w')
    for protein, tissue in low_expression:
        if sys.intern(protein) in protein_set:
            unique_triples.append((protein, 'PROTEIN_EXPRESSED_IN', tissue))
            level_output.write(f'{protein}\tPROTEIN_EXPRESSED_IN\t{tissue}\tLOW\thpa\n')
    for protein, tissue in medium_expression:
        if sys.intern(protein) in protein_set:
            unique_triples.append((protein, 'PROTEIN_EXPRESSED_IN', tissue))
            level_output.write(f'{protein}\tPROTEIN_EXPRESSED_IN\t{tissue}\tMEDIUM\thpa\n')
    for protein, tissue in high_expression:
        if sys.intern(protein) in protein_set:
            unique_triples.append((protein, 'PROTEIN_EXPRESSED_IN', tissue))
            level_output.write(f'{protein}\tPROTEIN_EXPRESSED_IN\t{tissue}\tHIGH\thpa\n')

    structure_triples = []
    for cell, tissue in tissue_structure:
        structure_triples.append((cell, 'PART_OF_TISSUE', tissue))

    return unique_triples, structure_triples

# def get_all_disease_genetic_disorder(disease_set, genetic_disorder_set):
#     dgd = set()

#     for disease in disease_set:
#         mapped_mim = mesh_linker.convert_disease_to_omim([disease])
#         for mim in mapped_mim[0]:
#             title = f'MIM:{mim}'
#             if title in genetic_disorder_set:
#                 dgd.add((disease, title))

#     unique_triples = []
#     for disease, disorder in dgd:
#         unique_triples.append((disease, 'DISEASE_GENETIC_DISORDER', disorder))

#     return unique_triples


def write_protein_cellline_expressions(protein_set):
    pcl = set()
    hpa_other = join(other_root, 'hpa')
    makedirs(hpa_other) if not isdir(hpa_other) else None
    with open(join(data_root, 'hpa', 'hpa_cellines_exp.txt'), 'r') as fd:
        for line in fd:
            pro, _, cl_tissue, exp = line.strip().split('\t')
            cl = cl_tissue.split('#')[1]
            if cl == 'NA':
                continue
            exp = exp.split(':')[1]
            pcl.add((pro, cl, exp))

    with gzip.open(join(hpa_other, 'protein_cellline_expression.txt.gz'), 'wt') as fd:
        for protein, cellline, expression in pcl:
            if protein in protein_set:
                fd.write(f'{protein}\t{cellline}\t{expression}\n')


def write_triples(triples, output_fp):
    if len(triples[0]) == 3:
        with open(output_fp, 'w') as output:
            for s, p, o in triples:
                output.write(f'{s}\t{p}\t{o}\n')
    elif len(triples[0]) == 4:
        with open(output_fp, 'w') as output:
            for s, p, o, src in triples:
                output.write(f'{s}\t{p}\t{o}\t{src}\n')


def filter_ctd_drug_protein(protein_set):
    ctd_drug_protein = join(data_root, 'ctd', 'ctd_drug_protein_interactions.txt')
    ctd_drug_protein_filtered = join(links_root, 'ctd_drug_protein_interactions_SWISSPORT.txt' )
    with open(ctd_drug_protein_filtered, 'w') as output_fd:
        with open(ctd_drug_protein, 'r') as fd:
            for line in fd:
                s, p, o, prov = line.strip().split('\t')[:4]
                if o in protein_set:
                    output_fd.write(line)


def write_uniprot_metadata():
    uniprot_meta_dp = join(meta_root, 'protein')
    makedirs(uniprot_meta_dp) if not isdir(uniprot_meta_dp) else None

    meta_output_files = {
        'NAME': open(join(uniprot_meta_dp, 'uniprot_name.txt'), 'w'),
        'SHORT_NAME': open(join(uniprot_meta_dp, 'uniprot_shortname.txt'), 'w'),
        'FULL_NAME': open(join(uniprot_meta_dp, 'uniprot_fullname.txt'), 'w'),
        'ORGANISM_CLASS': open(join(uniprot_meta_dp, 'uniprot_organism_class.txt'), 'w'),
        'OTHER_ID': open(join(uniprot_meta_dp, 'uniprot_other_ids.txt'), 'w'),
        'RELATED_KEYWORD': open(join(uniprot_meta_dp, 'uniprot_related_keywords.txt'), 'w'),
        'RELATED_PUBMED_ID': open(join(uniprot_meta_dp, 'uniprot_related_pubmed_ids.txt'), 'w'),
        'SPECIES': open(join(uniprot_meta_dp, 'uniprot_species.txt'), 'w'),
    }
    with open(join(data_root, 'uniprot', 'uniprot_metadata.txt'), 'r') as fd:
        for line in fd:
            s, p, o = line.strip().split('\t')

            # Fail if metadata type is not in the map
            if p not in meta_output_files:
                raise Exception(f'Predicate not recognized {p}')

            meta_output_files[p].write(line)
    for fd in meta_output_files.values():
        fd.close()


def write_drugbank_metadata():
    drugbank_meta_dp = join(meta_root, 'drug')
    makedirs(drugbank_meta_dp) if not isdir(drugbank_meta_dp) else None

    meta_output_files = {
        'NAME': open(join(drugbank_meta_dp, 'drugbank_name.txt'), 'w'),
        'SYNONYM': open(join(drugbank_meta_dp, 'drugbank_synonym.txt'), 'w'),
        'TYPE': open(join(drugbank_meta_dp, 'drugbank_type.txt'), 'w'),
        'PRODUCT': open(join(drugbank_meta_dp, 'drugbank_product.txt'), 'w'),
        'PUBMED_ARTICLE': open(join(drugbank_meta_dp, 'drugbank_related_pubmed_ids.txt'), 'w'),
        'DIRECT_PARENT': open(join(drugbank_meta_dp, 'drugbank_direct_parent.txt'), 'w'),
        'KINGDOM': open(join(drugbank_meta_dp, 'drugbank_kingdom.txt'), 'w'),
        'SUPERCLASS': open(join(drugbank_meta_dp, 'drugbank_superclass.txt'), 'w'),
        'CLASS': open(join(drugbank_meta_dp, 'drugbank_class.txt'), 'w'),
        'SUBCLASS': open(join(drugbank_meta_dp, 'drugbank_subclass.txt'), 'w'),
        'ALTERNATIVE_PARENT': open(join(drugbank_meta_dp, 'drugbank_alternative_parent.txt'), 'w'),
        'SUBSTITUENT': open(join(drugbank_meta_dp, 'drugbank_substituent.txt'), 'w'),
        'PRODUCT_STAGE': open(join(drugbank_meta_dp, 'drugbank_product_stage.txt'), 'w')
    }
    files = [
        join(data_root, 'drugbank', 'db_meta.txt'),
        join(data_root, 'drugbank', 'db_product_stage.txt'),
        join(data_root, 'drugbank', 'db_classification.txt'),
    ]
    for f in files:
        with open(f, 'r') as fd:
            for line in fd:
                s, p, o = line.strip().split('\t')

                # Fail if metadata type is not in the map
                if p not in meta_output_files:
                    raise Exception(f'Predicate not recognized {p}')

                meta_output_files[p].write(line)

    for fd in meta_output_files.values():
        fd.close()

    with open(join(data_root, 'drugbank', 'db_pathways.txt'), 'r') as fd:
        with open(join(pathway_properties_root, 'pathway_category.txt'), 'w') as output:
            for line in fd:
                s, p, o = line.strip().split('\t')
                if p == 'PATHWAY_CATEGORY':
                    output.write(line)


# def write_pathway_go_annotations():
#     pred_file_map = {
#         'GO_BP': open(join(pathway_properties_root, 'pathway_go_biological_processes.txt'), 'w'),
#         'GO_CC': open(join(pathway_properties_root, 'pathway_go_cellular_components.txt'), 'w'),
#         'GO_MF': open(join(pathway_properties_root, 'pathway_go_molecular_functions.txt'), 'w')
#     }

#     with open(join(data_root, 'reactome', 'reactome_go_rels.txt'), 'r') as fd:
#         for line in fd:
#             protein, map_type, goid, pathway_id, species = line.strip().split('\t')
#             pred_file_map[map_type].write(f'{pathway_id}\tPATHWAY_{map_type}\t{goid}\n')
    
#     for f in pred_file_map.values():
#         f.close()


def write_mesh_metadata():
    mesh_meta_dp = join(meta_root, 'disease')
    makedirs(mesh_meta_dp) if not isdir(mesh_meta_dp) else None

    meta_output_files = {
        'NAME': open(join(mesh_meta_dp, 'mesh_name.txt'), 'w'),
        'TYPE': open(join(mesh_meta_dp, 'mesh_type.txt'), 'w')
    }
    files = [
        join(data_root, 'mesh', 'mesh_disease_meta.txt'),
        join(data_root, 'mesh', 'mesh_scr_disease_meta.txt')
    ]
    for fp in files:
        with open(fp, 'r') as fd:
            for line in fd:
                s, p, o = line.strip().split('\t')

                # Fail if metadata type is not in the map
                if p not in meta_output_files:
                    raise Exception(f'Predicate not recognized {p}')

                meta_output_files[p].write(line)

    for fd in meta_output_files.values():
        fd.close()

def write_hpo_metadata():
    hpo_meta_dp = join(meta_root, 'phenotype')
    makedirs(hpo_meta_dp) if not isdir(hpo_meta_dp) else None

    meta_output_files = {
        'NAME': open(join(hpo_meta_dp, 'hpo_name.txt'), 'w'),
        'ASPECT': open(join(hpo_meta_dp, 'hpo_aspect.txt'), 'w'),
        'DEFINITION': open(join(hpo_meta_dp, 'hpo_definition.txt'), 'w'),
    }
    
    with open(join(data_root, 'hpo', 'hpo_metadata.txt'), 'r') as fd:
        for line in fd:
            row = line.strip().split('\t')
            if len(row) != 3: continue
            s, p, o = row
            if p in meta_output_files:
                meta_output_files[p].write(line)

    for fd in meta_output_files.values():
        fd.close()

def copy_folder(src_folder, dst_folder, included_files=[]):
    for f in listdir(src_folder):
        if f in included_files:
            copy(join(src_folder, f), dst_folder)


def compress_folder(folder):
    for root, dirs, files in walk(folder):
        for fp in files:
            if fp.endswith('.gz'):
                continue
            src_fp = join(root, fp)
            dst_fp = join(root, fp+'.gz')
            with open(src_fp, 'rb') as f_in, gzip.open(dst_fp, 'wb') as f_out:
                f_out.writelines(f_in)
            remove(src_fp)


def generate_core_links(file_ent_type):
    output = open(join(core_root, 'unibiomap.links.tsv'), 'w')
    for f in listdir(links_root):
        if f.endswith('.txt') and f != 'README.txt':
            h_type, t_type = file_ent_type[f]
            with open(join(links_root, f), 'r') as fd:
                for line in fd:
                    line = f"{h_type}\t{t_type}\t" + line.strip() + '\n'
                    output.write(line)
    output.close()

# TODO: SMPDB中存在一些（极少，3对儿）名称和描述完全一致但是 ID 不一致的通路，实际检索其通路图可能不同，暂时不知道如何处理
def find_duplicate_map(d):
    reverse_map = defaultdict(list)
    for key, value in d.items():
        reverse_map[value].append(key)
    replace_map_list = [
        {k: x[0] for k in x[1:]} for x in reverse_map.values() if len(x) > 1
    ]
    replace_map = {
        k: v for x in replace_map_list for k, v in x.items()
    }
    return replace_map

def filter_descs(data_root, protein_set):  
    """
    TODO
    """
    links_df = pd.read_csv(join(core_root, 'unibiomap.links.tsv'), sep='\t', header=None, names=['htype', 'ttype', 'h', 'rel', 't', 'src'])
    # disease
    dh = links_df[links_df['htype'] == 'disease']['h'].tolist()
    dt = links_df[links_df['ttype'] == 'disease']['t'].tolist()
    disease_set = set(dh + dt)
    with open(join(data_root, 'umls', 'umls_desc.json'), 'r') as fd:
        umls_desc = json.load(fd)
    disease_desc = {
        k: umls_desc[k.split(':')[1]] for k in disease_set
    }
    with open(join(core_root, 'disease_desc.json'), 'w') as fd:
        json.dump(disease_desc, fd, indent=4)
    
    # go
    goh = links_df[links_df['htype'] == 'go']['h'].tolist()
    got = links_df[links_df['ttype'] == 'go']['t'].tolist()
    go_set = set(goh + got)
    with open(join(data_root, 'go', 'go_desc.json'), 'r') as fd:
        go_desc_preprocessed = json.load(fd)

    go_desc = {}
    go_alt_id_map = {}

    n_new_gos = 0
    for go_id in go_set:
        try:
            go_desc[go_id] = go_desc_preprocessed[go_id]
        except:
            url = f"https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/{go_id}"
            retries = 3
            for attempt in range(retries):
                try:
                    resp = requests.get(url)
                    if resp.ok:
                        break
                except requests.exceptions.RequestException as e:
                    if attempt < retries - 1:
                        time.sleep(2 ** attempt)  # Exponential backoff
                    else:
                        raise e
            data = resp.json()
            term = data["results"][0]
            term_go_id = term['id']
            if term_go_id != go_id:
                go_alt_id_map[go_id] = term_go_id # check
            go_desc[term_go_id] = {
                "name": term["name"],
                "definition": term["definition"]["text"],
                "type": term["aspect"]
            }
            # else:
            #     go_desc[go_id] = {
            #         "name": term["name"],
            #         "definition": term["definition"]["text"],
            #         "type": term["aspect"]
            #     }
    
    links_df['h'] = links_df['h'].replace(go_alt_id_map)
    links_df['t'] = links_df['t'].replace(go_alt_id_map)

    with open(join(core_root, 'go_desc.json'), 'w') as fd:
        json.dump(go_desc, fd, indent=4)

    # pathway
    ph = links_df[links_df['htype'] == 'pathway']['h'].tolist()
    pt = links_df[links_df['ttype'] == 'pathway']['t'].tolist()
    pathway_set = set(sys.intern(d) for d in ph + pt)
    with open(join(data_root, 'kegg', 'kegg_desc.json'), 'r') as fd:
        kegg_desc = json.load(fd)
    with open(join(data_root, 'smpdb',  'smpdb_desc.json'), 'r') as fd:
        smpdb_desc = json.load(fd)
    with open(join(data_root, 'reactome',  'reactome_desc.json'), 'r') as fd:
        reactome_desc = json.load(fd)

    # id_changer = {} # 给每个source的pathway加上前缀
    pathway_desc = {}
    remove_pathways = set()
    # n_new_pws = 0
    kegg_desc = {sys.intern(k): v for k, v in kegg_desc.items()}
    smpdb_desc = {sys.intern(k): v for k, v in smpdb_desc.items()}
    reactome_desc = {sys.intern(k): v for k, v in reactome_desc.items()}

    for pathway in pathway_set:
        # smpdb
        if pathway.startswith('SMP'):
            desc = smpdb_desc.get(sys.intern(pathway), None)
            if desc is None:
                remove_pathways.add(pathway)
            else:
                pathway_desc[pathway] = desc
                # id_changer[pathway] = 'smpdb:' + pathway
        # reactome
        elif pathway.startswith('R-'):
            desc = reactome_desc.get(sys.intern(pathway), None)
            if desc is not None:
                pathway_desc[pathway] = desc
                # id_changer[pathway] = 'reactome:' + pathway
            else:
                remove_pathways.add(pathway) #多数都是404，暂且移除
                
                # print(f"searching {pathway}")
                # url = f'https://reactome.org/ContentService/data/discover/{pathway}'
                # retries = 3
                # for attempt in range(retries):
                #     try:
                #         resp = requests.get(url)
                #         if resp.ok:
                #             break
                #     except requests.exceptions.RequestException as e:
                #         if attempt < retries - 1:
                #             time.sleep(2)  # Exponential backoff
                #         else:
                #             raise e
                #             # remove_pathways.append(pathway)
                # pathway_data = resp.json()
                # # 如果有code字段且为404
                # # TODO: 现在采取的是扔掉已删除通路，后续考虑通过持久ID映射表找到对应新ID
                # if pathway_data.get('code', None) == 404:
                #     remove_pathways.append(pathway)
                # else:
                #     pathway_desc['reactome' + pathway] = {
                #         "name": pathway_data['name'],
                #         "definition": pathway_data['description']
                #     }
                #     id_changer[pathway] = 'reactome:' + pathway
                # n_new_pws += 1
                # print('\r' + f"Requested {n_new_pws} new GO terms",
                #               end='', flush=True)
        # kegg
        # elif re.match(r'^[a-z]{3}[0-9]', pathway): # 三个字母加数字开头（kegg id）
        else:
            desc = kegg_desc.get(sys.intern(pathway), None)
            if desc is None:
                remove_pathways.add(pathway)
            else:
                pathway_desc[pathway] = desc
                # id_changer[pathway] = 'kegg' + pathway
            
    with open(join(core_root, 'pathway_desc.json'), 'w') as fd:
        json.dump(pathway_desc, fd, indent=4)
    # remove remove_pathways from link_df
    # print('removing pathways')
    # links_df = links_df[~links_df['h'].isin(remove_pathways)]
    # links_df = links_df[~links_df['t'].isin(remove_pathways)]
    mask = ~links_df['h'].isin(remove_pathways) & ~links_df['t'].isin(remove_pathways)
    links_df = links_df[mask]

    # protein
    with open(join(data_root, 'uniprot', 'uniprot_desc.json'), 'r') as fd:
        uniprot_desc = json.load(fd)
    protein_desc = {
        k: uniprot_desc[k] for k in protein_set
    }
    with open(join(core_root, 'protein_desc.json'), 'w') as fd:
        json.dump(protein_desc, fd, indent=4)

    
    # chem
    ch = links_df[links_df['htype'] == 'compound']['h'].tolist()
    ct = links_df[links_df['ttype'] == 'compound']['t'].tolist()
    chem_set = set(sys.intern(x) for x in ch + ct)

    chem_desc = {}
    chem_df = pd.read_csv(join(data_root, 'all_compounds.tsv.gz'), sep='\t', header=0, compression='gzip')
    for i, row in chem_df.iterrows():
        uci = 'UCI:' + str(row['UCI'])
        if sys.intern(uci) in chem_set:
            if pd.isna(row['InChI']):
                inchi = row['InChI_unichem']
            else:
                inchi = row['InChI']
            # TODO: 可以直接 inchi = row['InChI_unichem']
            
            if pd.isna(row['SMILES']):
                smi = Chem.MolToSmiles(Chem.MolFromInchi(inchi))
            else:
                smi = row['SMILES']

            if pd.isna(row['Name']):
                name = None
            else:
                name = row['Name']

            chem_desc[uci] = {
                'inchikey': row['InChIKey'],
                'smiles': smi,
                'inchi': inchi,
                'name': name
            }
            # convert smiles to ipuac
    with open(join(core_root, 'compound_desc.json'), 'w') as fd:
        json.dump(chem_desc, fd, indent=4)

    # phenotype
    ph = links_df[links_df['htype'] == 'phenotype']['h'].tolist()
    pt = links_df[links_df['ttype'] == 'phenotype']['t'].tolist()
    phenotype_set = set(sys.intern(d) for d in ph + pt)
    with open(join(data_root, 'hpo', 'hpo_desc.json'), 'r') as fd:
        hpo_desc = json.load(fd)
    phenotype_desc = {
        k: hpo_desc[k] for k in phenotype_set
    }
    with open(join(core_root, 'phenotype_desc.json'), 'w') as fd:
        json.dump(phenotype_desc, fd, indent=4)

    # re-save links
    links_df.to_csv(join(core_root, 'unibiomap.links.tsv'), sep='\t', index=False, header=False)




def generate_props(name, folder):
    output = open(join(core_root, f'unibiomap.properties.{name}.tsv'), 'w')
    for root, dirs, files in walk(folder):
        for fp in files:
            if fp.endswith('.txt') and fp != 'README.txt':
                with open(join(root, fp)) as fd:
                    for line in fd:
                        output.write(line)
    output.close()

def generate_core_props():
    generate_props('protein', protein_properties_root)
    generate_props('drug', drug_properties_root)
    generate_props('cell', cell_properties_root)
    generate_props('pathway', pathway_properties_root)
    generate_props('disease', disease_properties_root)
    generate_props('genetic_disorder', join(properties_root, 'genetic_disorders'))



def generate_meta(name, folder):
    output = open(join(core_root, f'unibiomap.metadata.{name}.tsv'), 'w')
    for root, dirs, files in walk(folder):
        for fp in files:
            if fp.endswith('.txt') and fp != 'README.txt':
                with open(join(root, fp)) as fd:
                    for line in fd:
                        output.write(line)
    output.close()


def generate_core_metadata():
    generate_meta('protein', join(meta_root, 'protein'))
    generate_meta('drug', join(meta_root, 'drug'))
    generate_meta('disease', join(meta_root, 'disease'))
    generate_meta('pathway', join(meta_root, 'pathway'))
    generate_meta('phenotype', join(meta_root, 'phenotype'))


def write_cell_names():
    cellname_set = set()
    cell_meta = join(meta_root, 'cell')
    makedirs(cell_meta) if not isdir(cell_meta) else None

    with open(join(data_root, 'cellosaurus', 'cl_map.txt')) as fd:
        for line in fd:
            name, cell_line = line.strip().split('\t')

            cellname_set.add((cell_line, name))

    with open(join(cell_meta, 'cell_names.txt'), 'w') as output:
        for cell_line, name in cellname_set:
            output.write(f'{cell_line}\tNAME\t{name}\n')


def write_pathway_names():
    pathway_names = set()
    pathway_meta = join(meta_root, 'pathway')
    reactome_pathway_names = join(data_root, 'reactome', 'reactome_pathway_names.txt')
    smpdb_pathway_names = join(data_root, 'smpdb', 'smpdb_pathway_names.txt')
    makedirs(pathway_meta) if not isdir(pathway_meta) else None

    kegg_pathways = kegg_linker.pathway_ids

    for pathway in kegg_pathways:
        names = kegg_linker.convert_pathway_to_names([pathway])
        for name in names[0]:
            pathway_names.add((pathway, name))

    with open(reactome_pathway_names, 'r', encoding="utf-8") as fd:
        for line in fd:
            pathway, _, name = line.strip().split('\t')
            pathway_names.add((pathway, name))

    # with open(smpdb_pathway_names, 'r') as fd:
    #     for line in fd:
    #         pathway, _, name = line.strip().split('\t')
    #         pathway_names.add((pathway, name))

    with open(join(pathway_meta, 'pathway_names.txt'), 'w') as output:
        for pathway, name in pathway_names:
            output.write(f'{pathway}\tNAME\t{name}\n')





