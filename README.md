# Database recombination in GLASS2

A comprehensive database compilation pipeline for GPCR-ligand association data, integrating multiple public databases.

## Prerequisites

Before running the pipeline, complete the following preparation steps:

### 1. Prepare Manual Data Files

Refer to `manually_prepared/README.md` for instructions on downloading and preparing required data files that cannot be automatically retrieved. This includes:

- `KiDatabase.csv` from PDSP Ki Database
- `llm.jsonl` for text-mining data
- `BindingDB_All.zip` from BindingDB

### 2. Configure Database Credentials

Edit `src/sources.ini` and update the DrugBank credentials:

```ini
[default]
username = "your_drugbank_username_here"
password = "your_drugbank_password_here"
```

Replace the placeholder values with your actual DrugBank account credentials. Other configuration options in `sources.ini` can be modified as needed.

## Installation

Install required Python dependencies:

```bash
pip install -r requirements.txt
```

## Usage

Run the compilation pipeline:

```bash
python compile.py
python chem_cal.py
```

The pipeline will:

1. Download data from multiple public databases
2. Parse and standardize compound-protein interaction data
3. Merge and deduplicate entries across databases
4. Generate the final GLASS2 dataset

## Output Structure

The compiled data will be saved in `database/glass2/` with the following structure:

```
database/glass2/
├── glass2_cls.csv           # Classification dataset (active/inactive labels)
├── glass2_reg_act.csv       # Regression dataset (active compounds with quantitative values)
├── glass2_reg_inact.csv     # Regression dataset (inactive compounds with quantitative values)
├── glass2_full.tsv          # Complete dataset with all interactions
├── ligands.tsv              # Compound information (InChIKey, SMILES, names, etc.)
├── protein.json             # GPCR protein entries from UniProt
├── inchikey2uci.json        # InChIKey to UniChem identifier mapping
└── uci_xref.json            # Cross-references for compound identifiers
```

### Dataset Descriptions

- **Classification datasets** (`glass2_cls.csv`): Binary labels for compound-protein interactions
- **Regression datasets** (`glass2_reg_act.csv`, `glass2_reg_inact.csv`): Quantitative binding affinity values (Ki, Kd, IC50, EC50) for active and inactive compounds
- **Full dataset** (`glass2_full.tsv`): Complete interaction records with all metadata
- **Ligand information** (`ligands.tsv`): Chemical structures and identifiers
- **Protein information** (`protein.json`): GPCR target details and annotations
- **Cross-reference mappings**: Identifier mappings for data integration

## Notes

- The pipeline automatically checks MD5 hashes to avoid redundant processing of unchanged files
- Processing time varies depending on data size and system performance
- Ensure sufficient disk space for downloaded source files and processed outputs
