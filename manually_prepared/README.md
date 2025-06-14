# Manually Prepared Data

This folder contains data that needs to be manually prepared by users before running the pipeline.

## Required Files

### KiDatabase.csv

**Source**: https://pdsp.unc.edu/databases/kiDownload/

**Instructions**:
1. Visit the PDSP Ki Database download page: https://pdsp.unc.edu/databases/kiDownload/
2. Click the download button to obtain `KiDatabase.csv`
3. Place the downloaded `KiDatabase.csv` file in this folder (`manually_prepared/`)

**Description**: This file contains Ki (inhibition constant) data from the PDSP (Psychoactive Drug Screening Program) database, which is essential for the binding affinity analysis pipeline.

## File Structure

After manual preparation, this folder should contain:
```
manually_prepared/
├── readme.md (this file)
└── KiDatabase.csv (downloaded from PDSP)
```

## Notes

- Ensure the file is named exactly `KiDatabase.csv`
- The file should be placed directly in the `manually_prepared/` folder
- This data is required before running any analysis scripts