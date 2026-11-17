# Dataset provenance

- Official source: https://services.healthtech.dtu.dk/services/DeepLoc-2.0/4-Data.php
- CSV: https://services.healthtech.dtu.dk/services/DeepLoc-2.0/data/Swissprot_Train_Validation_dataset.csv
- Paper: Thumuluri et al. (2022), https://doi.org/10.1093/nar/gkac278
- Underlying UniProt data: CC BY 4.0, https://www.uniprot.org/help/license

The official source CSV is downloaded by `code/01_download.R` and ignored by Git. The generated `clean_proteins.csv` includes accession, sequence, original partition, class and length. It is a filtered derivative of the official data. Cite both DeepLoc and UniProt. The source checksum is recorded in `results/data_provenance.txt` on execution.
