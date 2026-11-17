# CellNav AI | Setup and submission guide

## 1. Run the complete analysis

- Install RStudio Desktop: https://posit.co/download/rstudio-desktop/ (R is separate and should already be installed).
- Extract the ZIP to a normal folder, e.g. Documents/CellNavAI. Double-click `CellNavAI.Rproj`.
- Edit the student field in `report/README_template.md`, not just in the generated README.
- In the RStudio Console, run:

```r
getwd()
source("run_all.R")
```

The working directory must be the folder that contains `run_all.R`. An internet connection is needed on the first run. The script installs one CRAN package (`ranger`) and its dependencies. The first run can take several minutes depending on your Mac and dataset size. Do not interrupt while a model is training. The pipeline uses two CPU threads and no GPU.

## 2. What gets generated

- `data/clean_proteins.csv`: cleaned dataset with source accessions and partition IDs.
- `results/quality_control.csv`, `split_counts.csv`, `length_summary.csv`: dataset summaries.
- `results/metrics.csv`, `per_class.csv`, `test_predictions.csv`: actual classifier performance.
- `results/confusion_*.csv`, `feature_importance.csv`: error and feature analysis.
- `results/session_info.txt`, `data_provenance.txt`: reproducibility records.
- `results/model_*.rds`: local trained models (ignored by Git, reproducible from scripts).
- Six PNG figures under `figures/`.
- Root `README.md`: full report with real numbers, embedded figures and references.

No R Markdown/Quarto notebook is used in this version, so there is no notebook rendering deliverable. The assessment explicitly permits `.R` files. The root README is the report.

## 3. Understand the experiment

The model predicts one of Cytoplasm, Mitochondrion, Nucleus for proteins from a restricted single-location dataset. AAC means amino acid composition; Bigrams means normalized adjacent amino acid pair frequencies; AAC_Nterm means AAC plus frequencies in the first 30 residues. All three use random forests, and all are compared with a majority-class baseline. The official homology partitions 0–3 train the models and partition 4 is held out. Macro F1 is the primary metric. Do not describe the results as a medical diagnostic or general-purpose protein locator.

## 4. GitHub workflow (real commits only)

The assignment requires a public repository named exactly with your registration number and at least five meaningful commits. The listed deadline is 15 November 2026, already past as of 17 November 2026. Check the lecturer's policy before submitting. Never backdate commits or claim that a generated project was your independent work; follow your institution's AI assistance rules.

You can use GitHub Desktop (https://desktop.github.com/) or the macOS Terminal. Install Git if necessary (`git --version` in Terminal; macOS may prompt you to install developer tools). For a genuine history, commit each stage as you review and execute it. In RStudio's Terminal, from the extracted project folder:

```sh
git init
git branch -M main
git add .gitignore CellNavAI.Rproj START_HERE.md run_all.R README.md report/README_template.md
git commit -m "Initialize project and define research question"
git add code/01_download.R data/README.md
git commit -m "Document and import DeepLoc dataset"
git add code/02_clean.R code/03_features.R
git commit -m "Implement sequence cleaning and feature extraction"
git add code/04_eda.R
git commit -m "Add exploratory sequence visualizations"
git add code/05_models.R code/06_report.R
git commit -m "Implement localization classifiers and report generation"
```

If Git requests an identity, configure your own name and GitHub email using `git config user.name "Your Name"` and `git config user.email "your-email"`. Before the first commit, open and review the report template; update the student field. Now execute `source("run_all.R")` if not already done and review the generated outputs. Commit actual results after checking them:

```sh
git add README.md data/clean_proteins.csv figures results
git commit -m "Record measured results and final report"
git status
```

On GitHub create a new public empty repository with your registration number as its name. Do not initialize it with a README. Copy its actual HTTPS URL and use it in place of `YOUR_REPOSITORY_URL` below:

```sh
git remote add origin YOUR_REPOSITORY_URL
git push -u origin main
```

Authenticate using the browser/credential manager or GitHub Desktop if required. Confirm that the public repository shows the report with six figures and your real commits. Submit its URL through LMS only if the lecturer permits submission.

## 5. Troubleshooting

- If the package installation fails, try `install.packages("ranger", repos="https://cloud.r-project.org")` and read the first error. On macOS choose the compatible binary when available. Update R if the package requires a newer version.
- If the data download fails, use the official download link in `data/README.md` manually and save the CSV with the exact filename expected by the script.
- If an unexpected column error appears, inspect `names(read.csv("data/Swissprot_Train_Validation_dataset.csv", check.names=FALSE))`; the pipeline intentionally stops rather than silently using incorrect data.
- If a run stops, fix the error, then rerun `source("run_all.R")` from the beginning.
- Review `results/session_info.txt` for version information. Do not submit until the real outputs exist and you can explain them.

## 6. Assessment checks

- Real-world problem and research question: generated README sections 1–2.
- Real public dataset and preprocessing: `code/01_download.R`, `code/02_clean.R`.
- EDA: `code/04_eda.R` and figures 1–3.
- Predictive method: `code/03_features.R`, `code/05_models.R` and figures 4–6.
- Results, limitations and interpretation: generated README sections 5–7.
- Reproducibility and citations: generated README sections 8–9, session info and provenance.
- Public registration-number repository and actual commit history: complete on your own account.

This package is an AI-assisted project scaffold with an automatic results report; the actual R pipeline has not been executed in the authoring environment. Verify the code, data, numbers, discussion and your course's AI-use policy before submitting. Modify the report in your own words as appropriate and acknowledge assistance where required.
