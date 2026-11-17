# CellNav AI: Protein Subcellular Localization from Sequence Features

**IT4133 — Bioinformatics and Computational Biology | ICA 3**  
**2021/ICT/78 B.H.D.N.M.Wijayasinghe**
## 1. Abstract

Protein localization is fundamental to cellular function. This project investigates whether three intracellular compartments can be predicted from amino acid sequences without structural information or large protein language models. Experimentally annotated proteins from the DeepLoc 2.0 SwissProt dataset were filtered to single-location records in the cytoplasm, mitochondrion or nucleus. Three representations were compared using random forest classifiers: whole-sequence amino acid composition, adjacent amino acid pair frequencies and whole-sequence composition augmented with N-terminal composition. The original homology-based partitions were retained to evaluate generalization. @@ABSTRACT@@ The study provides a reproducible R workflow for comparing interpretable sequence representations and demonstrates their potential and limitations for protein annotation.

## 2. Introduction

Proteins perform a wide range of functions in living cells. Nuclear proteins participate in activities such as transcription, DNA replication and repair, whereas mitochondrial proteins contribute to energy metabolism and other processes. Cytoplasmic proteins perform diverse metabolic, structural and regulatory roles. A protein delivered to an inappropriate compartment may fail to interact with its intended substrates or partners. 
Experimental approaches to localization include microscopy and biochemical fractionation. Although experimental evidence is essential for establishing biological function, performing experiments for every newly sequenced protein is not always practical. Computational prediction can help prioritize proteins for further investigation. In particular, the primary amino acid sequence offers a readily available input for predictive algorithms. The sequence may contain targeting signals and compositional patterns related to its cellular destination. Mitochondrial targeting sequences, for example, frequently occur near the amino terminus, although not all mitochondrial proteins use this route (Rapaport, 2003). Nuclear localization signals can occur at multiple positions and are not captured completely by simple amino acid counts.

Recent protein localization systems employ pretrained protein language models that produce contextual sequence representations. DeepLoc 2.0 uses such representations to address multiple possible destinations (Thumuluri et al., 2022). This project does not attempt to reproduce that system or compete directly with its reported performance. Instead, it investigates a smaller, interpretable machine learning experiment appropriate for an undergraduate bioinformatics module and ordinary computing hardware.

The research question is: **How effectively can three sequence representations predict cytoplasmic, mitochondrial and nuclear localization, and does adding information from the N-terminal region improve performance compared with whole-sequence composition alone?** The primary hypothesis is that N-terminal composition adds useful localization information. A secondary question is whether adjacent amino acid pairs offer a better representation than individual residue frequencies. Both comparisons are evaluated using macro F1 on a held-out dataset partition. The objective is to measure predictive performance rather than demonstrate a causal relationship between specific amino acids and localization.

## 3. Data

The data source is the public DeepLoc 2.0 training/validation dataset, derived from UniProtKB/Swiss-Prot release 2021_03. The dataset authors curated eukaryotic proteins with experimental localization annotations, excluding fragments and proteins of 40 residues or fewer. Their original problem includes ten localization classes and allows multiple labels per protein. They also supplied five partitions based on sequence homology, designed to reduce similarity across partitions. This is particularly important because a random train-test split can place closely related proteins in both sets and overestimate performance.

The official CSV is downloaded directly from the Technical University of Denmark service. The accession (`ACC`) identifies each protein, `Sequence` contains the amino acid string, `Partition` indicates the original split, and ten binary columns encode localization annotations. Other available metadata are not used as predictors. Only sequence-derived numerical features are passed to the classifiers. The raw source remains unchanged and its MD5 checksum is stored with the analysis provenance.

To define a single-label classification problem, the project retains proteins annotated to exactly one of the ten locations, provided that location belongs to the three selected compartments. Missing or nonbinary localization annotations are excluded. Protein sequences must use the twenty standard amino acid symbols and be longer than 40 residues. Records with invalid accessions or partition numbers are removed. Exact duplicate sequences and repeated accessions are removed in their entirety, avoiding conflicting labels or duplicates crossing the evaluation boundary. This deliberately conservative approach sacrifices some data to improve dataset integrity.

@@QC@@

@@SPLIT@@

This is a secondary analysis of publicly available molecular data, not identifiable patient records. UniProt makes copyrightable database content available under CC BY 4.0; the derivative dataset is attributed to both UniProt and the DeepLoc authors. The raw dataset can be retrieved through the download script rather than redistributed with the repository.

## 4. Methods

### 4.1 Workflow and exploratory analysis

The analysis was designed as a sequential R pipeline. The first scripts import and validate the source, remove records outside the stated inclusion criteria and produce a clean table. Partitions 0–3 form the training set; partition 4 is reserved for evaluation. The choice of partitions and feature definitions is fixed in advance. Training-only exploratory plots characterize class frequency, sequence length and N-terminal amino acid composition. Test labels are used only for the predefined evaluation, apart from checking that each class is represented and documenting split counts.

A class-frequency bar chart assesses imbalance. A boxplot displays the distribution of log-transformed protein lengths by compartment, since protein lengths vary substantially. A heatmap summarizes average residue percentages in the first 30 amino acids of each class. These visualizations help establish whether the biological groups differ in properties that might be captured by sequence features. Length is described but deliberately excluded from the model inputs, so the experiment concentrates on normalized sequence representations.

### 4.2 Feature extraction

The first representation, AAC, contains twenty amino acid frequencies. For residue type a in a protein of length L, the feature is the count of a divided by L. Each feature lies between zero and one, and the twenty frequencies sum to one. This representation is compact and insensitive to sequence length, but it does not retain the order or position of residues.

The second representation contains 400 overlapping bigram frequencies, one for every ordered pair of the twenty amino acids. For example, the sequence ACA contains the adjacent pairs AC and CA, each with frequency 0.5. Counts are divided by L minus one. This captures local sequence order, analogous to bigrams in natural language processing, but cannot represent long-distance relationships or preserve the absolute positions of motifs. The implementation uses integer indexing and tabulation rather than searching individually for every possible pair.

The third representation combines the twenty whole-sequence frequencies with twenty additional frequencies calculated over the first 30 residues. The window is chosen as a simple, biologically motivated approximation to the N-terminal region and is not optimized using evaluation data. It can capture local enrichment associated with targeting signals, but does not explicitly identify motifs or cleavage sites. All representations are deterministic and calculated independently for each protein; no feature scaling statistics or vocabulary are learned from the test set. A small unit check verifies the expected composition and bigram frequencies for the sequence ACA.

### 4.3 Classification and evaluation

The same random forest framework is used for all three feature sets, so the experiment primarily varies the representation rather than the algorithm. Random forests combine many decision trees built using resampling and randomized feature selection (Breiman, 2001). They can model nonlinear relationships and do not require standardization of these frequency inputs. The R package `ranger` implements the forests efficiently (Wright and Ziegler, 2017). Each probability forest uses 300 trees, the default feature-sampling rule for its dimensionality and a fixed seed. The predicted class has the largest estimated probability. The majority-class baseline always predicts the most common training class, providing a minimal performance reference. Hyperparameter optimization and class reweighting are deliberately omitted to keep the experiment focused.

Accuracy is the fraction of correctly classified proteins. Precision measures the proportion of predicted instances of a given class that are correct, while recall measures the proportion of actual instances recovered. F1 is the harmonic mean of precision and recall; undefined class metrics are assigned zero. Macro F1 averages the three class F1 scores equally, and balanced accuracy averages the class recalls. Macro F1 is the primary metric because overall accuracy alone can obscure poor performance on minority classes. A confusion matrix identifies which compartments are most frequently confused. The combined model is selected in advance for detailed error and permutation-importance analysis, regardless of its observed rank. No test-set score is used to tune any forest. Permutation importance is based on the fitted forest's out-of-bag observations and indicates how predictive performance changes when an individual feature is disrupted; it is not causal evidence.

## 5. Results

### 5.1 Exploratory findings

![Training class frequencies](figures/01_class_balance.png)

*Figure 1. Number of proteins in each localization class in training partitions 0–3.*

![Sequence length distributions](figures/02_sequence_length.png)

*Figure 2. Log10-transformed sequence lengths by class. Outliers are omitted from the boxplot display only, not from the analysis.*

@@LENGTH@@

![N-terminal composition](figures/03_nterminal_composition.png)

*Figure 3. Mean amino acid percentages in the first 30 residues of the training proteins. Values are percentages, not counts or significance tests.*

@@EDA@@

### 5.2 Predictive performance

@@METRICS@@

![Model performance](figures/04_model_comparison.png)

*Figure 4. Accuracy, macro F1 and balanced accuracy for the baseline and three random forests on partition 4.*

@@PERFORMANCE@@

### 5.3 Error patterns and feature importance

![Confusion matrix](figures/05_confusion_matrix.png)

*Figure 5. Counts for the predefined combined model. Rows show actual labels and columns show predictions.*

@@ERRORS@@

@@PERCLASS@@

![Permutation feature importance](figures/06_feature_importance.png)

*Figure 6. The 12 highest-ranked permutation importance features of the combined model; AAC denotes whole-sequence composition and NT denotes N-terminal composition.*

@@IMPORTANCE@@

## 6. Discussion

@@DISCUSSION@@

Whole-sequence composition represents a useful reference because it is easy to calculate and explain. Bigrams add local order but greatly increase feature dimensionality. Higher dimensionality can make learning more difficult, especially when many pairs are rare. The combined representation retains a compact input while including partial positional information. Its relationship to known targeting mechanisms makes it especially interesting, even if its observed performance is not the highest.

An important strength is the use of the dataset authors' homology-based partitions rather than a newly randomized split. This reduces one important source of optimistic performance estimates in protein classification. Exact duplicate checks provide an additional safeguard, although the original homology criterion is not independently recomputed here. The evaluation uses only one of the five partitions and is therefore sensitive to the particular proteins assigned to that partition. The reported metrics are point estimates; no formal significance test or confidence interval is used, so small differences should be interpreted cautiously. Class imbalance can affect the forests, which is why the results include class-level metrics and a majority baseline.

The study also has biological limitations. Protein targeting can depend on longer motifs, three-dimensional structure, interactions with other proteins, alternative isoforms and cell-specific conditions. A frequency-based representation does not model these mechanisms directly. The chosen N-terminal window may miss signals outside the first 30 residues and cannot distinguish different residue arrangements within the window. Feature importance is affected by correlated features and by the training distribution; a highly ranked residue should not automatically be interpreted as an experimentally verified sorting signal. In particular, the study has no experimental intervention that could establish causality.

Another limitation is the restriction to single-location proteins in three compartments. Many proteins genuinely occupy multiple locations, and the original DeepLoc system was designed to handle this complexity. A protein from an excluded compartment would still be assigned one of the three classes by the present model. Predicted probabilities are not validated as calibrated uncertainty estimates, and the system should not be used as a definitive biological annotation tool. The evaluation population consists of curated eukaryotic proteins, which may not represent all newly sequenced organisms or difficult-to-annotate proteins. Thus, the findings apply to the defined dataset and task rather than to universal cellular localization.

Future work could perform full five-fold homology-aware cross-validation, quantify the uncertainty of differences, incorporate multiple localization labels and investigate additional sequence features, such as physicochemical properties or longer k-mers. A comparison with pretrained protein embeddings could test whether more sophisticated representations improve performance under the same split and label definition. An independent experimental dataset would provide a stronger evaluation of generalization. Finally, a user-facing interface could show predicted compartments together with explicit out-of-scope warnings.

## 7. Conclusion

@@CONCLUSION@@

CellNav AI provides an accessible example of machine learning applied to protein sequences. Its contribution is a controlled comparison of inexpensive and interpretable feature representations, not a novel biological discovery or a replacement for experimental localization. 
## 8. Reproducibility

Open `CellNavAI.Rproj` in RStudio, ensure an internet connection is available for the initial dataset and package installation, and run `source("run_all.R")` from the project root. The script downloads the source if necessary, runs the numbered scripts in `code/`, creates all figures and result tables, and regenerates this report. Alternatively, use `Rscript run_all.R` from a terminal in the project directory. Run the scripts in order because later stages depend on objects created by earlier stages. The raw source is cached locally and its checksum is saved to `results/data_provenance.txt`. The cleaned dataset, test predictions, confusion matrices and model metrics are available as CSV files. The software environment and package versions are recorded in `results/session_info.txt`.

@@VERSIONS@@

## 9. References

- Breiman, L. (2001). Random forests. *Machine Learning*, 45, 5–32. https://doi.org/10.1023/A:1010933404324
- Thumuluri, V., Almagro Armenteros, J. J., Johansen, A. R., Nielsen, H., & Winther, O. (2022). DeepLoc 2.0: multi-label subcellular localization prediction using protein language models. *Nucleic Acids Research*, 50(W1), W228–W234. https://doi.org/10.1093/nar/gkac278
- DTU Health Tech. DeepLoc 2.0 training and validation dataset. https://services.healthtech.dtu.dk/services/DeepLoc-2.0/4-Data.php
- UniProt Consortium (2021). UniProt: the universal protein knowledgebase in 2021. *Nucleic Acids Research*, 49(D1), D480–D489. https://doi.org/10.1093/nar/gkaa1100
- UniProt. Data licence: https://www.uniprot.org/help/license
- Wright, M. N., & Ziegler, A. (2017). ranger: A fast implementation of random forests for high dimensional data in C++ and R. *Journal of Statistical Software*, 77(1), 1–17. https://doi.org/10.18637/jss.v077.i01
- Rapaport, D. (2003). Finding the right organelle: targeting signals in mitochondrial outer-membrane proteins. *EMBO Reports*, 4, 948–952. https://pmc.ncbi.nlm.nih.gov/articles/PMC1326395/
- R Core Team. R: A language and environment for statistical computing. https://www.R-project.org/