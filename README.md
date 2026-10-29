# AR-aggression-cross-population-analysis

Reproducible analysis code and supporting materials for the manuscript:

**Androgen Receptor CAG Repeat Length and Aggressive Behavior Across Armenian, Hadza, and Datoga Men: A Cross-Population Genetic Association Analysis**

Target journal: *Annals of Human Genetics*

## Scope

This repository documents the statistical workflow used to reanalyse aggression-related phenotypes and genetic variants in an Armenian student sample and to integrate the shared **AR(CAG)n**, **DRD4 exon III VNTR**, and **HTR2A rs6311** variants across Armenian, Hadza, and Datoga male samples.

The principal confirmatory result is the association of longer **AR(CAG)n** repeat length with lower physical aggression. Analyses of DRD4, HTR2A, population-context heterogeneity, phenotypic dispersion, age sensitivity, and multilocus models are also documented.

## Repository structure

- `scripts/` — reproducible Python analysis scripts
- `documentation/` — coding rules, workflow notes, QC documentation, and provenance
- `supporting_information/` — publication-ready supporting tables and metadata that can be shared publicly
- `requirements.txt` — Python dependencies
- `CITATION.cff` — citation metadata for the repository

## Data availability and ethical restrictions

This repository **does not contain individual-level Hadza or Datoga source data**. The archival field files include information subject to ethical, provenance, and Indigenous-data-governance restrictions and therefore are not deposited in a public repository.

The public repository contains code, coding dictionaries, analysis documentation, and non-identifying summary outputs sufficient to document the analytical workflow. Restricted source data may be considered for controlled access through the corresponding author and relevant data custodians, subject to the original consent and governance conditions.

No participant names, identifying field labels, or linkage tables should be uploaded to this repository.

## Main analysis families

1. Armenian seven-variant single-locus analysis with HC3-robust standard errors and FDR control
2. Armenian seven-variant joint models, LASSO stability analysis, and exploratory pairwise interactions
3. Hadza–Datoga nine-variant models and genotype-by-population interaction block
4. Three-population shared-locus analysis of AR, DRD4, and HTR2A
5. Continuous AR-repeat sensitivity analysis
6. Age sensitivity analyses
7. Phenotypic variance, dispersion, and skewness analyses
8. Power and confirmatory sample-size calculations

## Reproducibility

The final release will contain scripts corresponding to the analysis families above. Each script will read locally supplied input files from paths excluded by `.gitignore` and write non-identifying summary outputs.

A Zenodo release and DOI will be created from the frozen GitHub version used for journal submission.

## Authors

A.M. Kulikov; P.A. Proshakov; L.A. Revyakina; R.A. Mkrtchyan; M.L. Butovskaya; O.E. Lazebny

Corresponding author: **O.E. Lazebny** — o.e.lazebny@idbras.ru

## License

The analysis code in this repository is released under the MIT License.  
This software license does **not** grant rights to restricted participant-level data or third-party datasets.
