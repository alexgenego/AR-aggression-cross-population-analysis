# AR-aggression-cross-population-analysis

Reproducible analysis code and supporting materials for the manuscript:

**Androgen Receptor CAG Repeat Length and Aggressive Behavior Across Armenian, Hadza, and Datoga Men: A Cross-Population Genetic Association Analysis**

Target journal: *Annals of Human Genetics*

## Scope

This repository documents the statistical workflow used to reanalyse aggression-related phenotypes and genetic variants in an Armenian student sample and to integrate the shared **AR(CAG)n**, **DRD4 exon III VNTR**, and **HTR2A rs6311** variants across Armenian, Hadza, and Datoga male samples.

The principal confirmatory result is the association of longer **AR(CAG)n** repeat length with lower physical aggression. Analyses of DRD4, HTR2A, population-context heterogeneity, phenotypic dispersion, age sensitivity, multilocus models, and power are also documented.

## Repository structure

- `scripts/` — reproducible Python analysis scripts
- `documentation/` — coding rules, workflow notes, QC documentation, and provenance
- `supporting_information/` — publication-ready, non-identifying supporting files
- `requirements.txt` — frozen Python package versions used for the reconstructed workflow
- `CITATION.cff` — citation metadata

## Reproducibility status

The repository contains scripts `01` through `09`, covering the principal analytical blocks reported in Supporting Information S2-S9.

Most analytical blocks are deterministic reconstructions and contain internal checks against frozen numerical targets. Two stochastic components require a specific qualification:

- the historical LASSO random seed/CV split objects were not preserved;
- the historical bootstrap random state for variance/skewness intervals was not preserved.

For these components, the repository uses a new fixed seed (`2026`) and explicitly labels generated stochastic quantities as **reconstructions**. The frozen Supporting Information remains the archival record of the historical stochastic runs.

See `scripts/README.md` for the script-to-table map and execution order.

## Data availability and ethical restrictions

This repository **does not contain individual-level Hadza or Datoga source data**. The archival field files include information subject to ethical, provenance, and Indigenous-data-governance restrictions and therefore are not deposited in a public repository.

The public repository contains code, coding dictionaries, analysis documentation, and non-identifying summary outputs sufficient to document the analytical workflow. Restricted source data may be considered for controlled access through the corresponding author and relevant data custodians, subject to the original consent and governance conditions.

No participant names, identifying field labels, or linkage tables should be uploaded to this repository.

## Main analysis families

1. Armenian seven-variant single-locus analysis with HC3-robust standard errors and FDR control
2. Armenian joint seven-variant models, LASSO stability analysis, and exploratory pairwise interactions
3. Hadza-Datoga nine-variant models and genotype-by-population interaction block
4. Three-population shared-locus analysis of AR, DRD4, and HTR2A
5. Continuous AR-repeat sensitivity analysis
6. Age sensitivity and nonlinearity analyses
7. Phenotypic variance, dispersion, and skewness analyses
8. Power and confirmatory sample-size calculations

## Environment

The reconstructed workflow was checked with Python 3.13.5 and the package versions frozen in `requirements.txt`.

## Zenodo

A Zenodo release and DOI will be created from the frozen GitHub version used for journal submission. The DOI will be added here and to the manuscript/Supporting Information after release.

## Authors

A.M. Kulikov; P.A. Proshakov; L.A. Revyakina; R.A. Mkrtchyan; M.L. Butovskaya; O.E. Lazebny

Corresponding author: **O.E. Lazebny** — o.e.lazebny@idbras.ru

## License

The analysis code in this repository is released under the MIT License.  
This software license does **not** grant rights to restricted participant-level data or third-party datasets.
