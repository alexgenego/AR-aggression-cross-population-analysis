# Reproducibility workflow

## Public components

The public repository contains:

- analysis code;
- genotype/phenotype coding rules;
- variant nomenclature and reference-sequence mapping;
- quality-control and reconciliation notes;
- a one-command workflow runner;
- non-identifying Supporting Information metadata.

## Restricted components

Participant-level Hadza and Datoga source files remain outside GitHub because the archival datasets contain provenance-sensitive and potentially identifying information. They must be supplied locally in the gitignored `data/` directory for full numerical reproduction.

## Validation status

Deterministic analysis blocks contain internal checks against selected frozen analytical targets from the manuscript Supporting Information.

Two historical stochastic states were not preserved:

1. LASSO stability/CV random splits;
2. bootstrap state for selected variance/skewness intervals and tests.

These components are reproduced as explicitly labeled seeded reconstructions (default seed `2026`). The frozen Supporting Information remains the archival record of the historical stochastic results.

## Frozen-release principle

The Zenodo DOI will be assigned to the repository state used for journal submission only after:

1. scripts and metadata are synchronized;
2. manuscript and Supporting Information use the same title and repository URL;
3. no restricted participant-level data are present in the public repository;
4. the release metadata and citation file are checked;
5. the final release/tag is created.

After the DOI is reserved or minted, it will be inserted into the manuscript, Supporting Information, repository README, and citation metadata.
