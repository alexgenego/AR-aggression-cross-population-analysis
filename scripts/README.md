# Analysis scripts

The repository contains the reconstructed and validated analysis workflow used for the manuscript.

## Required local inputs

Restricted participant-level input files must be placed in a local `data/` directory. That directory is excluded from Git version control.

Expected filenames:

- `data/Armenian_full_analysis.xlsx`
- `data/Hadza_Datoga_2012_cleaned_males.xlsx`
- `data/Hadza_Datoga_adults_enriched_by_name_v2.xlsx`

The source files are **not** distributed in this public repository because they contain restricted, provenance-sensitive, or potentially identifying information.

## Scripts

1. `01_three_population_shared_loci.py` — Three-population AR/DRD4/HTR2A analysis; common effects, heterogeneity, three-locus block, continuous AR, age sensitivity, and Figure 1 source data. Reproduces S6a-S6d and S6f targets.
2. `02_armenian_single_locus_fdr.py` — Armenian single-locus models, genotype-category omnibus tests, 84-test BH-FDR family, and hierarchical gene-level FDR. Reproduces S2a-S2d.
3. `03_armenian_joint_seven_locus.py` — Sex- and trait-specific joint seven-variant models with age covariate and HC3 covariance. Reproduces S3a.
4. `04_armenian_pairwise_interactions.py` — All 252 exploratory pairwise variant-variant interaction tests. Reproduces S4.
5. `05_armenian_lasso_reconstruction.py` — Seeded reconstruction of LASSO stability selection and repeated 5-fold CV. Uses the exact frozen alpha_min and alpha_1se values, but not the lost historical random split objects. It reproduces the pattern and procedure rather than historical stochastic values bit-for-bit.
6. `06_hadza_datoga_nine_variant_block.py` — Pooled Hadza-Datoga nine-variant genetic block, individual coefficients, and genotype-by-population interaction block. Reproduces S5a-S5c.
7. `07_variance_dispersion_skewness.py` — Within- and between-population variance tests plus skewness analyses. Deterministic targets reproduce the frozen analysis; bootstrap quantities are seeded reconstructions.
8. `08_age_models.py` — Complete Armenian, Hadza, Datoga, and pooled Hadza-Datoga age models. Reproduces S8.
9. `09_power_analysis.py` — Reproduces S9 power and future-sample-size calculations from the observed HC3 signal-to-noise ratio and noncentral t distribution.

## Typical execution

From the repository root:

```bash
python scripts/01_three_population_shared_loci.py
python scripts/02_armenian_single_locus_fdr.py
python scripts/03_armenian_joint_seven_locus.py
python scripts/04_armenian_pairwise_interactions.py
python scripts/05_armenian_lasso_reconstruction.py
python scripts/06_hadza_datoga_nine_variant_block.py
python scripts/07_variance_dispersion_skewness.py
python scripts/08_age_models.py
python scripts/09_power_analysis.py
```

Each script writes non-identifying summary CSV files to `outputs/` by default and performs internal reconciliation checks against selected frozen analytical targets.

## Stochastic analyses

Two historical random states were not preserved: the LASSO stability/CV split state and the bootstrap state for variance/skewness intervals. For those components the repository uses a fixed default seed of `2026` and labels generated stochastic quantities as reconstructions. The frozen Supporting Information remains the archival record of the historical stochastic runs.
