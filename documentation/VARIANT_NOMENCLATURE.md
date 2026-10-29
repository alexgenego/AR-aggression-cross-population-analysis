# Variant nomenclature and analytical coding

Reference framework: **GRCh38.p14** for chromosome-level genomic coordinates, together with curated RefSeq/RefSeqGene accession.version identifiers.

dbSNP `rs` identifiers are assembly-independent. HGVS expressions are reference-sequence dependent. Archived assay-oriented allele labels are retained where necessary to preserve correspondence with the original genotyping files.

| Gene | Variant | Identifier | Reference sequence(s) | Nomenclature / note | Analytical coding |
|---|---|---|---|---|---|
| **AR** | CAG repeat-length variant | repeat-length variant | NG_009014.2; NM_000044.6 | CAG repeat count in exon 1; reported by repeat number | Men: continuous repeat count and ≥23 vs ≤22; Armenian women: prespecified mean/min/max sensitivity codings |
| **DRD4** | exon III 48-bp VNTR | repeat-length variant | NG_021241.1; NM_000797.4 | 48-bp repeat-length variant | Primary: SS vs L+ (at least one allele with ≥7 repeats); SS/SL/LL exploratory |
| **OXTR** | rs53576 | rs53576 | NM_000916.4; NC_000003.12 | archived A/G assay labels retained | A-allele dosage 0-2; genotype-category sensitivity analyses |
| **HTR2A** | rs6311 | rs6311 | NG_013011.1; NM_000621.5; NC_000013.11 | HTR2A is on the reverse strand. Archived gene-oriented A/G assay labels correspond to the same variant that is represented as C/T on the forward genomic strand | A-allele dosage in the archived assay orientation |
| **HTR2C** | rs6318 | rs6318 | NG_012082.3; NM_000868.4 | historically described as Cys23Ser; X-linked | Sex-specific C-allele dosage; men hemizygous |
| **BDNF** | rs6265 (Val66Met) | rs6265 | NG_011794.2; NM_001709.5 | Val66Met | A-allele dosage 0-2 |
| **ERCC5** | rs1047768 | rs1047768 | NG_007146.1; NM_000123.4 | synonymous variant | **T-allele dosage 0-2** in the Armenian analytical files |

## Reference sources

- AR: https://www.ncbi.nlm.nih.gov/gene/367
- DRD4: https://www.ncbi.nlm.nih.gov/gene/1815
- OXTR: https://www.ncbi.nlm.nih.gov/gene/5021
- HTR2A: https://www.ncbi.nlm.nih.gov/gene/3356
- HTR2C: https://www.ncbi.nlm.nih.gov/gene/3358
- BDNF: https://www.ncbi.nlm.nih.gov/gene/627
- ERCC5: https://www.ncbi.nlm.nih.gov/gene/2073

## Important coding distinction

Nomenclature and analytical coding are separate issues. The reference allele shown by a forward genomic sequence does not automatically determine the allele used as the dosage variable in an archived assay. The scripts therefore preserve the historical analytical coding while documenting the current reference sequence and strand context.
