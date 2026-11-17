# Pediatric AML NMF — Integrated CpG Molecular Neighborhood Atlas

A **standalone, static GitHub Pages website** for the Script 11D selected CpG suite, incorporating the original GRCh38 Script 12C annotations and Script 13C directional clinical context. **This expands the previously deployed promoter-only atlas without touching the separate clinical network repository.**

## Explore

- **Promoter overlap** is the starting view. Switch to **gene-body overlap**, **intergenic (nearest-gene label only)** or **all genomic relationships**. The Script 12C promoter interval was 2,000 bp upstream to 500 bp downstream of a transcript TSS.
- Search an annotated gene/transcript, or jump directly to any exact CpG ID (`cg...`). Switch between multi-CpG gene summaries and one-CpG mode.
- Original NMF-group circles display the **descriptive mean beta of available unique selected CpGs** in the original primary group, not differentially methylated vs controls.
- SOURCE→TARGET arrows display selected CpG association tests for a choice of contrasts: incoming bridge I − source nonbridge N (default), incoming bridge I − target primary P, or target primary P − source nonbridge N. Arrows are red for positive effects, blue for negative, and gray/dashed for mixed-sign qualifying CpGs.
- Filter by overlapping provenance categories (06, 06b, 07a, 11C, 11D), genomic context, CpG-island context, W target-loading percentile, FDR and effect threshold, original source and H-neighborhood target.
- Select circles/arrows/table CpGs to inspect observed beta means and SDs, observed CpG-specific patient counts, all Welch/TOST p/q values, W statistics, GRCh38 coordinates and all original suite labels.
- Inspect any Script 13C clinical direction, even when the current gene has no qualifying methylation arrow. The right-hand panel shows pooled, female, and male **clinical** summary rows. These are not sex-specific beta-value results.
- Download selected CpG evidence as a CSV or the redacted 13C clinical summary as a CSV.

## Audited data included

| Source / category | Actual selected data |
| --- | ---: |
| Script 11D / 12C association rows | **79,920** |
| Unique annotated CpGs | **50,013** |
| Unique promoter CpGs | **15,048** |
| Unique gene-body CpGs | **26,012** |
| Unique intergenic CpGs (nearest gene descriptive) | **8,953** |
| CpG-selected 11D directed comparisons | **119** |
| Annotated gene/transcript labels across all contexts | **20,277** |
| Original NMF-group CpG mean beta records | **134,073** |
| Original NMF groups plotted | **23** (19 with recovered means) |
| 13C directional clinical summary records | **1,254** (418 pair records × 3 populations) |
| Eligible 13C clinical directions | **153 pooled / 74 female / 65 male** |

The file `data/build-audit.json` contains the build-verified counts, all 22 label frequencies, and public redaction information. The 11D and 12C CpG/direction keys were verified to agree exactly in order. Original-group beta estimates were cross-checked when directly reported and reconstructed from source bridge/nonbridge observations.

This website includes **all selected CpGs available in the supplied 11D/12C suite**. It does **not** claim to contain the full EPIC assay, all DMRs, all 13B clinical covariates, gene-expression results, sex-specific methylation tests, or external validation.

## Deploy or update an existing standalone repository

See **[UPDATE_EXISTING_REPOSITORY.md](UPDATE_EXISTING_REPOSITORY.md)**. Keep `index.html`, `assets/` and `data/` at the repository root. GitHub Pages should publish from the `main` branch and `/(root)` (or the actual default branch).

## Rebuild after new analyses

A copy of the aggregate-only dataset is already included; you do **not** need to run Python just to host the page. If regenerating locally:

```bash
python -m pip install pandas
python scripts/build_full_atlas.py \
  --annotation /path/to/Table_12C_B_ALL_Directional_CpG_Annotations_GRCh38.csv \
  --clinical /path/to/Table_13C_02_Clinical_Directional_CpG_Integrated_Atlas.csv \
  --11d /path/to/Table_11D_C_ALL_Selected_CpG_Direction_Associations.csv \
  --output data/atlas.js
python tests/check_public_data.py
```

`--11d` is an integrity-check input: 12C already includes all the relevant 11D columns. No raw IDAT files, primary individual-level beta matrix, private patient scores, H matrices, or USIs are published or required for website rendering.

## Interpretation and ethics

See **[METHODS_AND_LIMITATIONS.md](METHODS_AND_LIMITATIONS.md)**. All figures are descriptive discovery-cohort summaries, with per-CpG significance, not gene-level or independent tests. All clinical counts/characteristics for insufficient-population 13C rows are withheld in the public payload, and small subgroup counts are not published for those rows. Nevertheless, before making the repository public, verify the cohort's sharing conditions and institutional approval to release aggregate rare-disease results, including top clinical associations.

All content loads locally in the browser with no external server/API requests. The bundled data are approximately 31 MB uncompressed, so mobile devices may load more slowly than the former promoter-only edition.
