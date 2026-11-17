# Methods, provenance, and limitations — Integrated CpG Molecular Atlas v2

## Source files and traceability

The public payload is built directly from the supplied:

- **11D:** `Table_11D_C_ALL_Selected_CpG_Direction_Associations.csv` — 79,920 CpG × SOURCE→TARGET selected associations and 50,013 unique probes.
- **12C:** `Table_12C_B_ALL_Directional_CpG_Annotations_GRCh38.csv` — 79,920 matching 11D associations with GRCh38 positions, class, CpG island context, promoter/gene-body overlap, nearest-gene distance and original 11D statistics/overlap flags.
- **13C:** `Table_13C_02_Clinical_Directional_CpG_Integrated_Atlas.csv` — 1,254 direction × clinical-population summary records; 418 SOURCE→TARGET combinations × All/Female/Male.

The build validates that 11D and 12C identifiers (`CpG`, `Source`, `Target`, `Direction`) match in **exact order** and that CpG metadata are identical across directed occurrences. It does not silently join or discard rows.

### Genomic assignment choices

Within Script 12C the promoter annotation is transcript TSS−2,000 through TSS+500 bp. The mutually exclusive **primary region label** is `Promoter`, `Gene body`, or `Intergenic`. A promoter CpG is indexed by its documented `Promoter_Overlap_Genes`, a gene-body CpG by `GeneBody_Overlap_Genes`, and an intergenic CpG by `Nearest_Gene`, explicitly labeled a distance-only association. Multiple overlapping genes are retained; one CpG can appear under several gene symbols, but is counted **once per group/site** within a selected gene. Gene-body overlaps can also exist for a probe labeled Promoter; the selected label is determined by the primary Script 12C `Region_GRCh38` class.

The original genome annotation contains 24,792 promoter, 41,397 gene-body, and 13,731 intergenic CpG × direction rows. These correspond to **15,048**, **26,012**, and **8,953** *unique CpGs* respectively; summing the unique counts gives 50,013. Island categories (Island, Shore, Shelf, Open sea) are directly reported by 12C; no genome annotations are inferred in browser code.

## Original NMF and H-defined direction terminology

For SOURCE→TARGET, target primary **P** means original primary assignment is TARGET; incoming **I** means primary assignment is SOURCE and TARGET is in the operational Script 03 H neighborhood; SOURCE nonbridge **N** means primary assignment is SOURCE and TARGET is **not** in that H neighborhood. P, I and N are disjoint, and the SOURCE primary assignment is not changed for incoming patients.

This is not movement through time, cell lineage progression, or a follow-up sample comparison. The original NMF model and patient memberships are held fixed.

## Circle methylation: original-primary absolute beta

For CpG j in group g, target primary means are directly available as `Mean_Beta_Primary` when g is TARGET. When g is SOURCE, source original-primary can be reconstructed from disjoint I and N populations using CpG-specific finite-observation n:

\[\bar\beta_{g,j}=\frac{n_{I,j}\bar\beta_{I,j}+n_{N,j}\bar\beta_{N,j}}{n_{I,j}+n_{N,j}}.\]

The builder checks the equality of all duplicate group×CpG estimates (including direct versus reconstructed means) at 1e−6 tolerance. The public data contain **134,073 unique group×CpG means across 19 eligible groups**, with 25,767 consistency cross-checks and no conflicts. Four other NMF labels are shown as neutral without fabricated values.

For the current gene, region, signature and CpG-island filters, each group node displays the **unweighted average of its available unique CpG mean betas**. It is not a statistical test; a node colored light red has higher *absolute β* (toward 1), and light blue lower *absolute β* (toward 0). No normal marrow comparator or group contrast is included in node color. CpGs can have unequal availability among groups; click the circle for site coverage. The CpG-wise q/effect arrow filters do not redefine which CpGs contribute to node means (the category/region/island/W filters do).

## Arrows: selected CpG contrasts, not methylation trajectories

For each 11D selected CpG j at SOURCE→TARGET:

- Default effect \(\Delta_{IN,j}=\bar\beta_{I,j}-\bar\beta_{N,j}\) (BH q: incoming vs source nonbridge).
- Alternative effect \(\Delta_{IP,j}=\bar\beta_{I,j}-\bar\beta_{P,j}\) (BH q: primary vs incoming; the effect sign is reversed from the source's `Delta_Primary_Minus_Incoming`).
- Alternative effect \(\Delta_{PN,j}=\bar\beta_{P,j}-\bar\beta_{N,j}\) (BH q: primary vs source nonbridge).

A CpG qualifies for an arrow if it meets the selected signed effect threshold (default |Δβ|≥0.15) and CpG-level adjusted q cutoff (default q≤0.05), plus chosen context/island/flag/W/group filters. **Only finite tested q-values qualify**; selecting "any tested q" never substitutes missing q-values. W values prioritize but do not set the sign. Different arrows may have different qualifying CpG counts.

If all qualifying effects for a direction are positive, the arrow is red; all negative, blue. If some are positive and some negative, the arrow is **gray dashed and labeled mixed**. Mean effect for an arrow is the descriptive unweighted mean across qualifying selected CpGs, not a gene-level test or a gene-specific FDR. Arrow line thickness encodes qualifying CpG count, not patient count. Hover/click gives pooled 13C directional patient counts and detailed CpG stats.

## Original statistical inference and overlapping flags

11D tested Welch P vs I, P vs N and I vs N, plus P/I TOST equivalence within ±0.10 beta, with **separate genome-wide BH adjustment within each direction and test family**. The website retains all four raw p and BH q columns. Script 11D's original group-identity and neighborhood-conservation criteria are retained; the website does not re-run the statistical tests.

Each CpG×direction retains the complete **22-bit original label mask** across Scripts 06, 06b, 07a, 11C, and 11D. Labels intentionally overlap. Some 06b formal shared flags have **zero matches in this selected 11D association subset**; that is not evidence that the original entire dataset contained no formally shared CpGs. Sensitivity flags `_010` remain preserved and are not substituted for the `_015` primary criterion.

`W_Target_Loading_Percentile`, `W_Target_Component_Specificity`, `W_Target_Column_Share` and `W_Target_Sensitivity_Weight` are displayed as prioritization annotations only. Nonnegative W coefficients do **not** tell whether methylation is higher or lower.

## 13C clinical association summaries

The site retains 13C's available 418 direction pairs for each of All, Female, and Male (`1254` rows). 13C marks **153**, **74**, and **65** eligible clinical directions respectively. The displayed **Top_Feature** and **Top_Level** are the **single reported top characteristic** per direction/population, not the full set of tests from 13B. Top `FDR_Within_Population` and `FDR_Overall` refer **only to that clinical association**, not CpG methylation. A top feature can have q>0.05; do not infer statistical significance merely because it was reported as the top feature.

For ineligible clinical rows, public counts and top clinical feature/effects/q are replaced by null, rather than disclosing small bridge/nonbridge categories. The output's default minimum subgroup count is 5. Of the 1,254 summary records, 962 are ineligible/withheld. Pooled 11D CpG-label counts from 13C are reported separately from clinical significance. These safeguards are not a complete release risk assessment.

The **Female/Male labels in 13C are sex-stratified CLINICAL analyses only**, not sex-specific Welch/TOST methylation comparisons, and do not justify female/male red/blue methylation arrows. To color by sex, re-run observed methylation contrasts within sex-specific P/I/N subsets and implement separate effective-sample-size/multiple-testing controls.

## Scope and caveats

1. This is the **selected 11D discovery suite** (50,013 unique CpGs), not all ~729,000 EPIC probes, independent hold-out analysis, or differentially methylated region (DMR) results.
2. CpGs were selected in data related to the same H-defined neighborhoods being visualized; observed agreement is **not independent validation**. Additional H-blind methylation/PCA/UMAP comparisons and held-out cohort verification remain separate analyses.
3. Many tested CpGs are correlated within a gene/promoter. A mean of per-CpG effects and the number of qualifying probes **cannot** be interpreted as a gene-level meta-analysis or q-value.
4. Promoter methylation, gene-body methylation, nearest-gene labels and CpG-island category do not prove transcriptional consequences, gene expression, causality, lineage evolution, or clinical actionability. Use matched RNA-seq to assess associations separately.
5. Potential confounding includes sex, batch, blast percentage, cellular mixture, molecular subtype and treatment selection. No beta-value models adjusted for those variables are exported here.
6. Publishing high-dimensional aggregate methylation and clinical results in pediatric rare disease may require institutional and cohort data-use approval even without USIs. All GitHub Pages data files are public if the repository is public.
7. The page cannot calculate DMRs from the CpG summary alone. Future DMR integration needs region boundaries, region-level tests/effects, the overlapping CpGs and gene contexts plus source-specific P/I/N summaries.

## Reproducibility / export

- `scripts/build_full_atlas.py` performs source-key checks, gene mapping, no-duplication checks for node means, count-based source reconstruction, 22-label retention and eligibility redaction.
- `data/atlas.js` is a browser-readable public aggregate dataset, not a secret file.
- `data/build-audit.json` lists counts and provenance; `tests/check_public_data.py` verifies schema, region sums, row counts and redaction.
- All computations and filters in the browser are descriptive reformatting of the reported CpG-wise comparisons. No new significance values are computed.
