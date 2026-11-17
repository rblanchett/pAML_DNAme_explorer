# Methods, interpretation, and audit — Promoter Methylation Network v1

## Verified input provenance

- `Table_11D_C_ALL_Selected_CpG_Direction_Associations.csv`: 79,920 selected CpG × SOURCE→TARGET association rows, covering 50,013 distinct CpGs and 119 directed pair comparisons across all genomic contexts.
- `Table_12C_B_ALL_Directional_CpG_Annotations_GRCh38.csv`: those same 79,920 association rows plus GRCh38 genomic annotation. All 79,920 CpG/source/target/direction keys matched the supplied 11D table **in order**; no independent inference or unverified merge was performed.
- `Table_13C_02_Clinical_Directional_CpG_Integrated_Atlas.csv`: pooled (`Population=All`) eligible-direction bridge/nonbridge counts and clinical metadata. No sex-specific methylation statistics are inferred from the Female/Male clinical rows.
- `Region_GRCh38 == "Promoter"` and an actual `Promoter_Overlap_Genes` annotation yields **24,792 CpG–direction rows, 15,048 distinct CpGs, 9,968 gene annotation labels** (after splitting semicolon-delimited multiple overlapping promoters), and **113 directions**. The genome annotation's promoter definition in Script 12C is 2 kb upstream / 500 bp downstream of the annotated TSS. Because gene labels can be noncoding transcripts, these are *promoter-overlap gene/transcript labels*, not all verified coding genes.

## Population definitions, preserving Script 11D exactly

For SOURCE→TARGET, the **P population** is original primary TARGET patients, **I** is original SOURCE patients including TARGET in their operational Script 03 H neighborhood, and **N** is original SOURCE patients excluding TARGET from their operational H neighborhood. These are mutually exclusive comparisons, not transitions in calendar time, disease evolution, or NMF reassignment.

## Original-primary group node means

At CpG `j`, TARGET's primary mean beta is directly given by `Mean_Beta_Primary`. For SOURCE, the incoming SOURCE-bridge and SOURCE-nonbridge subsets partition its original-primary patients. With CpG-specific observed counts `n_Ij` and `n_Nj`, we reconstruct:

```
mean_beta_SOURCE,j = (n_Ij * mean_beta_I,j + n_Nj * mean_beta_N,j) / (n_Ij + n_Nj)
```

The reconstructed source-group means match observed target-primary means wherever the CpG/group is available in both roles (at machine precision). **41,046 unique group–CpG mean values across 19 eligible primary groups were recovered with no consistency conflicts**. This is valid only because the incoming and nonbridge groups fully partition source original-primary membership and the counts refer to finite observed beta data at the *same CpG*.

When a gene is selected, the node color is the arithmetic mean of this group's available, selected **unique promoter CpG mean-beta values**. These are CpG-level means over original primary patients, not W weights, not fold changes, and not survival or gene-expression data. To avoid double counting, each CpG is counted only once per original group even when that group participates in multiple directed comparisons.

Nodes use a **continuous pale blue (β near zero) → white (β≈0.5) → pale red (β near one)** color scale. Red and blue on *nodes* therefore mean higher or lower **absolute methylation level**, **not differential hyper/hypomethylation against healthy marrow or a comparable reference group**. Values refer to the represented selected CpG panel; if a group has no value for those CpGs the node is neutral gray. Different groups may have different CpG coverage; the popup shows the number available so users can judge comparability. This is a descriptive average, not a validated gene-level methylation measure.

The graph displays all 23 original NMF labels. Only 19 were eligible in Script 11D (original-primary size ≥5, with an eligible directional analysis). The remaining 4 are shown neutrally without inventing original group sizes or methylation levels.

## Directed edges and statistics

The default per-CpG arrow effect is `Δβ_I−N = mean_beta_Incoming − mean_beta_Source_Nonbridge` for SOURCE→TARGET. **Red** means incoming bridge patients have **higher** observed methylation than SOURCE nonbridges. **Blue** means incoming patients have **lower** methylation. In the optional alternative contrast, `Δβ_I−P = mean_beta_Incoming − mean_beta_TargetPrimary` and the relevant reported BH-adjusted Welch q-value is used. Thus the red/blue reference changes when the user changes the contrast.

An individual CpG contributes to a displayed arrow if it has the selected suite-label category, meets the chosen minimum target W percentile (optional sensitivity prioritization), has |Δβ| ≥0.15 (changeable), and its **CpG-wise, per-direction BH** adjusted test q-value is ≤0.05 (changeable). A displayed arrow has at least one qualifying CpG. Direction colors are **not computed from W**.

At gene-level, multiple qualifying CpGs on the same arrow may show conflicting signed effects. If so the arrow is **gray dashed** and explicitly reports the number of higher/lower sites. If all qualifying sites agree, use red or blue based on their shared direction. Arrow line width scales weakly with the count of qualifying CpGs. The mean Δβ in the detail panel is the unweighted arithmetic mean *among qualifying CpGs*; no new gene-level p-value is calculated. A single gene can have multiple promoter CpGs (and a CpG can overlap multiple promoters); associations are counted separately by gene annotation label but not repeatedly within a given gene.

Original `Primary_Restricted` and `Neighborhood_Conserved` labels originate from Script 11D's Welch and TOST framework and are **not** recalculated by the site. The full 22 suite provenance flags (from 06, 06b, 07a, 11C, 11D) are retained per original CpG × direction row; overlapping labels remain overlapping. The display is a *selected CpG-suite* view. It **does not screen unselected CpGs or redo a genome-wide multiple-testing correction**.

## W weights

The `W_Target_Loading_Percentile`, `W_Target_Component_Specificity` and `W_Target_Column_Share` columns from 11D/12C are preserved. The viewer can *filter* on target W percentile to prioritize selected CpGs, but W is nonnegative and does **not** provide hyper/hypomethylation direction. Group and arrow colors come exclusively from observed beta means or their signed contrasts.

## Public export, provenance, and limitations

The page ships an in-browser `atlas.js` containing aggregate measurements and provenance; it never contacts your R environment or GitHub API. It does **not** export individual beta measurements, USIs, individual clinical records or H profiles. The public data are nonetheless high-dimensional summaries from real pediatric patient groups. Please check the source cohort's data use restrictions before public release.

The CpG panel is based on H-defined neighborhoods and selected within this same discovery cohort, so shared molecular patterns cannot be treated as independent confirmation of the H-based group definitions. Welch and TOST statistical adjustments in Script 11D were separate per directed comparison and test family, not familywise across every CpG *and* 113 directed comparisons. CpG-level q-values are not gene-level q-values. Cross-sectional differences do not establish gene transcription, promoter functional status, differentiation, causality or survival implications. Sex-, treatment-, batch-, cellular-composition- and risk-adjusted observed-beta summaries are **not** available from these source tables.

A future DMR map should use separately tested genomic *regions* (with region-level mean methylation, delta, uncertainty, q and mapping), not retroactively designate averaged individual CpG values as tested DMR results.
