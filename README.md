# Pediatric AML Promoter Methylation Explorer

Interactive, aggregate-only promoter CpG methylation explorer for pediatric AML NMF molecular neighborhoods. **This repository and website are independent of the clinical bridge explorer.**

## Website

Once GitHub Pages is enabled on the repository main branch at `/(root)`, your website will be available at:

`https://YOUR_GITHUB_USERNAME.github.io/YOUR_REPOSITORY_NAME/`

Replace the placeholders with your actual account and new repository name. Do not use the GitHub source-code view URL for the website.

## Interactive features

- Search by promoter-overlapping gene or individual CpG.
- Original NMF group nodes colored by observed original-primary group mean beta for the selected CpGs.
- Directional bridge arrows colored by the sign of observed beta difference (incoming vs source nonbridge, or incoming vs target primary).
- Retain overlapping methylation-signature labels from the CpG suite, adjustable BH FDR, delta-beta and W prioritization.
- Display underlying CpG evidence and download a CSV of shown aggregate results.

All methylation statistics shown are pooled-cohort, exploratory and post-selection; do not describe the red/blue as longitudinal changes or causal transitions. The website has no private patient-level data, though aggregate pediatric oncology findings still require a publication and data-governance review before making the repository public.

See [SETUP_NEW_REPOSITORY.md](SETUP_NEW_REPOSITORY.md) for step-by-step setup, [METHODS_AND_LIMITATIONS.md](METHODS_AND_LIMITATIONS.md) for methodology, and [scripts/build_promoter_explorer.py](scripts/build_promoter_explorer.py) for updating the aggregate data.
