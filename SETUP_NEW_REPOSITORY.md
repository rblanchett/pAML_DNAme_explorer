# Set up a **separate** GitHub repository

This edition is packaged for a new, independent repository. **Do not add these files to the clinical atlas repository.**

1. On GitHub, choose **New repository**. Give it a name such as `pAML_Promoter_Methylation_Atlas`. Do not confuse this with the existing clinical atlas repository. Choose public visibility **only after confirming institutional/data-sharing permission** for the aggregate results in `data/atlas.js`. GitHub Pages can require public visibility depending on plan and configuration.
2. Download/unzip `pAML_Promoter_Methylation_Explorer_Standalone_Repo.zip` locally.
3. In the new GitHub repository, upload **the contents of the unzipped folder**. That means `index.html`, `assets/`, `data/`, `README.md`, etc. must be directly at the repository root, **not** under a `promoter-methylation/` subfolder. Do not upload the ZIP itself.
4. Commit the files.
5. Open **Settings → Pages** in the new repository. Choose **Deploy from a branch**, select the default branch (usually `main`) and **`/(root)`**, and Save.
6. When GitHub Pages has published, open `https://YOUR_GITHUB_USERNAME.github.io/YOUR_REPOSITORY_NAME/` (substitute your actual names). The GitHub repository page showing the HTML source code is **not** the live explorer page.
7. Check that the demo search for `GATA3` loads 23 NMF circles, 3 arrows under the default filters, and a CpG evidence table. If there is a load error, verify that `data/atlas.js` is in the repository root's `data/` folder.

**Optional:** Add a link from the clinical atlas README *later*, if you decide you want cross-navigation. The repositories and deployments remain independent either way.

## Updating the source data locally

From the new repository's root folder, with Python 3.9+ and pandas:

```bash
python scripts/build_promoter_explorer.py \
  --annotation Table_12C_B_ALL_Directional_CpG_Annotations_GRCh38.csv \
  --clinical Table_13C_02_Clinical_Directional_CpG_Integrated_Atlas.csv \
  --11d Table_11D_C_ALL_Selected_CpG_Direction_Associations.csv \
  --output data/atlas.js
```

Do **not** upload the 11D/12C/13C raw tables or any individual-patient files to the public site. Check that annotations and statistics come from compatible runs. Commit the regenerated `data/atlas.js` and `data/build-audit.json` after review. For a quick smoke test with Node 18+: `node tests/smoke_test.js`.

### Data safety

This package includes aggregate beta statistics and sample counts, genomic coordinates, and statistical tests, but **not** patient-level beta matrices or USIs. Rare-group summary outputs may still require institutional permission before public dissemination. A private repository is not automatically a published GitHub Pages site on all plans.
