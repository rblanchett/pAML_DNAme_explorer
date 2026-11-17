# Update your existing standalone promoter repository to include the complete selected CpG suite

**Do not make a new repository.** This package is the next version of your existing standalone promoter explorer and has `index.html` at its top level. Your separate clinical molecular neighborhood atlas should remain untouched.

## Option A — GitHub website

1. Extract `pAML_Integrated_CpG_Atlas_Standalone_Update.zip` on your computer.
2. Open **the standalone methylation repository** on GitHub, not the clinical one. The website should currently have `index.html` at its root.
3. Make a new branch for this update if you want a rollback point (recommended).
4. Upload/replace the contents of the extracted folder at the **repository root**. Especially important: overwrite `index.html`, `assets/app.js`, `assets/style.css` and `data/atlas.js`. Also upload `.nojekyll`, the docs, `scripts/build_full_atlas.py`, `data/build-audit.json`, and `tests/check_public_data.py`.
5. If your repository still has the old root-level `build_promoter_explorer.py`, delete it or move it to an archive folder. It belongs to the old promoter-only schema, not the new v2 dataset. Do not use the old script to regenerate `data/atlas.js`.
6. Commit, review the GitHub Pages deployment, then visit your existing website URL. No GitHub Pages settings or repository name need to change.

**Do not upload the ZIP archive itself as a website file.** Upload its unpacked contents, with `index.html` directly at repository root (not inside an additional parent directory). GitHub's web upload interface may require uploading files from the relevant directories separately.

## Option B — Git in a terminal (recommended for many files)

Use this *inside a clean local clone of your existing methylation repository*, replacing both paths with your own:

```bash
cd /path/to/your/standalone/methylation/repo
git status
git switch -c expand-all-11d-cpgs
rsync -av --exclude='.git' /path/to/EXTRACTED_PACKAGE/ ./
# If present, remove the obsolete promoter-only builder:
# git rm build_promoter_explorer.py
python tests/check_public_data.py
git add -A
git commit -m "Expand methylation atlas to all 11D CpGs and 13C clinical summaries"
git push -u origin expand-all-11d-cpgs
```

Review the branch and merge it to the branch used by GitHub Pages. If you push directly to `main`, omit the branch creation step.

## What changed

- All selected 11D CpGs, not just promoter-overlapping sites; **promoter remains the default**.
- Search full annotated gene list (including transcript labels), or directly search by a CpG ID.
- CpG island and genomic-context filtering; 22 retained suite flags and four Welch/TOST p/q families.
- Original group means across selected CpGs, W statistics, per-site n/SDs, and an enlarged details panel.
- All 13C clinical direction/population summaries, with safety redaction for ineligible small groups; this is **not** sex-stratified methylation.
- Individual filtered CpG and public 13C CSV export.

## Prepublication safety checks

- Only publish data you are permitted to share. This package has **aggregate** selected CpG beta summaries and clinical top-association metadata; it is not anonymized by a controlled disclosure assessment.
- No private patient-level CSVs (`*_PRIVATE.csv`), USIs, individual methylation profiles, raw count matrices, RDS files, or original unredacted clinical tables should be committed.
- Verify gene annotation scope, current source consent/data-use restrictions, and the public-redaction policy with your team.
- GitHub Pages serves all files in the public repo; do not assume that hiding columns in the interface hides data.

## Verification after deployment

The heading is now **Integrated CpG molecular network**. The header should show `50,013 unique CpGs` and `79,920 CpG × direction rows`. GATA3 loads in **Promoter** mode; switching to Gene body or All changes the evidence. You can inspect WT1, select an arrow, see All/Female/Male clinical summaries, search for an exact CpG, and export filtered evidence.
