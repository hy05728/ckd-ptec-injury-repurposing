# Single-cell injury signature of proximal tubules guides drug repurposing in CKD

Code and derived result tables for the manuscript:

> **A single-cell injury signature of proximal tubules guides drug repurposing in chronic kidney disease**
> Yue Hu, Yangzhou University Medical College

This repository contains the analysis pipeline behind the reported results: construction of a
proximal tubule epithelial cell (PTEC) injury signature from human single-cell RNA sequencing,
reversal-based compound ranking, structure-based validation against ALK5 (TGFBR1), pathway
enrichment, independent-cohort reproducibility, and molecular dynamics / MM-PBSA.

---

## Data sources

| Resource | Accession / DOI | Used for |
|---|---|---|
| Human diabetic kidney scRNA-seq (Wilson et al., 2019) | GEO **GSE131882** | Discovery cohort; injury signature |
| Human kidney atlas (Lake et al., 2023) | CELLxGENE collection `bcb61471-2a44-4d00-a0af-ff085512674c` | Independent-cohort validation |
| DrugReflector model checkpoints | Zenodo **10.5281/zenodo.16912444** | Compound ranking |
| CMap / L1000 metadata (GSE92742) | GEO **GSE92742** | Compound ID to name mapping |
| ALK5 kinase domain structure | PDB **2WOU** chain A | Molecular docking and MD |
| LINCS compound SMILES | GSE92742 `pert_info` | Ligand construction |

All source data are public. No new experimental data were generated.

---

## Environment

Two environments are used, because single-cell analysis and MD have incompatible dependency sets.

**1. Single-cell / docking / figures (Python 3.13)**

```bash
python -m venv .venv && source .venv/bin/activate     # Windows: .venv\Scripts\activate
pip install scanpy anndata leidenalg pandas numpy scipy matplotlib \
            rdkit openbabel-wheel drugreflector
```

`drugreflector` is installed from source: `pip install git+https://github.com/Cellarity/drugreflector.git`

R 4.4+ with `SingleCellExperiment`, `Matrix`, `org.Hs.eg.db`, `clusterProfiler`, `enrichplot`.

Docking uses the AutoDock Vina 1.1.2 binary (`vina`), not the Python bindings.

**2. Molecular dynamics (conda, Python 3.11)**

```bash
conda env create -f md/environment.yml     # or: bash md/setup_env.sh on a Linux GPU host
```

See `md/README_CN.md` for a full walkthrough (Chinese) of the AutoDL deployment.

---

## Pipeline

Scripts are numbered in execution order. Each script defines its working directory at the top;
adjust `ROOT` / `D` / `RES` to your local paths before running.

| Step | Script | Output |
|---|---|---|
| 0 | `scripts/00_download.py` | Downloads GEO supplementary files and Zenodo checkpoints |
| 1 | `scripts/convert_rds.R` | Converts double-gzipped `dgecounts` objects to `matrix.mtx` + `genes.tsv` + `barcodes.tsv` |
| 2 | `scripts/01_load_merge.py` | Merges six samples into one AnnData object (121,788 cells before QC) |
| 3 | `scripts/02_annotate_vscore.py` | QC (20,703 cells), clustering, PTEC annotation (7,788 PTECs), v-score signature |
| 4 | `scripts/03_rank.py` | DrugReflector reversal ranking of 9,597 compounds, drug-class recall |
| 5 | `scripts/04_dock_alk5.py` | AutoDock Vina docking against ALK5 (2WOU chain A) |
| 6 | `scripts/05_benchmark.py` | 60-decoy enrichment benchmark (EF, ROC-AUC) |
| 7 | `scripts/07_enrich_signature.R` | GO / KEGG enrichment of the signature |
| 8 | `scripts/08b_lake_disease_validate.py` | Independent-cohort reproducibility in the Lake 2023 atlas |
| 9 | `scripts/10_make_figures.py` | Figures 1-5 (PNG 300 dpi, TIFF 600 dpi, PDF/EPS) |
| 10 | `scripts/12_verify_numbers.py` | Recomputes every ranking statistic quoted in the manuscript |
| 11 | `scripts/14_permutation_null.py` | Compound-level permutation test (100 shuffled signatures) |
| 12 | `scripts/14b_category_null.py` | Class-level permutation test (200 shuffled signatures) |
| 13 | `scripts/18_marker_vs_unsupervised.py` | Marker-free re-clustering; consistency with the marker-based injury label |
| 14 | `scripts/19_sensitivity_ranking.py` | Ranking sensitivity to the injury-state definition |
| 15 | `scripts/17_build_benchmark_set.py` | Builds the ChEMBL-scale docking benchmark (actives from CHEMBL4439 + property-matched decoys) |
| 16 | `scripts/20_dock_benchmark_v2.py` | Docks the 518-compound benchmark set |
| 17 | `scripts/21_recompute_metrics.py` | ROC-AUC / EF / BEDROC with bootstrap CI (note: score = -affinity) |
| 18 | `scripts/16_make_graphical_abstract.py` | Graphical abstract (`scripts/wf.py` draws the workflow) |
| 19 | `scripts/13_selfcheck.py`, `scripts/audit_prose.py` | Pre-submission numeric and prose checks |
| MD | `md/01_prep.sh` → `md/02_md.py` → `md/03_mmpbsa.sh` → `md/04_analyze.py` | 100 ns MD and MM/PBSA binding free energies |

---

## Key derived tables (`results/`)

| File | Content |
|---|---|
| `vscore_injury_normal2injured.csv` | PTEC injury signature (healthy → injured), one v-score per gene |
| `vscore_reversal_injured2normal.csv` | Negated signature used as the reversal query |
| `ranks_reversalinjured_normal.csv` | Full ranking of 9,597 compounds, position 1 = most likely to reverse |
| `docking_results.csv` | Vina affinities and SMILES of the 11 docked compounds |
| `benchmark_decoys.csv` | Decoy docking scores underlying the enrichment analysis |
| `marker_concordance.csv` | Injury-marker v-scores across the discovery and validation cohorts |
| `enrichment_summary.txt` | Top GO and KEGG terms for the upregulated and downregulated signature |

The full ranked compound list in `ranks_reversalinjured_normal.csv` is indexed by LINCS
perturbation identifiers (BRD-...); map them to compound names with the L1000 metadata
(`GSE92742_Broad_LINCS_sig_info.txt.gz`).

---

## Selected results

| Analysis | Result |
|---|---|
| QC | 121,788 → 20,703 cells; 7,788 PTECs (6,269 healthy, 1,519 injured) |
| Injury markers (v-score) | VCAM1 1.03, KRT8 0.77, HAVCR1 0.34, LCN2 0.08 |
| Reversal ranking | RAS inhibitors in the top 1.3% (irbesartan 54, telmisartan 66, enalapril 125) |
| Top reversal hits | Y-27632 (1), Y-39983 (2), LY-364947 (6), YC-1 (10), fasudil (11) |
| Docking benchmark | EF5% = 8.0, EF10% = 4.0, EF20% = 3.0, ROC-AUC 0.873 |
| Independent cohort | Pearson r = 0.73, Spearman ρ = 0.67 across 21,772 genes |
| Permutation test | ROCK P = 0.005, ALK5 P = 0.020; RAS inhibitors not significant (0.29 / 0.22) |
| Marker-free redefinition | ROCK P = 0.030, ALK5 P = 0.010; signature concordance r = 0.80 |
| Docking benchmark (ChEMBL-scale) | ROC-AUC 0.904 (95% CI 0.873-0.932), BEDROC 3.80, 120 actives + 398 decoys |

---

## Notes on reproducibility

- The injured-PTEC label is marker-based (any of HAVCR1, VCAM1, LCN2, KRT8), not derived from
  unsupervised state annotation. This choice is deliberate (it transfers across cohorts) and is
  discussed as a limitation in the manuscript.
- `GSE131882` sample `control.s2` is deposited as an unfiltered droplet matrix (~92,000 barcodes);
  the QC thresholds remove the empty droplets. Final cell numbers match the original publication.
- Docking was run at exhaustiveness 8 with random seed 42, so scores are reproducible on the
  same Vina build.
- The comparison with nephrotoxic agents rests on five compounds present in the L1000 library and
  is reported descriptively (Mann-Whitney one-sided P = 0.13).

---

## License

MIT (see `LICENSE`). The underlying datasets retain their original licenses.
