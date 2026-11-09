# PLSC 508

Coursework for **PLSC 508** (Penn State, Fall 2026): network visualisation, statistical network models and a replication study.

## Contents

| Part | Topic | Data | Language |
|---|---|---|---|
| [Homework 1](#homework-1) | Network visualisation and descriptive measures | Swiss climate policy network | Python |
| [Homework 2](#homework-2) | CUG, QAP, ERGM and latent space models, compared | Bonobo grooming network | R |
| [Final Project](#final-project) | Research proposal update and replication of Cranmer et al. (2017) | Swiss climate policy network | R, LaTeX |

## Repository layout

```
PLSC-508/
├── README.md
├── Homework 1/
│   ├── Kadir_Cihan_Duran_PLSC508_HW1.ipynb
│   └── Datasets/            nodes.csv, edges.csv
├── Homework 2/
│   ├── Kadir Cihan Duran - PLSC 508 HW2.pdf
│   ├── Kadir Cihan Duran - PLSC 508 HW2 code.R
│   └── data/                bonobo_nodes.csv, bonobo_grooming_edges.csv
└── Final Project/
    ├── replicate_Cranmer2017.R
    ├── check_lsm_transposed.R
    ├── data/  results/  compat/  report/
    └── README.txt
```

---

## Homework 1

Network visualisation and descriptive network measures for the Swiss climate policy collaboration network (Ingold 2008; 34 actors).

| File | Role |
|---|---|
| `Kadir_Cihan_Duran_PLSC508_HW1.ipynb` | Python notebook with code, figures and written interpretation. |
| `Datasets/nodes.csv` | Node list. |
| `Datasets/edges.csv` | Edge list. |

---

## Homework 2

Four methods applied to one network, and a comparison of whether they agree: conditional uniform graph (CUG) tests, QAP, ERGM and latent space models (LSM). Convergence and model fit are checked for the ERGM and the LSM.

| File | Role |
|---|---|
| `Kadir Cihan Duran - PLSC 508 HW2.pdf` | The write-up: results, interpretation and comparison across methods. |
| `Kadir Cihan Duran - PLSC 508 HW2 code.R` | The one script that produces every table and figure in the write-up. |
| `data/bonobo_nodes.csv` | Node attributes: `id`, `sex`, `age` (years), `rearing` (Parent = reared by the mother, Other = atypically reared). |
| `data/bonobo_grooming_edges.csv` | Undirected weighted edge list: `from`, `to`, `weight` (proportion of scans with grooming). |

**Data source:**
- Torfs, J. R. R., Stevens, J. M. G., Verspeek, J., Laméris, D. W., Guéry, J.-P., Eens, M., & Staes, N. (2023). Multi-group analysis of grooming network position in a highly social primate. *PLoS ONE, 18*(4), e0284361. https://doi.org/10.1371/journal.pone.0284361
- Network curated in the Animal Social Network Repository: Sah, P., Méndez, J. D., & Bansal, S. (2019). A multi-species repository of social networks. *Scientific Data, 6*, 44. File `Network_1027.graphml`, github.com/bansallab/asnr, commit 3678a65. The two CSV files are a direct conversion of that file.

---

## Final Project

Updated research proposal and replication of Cranmer, Leifeld, McClurg & Rolfe (2017), *American Journal of Political Science* 61(1): 237-251 (data: doi:10.7910/DVN/2XP8YF), using the Swiss climate policy network.

**How to run.** From inside `Final Project/`, run `Rscript replicate_Cranmer2017.R`. It produces Table 1 (logit, MRQAP, latent space model, ERGM) and Figures 2-5 (goodness of fit) of `report/Kadir Cihan Duran PLSC 508 Update.pdf`. Required R packages: `network`, `sna`, `ergm`, `latentnet`, `ggplot2`.

**Authors' package versions.** To install the versions used in the article (`ergm` 3.8.0, `latentnet` 2.7.1) into a separate folder, run `bash compat/install_authors_versions.sh`, then `AUTHORS_VERSIONS=1 Rscript replicate_Cranmer2017.R`.

**Result.** The logit, MRQAP and ERGM columns of the article's Table 1 reproduce (coefficients within 0.02, BIC within 1). The latent space model column reproduces only with `latentnet` 2.7.1; see `Final Project/README.txt`.

| File or folder | Role |
|---|---|
| `replicate_Cranmer2017.R` | The one script that produces every result. |
| `check_lsm_transposed.R` | Optional check of why the latent space model differs between package versions (`latentnet` 2.7.1 transposed asymmetric covariates; fixed in 2.8.0). |
| `data/` | The seven raw `climate0205-*.csv` files and the authors' `Codebook.pdf`. |
| `results/` | Every figure and table in the PDF (made with the authors' package versions), plus the fitted ERGM and LSM. |
| `results_current_versions/` | Created when the script runs with the packages installed on your computer. |
| `compat/` | Installer for the authors' package versions. |
| `report/` | LaTeX source and output PDF. |
| `README.txt` | Detailed description of each file in `Final Project/`. |
| `.gitignore` | Keeps LaTeX build files out of git; the output PDF is kept. |
