# PLSC 508

My work for **PLSC 508** (Penn State, Fall 2026). Both parts use the Swiss climate policy collaboration network (Ingold 2008; 34 actors), available from the replication archives of Cranmer et al. (2017) and Minhas, Hoff & Ward (2019).

## Folders

| Folder | Role |
|---|---|
| `Homework 1/` | Network visualisation and descriptive network measures. |
| `Homework 1/Kadir_Cihan_Duran_PLSC508_HW1.ipynb` | The Python notebook with the code, figures and written interpretation. Run it from inside `Homework 1/Datasets/`. |
| `Homework 1/Datasets/` | `nodes.csv` and `edges.csv`: the network as a node list and an edge list. |
| `Final Project/` | Updated research proposal and replication of Cranmer, Leifeld, McClurg & Rolfe (2017), *American Journal of Political Science* 61(1): 237-251 (data: doi:10.7910/DVN/2XP8YF). Details below and in `Final Project/README.txt`. |

## Final Project

Everything needed to produce the figures and tables of `Final Project/report/Kadir Cihan Duran PLSC 508 Update.pdf` is in this folder. Run `Rscript replicate_Cranmer2017.R` from inside `Final Project/` (R packages: `network`, `sna`, `ergm`, `latentnet`, `ggplot2`).

| File or folder | Role |
|---|---|
| `replicate_Cranmer2017.R` | The one script that produces every result: Table 1 (logit, MRQAP, latent space model, ERGM) and Figures 2-5 (goodness of fit). |
| `check_lsm_transposed.R` | Optional check of why the latent space model differs between package versions (`latentnet` 2.7.1 transposed asymmetric covariates; fixed in 2.8.0). |
| `data/` | The seven raw `climate0205-*.csv` files and the authors' `Codebook.pdf`. |
| `results/` | Every figure and table in the PDF (made with the authors' package versions), and the fitted ERGM and LSM. |
| `results_current_versions/` | Created when the script is run with the packages installed on your computer. |
| `compat/` | Installs the authors' package versions (`ergm` 3.8.0, `latentnet` 2.7.1) into a separate folder: `bash compat/install_authors_versions.sh`, then `AUTHORS_VERSIONS=1 Rscript replicate_Cranmer2017.R`. |
| `report/` | The LaTeX source and the output PDF. |
| `README.txt` | The detailed description of each file in `Final Project/`. |
| `.gitignore` | Keeps LaTeX build files out of git; the output PDF is kept. |

The logit, MRQAP and ERGM columns of the article's Table 1 reproduce (coefficients within 0.02, BIC within 1). The latent space model column only reproduces with `latentnet` 2.7.1; see `Final Project/README.txt`.
