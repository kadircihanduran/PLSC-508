PLSC 508 - Final project: replication of Cranmer, Leifeld, McClurg & Rolfe (2017)
================================================================================

Kadir Cihan Duran, PLSC 508, final paper update (October 2026)

Replicated article:
  Cranmer, S. J., P. Leifeld, S. D. McClurg, and M. Rolfe. 2017. "Navigating the Range of
  Statistical Tools for Inferential Network Analysis." American Journal of Political Science
  61(1): 237-251.
  Replication data: Harvard Dataverse, doi:10.7910/DVN/2XP8YF
  (Swiss climate policy network; originally collected by Karin Ingold, 2008)

This folder holds the data and the R code that produce every figure and table in
"report/Kadir Cihan Duran PLSC 508 Update.pdf":
  - a reproduction of the article's Table 1 (logit, MRQAP, latent space model, ERGM) and
    Figures 2-5 (goodness of fit of the four models) (in the PDF: Table 1, Figures 1-5).
Everything was run on a Mac (Apple Silicon, R 4.5.2).


HOW TO RUN
----------
1. Install the R packages once:
     install.packages(c("network", "sna", "ergm", "latentnet", "ggplot2"))
2. In a terminal, from this folder (Final Project):
     Rscript replicate_Cranmer2017.R
   This uses the packages installed on your computer and writes to results_current_versions/.
   The first run takes about 5-15 minutes; the fitted ergm and LSM are saved and reused after that.
3. To use the package versions the authors used (ergm 3.8.0, latentnet 2.7.1), install them once
   into a separate folder (your normal R library is not touched), then run with the switch:
     bash compat/install_authors_versions.sh
     AUTHORS_VERSIONS=1 Rscript replicate_Cranmer2017.R
   This writes to results/, which holds the results in the PDF.
4. Rebuild the PDF (needs LaTeX):
     cd report
     pdflatex "Kadir Cihan Duran PLSC 508 Update.tex"
     pdflatex "Kadir Cihan Duran PLSC 508 Update.tex"


WHY TWO VERSIONS (the latent space model, LSM)
----------------------------------------------
The logit, MRQAP and ERGM columns of Table 1 are the same with old and new package versions
(ERGM within 0.021). The LSM column is not: with latentnet 2.12 the same code and seed give an
intercept of -2.00 instead of 0.95 and a BIC of 918 instead of 856. latentnet 2.7.1 (the version
the authors used) had a bug that transposed asymmetric covariate matrices (edgecov, nodeicov, ...) in
directed networks; it was fixed in 2.8.0. Fitting the model with 2.7.1 on the transposed
covariates reproduces the current-version numbers to within 0.01. So the published LSM column
was estimated with transposed covariates. check_lsm_transposed.R tests this on your computer:
on a current latentnet its "transposed" fit should reproduce the PUBLISHED LSM column; on
latentnet 2.7.1 the two are swapped.
The script only reproduces the results. The comparison with the published values (Table 2 of
the PDF) is typed into the LaTeX source; those values come from the article's Table 1.


FILES
-----
README.txt                 This file.
.gitignore                 Keeps LaTeX build files (aux, log, ...) and the compat/rlib folder out of git;
                           the output PDF is kept.
replicate_Cranmer2017.R    The ONE script that produces every result:
                             - reads data/*.csv and builds the model inputs (as in the authors' analysis.R),
                             - Table 1: logit, MRQAP, LSM and ERGM (authors' specifications, MCMC
                               settings and seed 12345),
                             - Figures 2-5: goodness of fit of the four models.
check_lsm_transposed.R     Optional. Fits the LSM with transposed asymmetric covariates and prints it
                           next to the published column and your own fit (see above); writes
                           LSM_transposed_check.csv to the same results folder.

data/                      INPUT DATA (Dataverse archive; collected by Ingold 2008)
  climate0205-collab.csv     Dependent network: directed collaboration ties among the 34 actors
                             (row actor reports collaborating with column actor; 1/0).
  climate0205-rep.csv        Influence attribution: row actor names column actor as particularly
                             influential (1/0). Its column sums are the vertex attribute
                             "influence reputation".
  climate0205-allop.csv      Perceived alliance (+1), neutral (0) or opposition (-1) between row and
                             column actor (dyadic covariate).
  climate0205-prefdist.csv   Preference dissimilarity: distance between the two actors' positions on
                             four policy issues (dyadic covariate).
  climate0205-committee.csv  Actor x policy forum membership (34 x 3). Its one-mode projection is the
                             "joint forum participation" covariate.
  climate0205-type.csv       Actor type (gov, ngo, private, science, party).
  climate0205-dm.csv         Formal decision-making power (1/0); not used.
  Codebook.pdf               The authors' codebook for these files.

results/                   RESULTS used in the PDF (made with ergm 3.8.0 / latentnet 2.7.1)
  Table1_reproduced.csv/.tex   Table 1 of the PDF = the article's Table 1, reproduced.
  Figure2_gof_logit.pdf        Figure 1 of the PDF = the article's Figure 2 (logit).
  Figure3_gof_mrqap.pdf        Figure 2 of the PDF = the article's Figure 3 (MRQAP).
  Figure4_gof_lsm.pdf          Figure 3 of the PDF = the article's Figure 4, first part (LSM).
  Figure4_latent_positions.pdf Figure 4 of the PDF = the article's Figure 4, second part.
  Figure5_gof_ergm.pdf         Figure 5 of the PDF = the article's Figure 5 (ERGM).
  Figures2to5_gof_coverage.csv Share of the observed network's statistics inside the 95% range of the
                               simulated networks, per model.
  LSM_transposed_check.csv     Output of check_lsm_transposed.R (means and 95% intervals).
  ergmEst.rds, lsmEst.rds      The fitted ERGM and LSM (saved so that reruns are quick).
results_current_versions/  Created when you run the script with the packages on your computer (step 2).

compat/                    ONLY NEEDED TO USE THE AUTHORS' PACKAGE VERSIONS (step 3)
  install_authors_versions.sh   Installs ergm 3.8.0 and latentnet 2.7.1 into compat/rlib. Run on
                                macOS Apple Silicon (R 4.5.2). For R >= 4.5 it
                                renames Calloc/Free/Realloc to R_Calloc/R_Free/R_Realloc in ergm, and it
                                links latentnet without the Fortran runtime (no gfortran needed).
  pkgs/                         Source files of ergm 3.8.0 and latentnet 2.7.1 (from CRAN archive; GPL).
  latentnet_2.7.1_R4.patch      Lets latentnet 2.7.1 build under R >= 4.2 and recent compilers
                                (compiler flags and two scalar '&&' conditions; no estimator changed).
  trust/                        Small stand-in for the 'trust' package, used only if CRAN's cannot be installed.
  rlib/                         Created by the install script (not kept in git).

report/
  Kadir Cihan Duran PLSC 508 Update.tex   LaTeX source of the PDF (reads ../results/).
  Kadir Cihan Duran PLSC 508 Update.pdf   The PDF.
