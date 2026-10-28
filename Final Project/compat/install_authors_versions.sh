#!/bin/bash
###############################################################################
# Installs the package versions the authors of Cranmer et al. (2017) used for the
# latent space model, into a SEPARATE library folder (compat/rlib). Your normal R
# library is not touched.
#
#   ergm 3.8.0, latentnet 2.7.1   (+ statnet.common 4.9.0 and network 1.18.2, the versions
#                                   this combination was tested with)
#
# Run it once, from the "Final Project" folder:
#     bash compat/install_authors_versions.sh
# Then run the analysis with the old versions:
#     AUTHORS_VERSIONS=1 Rscript replicate_Cranmer2017.R
#
# Needs: R, a C compiler (macOS: Xcode command line tools, "xcode-select --install"), and
# an internet connection (CRAN). No Fortran compiler is needed. Takes about 5-10 minutes.
# latentnet 2.7.1 is patched so that it builds under R >= 4.2 and recent compilers
# (compat/latentnet_2.7.1_R4.patch: compiler flag -fcommon, FCONE definitions, and two scalar
# '&&' conditions). No estimator is changed.
###############################################################################
set -e
cd "$(dirname "$0")/.."                       # the Final Project folder
LIB="$PWD/compat/rlib"; mkdir -p "$LIB"
TMP=$(mktemp -d)
REPO=https://cloud.r-project.org
ARCHIVE=https://cran.r-project.org/src/contrib/Archive

# Compiler flags for old C code on new compilers (used only for this installation)
MK="$TMP/Makevars"
echo "CFLAGS = -g -O2 -std=gnu17 -fcommon -Wno-implicit-function-declaration -Wno-int-conversion -Wno-incompatible-pointer-types" > "$MK"
echo "FLIBS =" >> "$MK"      # these packages have no Fortran code; do not link the gfortran runtime
export R_MAKEVARS_USER="$MK"
export R_LIBS="$LIB"

if [ -z "$SKIP_DEPS" ]; then
  echo "== 1. Other packages the old versions need (only those you do not have) =="
  Rscript -e 'need <- setdiff(c("coda","Matrix","MASS","robustbase","lpSolve","rle","trust","sna","mvtnorm","abind"),
                              rownames(installed.packages()))
              if (length(need)) install.packages(need, lib = Sys.getenv("R_LIBS"), repos = "'$REPO'")'
  echo "== 2. statnet.common 4.9.0 and network 1.18.2 (pinned) =="
  Rscript -e 'for (u in c("'$ARCHIVE'/statnet.common/statnet.common_4.9.0.tar.gz",
                          "'$ARCHIVE'/network/network_1.18.2.tar.gz"))
                install.packages(u, repos = NULL, type = "source", lib = Sys.getenv("R_LIBS"))'
fi
# 'trust': CRAN's version if it is installed, else the small local re-implementation
Rscript -e 'quit(status = !requireNamespace("trust", quietly = TRUE))' || R CMD INSTALL -l "$LIB" compat/trust

echo "== 3. ergm 3.8.0 =="
# R >= 4.5 removed the old Calloc/Free/Realloc macros; ergm 3.8.0 uses them in one file.
# Rename them to R_Calloc/R_Free/R_Realloc (same functions) in a temporary copy.
mkdir -p "$TMP/ergm_src" && tar xzf compat/pkgs/ergm_3.8.0.tar.gz -C "$TMP/ergm_src"
perl -pi -e 's/\b(Calloc|Free|Realloc)\(/R_$1(/g' "$TMP/ergm_src/ergm/src/changestats.dgw_sp.c"
R CMD INSTALL --preclean -l "$LIB" "$TMP/ergm_src/ergm"

echo "== 4. latentnet 2.7.1 (patched for R >= 4.2) =="
tar xzf compat/pkgs/latentnet_2.7.1.tar.gz -C "$TMP"
PATCH="$PWD/compat/latentnet_2.7.1_R4.patch"
(cd "$TMP" && patch -p1 < "$PATCH")
# latentnet has no Fortran code: drop $(FLIBS) so it does not link a gfortran runtime that may be missing
printf 'PKG_LIBS = $(LAPACK_LIBS) $(BLAS_LIBS)\nPKG_CFLAGS = -fcommon\n' > "$TMP/latentnet/src/Makevars"
R CMD INSTALL --preclean -l "$LIB" "$TMP/latentnet"

echo "== Check =="
Rscript -e 'suppressPackageStartupMessages({library(ergm); library(latentnet)})
            cat("ergm", as.character(packageVersion("ergm")), " latentnet", as.character(packageVersion("latentnet")), "\n")
            stopifnot(as.character(packageVersion("latentnet")) == "2.7.1")'
echo "Done. Now run:  AUTHORS_VERSIONS=1 Rscript replicate_Cranmer2017.R"
