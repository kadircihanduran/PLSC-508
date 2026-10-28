.f <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
if (length(.f)) setwd(dirname(normalizePath(.f[1])))

dataPath <- "data/"
USE_AUTHORS_VERSIONS <- identical(Sys.getenv("AUTHORS_VERSIONS"), "1")
NSIM     <- 100

if (USE_AUTHORS_VERSIONS) {
  if (!dir.exists("compat/rlib")) stop("compat/rlib not found. Run first:  bash compat/install_authors_versions.sh")
  .libPaths(c(normalizePath("compat/rlib"), .libPaths()))
}
lnVersion <- tryCatch(as.character(packageVersion("latentnet")), error = function(e) "not installed")
if (USE_AUTHORS_VERSIONS && lnVersion != "2.7.1")
  stop("latentnet ", lnVersion, " is loaded, not 2.7.1. Start a fresh R session (packages from an earlier ",
       "session stay loaded), or re-run bash compat/install_authors_versions.sh")
outPath <- if (lnVersion == "2.7.1") "results/" else "results_current_versions/"
cat("latentnet", lnVersion, "-> results are written to", outPath, "\n")

seed <- 12345
set.seed(seed)
dir.create(outPath, showWarnings = FALSE, recursive = TRUE)
out <- function(f) file.path(outPath, f)

pkgs <- c("network", "sna", "ergm", "latentnet", "ggplot2")
miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
if (length(miss)) stop("Please install: ", paste(miss, collapse = ", "))
suppressPackageStartupMessages({
  library(network); library(sna); library(ergm); library(latentnet); library(ggplot2)
})
theme_set(theme_bw(base_size = 10))

cached <- function(file, expr) {
  f <- out(file)
  if (file.exists(f)) return(readRDS(f))
  t0 <- Sys.time(); res <- expr
  cat(sprintf("  %s estimated in %.1f min\n", file,
              as.numeric(difftime(Sys.time(), t0, units = "mins"))))
  saveRDS(res, f); res
}
ergmFormula <- function(lhs) as.formula(paste(lhs,
  "~ edges + edgecov(collab.t) + edgecov(gov.ifactor) + edgecov(ngo.ofactor) +",
  "nodematch('orgtype') + edgecov(priv.ngo) + edgecov(forum) + edgecov(infrep) +",
  "nodeicov('influence') + absdiff('influence') + edgecov(prefdist) + edgecov(allopp) +",
  "odegree1.5 + twopath + gwidegree(2, fixed = TRUE) + gwesp(1, fixed = TRUE) +",
  "gwodegree(0.5, fixed = TRUE)"), env = globalenv())
lsmFormula <- function(lhs) as.formula(paste(lhs,
  "~ euclidean(d = 2, G = 0) + edgecov(gov.ifactor) + edgecov(ngo.ofactor) +",
  "nodematch('orgtype') + edgecov(priv.ngo) + edgecov(forum) + edgecov(infrep) +",
  "nodeicov('influence') + absdiff('influence') + edgecov(prefdist) + edgecov(allopp)"),
  env = globalenv())


rd <- function(f) as.matrix(read.table(file.path(dataPath, f), header = TRUE,
                                       row.names = 1, sep = ";"))
forumAff <- rd("climate0205-committee.csv")  
infrep   <- rd("climate0205-rep.csv")        
collab   <- rd("climate0205-collab.csv")     
allopp   <- rd("climate0205-allop.csv")      
prefdist <- rd("climate0205-prefdist.csv")   
types <- as.character(read.table(file.path(dataPath, "climate0205-type.csv"),
                                 header = TRUE, row.names = 1, sep = ";")[, 2])
n <- nrow(collab)

forum <- forumAff %*% t(forumAff); diag(forum) <- 0   
influence <- sna::degree(infrep, cmode = "indegree")    
nw.collab <- network(collab)
set.vertex.attribute(nw.collab, "orgtype", types)
set.vertex.attribute(nw.collab, "influence", influence)

collab.t <- t(collab)
priv.ngo <- outer(types, types, function(a, b)
  1 * ((a == "private" & b == "ngo") | (a == "ngo" & b == "private")))
gov.ifactor <- matrix(1 * (types == "gov"), n, n, byrow = TRUE)
ngo.ofactor <- matrix(1 * (types == "ngo"), n, n, byrow = FALSE)
type.nodematch <- 1 * outer(types, types, "==")
influence.icov <- matrix(influence, n, n, byrow = TRUE)
influence.absdiff <- abs(outer(influence, influence, "-"))

logit.data <- data.frame(collab = c(collab), prefdist = c(prefdist),
  ngo.ofactor = c(ngo.ofactor), gov.ifactor = c(gov.ifactor), forum = c(forum),
  infrep = c(infrep), collab.t = c(collab.t), type.nodematch = c(type.nodematch),
  priv.ngo = c(priv.ngo), influence.icov = c(influence.icov),
  influence.absdiff = c(influence.absdiff), allopp = c(allopp))

typeLab <- c(gov = "Government", ngo = "NGO", private = "Private / business",
             science = "Scientific", party = "Political party")
typeCol <- c(gov = "#2a78d6", ngo = "#1baf7a", private = "#eb6834",
             science = "#6250d6", party = "#e34948")

covariates <- list(priv.ngo, allopp, prefdist, forum, infrep, influence.icov,
                   influence.absdiff, gov.ifactor, ngo.ofactor, type.nodematch, collab.t)
names(covariates) <- c("priv.ngo", "allopp", "prefdist", "forum", "infrep", "influence.icov",
                       "influence.absdiff", "gov.ifactor", "ngo.ofactor", "type.nodematch", "collab.t")

cat("Estimating models (first run: ERGM ~10 min, LSM ~4 min; logit and MRQAP seconds)...\n")
model.logit <- glm(collab ~ priv.ngo + allopp + prefdist + forum + infrep + influence.icov +
                     influence.absdiff + gov.ifactor + ngo.ofactor + type.nodematch + collab.t,
                   data = logit.data, family = binomial)
set.seed(seed)
model.qap <- netlogit(nw.collab, covariates, reps = 200, nullhyp = "qap")
model.ergm <- cached("ergmEst.rds", {
  set.seed(seed)
  ergm(ergmFormula("nw.collab"), eval.loglik = TRUE,
       control = control.ergm(seed = seed, MCMC.samplesize = 5000, MCMC.interval = 5000)) })
model.ls <- cached("lsmEst.rds", {
  set.seed(seed)
  ergmm(lsmFormula("nw.collab"), seed = seed,
        control = control.ergmm(sample.size = 10000, burnin = 50000, interval = 100)) })

vars <- c("Intercept/Edges", "Business vs. NGO", "Opposition/alliance",
  "Preference dissimilarity", "Joint forum participation", "Influence attribution",
  "Alter's influence indegree", "Influence absolute diff.", "Alter = Government actor",
  "Ego = Environmental NGO", "Same actor type", "Mutuality", "Outdegree popularity",
  "Twopaths", "GWIdegree (2.0)", "GWESP (1.0)", "GWOdegree (0.5)")
N <- n * (n - 1)                              

lc <- summary(model.logit)$coefficients
logitT <- data.frame(est = lc[, 1], u = lc[, 2], star = lc[, 4] < 0.05)
# MRQAP
qapT <- data.frame(est = model.qap$coefficients, u = NA, star = model.qap$pgreqabs < 0.05)
# ERGM
es <- summary(model.ergm)
ec <- if (!is.null(es$coefficients)) es$coefficients else es$coefs      # ergm 4 / ergm 3
ergmOrder <- c(1, 6, 12, 11, 7, 8, 9, 10, 3, 4, 5, 2, 13:17)
ergmT <- data.frame(est = ec[ergmOrder, 1], u = ec[ergmOrder, 2], star = ec[ergmOrder, ncol(ec)] < 0.05)
# LSM
bs <- model.ls$sample$beta
cn <- model.ls$model$coef.names
if (is.null(cn) || length(cn) != ncol(bs)) cn <- colnames(bs)
lsNames <- c("(Intercept)", "edgecov.priv.ngo", "edgecov.allopp", "edgecov.prefdist",
  "edgecov.forum", "edgecov.infrep", "nodeicov.influence", "absdiff.influence",
  "edgecov.gov.ifactor", "edgecov.ngo.ofactor", "nodematch.orgtype")
lsOrder <- match(lsNames, cn)
if (anyNA(lsOrder)) stop("LSM coefficient names not found. This fit has: ", paste(cn, collapse = ", "))
lo <- apply(bs, 2, quantile, 0.025)[lsOrder]; hi <- apply(bs, 2, quantile, 0.975)[lsOrder]
lsmT <- data.frame(est = colMeans(bs)[lsOrder], lo = lo, hi = hi, star = lo > 0 | hi < 0)

bic <- c(Logit = BIC(model.logit), MRQAP = model.qap$bic,
         LSM = tryCatch(summary(model.ls)$bic$overall, error = function(e) NA_real_),
         ERGM = tryCatch(-2 * as.numeric(logLik(model.ergm)) + length(coef(model.ergm)) * log(N),
                         error = function(e) NA_real_))

fmt <- function(d, interval = FALSE) {
  m <- nrow(d); e <- s <- rep("", length(vars))
  e[1:m] <- paste0(sprintf("%.2f", round(d$est, 2) + 0), ifelse(d$star, " *", ""))
  if (interval) s[1:m] <- sprintf("[%.2f; %.2f]", d$lo, d$hi)
  else if (!all(is.na(d$u))) s[1:m] <- sprintf("(%.2f)", d$u)
  c(rbind(e, s)) }
table1 <- data.frame(Variable = c(rbind(vars, "")), Logit = fmt(logitT), MRQAP = fmt(qapT),
                     LSM = fmt(lsmT, TRUE), ERGM = fmt(ergmT))
table1 <- rbind(table1, data.frame(Variable = "BIC", Logit = sprintf("%.2f", bic["Logit"]),
  MRQAP = sprintf("%.2f", bic["MRQAP"]), LSM = sprintf("%.2f", bic["LSM"]),
  ERGM = sprintf("%.2f", bic["ERGM"])))
cat("\n=== Table 1 (reproduced). Logit/ERGM: (s.e.); LSM: [95% interval]; * p < 0.05 or 0 outside the interval ===\n")
print(table1, row.names = FALSE, right = FALSE)
write.csv(table1, out("Table1_reproduced.csv"), row.names = FALSE)

esc <- function(x) gsub(" \\*", "$^{*}$", x)
tex1 <- c("\\begin{tabular}{lrrrr}", "\\toprule",
          "Variable & Logit & MRQAP & LSM & ERGM \\\\", "\\midrule",
          apply(table1, 1, function(r) paste0(ifelse(r[1] == "", "", ""), paste(esc(r), collapse = " & "), " \\\\")),
          "\\bottomrule", "\\end{tabular}")
tex1 <- gsub("\\[", "{[}", gsub("\\]", "{]}", tex1)); tex1 <- gsub("{{[}", "{[}", tex1, fixed = TRUE)
writeLines(tex1, out("Table1_reproduced.tex"))

KMAX <- 12
gofStats <- function(A) {                                  
  nw <- network(A, directed = TRUE)
  cnt <- function(x, len) { x <- x[seq_len(len)]; x[is.na(x)] <- 0; as.numeric(x) }
  s <- summary(nw ~ dsp(0:KMAX) + esp(0:KMAX) + idegree(0:KMAX) + odegree(0:KMAX) + istar(2:KMAX))
  g <- sna::geodist(A, inf.replace = NA)$gdist; g <- g[row(g) != col(g)]
  geo <- c(sapply(1:6, function(k) sum(g == k, na.rm = TRUE)), sum(is.na(g)))
  c(setNames(s[1:(KMAX + 1)], paste0("dsp|", 0:KMAX)),
    setNames(s[(KMAX + 2):(2 * KMAX + 2)], paste0("esp|", 0:KMAX)),
    setNames(geo, paste0("geodesic|", c(1:6, "Inf"))),
    setNames(s[(2 * KMAX + 3):(3 * KMAX + 3)], paste0("ideg|", 0:KMAX)),
    setNames(s[(3 * KMAX + 4):(4 * KMAX + 4)], paste0("odeg|", 0:KMAX)),
    setNames(s[(4 * KMAX + 5):length(s)], paste0("istar|", 2:KMAX))) }
obsStats <- gofStats(collab)
simBern <- function(P) { Y <- matrix(rbinom(length(P), 1, P), nrow(P)); diag(Y) <- 0; Y }
xbeta <- function(b, int = TRUE) {                          
  eta <- Reduce("+", Map(function(x, k) x * k, covariates, b[-1])) + b[1]; eta }
simulateModel <- function(model) {
  set.seed(seed)
  if (model == "Logit") lapply(1:NSIM, function(i) simBern(plogis(xbeta(coef(model.logit)))))
  else if (model == "MRQAP") lapply(1:NSIM, function(i) simBern(plogis(xbeta(model.qap$coefficients))))
  else if (model == "ERGM") {
    s <- simulate(ergmFormula("nw.collab"), coef = unname(coef(model.ergm)), nsim = NSIM, seed = seed,
                  control = control.simulate.formula(MCMC.burnin = 30000, MCMC.interval = 10000))
    lapply(s, as.sociomatrix)
  } else {                                                  
    draws <- round(seq(1, nrow(bs), length.out = NSIM))
    tryCatch(lapply(draws, function(i) {
      eta <- Reduce("+", Map(function(x, b) x * b, model.ls$model$X, bs[i, ]))
      simBern(plogis(eta - as.matrix(dist(model.ls$sample$Z[i, , ])))) }),
      error = function(e) {                                
        message("Manual LSM simulation failed (", conditionMessage(e), "); using simulate()")
        lapply(simulate(model.ls, nsim = NSIM, seed = seed)$networks, as.sociomatrix) }) } }
statTitle <- c(dsp = "Dyad-wise shared partners", esp = "Edge-wise shared partners",
               geodesic = "Geodesic distances", ideg = "Indegree", odeg = "Outdegree",
               istar = "Incoming k-star")
fig <- c(Logit = 2, MRQAP = 3, LSM = 4, ERGM = 5)
coverage <- c()
for (mod in names(fig)) {
  sims <- t(sapply(simulateModel(mod), gofStats))
  long <- function(M, type) do.call(rbind, lapply(colnames(M), function(nm) {
    p <- strsplit(nm, "|", fixed = TRUE)[[1]]
    data.frame(stat = factor(statTitle[p[1]], statTitle), x = factor(p[2], unique(p[2])),
               value = M[, nm]) }))
  d <- do.call(rbind, lapply(names(statTitle), function(st) {
    cols <- grep(paste0("^", st, "\\|"), colnames(sims)); x <- long(sims[, cols, drop = FALSE])
    x$x <- factor(x$x, unique(sub(".*\\|", "", colnames(sims)[cols]))); x }))
  o <- do.call(rbind, lapply(names(statTitle), function(st) {
    cols <- grep(paste0("^", st, "\\|"), names(obsStats))
    data.frame(stat = factor(statTitle[st], statTitle),
               x = factor(sub(".*\\|", "", names(obsStats)[cols]), sub(".*\\|", "", names(obsStats)[cols])),
               value = unname(obsStats[cols])) }))
  d <- d[d$stat != "Incoming k-star" | TRUE, ]

  keep <- unique(c(paste(d$stat, d$x)[d$value > 0], paste(o$stat, o$x)[o$value > 0]))
  d <- d[paste(d$stat, d$x) %in% keep, ]; o <- o[paste(o$stat, o$x) %in% keep, ]
  g <- ggplot(d, aes(x, value)) + geom_boxplot(outlier.size = 0.4, fill = "grey90", linewidth = 0.3) +
    geom_line(data = o, aes(group = 1), colour = "#2a78d6", linewidth = 0.8) +
    geom_point(data = o, colour = "#2a78d6", size = 1.2) +
    facet_wrap(~stat, scales = "free", nrow = 2) +
    labs(x = NULL, y = "Count", title = sprintf("Figure %d: %s model - goodness of fit", fig[mod],
         c(Logit = "Logit", MRQAP = "MRQAP", LSM = "Latent space", ERGM = "ERGM")[mod])) +
    theme(axis.text.x = element_text(size = 6))
  ggsave(out(sprintf("Figure%d_gof_%s.pdf", fig[mod], tolower(mod))), g, width = 9, height = 6)
  inside <- sapply(colnames(sims), function(nm) {
    q <- quantile(sims[, nm], c(0.025, 0.975)); obsStats[nm] >= q[1] & obsStats[nm] <= q[2] })
  coverage[mod] <- mean(inside)
}
cat("\n=== Figures 2-5: share of the observed network's statistics inside the 95% range of the simulated networks ===\n")
print(round(coverage, 2))
write.csv(data.frame(model = names(coverage), share_inside_95 = round(coverage, 3)),
          out("Figures2to5_gof_coverage.csv"), row.names = FALSE)

Zm <- tryCatch(model.ls$mkl$Z, error = function(e) NULL)
if (is.null(Zm)) Zm <- apply(model.ls$sample$Z, c(2, 3), mean)
lp <- data.frame(x = Zm[, 1], y = Zm[, 2], actor = rownames(collab),
                 type = factor(types, levels = names(typeLab)))
figZ <- ggplot(lp, aes(x, y, colour = type)) + geom_text(aes(label = actor), size = 3) +
  scale_colour_manual(values = typeCol, labels = typeLab, name = "Actor type") + coord_equal() +
  labs(x = NULL, y = NULL, title = "Figure 4 (positions): estimated latent positions (MKL)")
ggsave(out("Figure4_latent_positions.pdf"), figZ, width = 6, height = 5)
cat("\nDone. Figures and tables are in", normalizePath(outPath), "\n")

# For editing code to successfully run in my Mac, I got help from Claude