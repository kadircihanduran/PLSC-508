source("replicate_Cranmer2017.R")          

cat("\n\n#### Transposed-covariate LSM, latentnet", as.character(packageVersion("latentnet")), "####\n")
gov.t <- t(gov.ifactor); ngo.t <- t(ngo.ofactor); infrep.t <- t(infrep); allopp.t <- t(allopp)
fT <- as.formula(paste("nw.collab ~ euclidean(d = 2, G = 0) + edgecov(gov.t) + edgecov(ngo.t) +",
  "nodematch('orgtype') + edgecov(priv.ngo) + edgecov(forum) + edgecov(infrep.t) +",
  "nodeocov('influence') + absdiff('influence') + edgecov(prefdist) + edgecov(allopp.t)"),
  env = globalenv())
set.seed(seed)
mT <- ergmm(fT, seed = seed, control = control.ergmm(sample.size = 10000, burnin = 50000, interval = 100))
bT <- mT$sample$beta; cnT <- mT$model$coef.names
ord <- match(c("(Intercept)", "edgecov.priv.ngo", "edgecov.allopp.t", "edgecov.prefdist",
               "edgecov.forum", "edgecov.infrep.t", "nodeocov.influence", "absdiff.influence",
               "edgecov.gov.t", "edgecov.ngo.t", "nodematch.orgtype"), cnT)
cmp <- data.frame(variable = vars[1:11],
  published = c(0.95, -1.37, 0.00, -1.76, 1.52, 0.08, 0.01, 0.04, -0.46, -0.60, 1.17),
  your_original_fit = round(lsmT$est, 2),
  your_transposed_fit = round(colMeans(bT)[ord], 2),
  transposed_2.5 = round(apply(bT, 2, quantile, 0.025)[ord], 2),
  transposed_97.5 = round(apply(bT, 2, quantile, 0.975)[ord], 2))
print(cmp, row.names = FALSE)
cat("\nIf you run latentnet >= 2.8: 'your_transposed_fit' should be close to 'published'\n",
    "and 'your_original_fit' close to the current-version numbers (intercept about -2.0).\n",
    "If you run latentnet 2.7.1: the two columns are swapped.\n", sep = "")
write.csv(cmp, out("LSM_transposed_check.csv"), row.names = FALSE)

# For editing code to successfully run in my Mac, I got help from Claude