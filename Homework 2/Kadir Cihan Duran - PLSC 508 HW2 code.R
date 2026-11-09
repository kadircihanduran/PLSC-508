library(network)
library(sna)
library(ergm)
library(latentnet)
library(coda)
library(ggplot2)
set.seed(508)

nodes <- read.csv(file.path("data", "bonobo_nodes.csv"), stringsAsFactors = FALSE)
edges <- read.csv(file.path("data", "bonobo_grooming_edges.csv"), stringsAsFactors = FALSE)
n <- nrow(nodes)

W <- matrix(0, n, n, dimnames = list(nodes$id, nodes$id))
W[cbind(edges$from, edges$to)] <- edges$weight
W <- W + t(W)

threshold <- median(edges$weight)
A <- (W > threshold) * 1

groom.net <- network(A, directed = FALSE, loops = FALSE, matrix.type = "adjacency")
groom.net %v% "sex"     <- nodes$sex
groom.net %v% "age"     <- nodes$age
groom.net %v% "rearing" <- nodes$rearing

desc.table <- data.frame(
  Network      = c("Any grooming", "Strong grooming (analysed)"),
  Nodes        = n,
  Edges        = c(sum(W > 0) / 2, network.edgecount(groom.net)),
  Density      = round(c(gden(W > 0, mode = "graph"), gden(A, mode = "graph")), 3),
  Transitivity = round(c(gtrans(W > 0, mode = "graph"), gtrans(A, mode = "graph")), 3),
  Mean.degree  = round(c(mean(rowSums(W > 0)), mean(rowSums(A))), 2),
  Centralization = round(c(centralization(W > 0, degree, mode = "graph"),
                           centralization(A, degree, mode = "graph")), 3)
)
desc.table
mixingmatrix(groom.net, "sex")

set.seed(508)
groom.xy <- network.layout.fruchtermanreingold(groom.net, layout.par = list())
ew <- W[as.matrix(groom.net, matrix.type = "edgelist")]
op <- par(mar = c(0, 0, 0, 0))
plot(groom.net, coord = groom.xy, displaylabels = TRUE, label.cex = 0.7, label.pos = 5,
     vertex.col = ifelse(nodes$sex == "Female", "#E69F00", "#56B4E9"),
     vertex.sides = ifelse(nodes$sex == "Female", 50, 4),
     vertex.cex = 1.2 + nodes$age / 15,
     vertex.border = ifelse(nodes$rearing == "Other", "red", "grey20"),
     vertex.lwd = ifelse(nodes$rearing == "Other", 3, 1),
     edge.col = "grey55", edge.lwd = 1 + 25 * ew)
legend("bottomleft", c("Female", "Male", "Atypically reared"), pch = c(21, 23, 1),
       pt.bg = c("#E69F00", "#56B4E9", NA), col = c("grey20", "grey20", "red"),
       pt.cex = 1.8, bty = "n")
par(op)

set.seed(508)
cug.density <- cug.test(groom.net, FUN = gden, mode = "graph", cmode = "size", reps = 1000)
set.seed(508)
cug.trans <- cug.test(groom.net, FUN = gtrans, mode = "graph", cmode = "edges", reps = 1000)
set.seed(508)
cug.cent <- cug.test(groom.net, FUN = centralization, mode = "graph", cmode = "edges",
                     reps = 1000, FUN.args = list(FUN = degree))

cug.table <- data.frame(
  Statistic = c("Density", "Transitivity", "Degree centralization"),
  Conditioning = c("size", "edges", "edges"),
  Observed = round(c(cug.density$obs.stat, cug.trans$obs.stat, cug.cent$obs.stat), 3),
  Null.mean = round(c(mean(cug.density$rep.stat), mean(cug.trans$rep.stat),
                      mean(cug.cent$rep.stat)), 3),
  Pr.X.ge.Obs = round(c(cug.density$pgteobs, cug.trans$pgteobs, cug.cent$pgteobs), 3),
  Pr.X.le.Obs = round(c(cug.density$plteobs, cug.trans$plteobs, cug.cent$plteobs), 3)
)
cug.table

par(mfrow = c(1, 3))
plot(cug.density, main = "Density | size")
plot(cug.trans, main = "Transitivity | edges")
plot(cug.cent, main = "Degree centralization | edges")
par(mfrow = c(1, 1))

female <- as.numeric(nodes$sex == "Female")
both.female <- outer(female, female) * 1
both.male   <- outer(1 - female, 1 - female) * 1
age.diff    <- abs(outer(nodes$age, nodes$age, "-"))
age.sum     <- outer(nodes$age, nodes$age, "+")
atypical    <- as.numeric(nodes$rearing == "Other")
rear.other  <- (outer(atypical, atypical, "+") > 0) * 1
for (m in c("both.female", "both.male", "age.diff", "age.sum", "rear.other")) {
  tmp <- get(m); diag(tmp) <- 0; assign(m, tmp)
}

qap.cor <- function(y, x, reps = 1000) {
  set.seed(508)
  q <- qaptest(list(y, x), gcor, g1 = 1, g2 = 2, reps = reps, mode = "graph")
  c(r = gcor(y, x, mode = "graph"), p.ge = q$pgreq, p.le = q$pleeq)
}
cov.list <- list(`Female-female` = both.female, `Male-male` = both.male,
                 `Age difference` = age.diff, `Age sum` = age.sum,
                 `Atypically reared (>=1)` = rear.other)
cor.table <- t(sapply(cov.list, function(x) qap.cor(A, x)))
cor.table <- data.frame(Covariate = rownames(cor.table), round(cor.table, 3), row.names = NULL)
cor.table

par(mfrow = c(1, 3))
for (k in c("Female-female", "Male-male", "Age sum")) {
  set.seed(508)
  plot(qaptest(list(A, cov.list[[k]]), gcor, g1 = 1, g2 = 2, reps = 1000, mode = "graph"))
  title(sub = k, col.sub = "navy", font.sub = 2)
}
par(mfrow = c(1, 1))

X <- list(both.female, both.male, age.diff, age.sum, rear.other)
xnames <- c("Intercept", "Female-female", "Male-male", "Age difference", "Age sum",
            "Atypically reared (>=1)")

logit.classical <- netlogit(A, X, mode = "graph", diag = FALSE, nullhyp = "classical")
set.seed(508)
logit.qap <- netlogit(A, X, mode = "graph", diag = FALSE, nullhyp = "qap", reps = 2000)
set.seed(508)
lm.qap <- netlm(W * 100, X, mode = "graph", diag = FALSE, nullhyp = "qap", reps = 2000)

qap.table <- data.frame(
  Predictor = xnames,
  logit.coef = round(logit.qap$coefficients, 3),
  odds.ratio = round(exp(logit.qap$coefficients), 2),
  logit.p.classical = round(logit.classical$pgreqabs, 3),
  logit.p.QAP = round(logit.qap$pgreqabs, 3),
  lm.coef.rate = round(lm.qap$coefficients, 3),
  lm.p.QAP = round(lm.qap$pgreqabs, 3)
)
qap.table

ctrl <- control.ergm(seed = 508, MCMC.samplesize = 5000, MCMC.burnin = 20000,
                     MCMC.interval = 1000, MCMLE.maxit = 30)
ergm0 <- ergm(groom.net ~ edges)
ergm1 <- ergm(groom.net ~ edges + nodematch("sex", diff = TRUE) + absdiff("age") +
                nodecov("age") + nodefactor("rearing", levels = "Other"))
ergm2 <- ergm(groom.net ~ edges + nodematch("sex", diff = TRUE) + absdiff("age") +
                nodecov("age") + nodefactor("rearing", levels = "Other") +
                gwesp(0.5, fixed = TRUE) + gwdegree(0.5, fixed = TRUE),
              control = ctrl)

ergm.tab <- function(m, tag) {
  s <- summary(m)$coefficients
  out <- data.frame(term = rownames(s), est = round(s[, 1], 3), se = round(s[, 2], 3),
                    p = round(s[, 5], 3), row.names = NULL)
  names(out)[2:4] <- paste0(names(out)[2:4], ".", tag)
  out
}
ergm.est <- Reduce(function(a, b) merge(a, b, by = "term", all = TRUE, sort = FALSE),
                   list(ergm.tab(ergm0, "M0"), ergm.tab(ergm1, "M1"), ergm.tab(ergm2, "M2")))
ergm.est

fit.table <- data.frame(Model = c("Model 0: edges", "Model 1: covariates", "Model 2: + gwesp + gwdegree"),
                        AIC = round(sapply(list(ergm0, ergm1, ergm2), AIC), 1),
                        BIC = round(sapply(list(ergm0, ergm1, ergm2), BIC), 1))
fit.table

mcmc.diagnostics(ergm2, which = c("plots"))

ergm.chain <- ergm2$sample
geweke.z <- geweke.diag(ergm.chain)[[1]]$z
ess <- effectiveSize(ergm.chain)
data.frame(statistic = names(geweke.z), geweke.z = round(geweke.z, 2),
           geweke.p = round(2 * pnorm(-abs(geweke.z)), 3), ESS = round(ess), row.names = NULL)

set.seed(508)
gof1 <- gof(ergm1, GOF = ~ model + degree + espartners + distance,
            control = control.gof.ergm(nsim = 200, seed = 508))
gof2 <- gof(ergm2, GOF = ~ model + degree + espartners + distance,
            control = control.gof.ergm(nsim = 200, seed = 508))
par(mfrow = c(2, 4), mar = c(4, 4, 2, 1))
plot(gof1, main = "")
plot(gof2, main = "")
par(mfrow = c(1, 1))

gof.p <- function(g) {
  sapply(list(degree = g$pval.deg, espartners = g$pval.espart, distance = g$pval.dist),
         function(m) sum(m[, "MC p-value"] < 0.05, na.rm = TRUE))
}
rbind(`Model 1 (# of bins with p < .05)` = gof.p(gof1),
      `Model 2 (# of bins with p < .05)` = gof.p(gof2))

lsm.ctrl <- ergmm.control(burnin = 100000, interval = 50, sample.size = 4000)
lsm.f <- groom.net ~ nodematch("sex", diff = TRUE) + absdiff("age") + nodecov("age") +
  nodefactor("rearing", levels = "Other")
set.seed(508)
lsm0 <- ergmm(groom.net ~ euclidean(d = 2), control = lsm.ctrl, seed = 508)
lsm1 <- ergmm(update(lsm.f, . ~ . + euclidean(d = 2)), control = lsm.ctrl, seed = 508)
lsm2 <- ergmm(update(lsm.f, . ~ . + euclidean(d = 2) + rsociality), control = lsm.ctrl, seed = 508)

l1 <- summary(lsm1)$pmean$coef.table
l2 <- summary(lsm2)$pmean$coef.table
lsm.est <- data.frame(term = c("(Intercept)", names(coef(ergm1))[-1]),
                      round(l1[, 1:3], 3), round(l2[, 1:3], 3), row.names = NULL)
names(lsm.est)[2:7] <- c("LSM1.mean", "LSM1.lower", "LSM1.upper",
                         "LSM2.mean", "LSM2.lower", "LSM2.upper")
lsm.est
round(summary(lsm2)$pmean$sociality.var, 2)

lsm.bic <- sapply(list(`LSM 0: positions` = lsm0, `LSM 1: covariates + positions` = lsm1,
                       `LSM 2: + sociality` = lsm2),
                  function(m) round(summary(m)$bic$overall, 1))
lsm.bic

mcmc.diagnostics(lsm1, which.diags = "trace")

coef.cols <- function(m) {
  ch <- as.mcmc.list(m)
  ch[, grep("^beta", varnames(ch)), drop = FALSE]
}
ch1 <- coef.cols(lsm1)
geweke.lsm <- geweke.diag(ch1)[[1]]$z
chains <- lapply(c(1, 22, 333, 4444), function(s)
  coef.cols(ergmm(update(lsm.f, . ~ . + euclidean(d = 2)), control = lsm.ctrl, seed = s))[[1]])
rhat <- gelman.diag(as.mcmc.list(chains), autoburnin = FALSE)$psrf
data.frame(parameter = c("(Intercept)", names(coef(ergm1))[-1]),
           ESS = round(effectiveSize(ch1)), geweke.z = round(geweke.lsm, 2),
           Rhat = round(rhat[, 1], 3), row.names = NULL)

ut <- upper.tri(A)
lsm.prob.draws <- function(m, nsim = 500) {
  idx <- round(seq(1, m$control$sample.size, length.out = nsim))
  lapply(idx, function(i) plogis(latentnet:::ergmm.eta(m$model, m$sample[[i]])))
}
draw.undirected <- function(P) {
  Y <- matrix(0, n, n); Y[ut] <- rbinom(sum(ut), 1, P[ut]); Y + t(Y)
}
gof.stats <- function(m) {
  m <- as.matrix(m); d <- rowSums(m)
  sp <- (m %*% m)[ut & m == 1]
  gd <- geodist(m)$gdist[ut]; gd[is.infinite(gd)] <- n
  list(degree = tabulate(d + 1, n), esp = tabulate(sp + 1, n - 1), dist = tabulate(gd, n))
}
ppc <- function(sims, obs = gof.stats(A)) {
  sapply(c("degree", "esp", "dist"), function(s) {
    sm <- t(sapply(sims, function(z) gof.stats(z)[[s]])); o <- obs[[s]]
    p <- sapply(seq_along(o), function(k) min(1, 2 * min(mean(sm[, k] >= o[k]), mean(sm[, k] <= o[k]))))
    list(list(sim = sm, obs = o, p = p))
  })
}
plot.ppc <- function(pp, lab) {
  xl <- list(degree = 0:(n - 1), esp = 0:(n - 2), dist = c(1:(n - 1), "Inf"))
  ttl <- c(degree = "degree", esp = "edge-wise shared partners", dist = "minimum geodesic distance")
  for (s in names(pp)) {
    sm <- pp[[s]]$sim; o <- pp[[s]]$obs
    nz <- which(colSums(sm) > 0 | o > 0)
    keep <- if (s == "dist") c(seq_len(max(setdiff(nz, n))), n) else seq_len(max(nz))
    boxplot(sm[, keep], names = xl[[s]][keep], xlab = ttl[s], ylab = "count", main = lab,
            col = "grey90", outline = FALSE)
    q <- apply(sm[, keep, drop = FALSE], 2, quantile, c(.025, .975))
    lines(seq_along(keep), q[1, ], lty = 2, col = "grey40"); lines(seq_along(keep), q[2, ], lty = 2, col = "grey40")
    lines(seq_along(keep), o[keep], lwd = 2); points(seq_along(keep), o[keep], pch = 16)
  }
}

set.seed(508)
P.l1 <- lsm.prob.draws(lsm1); P.l2 <- lsm.prob.draws(lsm2)
sims.l1 <- lapply(P.l1, draw.undirected); sims.l2 <- lapply(P.l2, draw.undirected)
ppc.l1 <- ppc(sims.l1); ppc.l2 <- ppc(sims.l2)
par(mfrow = c(2, 3), mar = c(4, 4, 2, 1))
plot.ppc(ppc.l1, "LSM 1"); plot.ppc(ppc.l2, "LSM 2")
par(mfrow = c(1, 1))
n.bad <- function(pp) sapply(pp, function(z) sum(z$p < .05))
rbind(`LSM 1 (# bins p < .05)` = n.bad(ppc.l1), `LSM 2 (# bins p < .05)` = n.bad(ppc.l2))

c(expected.ties.per.draw = mean(sapply(P.l1, function(P) sum(P[ut]))),
  builtin.simulate.median.ties = median(sapply(simulate(lsm1, nsim = 200)$networks, network.edgecount)))

auc <- function(p, y) { r <- rank(p); n1 <- sum(y); n0 <- length(y) - n1
  (sum(r[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0) }
stat.fun <- function(x) { m <- as.matrix(x)
  c(edges = sum(m[ut]), transitivity = gtrans(m, mode = "graph")) }
sim.tab <- function(sims) apply(t(sapply(sims, stat.fun)), 2, quantile, c(.025, .5, .975))
sims.e1 <- simulate(ergm1, nsim = 500, seed = 508)
sims.e2 <- simulate(ergm2, nsim = 500, seed = 508,
                    control = control.simulate.formula(MCMC.burnin = 20000, MCMC.interval = 1000))
tab <- rbind(observed = stat.fun(groom.net), sim.tab(sims.e1), sim.tab(sims.e2),
             sim.tab(sims.l1), sim.tab(sims.l2))
rownames(tab) <- c("observed", paste(rep(c("ERGM 1", "ERGM 2", "LSM 1", "LSM 2"), each = 3),
                                     c("2.5%", "50%", "97.5%")))
round(tab, 3)
c(AUC.ERGM1 = auc(predict(ergm1, output = "matrix")[ut], A[ut]),
  AUC.LSM1 = auc(Reduce(`+`, P.l1)[ut] / length(P.l1), A[ut]))

Z <- summary(lsm1)$pmean$Z
if (is.null(Z)) Z <- lsm1$mkl$Z
pos <- data.frame(id = nodes$id, x = Z[, 1], y = Z[, 2], sex = nodes$sex,
                  age = nodes$age, rearing = nodes$rearing)
el <- as.matrix(groom.net, matrix.type = "edgelist")
seg <- data.frame(x = Z[el[, 1], 1], y = Z[el[, 1], 2], xend = Z[el[, 2], 1], yend = Z[el[, 2], 2])
ggplot(pos, aes(x, y)) +
  geom_segment(data = seg, aes(x = x, y = y, xend = xend, yend = yend), colour = "grey75") +
  geom_point(aes(fill = sex, shape = rearing, size = age), colour = "grey20") +
  geom_text(aes(label = id), size = 2.6, vjust = -1.3) +
  scale_fill_manual(values = c(Female = "#E69F00", Male = "#56B4E9")) +
  scale_shape_manual(values = c(Other = 24, Parent = 21)) +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 4))) +
  coord_equal() + theme_minimal(base_size = 11) + labs(x = "Latent dim 1", y = "Latent dim 2")

s1 <- summary(ergm1)$coefficients; s2 <- summary(ergm2)$coefficients
fmt <- function(est, p) sprintf("%.2f%s", est, ifelse(p < .01, "**", ifelse(p < .05, "*", ifelse(p < .10, "+", ""))))
fmt.ci <- function(est, lo, hi) sprintf("%.2f [%.2f, %.2f]%s", est, lo, hi, ifelse(lo > 0 | hi < 0, "*", ""))
q.idx <- 2:6
comp <- data.frame(
  Effect = c("Female-female (vs mixed)", "Male-male (vs mixed)", "Age difference",
             "Age (per year, node)", "Atypically reared (node)"),
  QAP.logit = fmt(logit.qap$coefficients[q.idx], logit.qap$pgreqabs[q.idx]),
  ERGM.1 = fmt(s1[2:6, 1], s1[2:6, 5]),
  ERGM.2 = fmt(s2[2:6, 1], s2[2:6, 5]),
  LSM.1 = fmt.ci(l1[2:6, 1], l1[2:6, 2], l1[2:6, 3]),
  LSM.2 = fmt.ci(l2[2:6, 1], l2[2:6, 2], l2[2:6, 3])
)
comp <- rbind(comp, data.frame(Effect = "Triadic closure",
  QAP.logit = "n/a", ERGM.1 = "n/a", ERGM.2 = fmt(s2["gwesp.fixed.0.5", 1], s2["gwesp.fixed.0.5", 5]),
  LSM.1 = "latent space", LSM.2 = "latent space"))
comp

dyad.df <- data.frame(y = A[ut], sapply(X, function(x) x[ut]))
q.se <- summary(glm(y ~ ., data = dyad.df, family = binomial))$coefficients[q.idx, 2]
eff <- c("Female-female", "Male-male", "Age difference", "Age (node)", "Atypically reared")
forest <- rbind(
  data.frame(method = "QAP logit", effect = eff, est = logit.qap$coefficients[q.idx],
             lo = logit.qap$coefficients[q.idx] - 1.96 * q.se, hi = logit.qap$coefficients[q.idx] + 1.96 * q.se),
  data.frame(method = "ERGM 1", effect = eff, est = s1[2:6, 1], lo = s1[2:6, 1] - 1.96 * s1[2:6, 2], hi = s1[2:6, 1] + 1.96 * s1[2:6, 2]),
  data.frame(method = "ERGM 2", effect = eff, est = s2[2:6, 1], lo = s2[2:6, 1] - 1.96 * s2[2:6, 2], hi = s2[2:6, 1] + 1.96 * s2[2:6, 2]),
  data.frame(method = "LSM 1", effect = eff, est = l1[2:6, 1], lo = l1[2:6, 2], hi = l1[2:6, 3]),
  data.frame(method = "LSM 2", effect = eff, est = l2[2:6, 1], lo = l2[2:6, 2], hi = l2[2:6, 3]))
forest$effect <- factor(forest$effect, levels = rev(eff))
forest$method <- factor(forest$method, levels = c("QAP logit", "ERGM 1", "ERGM 2", "LSM 1", "LSM 2"))
ggplot(forest, aes(est, effect, colour = method)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.6), size = 0.3) +
  facet_wrap(~ effect == "Age difference" | effect == "Age (node)", scales = "free",
             labeller = as_labeller(c(`FALSE` = "Categorical effects", `TRUE` = "Age effects (per year)"))) +
  scale_colour_manual(values = c("#000000", "#0072B2", "#56B4E9", "#D55E00", "#E69F00")) +
  labs(x = "Estimate (log-odds)", y = NULL, colour = NULL) + theme_minimal(base_size = 11) +
  theme(legend.position = "bottom")

thr.set <- c(`Any grooming (>0)` = 0, `Above 25th pct` = unname(quantile(edges$weight, .25)),
             `Above median (main)` = threshold, `Above 75th pct` = unname(quantile(edges$weight, .75)))
sens <- lapply(names(thr.set), function(nm) {
  Ak <- (W > thr.set[nm]) * 1
  set.seed(508)
  cg <- cug.test(Ak, gtrans, mode = "graph", cmode = "edges", reps = 1000)
  set.seed(508)
  lq <- netlogit(Ak, X, mode = "graph", nullhyp = "qap", reps = 1000)
  nk <- network(Ak, directed = FALSE); nk %v% "sex" <- nodes$sex; nk %v% "age" <- nodes$age
  nk %v% "rearing" <- nodes$rearing
  ek <- ergm(nk ~ edges + nodematch("sex", diff = TRUE) + absdiff("age") + nodecov("age") +
               nodefactor("rearing", levels = "Other"))
  se <- summary(ek)$coefficients
  data.frame(Threshold = nm, Edges = sum(Ak) / 2, Transitivity = round(gtrans(Ak, mode = "graph"), 3),
             CUG.p.trans = cg$pgteobs,
             QAP.FF = fmt(lq$coefficients[2], lq$pgreqabs[2]),
             QAP.MM = fmt(lq$coefficients[3], lq$pgreqabs[3]),
             QAP.age = fmt(lq$coefficients[5], lq$pgreqabs[5]),
             ERGM1.MM = fmt(se[3, 1], se[3, 5]), ERGM1.age = fmt(se[5, 1], se[5, 5]))
})
do.call(rbind, sens)
