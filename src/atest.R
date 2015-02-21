library(lineprof)
source("dp.R")
source("dpEM.R")

l <- lineprof(dp.EM.infer.multivariate(X, alpha = 5, prior.mean = colMeans(X$x), prior.scale.matrix = cov(X$x), prior.df = 3, prior.precision.factor = 1, posterior.predictive = T, random.init.cluster = T, cluster.size.threshold = 4, stepwise.plot = T, inner.tol = 1e-1, tol = 1e-2))
shine(l)
