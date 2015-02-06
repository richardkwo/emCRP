library(lineprof)
source("dpEM.R")
X <- sample.dp.normal.mixture(N = 1000, alpha = 2, prior.mean = c(0,0))
l <- lineprof(dpem.results<-dp.EM.infer(X, posterior.predictive = T, verbose = T, cluster.size.threshold = 2))