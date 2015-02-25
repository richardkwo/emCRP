source("myutils.R")
source("dpEM_bottomUp_1dGaussian.R")

N <- 500
X <- mixture.of.gaussian(n = N, prop.vec = c(0.3,0.5,0.2), mean.vec = c(-2,0,2.5), 
                         sd.vec = sqrt(c(0.4,0.3,0.3)))

# X <- mixture.of.gaussian(n = N, prop.vec = c(1), mean.vec = c(-2), 
#                          sd.vec = sqrt(c(0.4)))
X$alpha <- 1
X$prior.mean <- mean(X$x)
X$prior.kappa <- 1
X$prior.a <- 1
X$prior.b <- 1

curve(X$density.fun(x), from=min(X$x)*1.1, to=max(X$x)*1.1, col="red", 
      main=paste("True density (red), and KDE (violet),N=",N))
curve(approxfun(density(X$x))(x), from=min(X$x)*1.1, to=max(X$x)*1.1, col="violet", add=T)
points(X$x, rep(0.01, length(X$x)))

result.1dGaussian <- bottom.up.dp.EM.1dGaussian(X, estimate.type = "Bayes", 
                                                W.tol = 1, inner.tol=0.1,
                                                init.K=10, cluster.size.threshold = 5, max.iter=20,
                                                verbose = 2, stepwise.plot=T)

# result.1dGaussian <- dp.EM.infer.1D.Gaussian(X, estimate.type = "MLE", random.init.cluster = F, 
#                                              cluster.size.threshold = 3, verbose = 1, stepwise.plot=T)
