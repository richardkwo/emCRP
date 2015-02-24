source("myutils.R")
source("dpEM_1dGaussian.R")
N <- 300
X <- mixture.of.gaussian(n = N, prop.vec = c(1), mean.vec = c(0), 
                         sd.vec = sqrt(c(0.05)))
X$alpha <- 2
X$prior.mean <- mean(X$x)
X$prior.kappa <- 1
X$prior.a <- 1
X$prior.b <- 1

curve(X$density.fun(x), from=min(X$x)*1.1, to=max(X$x)*1.1, col="red", 
      main=paste("True density (red), and KDE (violet),N=",N))
curve(approxfun(density(X$x))(x), from=min(X$x)*1.1, to=max(X$x)*1.1, col="violet", add=T)
points(X$x, rep(0.01, length(X$x)))

result.1dGaussian <- dp.EM.infer.1D.Gaussian(X, estimate.type = "MLE", random.init.cluster = F, 
                                             cluster.size.threshold = 10, verbose = 1, stepwise.plot=T)
