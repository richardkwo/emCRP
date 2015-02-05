library(MASS)
source("myutils.R")

# fast EM-style inference for DP mixture of Gaussians
dp.EM.infer <- function(X, alpha=NULL, prior.mean=NULL, prior.sd=NULL, sd=NULL, 
                        posterior.predictive=FALSE, cluster.size.threshold=1, 
                        max.iter=100, max.inner.iter=10, max.stable.iter=5, inner.tol=1e-5, 
                        verbose=0) {
    N <- X$N
    if (is.null(alpha)) {
        alpha <- X$alpha
    }
    if (is.null(prior.mean)) {
        prior.mean <- X$prior.mean
    }
    if (is.null(prior.sd)) {
        prior.sd <- X$prior.sd
    }
    if (is.null(sd)) {
        sd <- X$sd
    }
    p <- length(prior.mean)
    # init
    W <- cbind(rep(0, N)) # of shape N x (K+1)
    # posterior of theta_k ~ N(cluster.means[k], cluster.sds[k] x I)
    cluster.means <- NULL # (K+1) x p, the MLE estimate of the mean of clusters
    cluster.sds <- NULL # (K+1) x 1, the posterior sd of the mean of cluster
    K <- 0
    iter <- 0
    # history loggers
    W.history = list()
    # EM
    while (iter<=max.iter) {
        iter <- iter + 1
        # inner loop to get stable estimates before creating new clusters
        change <- 1e4
        inner.iter <- 0
        while (change>inner.tol & inner.iter<max.inner.iter) {
            inner.iter <- inner.iter + 1
            # E step
            W.old <- W
            for (i in 1:N) {
                prob.vec <- rep(0, K+1)
                # prob of *
                prob.vec[K+1] <-  alpha * gaussian.prob.kernel(x=(X$x[i,]-prior.mean), 
                                                               s2=sd^2 + prior.sd^2)
                # prob of clusters
                if (K>0) {
                    for (k in 1:K) {
                        if (posterior.predictive) {
                            # integrated likelihood against posterior of the cluster
                            prob.vec[k] <- sum(W.old[,k]) * gaussian.prob.kernel(x=(X$x[i,]-cluster.means[k, ]), 
                                                                                 s2=sd^2 + cluster.sds[k]^2)
                        } else {
                            # likelihood under point estimate
                            prob.vec[k] <- sum(W.old[,k]) * gaussian.prob.kernel(x=(X$x[i,]-cluster.means[k,]), 
                                                                                 s2=sd^2)
                        }                    
                    }
                }
                prob.vec <- prob.vec / sum(prob.vec)
                W[i,] <- prob.vec
            }
            # M-step
            if (K>0) {
                # re-estimate for the mean of each cluster
                cluster.means.old <- cluster.means
                for (k in 1:K) {
                    H <- get.gaussian.posterior(prior.mean = prior.mean, 
                                                prior.sd = prior.sd, 
                                                sd = sd, 
                                                x = X$x, 
                                                weights = W[, k])
                    cluster.means[k, ] <- H$mean
                    cluster.sds[k] <- H$sd 
                }
                change <- sqrt(mean((cluster.means-cluster.means.old)^2) * p)
            } else {
                change <- 0
            }
        } # inner loop finish
        cluster.sizes <- apply(W, 2, sum)
        new.cluster.size <- cluster.sizes[K+1]
        if (verbose>0) {
            cat(sprintf("Iter %d (with %d inner iters): K = %d\n", iter, inner.iter, K))
            if (K>0) {
                display.df <- data.frame(cluster.means)
                display.df <- rbind(display.df, rep(NA, p))
                display.df$size <- cluster.sizes
                print(display.df)
            }
        }
        # add/remove cluster
        # add a new cluster
        if (new.cluster.size>cluster.size.threshold) {
            iter.since.stable <- 0
            # spin off a new cluster
            H.new <-  get.gaussian.posterior(prior.mean = prior.mean, 
                                             prior.sd = prior.sd, 
                                             sd = sd, 
                                             x = X$x, 
                                             weights = W[, K+1])
            # init mean with a random draw from posterior
            stopifnot(nrow(cluster.means)==K)
            stopifnot(length(cluster.sds)==K)
            cluster.means <- rbind(cluster.means, mvrnorm(n=1, 
                                                          mu = H.new$mean, 
                                                          Sigma = H.new$sd^2 * diag(p)))
            cluster.sds <- c(cluster.sds, H.new$sd)
            K <- K+1
            W <- cbind(W, rep(0, N))
            if (verbose>0) {
                cat("new cluster created at", cluster.means[K, ], "\n")
            }
        } else {
            iter.since.stable <- iter.since.stable + 1
        }
        # remove a cluster
        if (iter.since.stable>max.stable.iter) {
            k <- which.min(cluster.size[1:K])
            if (cluster.size[k] < cluster.size.threshold & K>1) {
                if (verbose>0) {
                    
                }
                # delete it and put weights back to *
                cluster.means <- cluster.means[-k, ]
                
            }
        }
    }
}
