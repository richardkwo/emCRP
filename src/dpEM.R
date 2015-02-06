library(MASS)
library(mvtnorm)
library(foreach)
library(doMC)
registerDoMC()
source("myutils.R")

# fast EM-style inference for DP mixture of Gaussians
dp.EM.infer <- function(X, alpha=NULL, prior.mean=NULL, prior.std=NULL, std=NULL, 
                        posterior.predictive=FALSE, cluster.size.threshold=1, 
                        max.iter=100, max.inner.iter=100, max.stable.iter=10, inner.tol=1e-3, 
                        verbose=0) {
    N <- X$N
    if (is.null(alpha)) {
        alpha <- X$alpha
    }
    if (is.null(prior.mean)) {
        prior.mean <- X$prior.mean
    }
    if (is.null(prior.std)) {
        prior.std <- X$prior.std
    }
    if (is.null(std)) {
        std <- X$std
    }
    p <- length(prior.mean)
    # init
    W <- cbind(rep(0, N)) # of shape N x (K+1)
    # posterior of theta_k ~ N(cluster.means[k], cluster.sds[k] x I)
    cluster.means <- NULL # (K+1) x p, the MLE estimate of the mean of clusters
    cluster.sds <- NULL # (K+1) x 1, the posterior std of the mean of cluster
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
            old.cluster.sizes <- colSums(W)
            W <- foreach (i = 1:N, .combine=rbind) %dopar% { 
                prob.vec <- rep(0, K+1)
                # prob of *
                prob.vec[K+1] <-  alpha * dmvnorm(X$x[i,], prior.mean, (std^2+prior.std^2)*diag(p))
                # prob of clusters
                if (K>0) {
                    for (k in 1:K) {
                        if (posterior.predictive) {
                            # integrated likelihood against posterior of the cluster
                            prob.vec[k] <- old.cluster.sizes[k] * dmvnorm(X$x[i,], cluster.means[k,], (std^2 + cluster.sds[k]^2)*diag(p))
                        } else {
                            # likelihood under point estimate
                            prob.vec[k] <- old.cluster.sizes[k] * dmvnorm(X$x[i,], cluster.means[k,], (std^2)*diag(p))
                        }                    
                    }
                }
                prob.vec <- prob.vec / sum(prob.vec)
                prob.vec
            }
            # M-step
            if (K>0) {
                # re-estimate for the mean of each cluster
                cluster.means.old <- cluster.means
                for (k in 1:K) {
                    H <- get.gaussian.posterior(prior.mean = prior.mean, 
                                                prior.std = prior.std, 
                                                std = std, 
                                                x = X$x, 
                                                weights = W[, k])
                    cluster.means[k, ] <- H$mean
                    cluster.sds[k] <- H$std 
                }
                change <- sqrt(mean((cluster.means-cluster.means.old)^2) * p)
            } else {
                change <- 0
            }
        } # inner loop finish
        cluster.sizes <- apply(W, 2, sum)
        new.cluster.size <- cluster.sizes[K+1]
        if (verbose>0) {
            cat(sprintf("\nIter %d (with %d inner iters): K = %d\n", iter, inner.iter, K))
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
                                             prior.std = prior.std, 
                                             std = std, 
                                             x = X$x, 
                                             weights = W[, K+1])
            # init mean with a random draw from posterior
            stopifnot(nrow(cluster.means)==K)
            stopifnot(length(cluster.sds)==K)
            cluster.means <- rbind(cluster.means, mvrnorm(n=1, 
                                                          mu = H.new$mean, 
                                                          Sigma = H.new$std^2 * diag(p)))
            cluster.sds <- c(cluster.sds, H.new$std)
            K <- K+1
            W <- cbind(W, rep(0, N))
            if (verbose>0) {
                cat("* new cluster created at", cluster.means[K, ], "\n")
            }
        } else {
            iter.since.stable <- iter.since.stable + 1
        }
        # remove a cluster
        if (iter.since.stable>=0) {
            k <- which.min(cluster.sizes[1:K])
            if (cluster.sizes[k] < cluster.size.threshold & K>1) {
                if (verbose>0) {
                    cat(sprintf("cluster %d (with size %f) removed \n", k, cluster.sizes[k]))
                }
                # delete it and put weights back to *
                cluster.means <- cluster.means[-k, ]
                cluster.sds <- cluster.sds[-k]
                W[, K+1] <- W[, K+1] + W[, k]
                W <- W[, -k]
                K <- K - 1
                iter.since.stable <- 0
            } else {
                if (iter.since.stable>max.stable.iter) {
                    break
                }
            }
        }
    }
    params <- list(alpha=alpha, prior.mean=prior.mean, prior.std=prior.std, std=std, 
                   posterior.predictive=posterior.predictive, cluster.size.threshold=cluster.size.threshold, 
                   max.iter=max.iter, max.inner.iter=max.inner.iter, max.stable.iter=max.stable.iter, 
                   inner.tol=inner.tol)
    # rerank the output by cluster size
    idx.new.to.old <- order(cluster.sizes[1:K], decreasing = T)
    cluster.sizes[1:K] <- cluster.sizes[idx.new.to.old]
    W[, (1:K)] <- W[, idx.new.to.old]
    cluster.means <- cluster.means[idx.new.to.old, ]
    cluster.sds <- cluster.sds[idx.new.to.old]
    display.df <- data.frame(cluster.means)
    display.df <- rbind(display.df, rep(NA, p))
    display.df$size <- cluster.sizes
    
    cat("\nResult \n")
    print(display.df)
    
    return(list(params = params, 
                K = K, N=N, p=p, 
                cluster.means = cluster.means, 
                cluster.sds = cluster.sds, 
                cluster.sizes = cluster.sizes[1:K], 
                W = W, 
                display.df = display.df))
}

analyze.dpem.result <- function(dpem.result, X) {
    K <- X$num.clusters
    K.est <- dpem.result$K
    cluster.means.true.df <- data.frame(X$cluster.means)
    cluster.means.true.df$z <- factor(1:K)
    est.df <- dpem.result$display.df[1:K.est, ]
    est.df$std <- dpem.result$cluster.sds
    est.df$z <- factor(1:K.est)
    if (!is.null(X)) {
        # estimated clustering
        plot.df <- data.frame(X1=X$x[,1], X2=X$x[,2], z=factor(X$z))
        plot.df$z.MAP.est <- factor(apply(dpem.result$W, 1, which.max))
        plot.df$z.confidence <- apply(dpem.result$W, 1, max)
        fig.1 <- ggplot(plot.df, aes(x=X1, y=X2)) + 
            geom_point(aes(color=z.MAP.est, alpha=z.confidence)) + 
            geom_point(aes(x=X1, y=X2, color=z), size=8, shape=8, data=est.df) + 
            geom_point(aes(x=X1, y=X2, shape="true"), size=3, data=cluster.means.true.df)
        print(fig.1)
        # compare centroid estimate to true means
        fig.2 <- ggplot(est.df, aes(x=X1, y=X2, color=z)) + 
            geom_point(aes(shape="est")) + 
            geom_errorbar(aes(ymin=X2-std, ymax=X2+std)) + 
            geom_errorbarh(aes(xmin=X1-std, xmax=X1+std)) + 
            geom_point(aes(x=X1, y=X2, color=z, shape="true"), size=3, data=cluster.means.true.df) 
        print(fig.2)
    }
}