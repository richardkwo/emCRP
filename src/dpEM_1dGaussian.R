library(MASS)
library(mvtnorm)
library(foreach)
library(doMC)
registerDoMC()
source("myutils.R")

# fast EM-style inference for DP mixture of 1-D Gaussians
# estimate parameters & density
dp.EM.infer <- function(X, alpha=NULL, prior.mean=NULL, prior.kappa=NULL, 
                        prior.a=NULL, prior.b=NULL, 
                        estimate.type="MLE",                  # c("MLE","MAP","Bayes")
                        random.init.cluster=FALSE,
                        cluster.size.threshold=2, 
                        max.iter=100, max.inner.iter=100, max.stable.iter=5, 
                        inner.tol=1e-2, tol=1e-3, 
                        verbose=1, stepwise.plot=FALSE) {
    
    # prior: mean | precision ~ N(prior.mean, (prior.kappa * precision)^-1)
    #        precision ~ Gamma(prior.a, rate=prior.b)
    
    N <- X$N
    if (is.null(alpha)) {
        alpha <- X$alpha
    }
    if (is.null(prior.mean)) {
        prior.mean <- X$prior.mean
    }
    if (is.null(prior.kappa)) {
        prior.kappa <- X$prior.kappa
    }
    if (is.null(prior.a)) {
        prior.a <- X$prior.a
    }
    if (is.null(prior.b)) {
        prior.b <- X$prior.b
    }
    stopifnot(length(prior.mean)==1)
    stopifnot(estimate.type="MLE" | estimate.type="MAP" | estimate.type=="Bayes")
    
    # init
    W <- cbind(rep(1, N)) # of shape N x (K+1), the last column corresponds to *
    # posterior of theta_k ~ N(cluster.means[k], cluster.sds[k] x I)
    cluster.x.bar <- NULL # (K+1), the MLE estimate of the mean of clusters (\bar{x}_c)
    cluster.var.bar <- NULL # (K+1), the MLE estimate of the var of clusters (\bar{v}_c)
    cluster.sizes <- NULL # (K+1), cluster sizes (n_c)
    cluster.m <- NULL     # (K+1), updated m for each cluster, also MAP estimate of mean
    cluster.kappa <- NULL # (K+1), updated kappa
    cluster.a <- NULL     # (K+1), updated a
    cluster.b <- NULL     # (K+1), updated b
    K <- 0
    iter <- 0
    iter.since.stable <- 0
    # history loggers
    history.W = list()
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
                # prob of *, a t distribution
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
        # stepwise visualization
        if (stepwise.plot & K>0 & iter.since.stable<2) {
            plot.df <- NULL
            z <- NULL
            draw.size <- 1000
            for (k in 1:K) {
                for (l in 1:draw.size){
                    # posterior predictive draw 
                    mu <- rmvnorm(n=1, cluster.means[k,], cluster.sds[k]^2 * diag(p))
                    sample.x <- rmvnorm(n=1, mu, std^2 * diag(p))
                    plot.df <- rbind(plot.df, sample.x)
                    z <- c(z, k)
                }
            }
            plot.df <- data.frame(plot.df)
            plot.df$z <- factor(z)
            if (!exists("data.df")) {
                data.df <- data.frame(X$x)
            }
            data.df$w.unexplained <- W[,K+1]
            fig <- ggplot(plot.df, aes(x=X1, y=X2)) + 
                stat_ellipse(aes(color=z)) + 
                geom_point(aes(alpha=w.unexplained), data=data.df) + 
                ggtitle(paste("Iter", iter, ", alpha=", alpha, ", gamma=", cluster.size.threshold))
            print(fig)
        }
        # add/remove cluster
        iter.since.stable <- iter.since.stable + 1
        # add a new cluster?
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
        }
        # remove a cluster ?
        k <- which.min(cluster.sizes[1:K]) 
        if (cluster.sizes[k] < cluster.size.threshold & K>1) {
            # if new cluster created, its size must be > the threshold, 
            # so it won't be removed
            if (verbose>0) {
                cat(sprintf("- cluster %d (with size %f) removed \n", k, cluster.sizes[k]))
            }
            # delete it and put weights back to *
            cluster.means <- cluster.means[-k, ]
            cluster.sds <- cluster.sds[-k]
            W[, K+1] <- W[, K+1] + W[, k]
            W <- W[, -k]
            K <- K - 1
            iter.since.stable <- 0
        }
        if (iter.since.stable>=max.stable.iter & change<tol) {
            break
        }
        
    }
    if (iter>max.iter) {
        warning(sprintf("Algorithm terminated when reaching max.iter=%d", max.iter))
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
