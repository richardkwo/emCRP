library(MASS)
library(mvtnorm)
library(foreach)
library(doMC)
registerDoMC()
source("myutils.R")

# fast EM-style inference for DP mixture of Gaussians
dp.EM.infer <- function(X, alpha=NULL, prior.mean=NULL, prior.std=NULL, std=NULL, 
                        posterior.predictive=TRUE, cluster.size.threshold=2, 
                        max.iter=100, max.inner.iter=100, max.stable.iter=5, 
                        inner.tol=1e-2, tol=1e-3, 
                        verbose=1, stepwise.plot=FALSE) {
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
    iter.since.stable <- 0
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
        # stepwise visualization
        if (stepwise.plot & K>0 & iter.since.stable<3) {
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
            plot.df$cluster.size <- c.size
            if (!exists("data.df")) {
                data.df <- data.frame(X$x)
            }
            fig <- ggplot(plot.df, aes(x=X1, y=X2)) + 
                geom_density2d(aes(color=z)) + 
                geom_point(data=data.df, alpha=.5) 
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

# fast EM-style inference for DP mixture of Gaussians
dp.EM.infer.multivariate <- function(X, alpha=NULL, 
                                     prior.mean=NULL, prior.scale.matrix=NULL, 
                                     prior.df=NULL, prior.precision.factor=NULL,
                                     posterior.predictive=TRUE, random.init.cluster=TRUE, 
                                     cluster.size.threshold=3, 
                                     max.iter=100, max.inner.iter=100, max.stable.iter=5, 
                                     inner.tol=1e-2, tol=1e-3, 
                                     verbose=1, stepwise.plot=FALSE) {
    N <- X$N
    p <- length(X$x[1, ])
    if (is.null(alpha)) {
        alpha <- X$alpha
        if (verbose>0) {
            cat("alpha = ", alpha, "\n")
        }
    }

    use.reference.prior <- FALSE
    if (is.null(prior.mean) & is.null(prior.scale.matrix) & is.null(prior.df) & is.null(prior.precision.factor)) {
        use.reference.prior <- TRUE
        prior.mean <- rep(0, p)
        prior.scale.matrix <- matrix(rep(0, p*p), nrow=p)
        prior.df <- -1
        prior.precision.factor <- 0
        cat("Using reference prior\n")
    }

    # init
    W <- cbind(rep(0, N)) # of shape N x (K+1)
    
    # posterior of mu_k | cluster.scale.matrices[k] ~ 
    # N(cluster.means[k], (cluster.precision.factors[k] * Lambda_k)^{-1})
    cluster.means <- NULL # (K+1) x p, the MAP estimate of the mean of clusters
    cluster.precision.factors <- NULL # (K+1), kappa for each cluster
    # posterior of Prec ~ Wishart((cluster.scale.matrices[k])^{-1}, cluster.dfs[k])
    cluster.scale.matrices <- list() # each is a p x p matrix, prior mean of Prec = (cluste.scale.matrix)^{-1} x cluster.df
    cluster.dfs <- NULL # (K+1), degree of freedom
    cluster.cov.MAPs <- list() # each p x p, the MAP estimate for the cov of each cluster
    
    K <- 0
    iter <- 0
    iter.since.stable <- 0
    # history loggers
    W.history = list()

    # EM
    v <- 1:p
    s1 <- lgamma((prior.df+2-v)/2)
    s2 <- lgamma((prior.df+1-v)/2)
    log.prob.fixed.part <- -p/2 * log(pi) + p/2 * (log(prior.precision.factor) - log(prior.precision.factor+1)) + sum(s1-s2)
    while (iter<=max.iter) {
        iter <- iter + 1
        # inner loop to get stable estimates before creating new clusters
        change <- 1e4
        inner.iter <- 0
        while (change>inner.tol & inner.iter<max.inner.iter) {
            inner.iter <- inner.iter + 1
            # E step - estimate sample weights to each cluster
            old.cluster.sizes <- colSums(W)
            W <- foreach (i = 1:N, .combine=rbind) %dopar% { 
                prob.vec <- rep(0, K+1)
                # prob of * (denoted by K+1)
                # alpha x marginal likelihood of X integrated against prior
                ## TODO: what if under ref prior ??
                one.sample.scale.mat <- prior.scale.matrix + prior.precision.factor/(prior.precision.factor+1) * (prior.mean-X$x[i,]) %*% t(prior.mean-X$x[i,])
                log.prob.asterisk <- log.prob.fixed.part + prior.df/2*log(det(prior.scale.matrix)) - (prior.df+1)/2*log(det(one.sample.scale.mat))
                prob.vec[K+1] <-  alpha * exp(log.prob.asterisk)
                stopifnot(prob.vec[K+1]>=0)
                # prob of clusters
                if (K>0) {
                    for (k in 1:K) {
                        if (posterior.predictive) {
                            # TODO: how to compute the marginal likelihood in the reference prior case ??
                            # evaluated with posterior predictive of each cluster -- multivariate t
                            # x (est. of cluster size)
                            t.df <- cluster.dfs[k] - p + 1
                            t.scale.mat <- (cluster.precision.factors[k]+1)/(cluster.precision.factors[k] * t.df) * cluster.scale.matrices[[k]]
                            prob.vec[k] <- old.cluster.sizes[k] * dmvt(X$x[i,], delta=cluster.means[k,], sigma=t.scale.mat, df=t.df, log=FALSE)
                        } else {
                            # likelihood under point estimate
                            prob.vec[k] <- old.cluster.sizes[k] * dmvnorm(X$x[i,], mean=cluster.means[k,], sigma=cluster.cov.MAPs[[k]], log=FALSE)
                        }
                        stopifnot(prob.vec[k]>=0)                    
                    }
                }
                prob.vec <- prob.vec / sum(prob.vec)
                prob.vec # the return value to be combined by "foreach"
                
            } # foreach ends
            print(head(W))
            print(cluster.cov.MAPs)
            # M-step - re-estimate the params of each cluster
            if (K>0) {
                cluster.means.old <- cluster.means
                for (k in 1:K) {
                    # sufficient stats
                    x.sum <- colSums(W[,k] * X$x)
                    nc <- sum(W[,k])
                    x.bar <- x.sum / nc
                    scatter <- t(X$x - x.bar) %*% diag(W[,k]) %*% (X$x - x.bar)  # p x p scatter matrix
                    # hyper-param updates
                    cluster.means[k,] <- (prior.precision.factor*prior.mean + x.sum)/(prior.precision.factor + nc)
                    cluster.precision.factors[k] <- prior.precision.factor + nc
                    cluster.dfs[k] <- prior.df + nc
                    mean.difference.term <- prior.precision.factor * nc / (prior.precision.factor + nc) * (prior.mean - x.bar) %*% t((prior.mean - x.bar))
                    cluster.scale.matrices[[k]] <- prior.scale.matrix + scatter + mean.difference.term
                    if (use.reference.prior) {
                        if (nc-d-1>0) {
                            cluster.cov.MAPs[[k]] <- S / (nc-d-1)
                        } else {
                            cluster.cov.MAPs[[k]] <- S / nc # rare!
                        }        
                    } else {
                        # TODO: double-check the formula
                        cluster.cov.MAPs[[k]] <- cluster.scale.matrices[[k]] / (cluster.dfs[k] + p + 1)
                    }
                }
                change <- sqrt(mean((cluster.means[1:K,]-cluster.means.old[1:K,])^2) * p)
            } else {
                change <- 0
            }
        } # inner EM loop finish
        # print(cluster.means)
        # print(cluster.scale.matrices)
        # print(cluster.cov.MAPs)
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
        
        if (stepwise.plot & K>0 & iter.since.stable<2) {
            plot.df <- NULL
            z <- NULL
            for (k in 1:(K+1)) {
                for (l in 1:1000) {
                    # density estimate from draws ~ posterior predictive
                    if (k==K+1) {
                        R.draw <- rWishart(1, prior.df, solve(prior.scale.matrix))[,,1]
                        scale.mat.draw <- solve(R.draw)
                        mu.draw <- rmvnorm(1, prior.mean, 1/prior.precision.factor * scale.mat.draw)
                    } else {
                        R.draw <- rWishart(1, cluster.dfs[k], solve(cluster.scale.matrices[[k]]))[,,1]
                        scale.mat.draw <- solve(R.draw)
                        mu.draw <- rmvnorm(1, cluster.means[k,], 1/cluster.precision.factors[k] * scale.mat.draw)
                    }
                    plot.df <- rbind(plot.df, rmvnorm(1, mu.draw, scale.mat.draw))
                    if (k==K+1) {
                        z <- c(z, "prior")
                    } else {
                        z <- c(z, paste(k))
                    }
                }
            }
            colnames(plot.df) <- NULL
            plot.df <- data.frame(plot.df)
            plot.df$z <- factor(z)
            if (!exists("data.df")) {
                data.df <- data.frame(X$x)
            }
            fig <- ggplot(plot.df, aes(x=X1, y=X2)) + 
                geom_density2d(aes(color=z)) + 
                geom_point(data=data.df, alpha=.5) 
            print(fig)
        }
        
        # add/remove cluster
        iter.since.stable <- iter.since.stable + 1
        # add a new cluster?
        if (new.cluster.size>cluster.size.threshold) {
            iter.since.stable <- 0
            # spin off a new cluster
            stopifnot(nrow(cluster.means)==K) # TODO: remove those asserts
            stopifnot(length(cluster.precision.factors)==K)
            stopifnot(length(cluster.scale.matrices)==K)
            stopifnot(length(cluster.dfs)==K)
            stopifnot(length(cluster.cov.MAPs)==K)
            # init mean with a random draw from posterior
            ## sufficient stats
            x.sum <- colSums(W[,K+1] * X$x)
            nc <- sum(W[,K+1])
            x.bar <- x.sum / nc
            scatter <- t(X$x - x.bar) %*% diag(W[,K+1]) %*% (X$x - x.bar)  # p x p scatter matrix
            ## hyper-param updates
            mu.hat <- (prior.precision.factor*prior.mean + x.sum)/(prior.precision.factor + nc)
            r.hat <- prior.precision.factor + nc
            nu.hat <- prior.df + nc
            diff.term <- prior.precision.factor * nc / (prior.precision.factor + nc) * (prior.mean - x.bar) %*% t((prior.mean - x.bar))
            scale.hat <- prior.scale.matrix + scatter + diff.term         
            ## draw
            if (random.init.cluster) {
                # randomized hyper-params for the new cluster 
                R.draw <- rWishart(1, nu.hat, solve(scale.hat))[,,1]
                cluster.scale.matrices[[K+1]] <- solve(R.draw)
                cluster.means <- rbind(cluster.means, rmvnorm(1, mu.hat, 1/r.hat * cluster.scale.matrices[[K+1]]))
            } else {
                # posterior estimated from w(*)
                cluster.scale.matrices[[K+1]] <- scale.hat
                cluster.means <- rbind(cluster.means, mu.hat)
            }
            cluster.dfs[K+1] <- nu.hat
            cluster.precision.factors[K+1] <- r.hat
            cluster.cov.MAPs[[K+1]] <- scale.hat / (nu.hat + p + 1) # TODO: double-check & fix replicated code
            # increment
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
            # so it won't be removed immediately
            if (verbose>0) {
                cat(sprintf("- cluster %d (with size %f) removed \n", k, cluster.sizes[k]))
            }
            # delete it and put weights back to *
            cluster.means <- cluster.means[-k, ]
            cluster.scale.matrices <- cluster.scale.matrices[-k]
            cluster.dfs <- cluster.dfs[-k]
            cluster.precision.factors <- cluster.precision.factors[-k]
            cluster.cov.MAPs <- cluster.cov.MAPs[-k]
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
    params <- list(alpha=NULL, 
        prior.mean=prior.mean, prior.scale.matrix=prior.scale.matrix, 
        prior.df=prior.df, prior.precision.factor=prior.precision.factor,
        posterior.predictive=posterior.predictive, random.init.cluster=random.init.cluster, 
        cluster.size.threshold=cluster.size.threshold, 
        max.iter=max.iter, max.inner.iter=max.inner.iter, max.stable.iter=max.stable.iter, 
        inner.tol=inner.tol, tol=tol)
    
    # rerank the output by cluster size
    idx.new.to.old <- order(cluster.sizes[1:K], decreasing = T)
    cluster.sizes[1:K] <- cluster.sizes[idx.new.to.old]
    W[, (1:K)] <- W[, idx.new.to.old]
    cluster.means <- cluster.means[idx.new.to.old, ]
    cluster.precision.factors <- cluster.precision.factors[idx.new.to.old]
    cluster.scale.matrices <- cluster.scale.matrices[idx.new.to.old]
    cluster.dfs <- cluster.dfs[idx.new.to.old]
    cluster.cov.MAPs <- cluster.cov.MAPs[idx.new.to.old]
    display.df <- data.frame(cluster.means)
    display.df <- rbind(display.df, rep(NA, p))
    display.df$size <- cluster.sizes
    
    cat("\nResult \n")
    print(display.df)
    
    return(list(params = params, 
                K = K, N=N, p=p, 
                cluster.means = cluster.means, 
                cluster.precision.factors = cluster.precision.factors,
                cluster.scale.matrices = cluster.scale.matrices,
                cluster.dfs = cluster.dfs, 
                cluster.cov.MAPs = cluster.cov.MAPs,
                cluster.sizes = cluster.sizes[1:K], 
                W = W, 
                display.df = display.df))
}


analyze.dpem.result <- function(dpem.result, X) {
    K <- X$num.clusters
    K.est <- dpem.result$K
    cluster.means.true.df <- data.frame(X$cluster.means)
    cluster.means.true.df$z <- factor(1:K)
    cluster.means.true.df$size <- X$cluster.sizes
    est.df <- dpem.result$display.df[1:K.est, ]
    est.df$std <- dpem.result$cluster.sds
    est.df$z <- factor(1:K.est)
    est.df$size <- dpem.result$cluster.sizes
    if (!is.null(X)) {
        # estimated clustering
        plot.df <- data.frame(X1=X$x[,1], X2=X$x[,2], z=factor(X$z))
        plot.df$z.MAP.est <- factor(apply(dpem.result$W, 1, which.max))
        plot.df$z.confidence <- apply(dpem.result$W, 1, max)
        fig.1 <- ggplot(plot.df, aes(x=X1, y=X2)) + 
            geom_point(aes(color=z.MAP.est, alpha=z.confidence)) + 
            geom_point(aes(x=X1, y=X2, color=z), size=8, shape=8, data=est.df) + 
            geom_point(aes(x=X1, y=X2, shape="true", size=size), data=cluster.means.true.df) + 
            scale_size_continuous(trans="log10")
        print(fig.1)
        # compare centroid estimate to true means
        fig.2 <- ggplot(est.df, aes(x=X1, y=X2, color=z)) + 
            geom_point(aes(shape="est", size=size)) + 
            geom_errorbar(aes(ymin=X2-std, ymax=X2+std)) + 
            geom_errorbarh(aes(xmin=X1-std, xmax=X1+std)) + 
            geom_point(aes(x=X1, y=X2, color=z, shape="true", size=size), data=cluster.means.true.df) +
            scale_size_continuous(trans="log10")
        print(fig.2)
    }
}