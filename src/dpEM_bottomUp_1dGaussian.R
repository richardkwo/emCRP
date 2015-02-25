library(MASS)
library(foreach)
library(doMC)
registerDoMC()
source("myutils.R")

# fast EM-style inference for DP mixture of 1-D Gaussians
# estimate parameters & density
bottom.up.dp.EM.1dGaussian <- function(X, alpha=NULL, 
                                    prior.mean=NULL,   # prior.mean: m_0
                                    prior.kappa=NULL,  # prior.kappa: # of pseudo-points on mean
                                    prior.a=NULL, prior.b=NULL, 
                                    estimate.type="MLE",  # c("MLE","MAP","Bayes")
                                    init.K=30,
                                    cluster.size.threshold=2, 
                                    max.iter=100, max.inner.iter=500, max.stable.iter=20, 
                                    inner.tol=1e-2, tol=1e-2, W.tol=0.05,
                                    verbose=1, stepwise.plot=FALSE) {
    
    # prior: mean | precision ~ N(prior.mean, (prior.kappa * precision)^-1)
    #        precision ~ Gamma(prior.a, rate=prior.b)
    
    N <- length(X$x)
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
    stopifnot(estimate.type=="MLE" | estimate.type=="MAP" | estimate.type=="Bayes")
    # t distribution for the *
    ast.mean <- prior.mean
    ast.sigma <- sqrt(prior.b*(prior.kappa+1)/(prior.a * prior.kappa))
    ast.df <- 2 * prior.a
    # random init with K samples and nearest assignments
    K <- init.K
    W <- matrix(rep(0,(K+1)*N), nrow=N)
    cluster.x.bar <- sample(X$x, size=K)
    for (i in 1:N) {
        nearest.k <- which.min(abs(cluster.x.bar-X$x[i]))
        W[i, nearest.k] <- 0.9
        W[i, K+1] <- 0.1 # some init value put on *
    }
    cat(sprintf("W initialized with %d random samples with nearest assignment.\n", K))
    # init parameters
    cluster.x.bar <- NULL # (K+1), the MLE estimate of the mean of clusters (\bar{x}_c)
    cluster.var.bar <- NULL # (K+1), the MLE estimate of the var of clusters (\bar{v}_c)
    cluster.sizes <- NULL # (K+1), cluster sizes (n_c)
    cluster.m <- NULL     # (K+1), updated m for each cluster
    cluster.mean.map <- NULL # (K+1), MAP estimate for mean of each cluster
    cluster.var.map <- NULL  # (K+1), MAP estimate for var of each cluster
    cluster.kappa <- NULL # (K+1), updated kappa
    cluster.a <- NULL     # (K+1), updated a
    cluster.b <- NULL     # (K+1), updated b
    iter <- 0
    iter.since.stable <- 0
    to.exit = FALSE
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
            # M-step
            # re-estimate the parameter for each cluster
            cluster.sizes <- colSums(W)
            cluster.x.bar.old <- cluster.x.bar
            cluster.var.bar.old <- cluster.var.bar
            for (k in 1:(K+1)) { # parameter of * also estimated
                nc <- cluster.sizes[k]
                cluster.x.bar[k] <- sum(W[,k] * X$x) / nc
                cluster.var.bar[k] <- sum(W[,k] * (X$x - cluster.x.bar[k])^2) / nc
                cluster.m[k] <- (cluster.x.bar[k] + prior.kappa * prior.mean / nc) / (1 + prior.kappa / nc)
                cluster.a[k] <- prior.a + nc / 2
                cluster.kappa[k] <- prior.kappa + nc
                cluster.b[k] <- prior.b + nc * cluster.var.bar[k] / 2 + 
                    prior.kappa * nc * (cluster.x.bar[k] - prior.mean)^2 / (2 * (prior.kappa + nc))
                cluster.mean.map[k] <- cluster.m[k]
                cluster.var.map[k] <- cluster.b[k] / (cluster.a[k]) #TODO: fix me
            }
            ## track the change
            if (!is.null(cluster.x.bar.old)) {
                if (length(cluster.x.bar.old)==K+1) {
                    change <- max(sqrt(max((cluster.x.bar-cluster.x.bar.old)^2)), 
                                  sqrt(max((cluster.var.bar.old-cluster.var.bar)^2)))  
                } else {
                    change <- max(sqrt(max((cluster.x.bar[1:K]-cluster.x.bar.old)^2)), 
                                  sqrt(max((cluster.var.bar[1:K]-cluster.var.bar.old)^2))) 
                }
            } else {
                change <- 1e4
            }

            # E step
            # multiple iterations to solve the fixed point
            W.iter = 0
            while (K>0) {
                W.iter <- W.iter + 1
                cluster.sizes <- colSums(W)
                #                 print(cluster.sizes)
                W.old <- W
                W <- foreach (i = 1:N, .combine=rbind) %dopar% { 
                    prob.vec <- rep(0, K+1)
                    # prob of *, a non-standardized t distribution
                    # prior factor: (alpha + size of * excluding i)
                    prob.vec[K+1] <-  (alpha + cluster.sizes[K+1] - W.old[i,K+1])* dt((X$x[i] - ast.mean)/ast.sigma, df=ast.df) / ast.sigma
                    # prob of clusters
                    # prior factor: (size of k excluding i)
                    if (K>0) {
                        for (k in 1:K) {
                            if (estimate.type=="MLE") {
                                # likelihood under MLE
                                prob.vec[k] <- (cluster.sizes[k] - W.old[i,k]) * 
                                    dnorm(X$x[i], mean=cluster.x.bar[k], sd=sqrt(cluster.var.bar[k]))
                            } else if (estimate.type=="MAP") {
                                # likelihood under MAP
                                prob.vec[k] <- (cluster.sizes[k] - W.old[i,k]) * 
                                    dnorm(X$x[i], mean=cluster.mean.map[k], sd=sqrt(cluster.var.map[k]))
                            } else {
                                # Bayesian posterior predictive 
                                t.sigma <- sqrt(cluster.b[k]*(cluster.kappa[k]+1)/(cluster.a[k]*cluster.kappa[k]))
                                prob.vec[k] <- (cluster.sizes[k] - W.old[i,k]) * 
                                    dt((X$x[i]-cluster.m[k])/t.sigma, df=2*cluster.a[k]) / t.sigma
                                # don't forget Jacobian!
                            }                  
                        }
                    }
                    prob.vec <- prob.vec / sum(prob.vec)
                    prob.vec # returning value of dopar
                } # %dopar% finished
                cluster.sizes.new <- colSums(W)
                if (max(abs(cluster.sizes.new-cluster.sizes)) < W.tol) {
                    if (verbose>1) {
                        cat(sprintf("W computed with %d iters\n", W.iter))
                    }
                    break
                }
            } # E-step finished
            
            # display
            if (verbose>1) {
                if (K==0) {
                    clusters <- c("*")
                } else {
                    clusters <- c(1:K, "*")
                }
                display.df <- data.frame(clusters=clusters, 
                                         cluster.mean=cluster.x.bar, 
                                         cluster.var=cluster.var.bar, 
                                         size=cluster.sizes)
                cat(sprintf("--Inner iter %d--\n", inner.iter))
                print(display.df)
            }
            
            # estimate density
            density.est.mle <- function(y) {
                d <- rep(0, length(y))
                if (K>0) {
                    for (k in 1:K) {
                        d <- d + cluster.sizes[k]/(N+alpha) * dnorm(x=y, mean=cluster.x.bar[k], sd=sqrt(cluster.var.bar[k]))
                    }
                }
                d <- d + (cluster.sizes[K+1]+alpha)/(N+alpha) * dt(x=(y-ast.mean)/ast.sigma, df=ast.df) / ast.sigma
                return(d)
            }
            density.est.map <- function(y) {
                d <- rep(0, length(y))
                if (K>0) {
                    for (k in 1:K) {
                        d <- d + cluster.sizes[k]/(N+alpha) * dnorm(x=y, mean=cluster.mean.map[k], sd=sqrt(cluster.var.map[k]))
                    }
                }
                d <- d + (cluster.sizes[K+1]+alpha)/(N+alpha) * dt(x=(y-ast.mean)/ast.sigma, df=ast.df) / ast.sigma
                return (d)
            }
            density.est.Bayes <- function(y) {
                d <- rep(0, length(y))
                if (K>0) {
                    for (k in 1:K) {
                        t.sigma <- sqrt(cluster.b[k]*(cluster.kappa[k]+1)/(cluster.a[k]*cluster.kappa[k]))
                        d <- d + cluster.sizes[k]/(N+alpha) * dt(x=(y-cluster.m[k])/t.sigma, df=2*cluster.a[k]) / t.sigma
                    }
                }
                d <- d + (cluster.sizes[K+1]+alpha)/(N+alpha) * dt(x=(y-ast.mean)/ast.sigma, df=ast.df) / ast.sigma
                return (d)
            }
            
            # stepwise visualization
            if (stepwise.plot & K>=0 & iter.since.stable<2) {
                if (min(X$x)<0) {
                    x.min <- min(X$x) * 1.1
                } else {
                    x.min <- min(X$x) / 1.1
                }
                if (max(X$x)>0) {
                    x.max <- max(X$x) * 1.1
                } else {
                    x.max <- max(X$x) / 1.1
                }
                x.seq <- seq(from=x.min, to=x.max, length.out = 200)
                d.true <- X$density.fun(x.seq)
                d.mle <- density.est.mle(x.seq)
                d.map <- density.est.mle(x.seq)
                d.Bayes <- density.est.Bayes(x.seq)
                max.d <- max(max(d.mle), max(d.map), max(d.Bayes), max(d.true)) * 1.1
                plot(x.seq, d.mle, t="l", col="blue", main=sprintf("iter %d, K=%d (round %d)", iter, K, inner.iter), 
                     ylim=c(0,max.d))
                #             lines(x.seq, d.map, t="l", col="brown")
                lines(x.seq, d.Bayes, t="l", col="darkgreen")
                lines(x.seq, d.true, t="l", col="red")
                points(X$x, rep(0.05, N))
                if (K>0) {
                    points(cluster.x.bar[1:K], rep(0.05, K), col="cyan", pch=4, cex=3, lwd=2)
                }
            }
            
        } # inner EM loop finish
        
        if (to.exit) {
            break
        }
        
        # remove clusters below threshold
        cluster.sizes <- colSums(W)
        tiny.clusters <- (1:K)[cluster.sizes[1:K]<cluster.size.threshold]
        num.tiny <- length(tiny.clusters)
        if (verbose>0) {
            cat(sprintf("- %d tiny clusters removed, remaining %d clusters \n", num.tiny, K-num.tiny))
        }
        # delete them and put weights back to *
        if (num.tiny>0) {
            cluster.x.bar <- cluster.x.bar[-tiny.clusters]
            cluster.var.bar <- cluster.var.bar[-tiny.clusters]
            cluster.m <- cluster.m[-tiny.clusters]
            cluster.a <- cluster.a[-tiny.clusters]
            cluster.kappa <- cluster.kappa[-tiny.clusters]
            cluster.b <- cluster.b[-tiny.clusters]
            cluster.mean.map <- cluster.mean.map[-tiny.clusters]
            cluster.var.map <- cluster.var.map[-tiny.clusters]
            if (num.tiny>1) {
                W[, K+1] <- W[, K+1] + rowSums(W[, tiny.clusters])
            } else {
                W[, K+1] <- W[, K+1] + W[, tiny.clusters] # rowSums breaks down in 1 col case, yike!
            }            
            W <- W[, -tiny.clusters]
            K <- K - num.tiny
            
        }
        W.tol <- W.tol/2
        inner.tol <- inner.tol/2
        stopifnot(length(cluster.x.bar)==K+1)
        stopifnot(K>0)
        
        if (W.tol<1e-2) {
            break
        }
        
    }
    if (iter>max.iter) {
        warning(sprintf("Algorithm terminated when reaching max.iter=%d", max.iter))
    }
    params <- list(alpha=alpha, prior.mean=prior.mean, prior.kappa=prior.kappa, 
                   prior.a=prior.a, prior.b=prior.b, 
                   estimate.type=estimate.type,
                   random.init.cluster=random.init.cluster,
                   cluster.size.threshold=cluster.size.threshold)
    
    # rerank the output by cluster size
    idx.new.to.old <- order(cluster.sizes[1:K], decreasing = T)
    cluster.sizes[1:K] <- cluster.sizes[idx.new.to.old]
    W[, (1:K)] <- W[, idx.new.to.old]
    cluster.x.bar[1:K] <- cluster.x.bar[idx.new.to.old]
    cluster.var.bar[1:K] <- cluster.var.bar[idx.new.to.old]
    cluster.m[1:K] <- cluster.m[idx.new.to.old]
    cluster.a[1:K] <- cluster.a[idx.new.to.old]
    cluster.kappa[1:K] <- cluster.kappa[idx.new.to.old]
    cluster.b[1:K] <- cluster.b[idx.new.to.old]
    cluster.mean.map[1:K] <- cluster.mean.map[idx.new.to.old]
    cluster.var.map[1:K] <- cluster.var.map[idx.new.to.old]
    display.df <- data.frame(cluster=c(1:K,"*"),
                             cluster.x.bar=cluster.x.bar, 
                             cluster.var.bar=cluster.var.bar, 
                             cluster.sizes=cluster.sizes)
    cat("\nResult\n")
    print(display.df)
    
    return(list(params = params, 
                K = K, N=N,
                cluster.x.bar=cluster.x.bar,
                cluster.var.bar=cluster.var.bar,
                cluster.m=cluster.m,
                cluster.a=cluster.a,
                cluster.kappa=cluster.kappa,
                cluster.b=cluster.b,
                cluster.mean.map=cluster.mean.map,
                cluster.var.map=cluster.var.map,
                W = W, 
                display.df = display.df))
}
