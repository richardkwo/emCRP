require(MASS)
library(grid)
library(gridExtra)

source("myutils.R")

# Neal's algorithm 3 
dp.gibbs.sample <- function(X, burn.in=500, sample.size=1000, 
                            alpha=NULL, prior.mean=NULL, prior.sd=NULL, sd=NULL){    
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
    N <- X$N
    p <- length(prior.mean)
    # init
    z <- rep(1, N)
    theta <- rbind(prior.mean)
    K <- 1
    cluster.sizes <- c(N)
    z.samples <- c()
    K.samples <- c()
    cluster.size.samples <- list()
    theta.samples <- list()
    
    prgBar <- txtProgressBar(min=0, max=burn.in+sample.size, style = 3)
    for (iter in 1:(burn.in+sample.size)) {
        if (iter%%10==0) {setTxtProgressBar(prgBar, iter)}
        z.new <- rep(0, N)
        # sample assignments
        for (i in 1:N) {
            if (sum(z[-i]==z[i])==0) {
                # remove cluster
                cluster.sizes[z[i]] <- 0           
                # theta[z[i]] will be removed before next round
            }
            # compute probs of cluster assignments
            V <- length(cluster.sizes)
            prob.vec <- rep(0, V+1)
            prob.vec[V+1] <- alpha * gaussian.prob.kernel(x=(X$x[i,]-prior.mean), 
                                                          s2=sd^2 + prior.sd^2)
            for (k in (1:V)[cluster.sizes>0]) {
                nc.no.i <- sum(z[-i]==k)
                if (nc.no.i>0) {
                    prob.vec[k] <- nc.no.i * gaussian.prob.kernel(x=(X$x[i,]-theta[k,]), 
                                                                  s2=sd^2)
                }
            }
            prob.vec <- prob.vec / sum(prob.vec)
            # draw cluster assignment
            zi.old <- z[i]
            z[i] <- sample(1:(V+1), size=1, prob=prob.vec)
            cluster.sizes[zi.old] <- cluster.sizes[zi.old] - 1
            # create a new cluster
            if (z[i]==V+1) {
                # draw a new theta from posterior given x[i]
                H <- get.gaussian.posterior(prior.mean, prior.sd, sd, x=X$x[i,])
                new.theta <- mvrnorm(n=1, mu = H$mean, Sigma = (H$sd^2)*diag(p))
                theta <- rbind(theta, new.theta)
                row.names(theta) <- NULL
#                 stopifnot(nrow(theta)==V+1)
                cluster.sizes[V+1] <- 1
            } else {
                cluster.sizes[z[i]] <- cluster.sizes[z[i]] + 1
            }
        }
        # get non-empty clusters and reorder them by size
        idx.new.to.old <- order(cluster.sizes, decreasing = TRUE)
        idx.old.to.new <- rep(0, length(idx.new.to.old))
        for (k in 1:length(idx.old.to.new)) {
            idx.old.to.new[idx.new.to.old[k]] <- k
        }
        cluster.sizes <- cluster.sizes[idx.new.to.old]
        cluster.sizes <- cluster.sizes[cluster.sizes>0]
        K <- length(cluster.sizes)
        z <- sapply(z, function(x) idx.old.to.new[x])
#         # self-test
#         for (i in 1:N) {
#             stopifnot(0<z[i] & z[i]<=K)
#         }
#         for (k in 1:K) {
#             stopifnot(sum(z==k)==cluster.sizes[k])
#         }
#         stopifnot(sum(cluster.sizes)==N)
        theta <- theta[idx.new.to.old[1:K], ]
        # sample parameters for each cluster
        for (k in 1:K) {
            # draw a new theta from posterior conditioned on this cluster
            G <- get.gaussian.posterior(prior.mean, prior.sd, sd, x=X$x[z==k,])
            theta[k, ] <- mvrnorm(n=1, mu=G$mean, Sigma=(G$sd^2)*diag(p))
        }
        if (iter>burn.in) {
            z.samples <- rbind(z.samples, z)
            K.samples <- c(K.samples, K)
            j <- iter-burn.in
            if (j==1) {
                theta.samples <- list(theta)
                cluster.size.samples <- list(table(z))
            } else {
                theta.samples[[j]] <- theta
                cluster.size.samples[[j]] <- table(z)
            }
        }
    }
    return(list(sample.size=sample.size, 
                theta.samples=theta.samples, 
                z.samples=z.samples, 
                K.samples=K.samples, 
                cluster.size.samples=cluster.size.samples))
}


estimate.from.gibbs <- function(X, gibbs.samples) {
    # get estimate from gibbs samples and compare to the truth
    K <- X$num.clusters
    x.plot.df <- data.frame(X$x)
    # real clusters
    for (k in 1:mean(gibbs.samples$K.samples)) {
        # estimate of cluster mean
        cluster.mean.samples <- t(sapply(gibbs.samples$theta.samples[gibbs.samples$K.samples>=k], 
                                         function(q) q[k, ]))
        cluster.mean.samples <- data.frame(cluster.mean.samples)
        est.mean <- apply(cluster.mean.samples, 2, mean)
        est.sd <- apply(cluster.mean.samples, 2, sd)
        est.df <- data.frame(x=est.mean[1],
                             y=est.mean[2], 
                             xmin=est.mean[1]-est.sd[1], 
                             xmax=est.mean[1]+est.sd[1],
                             ymin=est.mean[2]-est.sd[2],
                             ymax=est.mean[2]+est.sd[2])
        # estimate of sample weights
        est.weights <- apply(gibbs.samples$z.samples, 2, function(q) mean(q==k))
        est.cluster.size <- sum(est.weights)
        x.plot.df$est.weights <- est.weights
        x.plot.df$z <- X$z
        # plot
        fig.est <- ggplot(cluster.mean.samples, aes(x=X1, y=X2)) + 
            geom_density2d() + 
            geom_point(aes(x=x, y=y), data=est.df, size=4, color="green") + 
            geom_errorbar(aes(x=x,y=y,ymin=ymin,ymax=ymax), data=est.df, color="green") + 
            geom_errorbarh(aes(x=x,y=y,xmin=xmin,xmax=xmax), data=est.df, color="green") + 
            geom_point(aes(x=X1, y=X2, alpha=est.weights), data=x.plot.df) 
        if (k<=K) {
            real.mean <- X$cluster.means[k, ]
            fig.est <- fig.est + 
                geom_point(x=real.mean[1], y=real.mean[2], size=4, shape=8, color="red") +
                geom_point(aes(x=X1, y=X2), data=x.plot.df[x.plot.df$z==k, ], alpha=.5, size=3, shape=5, color="purple") + 
                ggtitle(sprintf("cluster %d: size=%d, est.size=%.1f", k, X$cluster.sizes[k], est.cluster.size))
        } else {
            fig.est <- fig.est + ggtitle(sprintf("spurious cluster %d: est.size=%.1f", k, est.cluster.size))
        }
        print(fig.est)
    }
}



