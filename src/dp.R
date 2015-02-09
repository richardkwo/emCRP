library(MASS)
library(ggplot2)

sample.dp.normal.mixture <- function(N, alpha, std=1, prior.mean=0, prior.std=1) {
    # DP mixture of MVN Gaussians with isotropic cov
    i=1
    p <- length(prior.mean)
    prior.Sigma <- prior.std^2 * diag(p)
    Sigma <- std^2 * diag(p)
    while(i<=N) {
        if (i==1) {
            # creating the first cluster 
            num.clusters <- 1
            cluster.assignments <- 1
            cluster.sizes <- c(1)
            cluster.means <- rbind(mvrnorm(n=1, mu=prior.mean, Sigma=prior.Sigma))
            x <- mvrnorm(n=1, mu=cluster.means, Sigma=Sigma)
            i <- i+1
        }
        else {
            # CRP
            table.probs <- c(alpha, cluster.sizes)
            table.probs <- table.probs / sum(table.probs)
            z <- sample(0:num.clusters, size=1, prob=table.probs)
            if (z==0) {
                # new cluster
                num.clusters <- num.clusters + 1
                z <- num.clusters
                cluster.sizes <- c(cluster.sizes, 1)                
                m <- mvrnorm(n=1, mu=prior.mean, Sigma=prior.Sigma)
                cluster.means <- rbind(cluster.means, m)
                row.names(cluster.means) <- NULL
            } else {
                cluster.sizes[z] <- cluster.sizes[z] + 1
                m <- cluster.means[z, ]
            }
            cluster.assignments <- c(cluster.assignments, z)
            x <- rbind(x, mvrnorm(n=1, mu=m, Sigma=Sigma))
            i <- i+1
        }
    }
    colnames(x) <- NULL
    row.names(x)<- NULL
    
    # re-ordering clusters in decreasing size
    old.to.new <- rank(-cluster.sizes, ties.method = "first")
    new.to.old <- rep(0, num.clusters)
    for (i in 1:num.clusters) {
        new.to.old[old.to.new[i]] <- i
    }
    cluster.assignments <- old.to.new[cluster.assignments]
    cluster.sizes <- cluster.sizes[new.to.old]
    cluster.means <- cluster.means[new.to.old, ]
    
    return(list(x=x, 
                N=N, 
                alpha = alpha, 
                prior.mean = prior.mean, 
                prior.std = prior.std, 
                std = std,
                num.clusters=num.clusters,
                z=cluster.assignments, 
                cluster.sizes=cluster.sizes, 
                cluster.means=cluster.means))
}


plot.dp.mixture <- function(X) {
    # plot the data clusters in 2D
    plot.df <- data.frame(x.1=X$x[,1], x.2=X$x[,2], z=factor(X$z))
    mean.df <- data.frame(x.1=X$cluster.means[,1], x.2=X$cluster.means[,2], z=factor(1:(X$num.clusters)))
    fig <- ggplot(plot.df) + 
        geom_point(aes(x=x.1, y=x.2, color=z)) + 
        geom_point(aes(x=x.1, y=x.2, color=z), size=8, shape=8, data=mean.df)
    return(fig)
}

plot.dp.cluster.size <- function(X) {
    # plot the relative size distribution
    relative.sizes <- X$cluster.sizes / sum(X$cluster.sizes)
    v <- 1:X$num.clusters
    plot(v, relative.sizes, type="p", log="y")
    predicted.sizes <- exp(-v/X$alpha)
    lines(v, predicted.sizes, type="l")
}