require(ggplot2)
require(proto)

StatEllipse <- proto(ggplot2:::Stat,
{
    required_aes <- c("x", "y")
    default_geom <- function(.) GeomPath
    objname <- "ellipse"
    
    calculate_groups <- function(., data, scales, ...){
        .super$calculate_groups(., data, scales,...)
    }
    calculate <- function(., data, scales, level = 0.95, segments = 51,...){
        dfn <- 2
        dfd <- length(data$x) - 1
        if (dfd < 3){
            ellipse <- rbind(c(NA,NA))	
        } else {
            require(MASS)
            v <- cov.trob(cbind(data$x, data$y))
            shape <- v$cov
            center <- v$center
            radius <- sqrt(dfn * qf(level, dfn, dfd))
            angles <- (0:segments) * 2 * pi/segments
            unit.circle <- cbind(cos(angles), sin(angles))
            ellipse <- t(center + radius * t(unit.circle %*% chol(shape)))
        }
        
        ellipse <- as.data.frame(ellipse)
        colnames(ellipse) <- c("x","y")
        return(ellipse)
    }
}
)

stat_ellipse <- function(mapping=NULL, data=NULL, geom="path", position="identity", ...) {
    StatEllipse$new(mapping=mapping, data=data, geom=geom, position=position, ...)
}

get.gaussian.posterior <- function(prior.mean, prior.std, std, x, weights=1) {
    # x must take the form of N x p
    if (is.null(nrow(x))) {
        x <- rbind(x)
    }
    n <- nrow(x)
    if (length(weights)==1) {
        weights <- rep(1, weights)
    } else {
        stopifnot(length(weights)==n)
    }
    nc <- sum(weights)
    posterior.std <- 1/sqrt(1/prior.std^2 + nc/std^2)
    posterior.mean <- posterior.std^2 * (prior.mean/prior.std^2 + colSums(weights * x)/std^2)
    return(list(mean=posterior.mean, std=posterior.std))
}

gaussian.prob.kernel <- function(x, s2) {
    p <- length(x)
    return ((s2)^(-p/2) * exp(-1/(2*s2) * sum(x^2)))
}

multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
    require(grid)
    
    # Make a list from the ... arguments and plotlist
    plots <- c(list(...), plotlist)
    
    numPlots = length(plots)
    
    # If layout is NULL, then use 'cols' to determine layout
    if (is.null(layout)) {
        # Make the panel
        # ncol: Number of columns of plots
        # nrow: Number of rows needed, calculated from # of cols
        layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                         ncol = cols, nrow = ceiling(numPlots/cols))
    }
    
    if (numPlots==1) {
        print(plots[[1]])
        
    } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
        
        # Make each plot, in the correct location
        for (i in 1:numPlots) {
            # Get the i,j matrix positions of the regions that contain this subplot
            matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
            
            print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                            layout.pos.col = matchidx$col))
        }
    }
}