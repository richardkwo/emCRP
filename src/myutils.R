require(ggplot2)

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