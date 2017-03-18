## The Hansen-Racine nonparametric bootstrap model average unit root test

hr.test <- function(x=NULL,
                    K.vec=NULL,
                    random.seed=42,
                    B=399,
                    alpha=0.05,
                    trend=TRUE,
                    type=c("ct","c","nc","all"),
                    verbose=TRUE,
                    S=12,
                    method=c("mma","jma"),
                    q=c(0.005,0.01,0.025,0.05,0.95,0.975,0.99,0.995)) {

    ## Some basic input checking
    
    method <- match.arg(method)
    type <- match.arg(type)
    if(type=="all") type <- c("nc", "c", "ct")
    
    if(is.null(x)) stop("You must provide data")
    if(!is.ts(x)) x <- ts(x)
    if(is.unsorted(K.vec)) stop("Lag vector K.vec must be sorted")
    if(any(K.vec<1)) stop("Lag vector K.vec must contain positive integers")
    if(alpha <= 0 | alpha >= 0.5) stop("Size (alpha) must lie in (0,0.5)")
    if(B < 1) stop("Number of bootstrap replications (B) must be a positive integer (e.g. 399)")
    if(any(q<=0) | any(q>=1)) stop("The quantile vector entries must lie in (0,1)")
    
    ## Save any existing random seed and restore upon exit

    if(exists(".Random.seed", .GlobalEnv)) {
      save.seed <- get(".Random.seed", .GlobalEnv)
      exists.seed = TRUE
    } else {
      exists.seed = FALSE
    }
    
    ## Set the random seed
    
    set.seed(random.seed)

    ## Set the length of the vector
    
    n <- length(x)
    
    ## Use Schwert's ad-hoc rule for the maximum lag for the candidate models 
    ## if none is provided
    
    if(is.null(K.vec)) K.vec <- 0:round(S*(n/100)^0.25)
    K <- length(K.vec)

    ## A simple function that returns its argument for the tsboot() call
    
    stat <- function(s) {s}

    if(verbose) cat("\rComputing statistics and model averaging weights")
    
    t.stat <- NULL
    rank.vec <- NULL
    
    first <- TRUE
    for(k in 1:K) {
        for(t in type) {
            out <- suppressWarnings(adfTest(x,lags=K.vec[k],type=t))
            if(method=="mma") {
                r <- residuals(out@test$lm)
            } else {
                r <- jackknife.prediction(out@test$lm) 
            }
            if(first) {
                ma.mat <- as.matrix(r)
            } else {
                n.r <- length(r)
                n.rm <- nrow(ma.mat)
                ma.mat <- cbind(ma.mat[(n.rm-n.r+1):n.rm,],r)
            }
            rank.vec <- c(rank.vec,out@test$lm$rank)
            t.stat <- c(t.stat,out@test$statistic)
            first <- FALSE
        }
    }
    
    print("Here we are")

    ## Model average weights (solve a simple quadratic program)

    M.dim <- ncol(ma.mat)
    sigsq.largest <- summary(out@test$lm)$sigma**2
    Dmat <- t(ma.mat)%*%ma.mat
    if(qr(Dmat)$rank<M.dim) Dmat <- Dmat + diag(1e-10,M.dim,M.dim)
    Amat <- cbind(rep(1,M.dim),diag(1,M.dim,M.dim))
    bvec <- c(1,rep(0,M.dim))
    if(method=="mma") {
        dvec <- -rank.vec*sigsq.largest
    } else {
        dvec <- t(as.matrix(x[(n-nrow(ma.mat)+1):n]))%*%ma.mat
    }
    print(K.vec)
    print(dim(ma.mat))
    print(dim(Dmat))
    print(M.dim)
    print(length(t.stat))
    print(length(rank.vec))
    w.hat.ma <- solve.QP(Dmat,dvec,Amat,bvec,1)$solution
    
    #print(dim(ma.mat))
    #print(length(t.stat))
    #print(length(K.vec))

    ## The model average test statistic is a weighted average of each of the above 
    ## candidate model's test statistics (all t-statistics for the coefficient on 
    ## the first lag of the series)

    t.stat.ma <- sum(t.stat*w.hat.ma)

    ## Impose the null with a model-free difference (See Swensen (2003) for a 
    ## similar procedure)

    e <- diff(x,1)
    l <- b.star(e,round=TRUE)[1,1]
    
    ## Vector to hold the bootstrap statistics

    t.stat.boot.ma <- numeric(length=B)
    
    if(verbose) cat("\r                                                ")

    for(b in 1:B) {
        
        if(verbose) cat(paste("\rBootstrap replication",b,"of",B))
        
        ## Generate a bootstrap resample under the null
        
        x.boot <- ts(c(x[1],cumsum(tsboot(e,stat,R=1,l=l,sim="geom")$t)),
                     frequency=frequency(x),
                     start=start(x))
        
        ## Recompute all candidate models and their test statistics
        
        t.stat.boot <- NULL
        
        for(k in 1:K) {
            for(t in type) {
                t.stat.boot <- c(t.stat.boot,suppressWarnings(adfTest(x.boot,lags=K.vec[k],type=t)@test$statistic))
            }
        }
        
        ## Compute the model average bootstrap statistics

        t.stat.boot.ma[b] <- sum(t.stat.boot*w.hat.ma)
        
    }
    
    ## Set return objects

    if(verbose) cat("\r                               ")

    decision <- paste("Fail to reject at the ",100*alpha,"% level (unit root)",sep="")
    if(t.stat.ma < quantile(t.stat.boot.ma,probs=alpha/2,type=1)) decision <- paste("Reject at the ",100*alpha,"% level (stationary)",sep="")
    if(t.stat.ma > quantile(t.stat.boot.ma,probs=1-alpha/2,type=1)) decision <- paste("Reject at the ",100*alpha,"% level (explosive)",sep="")

    reject <- as.numeric(ifelse(t.stat.ma < quantile(t.stat.boot.ma,probs=alpha/2,type=1) |
                                t.stat.ma > quantile(t.stat.boot.ma,probs=1-alpha/2,type=1),1,0))
                        

    if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)
    
    hrtest(tau = t.stat.ma,
           tau.alpha.low = quantile(t.stat.boot.ma,probs=alpha/2,type=1),
           tau.alpha.up = quantile(t.stat.boot.ma,probs=1-alpha/2,type=1),
           decision = decision,
           reject = reject,
           quantiles = quantile(t.stat.boot.ma,q,type=1),
           alpha = alpha,
           trend = trend,
           ma.weights = w.hat.ma,
           tau.boot = sort(t.stat.boot.ma),
           e.block.length = l,
           boot.num = B,
           adf.lags = K.vec)

}

## S3 functions for summary() and print()

hrtest <- function(tau,
                   tau.alpha.low,
                   tau.alpha.up,
                   decision,
                   reject,
                   quantiles,
                   alpha,
                   trend,
                   ma.weights,
                   tau.boot,
                   e.block.length,
                   boot.num,
                   adf.lags) {
   
    thr <- list(tau = tau,
                tau.alpha.low = tau.alpha.low,
                tau.alpha.up = tau.alpha.up,
                decision = decision,
                reject = reject,
                quantiles = quantiles,
                alpha = alpha,
                trend = trend,
                ma.weights = ma.weights,
                tau.boot = tau.boot,
                e.block.length = e.block.length,
                boot.num = boot.num,
                adf.lags = adf.lags)

    class(thr) = "hrtest"
    
    thr
}

print.hrtest <- function(x, ...){
    cat("\n        Bootstrap Model Averaged Unit Root Test\n",
        "\nTest statistic: ",x$tau,
        "\n",100*x$alpha,"% critical values: (",x$tau.alpha.low,",",x$tau.alpha.up,")",
        "\n",x$decision,
        "\nBootstrap replications: ",x$boot.num,
        "\nAutomatic expected block length: ",x$e.block.length,"\n",sep="")
}

summary.hrtest <- function(object, ...) {
    print(object)
}

jackknife.prediction <- function(model) {
    
    htt <- hatvalues(model)
    return(fitted(model) - htt*residuals(model)/(1-htt))
    
}
