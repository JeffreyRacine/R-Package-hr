## The Hansen-Racine nonparametric bootstrap model average unit root test

hr.test <- function(x=NULL,
                    adf.type=c("c","ct","nc","all"),
                    alpha=0.05,
                    alternative=c("both","stationary","explosive"),
                    B=399,
                    boot.method=c("geom","fixed","iid"),
                    df.type=c("nc","nccct","ncc","c","ct","none"),
                    group.start=4,
                    group.by=4,                    
                    lag.vec=NULL,
                    method=c("jma","mma"),
                    quantile.vec=c(0.005,0.01,0.025,0.05,0.95,0.975,0.99,0.995),
                    random.seed=42,
                    S=12,
                    verbose=TRUE) {

    ## Some basic input checking and conversion
    
    alternative <- match.arg(alternative)
    boot.method <- match.arg(boot.method)
    method <- match.arg(method)
    adf.type <- match.arg(adf.type)
    if(adf.type=="all") adf.type <- c("nc", "c", "ct")
    df.type <- match.arg(df.type)
    if(df.type=="ncc") {
        df.type <- c("nc", "c")
    } else if(df.type=="nccct") {
        df.type <- c("nc", "c", "ct")        
    } else if(df.type=="none") {
        df.type <- NULL
    }
    
    if(is.null(x)) stop("You must provide data")
    if(!is.ts(x)) x <- ts(x)
    if(is.unsorted(lag.vec)) sort(lag.vec)
    if(any(lag.vec<1)) stop("Lag vector lag.vec must contain positive integers")
    if(alpha <= 0 | alpha >= 0.5) stop("Size (alpha) must lie in (0,0.5)")
    if(B<0) stop("Number of bootstrap replications (B) must be positive")
    if(alpha*(B+1)!=floor(alpha*(B+1))) stop("alpha*(B+1) must be an integer")
    if(any(quantile.vec<=0) | any(quantile.vec>=1)) stop("The quantile vector entries must lie in (0,1)")
    
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
    
    ## Variation on Schwert's ad-hoc rule for the maximum lag for the 
    ## candidate models if none is provided (constant may differ from S=12 
    ## that he uses). Hansen (2014) suggests blocks of 4 or more regressors
    ## improves MSE.
    
    if(is.null(lag.vec)) lag.vec <- seq(group.start,ceiling(S*(n/100)^.25),by=group.by)
    K <- length(lag.vec)

    ## A simple function that returns its argument for the tsboot() call
    
    stat <- function(s) {s}

    if(verbose) cat("\rComputing statistics and model averaging weights")
    
    ## Dickey-Fuller models (no lags)
    
    for(t in df.type) {
        out <- suppressWarnings(adfTest(x,lags=0,type=t))
        if(!exists("ma.mat")) {
            ma.mat <- as.matrix(Dmat.func(out@test$lm,method=method))
            rank.vec <- out@test$lm$rank
            t.stat <- out@test$statistic
        } else {
            ma.mat <- cbind(ma.mat,Dmat.func(out@test$lm,method=method))
            rank.vec <- c(rank.vec,out@test$lm$rank)
            t.stat <- c(t.stat,out@test$statistic)
        }
    }
    
    ## Augmented Dickey-Fuller models (lagged first differences)

    for(k in 1:K) {
        for(t in adf.type) {
            out <- suppressWarnings(adfTest(x,lags=lag.vec[k],type=t))
            r <- Dmat.func(out@test$lm,method=method) 
            if(!exists("ma.mat")) {
                ma.mat <- as.matrix(r)
                rank.vec <- out@test$lm$rank
                t.stat <- out@test$statistic
            } else {
                n.r <- length(r)
                n.rm <- nrow(ma.mat)
                ma.mat <- cbind(ma.mat[(n.rm-n.r+1):n.rm,],r)
                rank.vec <- c(rank.vec,out@test$lm$rank)
                t.stat <- c(t.stat,out@test$statistic)
            }
        }
    }
    
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
    
    w.hat.ma <- solve.QP(Dmat,dvec,Amat,bvec,1)$solution

    ## The model average test statistic is a weighted average of each of the above 
    ## candidate model's test statistics (all t-statistics for the coefficient on 
    ## the first lag of the series)

    t.stat.ma <- sum(t.stat*w.hat.ma)

    ## Impose the null with a model-free difference (See Swensen (2003) for a 
    ## similar procedure)

    e <- diff(x,1)
    if(boot.method == "iid") {
        ## An IID bootstrap is obtained via a a fixed block bootstrap with a
        ## block length of 1
        l <- 1
        boot.method <- "fixed"
    } else {
        l <- b.star(e,round=TRUE)[1,1]
    }
    
    ## Vector to hold the bootstrap statistics

    t.stat.boot.ma <- numeric(length=B)
    
    if(verbose) cat("\r                                                ")

    for(b in 1:B) {
        
        if(verbose) cat(paste("\rBootstrap replication",b,"of",B))
        
        ## Generate a bootstrap resample under the null
        
        x.boot <- ts(c(x[1],cumsum(tsboot(e,stat,R=1,l=l,sim=boot.method)$t)),
                     frequency=frequency(x),
                     start=start(x))
        
        ## Recompute all candidate models and their test statistics
        
        t.stat.boot <- NULL
        for(t in df.type) {
            t.stat.boot <- c(t.stat.boot,suppressWarnings(adfTest(x.boot,lags=0,type=t)@test$statistic))
        }        
        
        for(k in 1:K) {
            for(t in adf.type) {
                t.stat.boot <- c(t.stat.boot,suppressWarnings(adfTest(x.boot,lags=lag.vec[k],type=t)@test$statistic))
            }
        }
        
        ## Compute the model average bootstrap statistics

        t.stat.boot.ma[b] <- sum(t.stat.boot*w.hat.ma)
        
    }
    
    ## Set return objects

    if(verbose) cat("\r                               ")

    decision <- paste("Fail to reject the null hypothesis at the ",100*alpha,"% level (unit root)",sep="")
    if(alternative=="both") {
        if(t.stat.ma < quantile(t.stat.boot.ma,probs=alpha/2,type=1)) decision <- paste("Reject the null hypothesis at the ",100*alpha,"% level (stationary)",sep="")
        if(t.stat.ma > quantile(t.stat.boot.ma,probs=1-alpha/2,type=1)) decision <- paste("Reject the null hypothesis at the ",100*alpha,"% level (explosive)",sep="")
        reject <- as.numeric(ifelse(t.stat.ma < quantile(t.stat.boot.ma,probs=alpha/2,type=1) |
                                    t.stat.ma > quantile(t.stat.boot.ma,probs=1-alpha/2,type=1),1,0))
        tau.alpha.low <- quantile(t.stat.boot.ma,probs=alpha/2,type=1)
        tau.alpha.up <- quantile(t.stat.boot.ma,probs=1-alpha/2,type=1)
        alternative <- "stationary or explosive"
    } else if(alternative=="stationary") {
        if(t.stat.ma < quantile(t.stat.boot.ma,probs=alpha,type=1)) decision <- paste("Reject the null hypothesis at the ",100*alpha,"% level (stationary)",sep="")
        reject <- as.numeric(ifelse(t.stat.ma < quantile(t.stat.boot.ma,probs=alpha,type=1),1,0))
        tau.alpha.low <- quantile(t.stat.boot.ma,probs=alpha,type=1)
        tau.alpha.up <- NA
    } else if(alternative=="explosive") {
        if(t.stat.ma > quantile(t.stat.boot.ma,probs=1-alpha,type=1)) decision <- paste("Reject the null hypothesis at the ",100*alpha,"% level (explosive)",sep="")
        reject <- as.numeric(ifelse(t.stat.ma > quantile(t.stat.boot.ma,probs=1-alpha,type=1),1,0))        
        tau.alpha.low <- NA
        tau.alpha.up <- quantile(t.stat.boot.ma,probs=1-alpha,type=1)
    }

    if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)
    
    hrtest(tau = t.stat.ma,
           tau.alpha.low = tau.alpha.low,
           tau.alpha.up = tau.alpha.up,
           decision = decision,
           reject = reject,
           alternative = alternative,
           quantiles = quantile(t.stat.boot.ma,quantile.vec,type=1),
           alpha = alpha,
           ma.weights = w.hat.ma,
           tau.boot = sort(t.stat.boot.ma),
           e.block.length = l,
           boot.num = B,
           adf.lags = lag.vec,
           varname = deparse(substitute(x)))

}

## S3 functions for summary() and print()

hrtest <- function(tau,
                   tau.alpha.low,
                   tau.alpha.up,
                   decision,
                   reject,
                   alternative,
                   quantiles,
                   alpha,
                   ma.weights,
                   tau.boot,
                   e.block.length,
                   boot.num,
                   adf.lags,
                   varname) {
   
    thr <- list(tau = tau,
                tau.alpha.low = tau.alpha.low,
                tau.alpha.up = tau.alpha.up,
                decision = decision,
                reject = reject,
                alternative = alternative,
                quantiles = quantiles,
                alpha = alpha,
                ma.weights = ma.weights,
                tau.boot = tau.boot,
                e.block.length = e.block.length,
                boot.num = boot.num,
                adf.lags = adf.lags,
                varname = varname)

    class(thr) <- "hrtest"
    
    thr
}

print.hrtest <- function(x, ...){
    cat("\n        Hansen-Racine Bootstrap Model Averaged Unit Root Test\n",
        "\nData: ",x$varname,
        "\nAlternative hypothesis: ",x$alternative,
        "\nTest statistic: ",x$tau,
        "\n",100*x$alpha,"% critical value(s): (",x$tau.alpha.low,",",x$tau.alpha.up,")",
        "\n",x$decision,
        "\nBootstrap replications: ",x$boot.num,
        "\nAutomatic expected block length: ",x$e.block.length,"\n",sep="")
}

summary.hrtest <- function(object, ...) {
    print(object)
}

Dmat.func <- function(model,method=c("mma","jma")) {
    method <- match.arg(method)
    ## Residuals (Mallows model averaging)
    if(method=="mma") return(residuals(model))
    ## Jackknife fitted values (jackknife model averaging)
    htt <- hatvalues(model)
    if(method=="jma") return(fitted(model) - htt*residuals(model)/(1-htt))    
}
