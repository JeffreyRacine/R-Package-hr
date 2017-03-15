hr.test <- function(x=NULL,
                    K.vec=NULL,
                    random.seed=42,
                    B=399,
                    alpha=0.05,
                    trend=TRUE,
                    verbose=TRUE,
                    q=c(0.005,0.01,0.025,0.05,0.95,0.975,0.99,0.995)) {

    if(is.null(x)) stop("You must provide data")
    if(!is.ts(x)) stop("You must provide a time series data object")
    if(is.unsorted(K.vec)) stop("Lag vector K.vec must be sorted")
    if(alpha <= 0 | alpha >= 0.5) stop("Size (alpha) must lie in (0,0.5)")
    if(B < 1) stop("Number of bootstrap replications (B) must be a positive integer (e.g. 399)")
    if(any(q<=0) | any(q>=1)) stop("The quantile vector entries must lie in (0,1)")

    if(verbose) cat("\rLoading required packages")
    
    suppressMessages(require(fUnitRoots))
    suppressMessages(require(boot))
    suppressMessages(require(quadprog))
    suppressMessages(require(np))

    if(exists(".Random.seed", .GlobalEnv)) {
      save.seed <- get(".Random.seed", .GlobalEnv)
      exists.seed = TRUE
    } else {
      exists.seed = FALSE
    }
    
    set.seed(random.seed)

    n <- length(x)
    if(is.null(K.vec)) K.vec <- 1:round(12*(n/100)^0.25)
    K <- length(K.vec)

    stat <- function(s) {s}

    if(verbose) cat("\rComputing statistics and model averaging weights")

    out.nc <- suppressWarnings(adfTest(x,lags=0,type="nc"))
    out.c <- suppressWarnings(adfTest(x,lags=0,type="c"))

    residual.mat <- cbind(residuals(out.nc@test$lm),residuals(out.c@test$lm))
    rank.vec <- c(out.nc@test$lm$rank,out.c@test$lm$rank)
    t.stat <- c(out.nc@test$statistic,out.c@test$statistic)

    if(trend) {
        out.ct <- suppressWarnings(adfTest(x,lags=0,type="ct"))
        residual.mat <- cbind(residual.mat,residuals(out.ct@test$lm))
        rank.vec <- c(rank.vec,out.ct@test$lm$rank)
        t.stat <- c(t.stat,out.ct@test$statistic)
    }

    for(k in 1:K) {
        if(trend) {
            out <- suppressWarnings(adfTest(x,lags=K.vec[k],type="ct"))
        } else {
            out <- suppressWarnings(adfTest(x,lags=K.vec[k],type="c"))
        }
        ## Residual vector shorter with lags, need to line up properly
        if(k==1) {
            residual.mat <- cbind(residual.mat[(1+K.vec[k]):nrow(residual.mat),],residuals(out@test$lm))
        } else {
            residual.mat <- cbind(residual.mat[(1+K.vec[k]-K.vec[k-1]):nrow(residual.mat),],residuals(out@test$lm))
        }
        rank.vec <- c(rank.vec,out@test$lm$rank)
        t.stat <- c(t.stat,out@test$statistic)
    }

    M.dim <- ncol(residual.mat)
    sigsq.largest <- summary(out@test$lm)$sigma**2
    Dmat <- t(residual.mat)%*%residual.mat
    if(qr(Dmat)$rank<M.dim) Dmat <- Dmat + diag(1e-10,M.dim,M.dim)
    Amat <- cbind(rep(1,M.dim),diag(1,M.dim,M.dim))
    bvec <- c(1,rep(0,M.dim))
    dvec <- -rank.vec*sigsq.largest
    w.hat.mma <- solve.QP(Dmat,dvec,Amat,bvec,1)$solution

    t.stat.mma <- sum(t.stat*w.hat.mma)

    e <- diff(x,1)
    l <- b.star(e,round=TRUE)[1,1]

    t.stat.boot.mma <- numeric(length=B)
    
    if(verbose) cat("\r                                                ")

    for(b in 1:B) {
        if(verbose) cat(paste("\rBootstrap replication",b,"of",B))
        x.boot <- ts(c(x[1],cumsum(tsboot(e,stat,R=1,l=l,sim="geom")$t)),
                     frequency=frequency(x),
                     start=start(x))
        t.stat.boot <- c(suppressWarnings(adfTest(x.boot,lags=0,type="nc")@test$statistic),
                         suppressWarnings(adfTest(x.boot,lags=0,type="c")@test$statistic))
        if(trend) t.stat.boot <- c(t.stat.boot,suppressWarnings(adfTest(x.boot,lags=0,type="ct")@test$statistic))
        for(k in 1:K) {
            if(trend) {
                t.stat.boot <- c(t.stat.boot,suppressWarnings(adfTest(x.boot,lags=K.vec[k],type="ct")@test$statistic))
            } else {
                t.stat.boot <- c(t.stat.boot,suppressWarnings(adfTest(x.boot,lags=K.vec[k],type="c")@test$statistic))
            }
        }

        t.stat.boot.mma[b] <- sum(t.stat.boot*w.hat.mma)
        
    }
    if(verbose) cat("\r                               ")

    decision <- paste("Fail to reject at the",alpha,"level (unit root)")
    if(t.stat.mma < quantile(t.stat.boot.mma,probs=alpha/2,type=1)) decision <- paste("Reject at the",alpha,"level (stationary)")
    if(t.stat.mma > quantile(t.stat.boot.mma,probs=1-alpha/2,type=1)) decision <- paste("Reject at the",alpha,"level (explosive)")

    reject <- as.numeric(ifelse(t.stat.mma < quantile(t.stat.boot.mma,probs=alpha/2,type=1) |
                                t.stat.mma > quantile(t.stat.boot.mma,probs=1-alpha/2,type=1),1,0))
                        

    if(exists.seed) assign(".Random.seed", save.seed, .GlobalEnv)
    
    sigtest(tau=t.stat.mma,
                decision=decision,
                reject=reject,
                quantiles=quantile(t.stat.boot.mma,q,type=1),
                alpha=alpha,
                trend=trend,
                mma.weights=w.hat.mma,
                tau.boot=t.stat.boot.mma,
                e.block.length=l,
                boot.num=B,
                adf.lags=K.vec)

}

sigtest <- function(tau,
                    decision,
                    reject,
                    quantiles,
                    alpha,
                    trend,
                    mma.weights,
                    tau.boot,
                    e.block.length,
                    boot.num,
                    adf.lags) {
    
    tsig <- list(tau=tau,
                 decision=decision,
                 reject=reject,
                 quantiles=quantiles,
                 alpha=alpha,
                 trend=trend,
                 mma.weights=mma.weights,
                 tau.boot=tau.boot,
                 e.block.length=e.block.length,
                 boot.num=boot.num,
                 adf.lags=adf.lags)

     class(tsig) = "hrtest"
    
    tsig
}

print.hrtest <- function(x, ...){
    cat("\nHansen-Racine Nonparametric Unit Root Test",
        "\nBootstrap (",x$boot.num," replications,")
}

summary.hrtest <- function(object, ...) {
    print(object)
}
