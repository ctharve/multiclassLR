##' @title MungeResults
##'

MungeResults <-
function(optobj, parmnames, classes, sig=0.05){
  require(MASS)
  n.c <- classes-1
  n.p <- length(parmnames)/n.c

  parms <- optobj$par
  #parms.se <- diag(solve(optobj$hessian))^(.5)
  parms.se <- diag(ginv(optobj$hessian))^(.5)
  parms.z <- parms * parms.se^(-1)
  parms.pval <- 2*pnorm(-abs(parms.z))
  betas <- exp(parms)
  dat.parms <- data.frame(betas, parms, parms.se, parms.z, parms.pval)
  rownames(dat.parms) <- parmnames
  
  index <- which(dat.parms$parms.pval<(sig))
  dat.sig <- dat.parms[index, ]

  ## TODO: these repetative tasks should be done over a list
  pmat.betas <- matrix(betas, ncol=n.c)
  pmat.parms <- matrix(parms, ncol=n.c)
  pmat.se <- matrix(parms.se, ncol=n.c)
  pmat.z <- matrix(parms.z, ncol=n.c)
  rownames(pmat.betas) <- parmnames[1:n.p]
  rownames(pmat.parms) <- parmnames[1:n.p]
  rownames(pmat.se) <- parmnames[1:n.p]
  rownames(pmat.z) <- parmnames[1:n.p]
  dat.mat <- list(pmat.betas=pmat.betas, pmat.parms=pmat.parms, pmat.z=pmat.z, pmat.se=pmat.se)
  
  retval <- list(dat.all=dat.parms, dat.sig=dat.sig, dat.mat=dat.mat)
  retval   
}
