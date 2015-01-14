##' @title MyMultinomialParallel
##'
#

MyMultinomialParallel <- function(parms, X, Y, var.names, hessian=TRUE, maxit=1E5){

  ## check argument validity ##
  ## TODO ##
  cat('=== parameters are valid ===', '\n\n')
  
  n.obs <- dim(X)[1]
  
  ## n.part @param number of partitions
  ## TODO: need to make a function for determining partitioning
  n.part <- floor(n.obs/5E4)+1  

  part.dat <- vector("list", n.part)

  ## outer loop
  ## 0.) initialize n.part parameter vectors
  ## 1.) in parallel, run optimization on each of n.part partitioned data sets for
  ##     n.iters, resulting in n.part parm sets
  ## 2.) for each of n.part parm sets evaluate the nllk on the entire data set 
  ## 3.) return to 1.) and start at the vector from 2.) with minimum nllk    
  ## 
  ##
  cat('=== begin optimization ===', '\n\n')
  st <- proc.time()
  for(ii in 1:max.iters){

    if(ii==1){
      final.vector <- rep(NA, n.parms)
    }
    
    parm.update <- ParallelLapply(part.dat, functions(dat){
      model <- optim(parms, fn=MultiNegLlk, X=as.matrix(X), Y=as.matrix(Y), method="BFGS", hessian=TRUE, control=list(maxit=maxit));

    })

    ## find optimal parm vector

    
  }
  en <- proc.time()
  cat('Run time: ', en-st, '\n\n');
  
  ## using the optimal vector, evaluate the hessian


  ## results: model fit and parameter estimates
  cat('=== processing ouput ===', '\n\n')
  nllk <- model$value 
  converged <- model$convergence
  hessian <- model$hessian
  model_res <- MungeResults(model, var.names, classes=4)

  retval <- list(nllk=nllk, converged=converged, hessian=hessian, results=model_res)
  retval
}
