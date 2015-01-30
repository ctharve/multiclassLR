## given a design matrix & class counts, returns a vector of zero weights
InitialWeights <- function(design.mat, n.classes)
{
  n.params.ni <- dim(design.mat)[2]
  initweights.ni <- rep(0, n.params.ni*n.classes)
  initweights.ni
}

## given a design matrix & class labels, returns a vector of complete variable names
VarNames <- function(design.mat, configurations)
{
  indvars <- colnames(design.mat)
  if(is.null(indvars)){stop("design matrix must have column names")}
  var.names <- matrix(sapply(configurations, function(config){t(paste0(indvars, config))}), ncol=1)[, 1]
  var.names
}
 
## multinomial data will be transformed for ggplot2
CreatePlotDat <- function(multiDat)
{
  if(!("results" %in% names(multiDat))){stop("\"results\" absent from input data")}
  if(!("dat.all" %in% names(multiDat$results))){stop("\"dat.all\" absent from results data")}
  
  var.names <-  rownames(multiDat$results$dat.all)
  all.data <- data.frame(var.names, multiDat$results$dat.all, stringsAsFactors=FALSE)
  intercept_exclude <- !(all.data$var.names %in% c('intercept_cf1', 'intercept_cf2', 'intercept_cf3'))
  feature <- sapply(all.data$var.names, function(id){strsplit(id, '_')[[1]][1]}, USE.NAMES=FALSE)
  config <- as.factor(sapply(all.data$var.names, function(id){strsplit(id, '_')[[1]][2]}, USE.NAMES=FALSE))
  plot_dat <- data.frame(feature, config, parms=all.data$parms, se=all.data$parms.se, stringsAsFactors=FALSE)[intercept_exclude, ] 
  plot_dat
}

## log-likelihood closure, bundling the intercept & number of classes into a likelihood function
## for optim
MultiNegLlk.ni <- function(intercepts, n.classes)
{  
  MyMultiNegLlk <- function(w, X, Y)
  {
    # add an intercept term & intercept weights  
    X.int <- cbind(rep(1, dim(X)[1]), X)
    #w.int <- c(intercepts[1], w[1:37], intercepts[2], w[38:74], intercepts[3], w[75:111])
    parm.classes <- n.classes-1
    n.all.parms <- length(w)
    if((n.all.parms %% parm.classes) != 0){stop("parms vector length is not a multiple of n.classes-1")}
    n.parms <- n.all.parms / parm.classes
    w.int <- c(intercepts[1], w[1:n.parms], intercepts[2], w[(n.parms+1):(2*n.parms)], intercepts[3], w[(2*n.parms+1):(3*n.parms)])

    Y_hat <- YHatMulti(w.int, X.int)  
    llk <- sum(diag(t(Y) %*% log(Y_hat)))  # trace is invariant to cyclic permutations
    -llk
  }  
  MyMultiNegLlk
}

## optimization for a multinomial model. intercept fit from an independent
## run, avoiding the use of a baseline class
MyMultinomial.ni <- function(intercepts, parms, X, Y, var.names, hessian=TRUE, maxit=1E5)
{
  ## check argument validity ##
  ## TODO ##
  n.classes <- dim(Y)[2]
  cat('=== parameters are valid ===', '\n\n')

  ## create gradient closure with intercepts
  MyMultiNegLlk <- MultiNegLlk.ni(intercepts, n.classes)
  
  ## model fitting
  cat('=== begin optimization ===', '\n\n')
  st <- proc.time()
  model <- optim(parms, fn=MyMultiNegLlk, X=as.matrix(X), Y=as.matrix(Y), method="BFGS", hessian=TRUE, control=list(maxit=maxit));
  en <- proc.time()
  cat('Run time: ', en-st, '\n\n');
  
  ## results: model fit and parameter estimates
  cat('=== processing ouput ===', '\n\n')
  nllk <- model$value 
  converged <- model$convergence
  hessian <- model$hessian
  model_res <- MungeResults(model, var.names, classes=n.classes)

  retval <- list(nllk=nllk, converged=converged, hessian=hessian, results=model_res)
  retval
}
