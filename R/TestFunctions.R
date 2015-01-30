InitialWeights <- function(design.mat, n.classes)
{
  n.params.ni <- dim(design.mat)[2]
  initweights.ni <- rep(0, n.params.ni*n.classes)
  initweights.ni
}

VarNames <- function(design.mat, configurations)
{
  indvars <- colnames(design.mat)
  var.names <- matrix(sapply(configurations, function(config){t(paste0(indvars, config))}), ncol=1)[, 1]
  var.names
}
 

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

###############################
## 1.5.) update functions to accept intercepts fit individually
###############################

MultiNegLlk.ni <- function(intercepts){
  
  MyMultiNegLlk <- function(w, X, Y){

  # add an intercept term & intercept weights  
  X.int <- cbind(rep(1, dim(X)[1]), X)
  ## TODO: classes=3 is hardcoded, and we need to change that
  #w.int <- c(intercepts[1], w[1:37], intercepts[2], w[38:74], intercepts[3], w[75:111])
  ## this is so bad and needs to be fixed, pronto
  w.int <- c(intercepts[1], w[1:36], intercepts[2], w[37:72], intercepts[3], w[73:108])

  Y_hat <- YHatMulti(w.int, X.int)  
  llk <- sum(diag(t(Y) %*% log(Y_hat)))  # trace is invariant to cyclic permutations
  -llk
  }  
  MyMultiNegLlk
}

MyMultinomial.ni <- function(intercepts, parms, X, Y, var.names, hessian=TRUE, maxit=1E5)
{
  ## check argument validity ##
  ## TODO ##
  cat('=== parameters are valid ===', '\n\n')

  ## create gradient closure with intercepts
  MyMultiNegLlk <- MultiNegLlk.ni(intercepts)
  
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
  model_res <- MungeResults(model, var.names, classes=4)

  retval <- list(nllk=nllk, converged=converged, hessian=hessian, results=model_res)
  retval
}

YHatMulti <- function (w, X, baseline = TRUE) 
{
    p <- dim(X)[2]
    #cat('parm vect length: ', length(w), '\n')
    #cat('design features: ', p, '\n')
    W <- matrix(w, nrow = p)
    class_inner_prod <- exp(X %*% W)
    denom <- (1 + rowSums(class_inner_prod))^-1
    if (baseline) {
        cbind(class_inner_prod * denom, denom)
    }
    else {
        class_inner_prod * denom
    }
}
