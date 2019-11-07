


#' Pseudo-observations for predictions for an illness-death model
#'
#' Computes pseudo-observations for modeling covariates effects on differents quantites of an
#' illness-death model.
#'
#' @param object an idmloo class object returned by a call to the idmloo function
#' @param s time point at which prediction is made
#' @param t time horizon for prediction
#' @param lifeExpect logical. if TRUE conmpute life expectancies, i.e., t=Inf
#' @param maxtime The upper limit of integration for calculations of life expectancies from Weibull parametrizations.
#' @param parallel Logical determinating parallel (=TRUE) or sequential (=FALSE) execution
#' @param cpus Numerical amount of CPUs requested for the cluster. If not set, values are total cores of the machines-1.
#' @param type Type of cluster. Can be 'SOCK', 'MPI', 'PVM' or 'NWS'. Default is 'SOCK' if parallel=TRUE.
#' @param ... others parameters
#'
#' @return pseudo
#' a matrix with each row corresponding to the pseudo-observations of one individual (ordered as in the
#' the original data set), each column presents the predicted quantities
#' @export
#'
#' @examples
#' pseudo.t10=pseudo(fit.loo.w, t=10,lifeExpect=FALSE, parallel = FALSE, cpus=NULL, type=NULL)
pseudo<-function(object, s=NA, t, lifeExpect,
                 parallel=FALSE, cpus=NULL, type=NULL,...){
  if (class(object)!="idmloo")stop("Argument 'object' must be an idmloo class object from idmloo function")
  if (missing(t) && lifeExpect==FALSE) stop("Argument t is missing.")
  if (lifeExpect==TRUE) t <- Inf
  # if (missing(s)) stop("Argument s is missing.")
  if (parallel==TRUE){
    if(is.null(type))type="SOCK"
    if(is.null(cpus)){
      cpus=parallel::detectCores()-1
    }
    else{
      if(is.character(cpus))stop("Argument 'cpus' has to be an integer")
      if(cpus<=1)stop("Argument 'cpus' has to be greater than 1")
      if(cpus>parallel::detectCores() & type!="MPI")stop("Argument 'cpus' must be at least the maximum of cores available (if type is not 'MPI')")
      if(length(cpus)>1)cpus=cpus[1]
    }
  }
  N=length(object)-1

  for (i in 1:N){
    if(class(object[[paste0("loo.",i)]])!="idm")stop("All leave-one-out illness-death model need to have been converged")
  }
  idm.N=object$loo.all
  object$loo.all<-NULL
  s<-ifelse(missing(s),0,s)
  s<-ifelse(idm.N$method=="splines", max(s,min(idm.N$knots01,idm.N$knots02,idm.N$knots12)),s)
  THETA=predict(object=idm.N, s=s, t=t, lifeExpect=lifeExpect,
                conf.int=FALSE)$transprob$Estimate

  predict.loo=function(object, s, t, lifeExpect){
    predict(object=object, s=s, t=t, lifeExpect=lifeExpect, conf.int=FALSE)$transprob$Estimate
  }

  if(parallel==FALSE){
    THETA.loo=lapply(object, predict.loo, s=s,t=t,lifeExpect=lifeExpect)
  }

  if(parallel==TRUE){
    if(sfIsRunning()==FALSE){
      sfInit(parallel = TRUE, cpus=cpus, type=type)
      sfLibrary(SmoothHazard)
      THETA.loo=sfLapply(object, predict.loo, s=s,t=t,lifeExpect=lifeExpect)
      sfStop()
    }
    else{
      sfLibrary(SmoothHazard)
      THETA.loo=sfLapply(object, predict.loo, s=s,t=t,lifeExpect=lifeExpect)

    }
  }
  ncol=ifelse(lifeExpect==FALSE, 10, 5)
  pseudo=matrix(ncol = ncol, nrow=N )
  for (i in 1:N){
    pseudo[i,]=N%*%THETA-(N-1)%*%THETA.loo[[i]]
  }
  if(lifeExpect==FALSE){
    colnames(pseudo)=c("p00","p01","p11","p12","p02_0","p02_1","p02","F01","F0.","RM")
  }
  else{
    colnames(pseudo)=c("LE.00","LE.0.","LE.01","LE.11","LTR")
  }
  return(pseudo)
}
