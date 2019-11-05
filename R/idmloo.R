

#' R function to estimate n illness-death models for leave-one-out data
#'
#' @param object an illness-death model estimate from SmoothHazard package
#' @param data data.frame used to estimate object
#' @param parallel Logical determinating parallel (=TRUE) or sequential (=FALSE) execution
#' @param cpus Numerical amount of CPUs requested for the cluster. If not set, values are total cores of the machines-1.
#' @param type Type of cluster. Can be 'SOCK', 'MPI', 'PVM' or 'NWS'. Default is 'SOCK' if parallel=TRUE.
#' @param ... others parameters
#'
#' @return a list of length N+1 with one indice for leave-one-out of the data and the last
#' the original estimate of object
#' @export
#'
#' @examples
#' data("Paq1000")
#' fit.splines <-  idm(formula02=Hist(time=t,event=death)~1,
#' formula01=Hist(time=list(l,r),event=dementia)~1,
#' formula12=Hist(time=t,event=death)~1,
#' method="Splines", kappa=c(200000,200000,20000),
#' data=Paq1000)
#' fit.loo.sp=idmloo(object=fit.splines, data=Paq1000, parallel = FALSE, cpus = NULL, type = NULL)
#'
#'
idmloo<-function(object, data,parallel=FALSE, cpus=NULL, type=NULL,...){
  #test for good utilisation
  if (class(object)!="idm")stop("Argument 'object' must be an illness-death  model from SmoothHazard")
  if(sum(object$NC)!=0)stop("Argument 'object' must have no covariate")
  if(!is.logical(parallel))stop("Argument 'parallel' must be a logical")

  if(missing(data)) stop("Need a data frame.")
  if(class(data)!="data.frame")stop("Argument 'data' must be a data.frame")
  if(dim(data)[1]!=object$N)stop("Argument 'data' is not the same used to fit 'object'" )

  truncated <- nchar(attr(object$responseAbs,"entry.type"))>1
  if(truncated==TRUE)stop("Left truncation is not available yet")


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
  N=object$N
  DATA_LOO=vector("list",N)
  for (i in 1:N){
    DATA_LOO[[i]]=data[-i,]
  }

  if(object$method=="Splines"){
    fit_loo=function(data, object){
      transition01=formula(object$terms$Formula01)
      transition02=formula(object$terms$Formula02)
      transition12=formula(object$terms$Formula12)

      fit<-try(idm(formula01=transition01,
                   formula02=transition02,
                   formula12=transition12,
                   data=data, conf.int=FALSE, method="Splines",CV=FALSE,
                   kappa=object$kappa,
                   n.knots=c(object$nknots01,object$nknots02,object$nknots12)))
      if (class(fit)=="try-error"){
        fit<-try(idm(formula01=transition01,
                     formula02=transition02,
                     formula12=transition12,
                     data=data, conf.int=FALSE, method="Splines",CV=TRUE,
                     kappa=object$kappa,
                     n.knots=c(object$nknots01,object$nknots02,object$nknots12)))
        if(class(fit)=="try-error"){
          fit<-try(idm(formula01=transition01,
                       formula02=transition02,
                       formula12=transition12,
                       data=data, conf.int=FALSE, method="Splines",CV=TRUE,
                       kappa=object$kappa,
                       n.knots=c(object$nknots01+2,object$nknots02+2,object$nknots12+2)))
          if(class(fit)=="try-error"){
            fit<-"pb"
          }
        }
      }
      return(fit)

    }

    if(parallel==FALSE){
      LA_LISTE_LOO=lapply(DATA_LOO,fit_loo, object=object)
    }
    if(parallel==TRUE){
      if(sfIsRunning()==FALSE){
        sfInit(parallel = TRUE, cpus=cpus, type=type)
        sfLibrary(SmoothHazard)
        LA_LISTE_LOO=sfLapply(DATA_LOO,fit_loo, object=object)
        sfStop()
      }
      else{
        sfLibrary(SmoothHazard)
        LA_LISTE_LOO=sfLapply(DATA_LOO,fit_loo, object=object)
      }
    }
  }
  if(object$method=="Weib"){
    fit_loo=function(data, object){
      data=as.data.frame(data)
      transition01=formula(object$terms$Formula01)
      transition02=formula(object$terms$Formula02)
      transition12=formula(object$terms$Formula12)

      fit<-try(idm(formula01=transition01,
                   formula02=transition02,
                   formula12=transition12,
                   data=data, conf.int=FALSE, method="Weib"))
      if(class(fit)=="try-error"){
        fit<-"pb"
      }
      return(fit)
    }
    if(parallel==FALSE){
      LA_LISTE_LOO=lapply(DATA_LOO,fit_loo, object=object)
    }
    if(parallel==TRUE){
      if(sfIsRunning()==FALSE){
        sfInit(parallel = TRUE, cpus=cpus, type=type)
        sfLibrary(SmoothHazard)
        LA_LISTE_LOO=sfLapply(DATA_LOO,fit_loo, object=object)
        sfStop()
      }
      else{
        sfLibrary(SmoothHazard)
        LA_LISTE_LOO=sfLapply(DATA_LOO,fit_loo, object=object)
      }
    }
  }
  names(LA_LISTE_LOO)=paste0("loo.",seq(1:N))
  for (i in 1:N){
    if(class(LA_LISTE_LOO[[paste0("loo.",i)]])=="try-error"){
      print(i)
      # warning(paste0("Model for leave-one-out number ",i," did not converged.
      # Please change it mannualy because the package is under development."))
    }
  }
  LA_LISTE_LOO$loo.all=object
  class(LA_LISTE_LOO)="idmloo"
  return(LA_LISTE_LOO)
}


