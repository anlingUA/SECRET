###########################################################################################
#' Normalization
#'
#' @param count the dataset need to be normalized
#' @param method select which normalization method to use
#'
#' @return normalized dataset
#' @export
#'
normcount <- function(count,method){
  if (method=='cpm'){
    count1<-t(t(count)/colSums(count)*1e6)
  } else if (method=='logcpm') {
    count1<-log2(t(t(count+0.5)/(colSums(count)+1)*1e6))
  } else if (method=='countNorm') {
    count1 <- t(t(count)/colSums(count))*median(colSums(count))
  } else if (method=='quantile') {
    count1 <- NOISeq::uqua(count)
  } else if (method=='tmm'){
    count1=NOISeq::tmm(count)
  } else if (method=='tpm'){
    DGEobj.utils::convertCounts(as.matrix(count),unit='TPM',
                                geneLength= rowMeans(count),log = FALSE,normalize = "none")
  } else { # no normalize
    count1=as.matrix(count)
  }
}

###########################################################################################
#' Semi-reference based cell type deconvolution method
#'
#' @param bulk raw filtered bulk data
#' @param reference raw filtered cell type reference data
#' @param withUnknown whether include unknown cell type
#' @param w weight for each genes
#' @param yNorm normalization method
#' @param bNorm normalization method
#' @param exprs_mode whether to calculate gene expression, this step takes long time so the default is F
#' @param bulkfull raw full bulk data
#'
#' @return estimated cell type proportion or cell type gene expression if exprs_mode=T
#' @export
#'

SECRET <-function(bulk,reference,withUnknown=T,w=w,yNorm,bNorm,exprs_mode=F,bulkfull=NULL){
  bulk1=normcount(bulk,yNorm)
  reference1=normcount(reference,bNorm)
  nct <- ncol(reference1)
  func <- function(P){
    return(sum(w*abs(reference1%*%P  - y), na.rm = TRUE))
  }
  Pinit <- rep(0,nct)
  if (withUnknown){
    estp=matrix(0,nrow=ncol(bulk1),ncol=ncol(reference1)+1)
    rownames(estp)=colnames(bulk1);colnames(estp)=c(colnames(reference1),'unknown')
    for (i in 1:ncol(bulk1)){
      y=bulk1[,i]
      hin <- function(P){# Creating inequality constraint function
        P_mat <- matrix(P, nrow = nct)
        c(1-colSums(P_mat),P)}
      aug_res <- alabama::auglag(Pinit, func, hin = hin,control.outer = list(NMinit=F,trace = F,itmax = 100))
      estp[i,1:ncol(reference1)]=aug_res$par
      estp[i,ncol(estp)]=1-sum(aug_res$par)}
  } else { #no missing
    estp=matrix(0,nrow=ncol(bulk1),ncol=ncol(reference1))
    rownames(estp)=colnames(bulk1);colnames(estp)=colnames(reference1)
    for (i in 1:ncol(bulk1)){
      y=bulk1[,i]
      hin <- function(P){
        P_mat <- matrix(P, nrow = nct)
        c(P)}
      heq <- function(P){ # equality constraint
        P_mat <- matrix(P, nrow = nct)
        c(1-colSums(P_mat))}
      aug_res <- alabama::auglag(Pinit, func, hin = hin,heq=heq,control.outer = list(NMinit=F,trace = F,itmax = 100))
      estp[i,]=aug_res$par
    }
  }
  if (exprs_mode){
    ##calculate gene expression for each gene
    bulkfull1=normcount(bulkfull,yNorm)
    Einit <- rep(0,ncol(estp))
    # calulate gene expression
    hinE <- function(E){# inequality constraint for proportion
      E_mat <- matrix(E, ncol = ncol(estp))
      c(E)}
    aug_E_list=list()
    library(parallel)
    ncores <- detectCores()
    aug_E_list <- mclapply(1:nrow(bulkfull1), function(i) {
      y <- bulkfull1[i,]
      func2 <- function(E){
        return(sum(abs(E%*% t(estp)  - y), na.rm = TRUE))
      }
      alabama::auglag(Einit, func2, hin = hinE, control.outer = list(NMinit = T, trace = F))
    }, mc.cores = ncores)
    gexp<- sapply(aug_E_list, function(x) x$par)
    rownames(gexp)=colnames(estp);colnames(gexp)=rownames(bulkfull1)
    ## gether information
    list(abs(estp),abs(gexp))
  } else {list(abs(estp))}
}
