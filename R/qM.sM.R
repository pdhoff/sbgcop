#' Matrix Quantiles
#' 
#' Computes quantiles along the third dimension of a 3-d array.
#' 
#' 
#' @param sM an m x n x s array
#' @param quantiles quantiles to be computed
#' @return an array of dimension m x n x l, where l is the length of
#' \code{quantiles}
#' @author Peter Hoff
#' @keywords array multivariate
#'
#' fit<-sbgcop.mcmc(swiss)
#' qC<-qM.sM(fit$C.psamp)
#' qR<-qM.sM(sR.sC(fit$C.psamp))
#'
"qM.sM" <-
function(sM,quantiles=c(0.025,.5,.975)) {

p1<-dim(sM)[1]
p2<-dim(sM)[2]
s <-dim(sM)[3]

qM<-array(dim=c(p1,p2,length(quantiles)))
dimnames(qM)<-list(dimnames(sM)[[1]],dimnames(sM)[[2]], 
              paste(quantiles*100,rep("% quantile",length(quantiles)),sep=""))

for(l in 1:length(quantiles)) {
qM[,,l]<-apply(sM,c(1,2),quantile,prob=quantiles[l],na.rm=TRUE)
                               }
qM
}

