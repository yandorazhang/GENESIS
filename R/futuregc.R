#' Illustration of futuregc()
#'
#' This function allows to predict future genomic control (GC) factor through simulations given the GWAS study sample size.
#' @param est parameter estimates by fitting either 2-component model, i.e., (pic, sigmasq, a); or 3-component model, i.e., (pic, p1, sigmasq1, sigmasq2, a).
#' @param n specifided future GWAS sample size.
#' @param nsim total number of simulations; by default, it is 1. 
#' @param M total number of SNPs in the reference panel; by default, it is the total number of common SNPs in Hapmap3 reference panel, which is equal to 1070777. 
#' @param seeds numeric random seeds used in simulation; by default, it is 123. 
#' @keywords 
#' @export
#' @examples futuregc(est,n,nsim=1,M=1070777,seeds=123)

futuregc <- function(est,n,nsim=1,M=1070777,seeds=123){
  
  lambdaGC = rep(0,nsim)
  
  # load the LD structure
  data("TaggingSNPinx"); data("pairwiseLD"); data(LDwindow1MB_cutoff0.1); SNPrsID = dataLD$SNPname; K=M;
  
  if(length(est)==3) components=2
  if(length(est)==5) components=3
  
  if(components==2){
    pic = est[1]; sigmasq = est[2]; a= est[3]; 
    if(a<0) a = 0
    
    for(iter in 1:nsim){
      # generate the joint effect size according to the fitted distribution
      set.seed(iter*seeds)
      z = rbinom(K, size=1, prob=pic); nonzero = sum(z)
      betajoint = rep(0, K)
      betajoint[which(z==1)] = rnorm(nonzero,mean=0,sd=sqrt(sigmasq))
      
      betamarginal = rep(0,K);
      for(k in 1:K){
        betamarginal[k] = crossprod(betajoint[unlist(TaggingSNPinx[k])], sqrt(unlist(pairwiseLD[k])))+rnorm(1,mean=0,sd=sqrt(a))
      }
      
      data(list=paste0("error_iter",iter))
      dftem = data.frame(cbind(error,SNP)); dftem0= data.frame(SNPrsID)
      dfmerge = merge(dftem0, dftem, by.x="SNPrsID", by.y="SNP",sort=F)
      betahat = betamarginal + as.numeric(as.character(dfmerge$error))/sqrt(n)
      varbetahat = rep(1/n,K)
      lambdaGC[iter] = median( (betahat/sqrt(varbetahat))^2)/qchisq(0.5,1)
    }
  }
  
  
  if(components==3){
    pic = est[1]; p1 = est[2]; sig1 = est[3]; sig2 = est[4];a=est[5]
    if(a<0)a=0
    
    for(iter in 1:nsim){
      # generate the joint effect size according to the fitted distribution
      set.seed(iter*123)
      z = sample(x=c(1,2,0), size=K,replace=T,prob=c(pic*p1, pic*(1-p1), 1-pic))
      nz1 = sum(z==1)
      nz2 = sum(z==2)
      betajoint = rep(0, K)
      betajoint[which(z==1)] = rnorm(nz1,mean=0,sd=sqrt(sig1))
      betajoint[which(z==2)] = rnorm(nz2,mean=0,sd=sqrt(sig2))
      
      betamarginal = rep(0,K);
      for(k in 1:K){
        betamarginal[k] = crossprod(betajoint[unlist(TaggingSNPinx[k])], sqrt(unlist(pairwiseLD[k])))+rnorm(1,mean=0,sd=sqrt(a))
      }
      
      data(list=paste0("error_iter",iter))
      dftem = data.frame(cbind(error,SNP)); dftem0= data.frame(SNPrsID)
      dfmerge = merge(dftem0, dftem, by.x="SNPrsID", by.y="SNP",sort=F)
      betahat = betamarginal + as.numeric(as.character(dfmerge$error))/sqrt(n)
      varbetahat = rep(1/n,K)
      lambdaGC[iter] = median( (betahat/sqrt(varbetahat))^2)/qchisq(0.5,1)
    }
  }
  return(mean(lambdaGC))
}
