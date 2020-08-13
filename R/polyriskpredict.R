#' Illustration of polyriskpredict()
#'
#' This function allows to make polygenic risk prediction at given sample size with SNPs included at optimum p-value threshold or genome-wide significance level (5e-8). 
#' 
#' @param N sample size, at which the predictive performance is calculated. It must be a scalar.
#' @param Ps a vector of two mixture weights in the effect-size distribution for susceptibility SNPs. The sum should be equal to one.
#' @param Sig2s a vector of two component-specific variances for susceptibility SNPs. In a case where the one-component 
#' normal distribution is assumed for the effect-size distribution, then the two elements of the vector should be set to the same value.
#' @param M the number of independent set of susceptibility SNPs with the default of 1070777. 
#' @param M1 the estimated number of susceptibility SNPs. 
#' @param type which threshold should be applied. Either "optimum" (default) or "GWAS" can be chosen.
#' @param alp.GWAS the (scalar) genome-wide significance level, if type=="GWAS". The default value is 5*10^(-8).
#' @param k.fold a vector of multiplicative constants, at which or higher risk then the average population risk, 
#' proportions of population and cases are calculated. The default value is set at 3:5.
#' @keywords 
#' @export
#' @examples polyriskpredict(N, Ps, Sig2s, M=1070777, M1, type="optimum", alp.GWAS=5*10^(-8), k.fold=3:5)

polyriskpredict <- function(N, Ps, Sig2s, M=1070777, M1, type="optimum", alp.GWAS=5*10^(-8), k.fold=3:5){
 
  ### functions used within polyriskpredict
  mu.func<-function(N, c.alp.half, Ps, Sig2s,M, M1){
    # this function calculates the expeced value of RN, defined in Chatterjee et al (2013, NG)
    # by assuming a mixture of two normal distributions with mean zero and two different variances
    # for the non-null regression coefficients, betas.
    # In a case where a single normal distribution is appropriate, then two variances are set to be the same,
    # and the sum of weights, Ps, is set to be 1. 
    cmp1<-M1*(Ps[1]*comp1(N, c.alp.half, Sig2s[1])+Ps[2]*comp1(N, c.alp.half, Sig2s[2]))
    cmp2<-M1*(Ps[1]*comp2(N, c.alp.half, Sig2s[1])+Ps[2]*comp2(N, c.alp.half, Sig2s[2]))
    cmp3<-M1*c.alp.half/N*(Ps[1]*comp3(N, c.alp.half, Sig2s[1])+Ps[2]*comp3(N, c.alp.half, Sig2s[2]))
    cmp4<-M1/sqrt(N)*(Ps[1]*comp4(N, c.alp.half, Sig2s[1])+Ps[2]*comp4(N, c.alp.half, Sig2s[2]))
    cmp5<-(M-M1)*(1-pnorm(c.alp.half))*2/N*(1+c.alp.half*2*dnorm(c.alp.half)/((1-pnorm(c.alp.half))*2))
    mus<-(cmp1+cmp4)/sqrt(cmp1+cmp2+cmp3+cmp4+cmp5)
    return(mus)
  }
  
  comp1<-function(N, c.alp.half, sig2.0){
    sig2.0*(1-pnorm(c.alp.half/sqrt(N),mean=0, sd=sqrt((1+N*sig2.0)/N))+pnorm(-c.alp.half/sqrt(N),mean=0, sd=sqrt((1+N*sig2.0)/N)))-
      N*sig2.0^2/(1+N*sig2.0)*(-c.alp.half/sqrt(N)*dnorm(-c.alp.half/sqrt(N),mean=0, sd=sqrt((1+N*sig2.0)/N))*2)
  }
  
  comp2<-function(N, c.alp.half, sig2.0){
    1/N*(1-pnorm(c.alp.half/sqrt(N),mean=0, sd=sqrt((1+N*sig2.0)/N))+pnorm(-c.alp.half/sqrt(N),mean=0, sd=sqrt((1+N*sig2.0)/N)))
  }
  
  comp3<-function(N, c.alp.half, sig2.0){
    (dnorm(-c.alp.half, mean=0, sd=sqrt(1+N*sig2.0))+dnorm(-c.alp.half, mean=0, sd=sqrt(1+N*sig2.0)))
  }
  
  comp4<-function(N, c.alp.half, sig2.0){
    sig2.0*sqrt(N)*c.alp.half/(1+N*sig2.0)*dnorm(c.alp.half, mean=0, sd=sqrt(1+N*sig2.0))*2
  }
  ### 

  ### The function definition begins here
  if (length(N)!=1) {stop("The length of N must be one.")}
  if (length(Sig2s)!=2) {stop("The length of Sig2s is not equal to two.")}
  if (length(Ps)!=2) {stop("The length of Ps is not equal to two.")}
  if (missing(M1)) {stop("M1 is missing!")}
  
  if (type=="optimum"){
    # to vectorize function mu.func over the argument c.alp.half
    mu.vec.func<-Vectorize(mu.func, "c.alp.half")
    c.seq<-exp(seq(log(10^(-8)), log(7), length=1000))	# generate candidate grid points for c.alp.half
    
    ## optimized threshold
    tmp.mus<-mu.vec.func(N=N, c.alp.half=c.seq, Ps=Ps, Sig2s=Sig2s, M=M, M1=M1)
    loc.max<-which.max(tmp.mus)
    
    if (loc.max==1){
      result.vec<-c(c.seq[loc.max], tmp.mus[loc.max])
    }else{
      result.vec<-unlist(optimize(mu.func, lower=c.seq[loc.max-1], upper=c.seq[loc.max+1], maximum=T, N=N, Sig2s=Sig2s, Ps=Ps, M=M, M1=M1, tol=10^(-5)))
    }
    result.vec[1]<-(1-pnorm(result.vec[1]))*2		# the first element is replaced with the optimal alpha level.
  }else{
    ## at the GWAS significance
    c.GWAS<-qnorm(1-alp.GWAS/2)
    result.vec<-c(alp.GWAS, mu.func(N=N, c.alp.half=c.GWAS, Ps=Ps, Sig2s=Sig2s, M=M, M1=M1))
  }
  # AUC
  AUC<-pnorm(result.vec[2]/sqrt(2)) 
  
  xi.k<-log(k.fold)/sqrt(result.vec[2])
  PPI<-1-pnorm(xi.k+sqrt(result.vec[2])/2)
  PCI<-1-pnorm(xi.k-sqrt(result.vec[2])/2)
  names(PPI)<-names(PCI)<-paste("at ", k.fold, "-fold or higher", sep="")
  
  list(alpha=result.vec[1], AUC=AUC, PPI=PPI, PCI=PCI)
}
