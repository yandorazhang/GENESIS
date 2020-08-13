#' Illustration of projection()
#'
#' This function allows to make future projections according to the fitted model.
#' @param est parameter estimates by fitting either 2-component model, i.e., (pic, sigmasq, a); or 3-component model, i.e., (pic, p1, sigmasq1, sigmasq2, a).
#' @param v covariance matrix of parameter estimates by fitting the 2- or 3-component model. 
#' @param n specifided future GWAS study sample size.
#' @param gwas.significance genome-wide significance level, by default it is 5e-8. 
#' @param tol tolerance accuracy vector for intgrate() function.
#' @param M total number of SNPs in the reference panel; by default, it is the total number of common SNPs in Hapmap3 reference panel, which is equal to 1070777. 
#' @param CI whether to calculate CI or not; by default, CI=FALSE.
#' @param nsim the total number of bootstrap samplers in order to calculate CI; by default, it is 1000.
#' @param CI.coverage coverage level of confidence interval; by default, it is 0.95, i.e., 95% CI. 
#' @param seeds numeric; random seeds used in simulation; by default, it is 123.
#' @keywords 
#' @export
#' @examples projection(est,v=NULL,n,gwas.significance=5e-8,tol=c(1e-12,1e-15),M=1070777,CI=FALSE,nsim=1000,CI.coverage=0.95,seeds=123)

projection <- function(est,v=NULL,n,gwas.significance=5e-8,tol=c(1e-12,1e-15),M=1070777,CI=FALSE,nsim=1000,CI.coverage=0.95,seeds=123){
  # within function
  pp <- function(est,n,gwas.significance=5e-8,tol=c(1e-12,1e-15),M=1070777){
    
    if(length(est)==3)components=2
    if(length(est)==5)components=3
    
    if(components==2){
      pic = est[1]
      sig = sqrt(est[2])
      den <- function(x){return(dnorm(x/sig)/sig )}
      herit0 <- pic*M*sig^2
    }
    
    if(components==3){
      pic = est[1]
      p1 = est[2]
      s1 = sqrt(est[3])
      s2 = sqrt(est[4])
      den <- function(x){return(p1 * dnorm(x/s1)/s1 + (1-p1)*dnorm(x/s2) /s2)}
      herit0 <- pic*M*(p1*est[3] + (1-p1)*est[4])
    }
    
    tem0 <- function(x){return(x^2*den(x))}

    c_gwsignificance = abs(qnorm(gwas.significance/2))
    pow <- function(x){return(1 - pnorm(c_gwsignificance - sqrt(n)*x) + pnorm(-c_gwsignificance - sqrt(n)*x) )}
    tem <- function(x){return(pow(x)*den(x))}
    Numdiscoveries = M*pic * integrate(tem, -Inf, Inf,rel.tol=tol[1], abs.tol=tol[2])[[1]]
    
    tem1 <- function(x){return(pow(x)*den(x)*x^2)}
    GVpercentage = M*pic * integrate(tem1, -Inf, Inf,rel.tol=tol[1], abs.tol=tol[2])[[1]]*100/herit0
    pheno.variance = GVpercentage*herit0/100
    
    return(list(Numdicoveries=Numdiscoveries, GVpercentage=GVpercentage, pheno.variance=pheno.variance,herit=herit0))
  }
  
  
  library(MASS)
  logest = log(est)
  
  if(CI==TRUE){
    set.seed((seeds))
    alpha = (1-CI.coverage)/2
    logv = diag(1/est)%*%v%*%diag(1/est)
    estmat = exp(mvrnorm(nsim,mu=logest,Sigma=logv))
    
    tem = pp(est,n,gwas.significance,tol,M)
    tem1 = apply(estmat, 1, function(t) {pp(t,n,gwas.significance,tol,M)})
    
    pest = tem$Numdicoveries;
    gvest = tem$GVpercentage;
    pheno.variance = tem$pheno.variance;
    herit = tem$herit
    
    re = unlist(lapply(tem1,function(t) t[1]))
    rere = apply(matrix(re,ncol=1),1,function(t) rbinom(1,size=M,prob=t/M))
    regv = unlist(lapply(tem1,function(t)t[2]))
    
    return(list(heritability = herit, 
                Numdiscoveries = c(pest,quantile(rere,alpha),quantile(rere,1-alpha)),
                pheno.variance = c(gvest,quantile(regv,alpha),quantile(regv,1-alpha))*herit/100, 
                GVpercentage = c(gvest,quantile(regv,alpha),quantile(regv,1-alpha))
                ))
  }
  
  if(CI==FALSE){
    tem = pp(est,n,gwas.significance,tol,M)
    pest = tem$Numdicoveries;
    gvest = tem$GVpercentage;
    pheno.variance = tem$pheno.variance; 
    herit = tem$herit;
    
    return(list(heritability = herit, Numdiscoveries = pest, pheno.variance=pheno.variance, GVpercentage = gvest))
  }
  
}


