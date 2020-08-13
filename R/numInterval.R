#' Illustration of numInterval()
#'
#' This function allows to calculate the number of SNPs with the absolute value of effect size falling into some interval according to the fitted mixture model.
#' @param lower the lower bound of the interval.
#' @param upper the upper bound of the interval.
#' @param est parameter estimates by fitting either 2-component model, i.e., (pic, sigmasq, a); or 3-component model, i.e., (pic, p1, sigmasq1, sigmasq2, a).
#' @param M total number of SNPs in the reference panel; by default, it is the total number of common SNPs in Hapmap3 reference panel, which is equal to 1070777. 
#' @keywords 
#' @export
#' @examples numInterval(lower,upper,est,M=1070777)

numInterval <- function(lower,upper,est,M=1070777){
  
  if(length(est)==3) components=2
  if(length(est)==5) components=3
  
  if(components==2){
    pic = est[1]
    sig = sqrt(est[2])
    den <- function(x){return(dnorm(x/sig)/sig )}
    cdf <- function(x){
      if(x<0) {res = pic*pnorm(x/sig)}
      if(x>=0) {res = pic*pnorm(x/sig) + 1-pic}
      return(res)
    }
  }
  
  if(components==3){
    pic = est[1]
    p0 = est[2]
    s1 = sqrt(est[3])
    s2 = sqrt(est[4])
    den <- function(x){return(p0 * dnorm(x/s1)/s1 + (1-p0)*dnorm(x/s2) /s2)}
    cdf <- function(x){
      if(x<0) {res = pic*(p0*pnorm(x/s1) + (1-p0)*pnorm(x/s2))}
      if(x>=0) {res = pic*(p0*pnorm(x/s1) + (1-p0)*pnorm(x/s2))  + 1-pic}
      return(res)
    }
  }
  
  return(M*2*(cdf(upper) - cdf(lower)))
  
}
