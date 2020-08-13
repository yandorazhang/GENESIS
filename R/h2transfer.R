#' Illustration of h2transfer
#'
#' This function allows to transfer the heritability in log-odds-ratio (frailty) scale to other scale (observed and liability scale).
#' @param h2.log heritability in log-odds-ratio (frailty) scale. 
#' @param se.h2.log standard error of heritability in log-odds-ratio scale.
#' @param population.prevalence population prevalence.
#' @param sample.prevalence sample prevalence. 
#' @keywords 
#' @export
#' @examples h2transfer(h2.log, se.h2.log, population.prevalence, sample.prevalence)

h2transfer <- function(h2.log, se.h2.log, population.prevalence, sample.prevalence){
  
  P <- sample.prevalence
  K <- population.prevalence
  z <- dnorm(qnorm(1-K))
  
  h2.observed <- h2.log*P*(1-P)
  se.h2.observed <- se.h2.log*P*(1-P)

  h2.liability <- h2.observed*(K^2)*((1-K)^2)/((z^2)*P*(1-P))
  se.h2.liability <- se.h2.observed*(K^2)*((1-K)^2)/((z^2)*P*(1-P))
  
  return(list(h2.observed=h2.observed, se.h2.observed=se.h2.observed,
              h2.liability=h2.liability,se.h2.liability=se.h2.liability))
}