\name{h2transfer}
\alias{h2transfer}
\title{
Function to  Transfer Heritability in Log-odds-ratio Scale to Other Scale
}
\description{
This function allows to transfer the heritability in log-odds-ratio (frailty) scale to other scale (observed and liability scale). 
}
\usage{
h2transfer(h2.log, se.h2.log, population.prevalence, sample.prevalence)
}
\arguments{
  \item{\code{h2.log}}{heritability in log-odds-ratio scale.}
  \item{\code{se.h2.log}}{standard error of heritability in log-odds-ratio scale.}
  \item{\code{population.prevalence}}{population prevalence.}
  \item{\code{sample.prevalence}}{sample prevalence.}
}
\value{
Heritability and its standard error in observed scale and liability scale. 
}
\author{
Yan Zhang <yandorazhang@gmail.com>
}
\examples{
h2transfer(.351, .075, 0.044, 17/(17+37))
}