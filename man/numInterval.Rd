\name{numInterval}
\alias{numInterval}
\title{
Function to Calculate Number of SNPs with Effect Size Absolute Value Falling in an Interval
}
\description{
This function calculates the number of SNPs with the absolute value of effect size falling into some interval according to the fitted mixture model.
}
\usage{
numInterval(lower,upper,est,M=1070777)
}
\arguments{
  \item{\code{lower}}{the lower bound of the interval.}
  \item{\code{upper}}{the upper bound of the interval.}
  \item{\code{est}}{the parameter estimates by fitting either 2-component model, i.e., (pic, sigmasq, a); or 3-component model, i.e., (pic, p1, sigmasq1, sigmasq2, a).}
  \item{\code{M}}{total number of SNPs in the reference panel. By default, it is the total number of common SNPs in Hapmap3 reference panel, which is equal to 1070777. }
}
\details{
This function calculates the number of SNPs with the absolute value of effect size falling into some interval according to the fitted mixture model. (corresponding with Supplementary Table 7 in Zhang et. al 2018)
}
\value{
The number of SNPs with effect size falling into some interval (lower, upper) according to the fitted mixture model.
}
\references{
Zhang, Yan, et al. "Estimation of complex effect-size distributions using summary-level statistics from genome-wide association studies across 32 complex traits" Nature Genetics (2018).
}
\author{
Yan Zhang <yandorazhang@gmail.com>
}
\examples{
numInterval(0.005,Inf,est=c(8.899809e-03, 9.476025e-02, 1.458650e-04, 2.227118e-05, 1.567643e-06))
}
