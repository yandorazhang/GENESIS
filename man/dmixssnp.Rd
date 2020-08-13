\name{dmixssnp}
\alias{dmixssnp}
\title{
Mixture Normal Probability Density for the Susceptibility SNP (sSNP)
}
\description{
This function allows to get the (mixture) normal probability density function for the effect-size distribution of susceptibility SNP (sSNP) based on the fitted 2- or 3-component model.
}
\usage{
dmixssnp(x,est)
}
\arguments{
  \item{\code{x}}{value of quantile.}
  \item{\code{est}}{parameter estimates by fitting either 2-component model, i.e., (pic, sigmasq, a); or 3-component model, i.e., (pic, p1, sigmasq1, sigmasq2, a).}
}
\details{
This function allows to get the (mixture) normal probability density function for the effect-size distribution of susceptibility SNP (sSNP) based on the fitted 2- or 3-component model.  It is a single normal distribution for 2-component model; two-component mixture normal distribution for 3-component model. Distributions with fatter tails imply that the underlying traits have relatively greater numbers of susceptibility SNPs with larger effects. (corresponding with Fig. 2 in Zhang et. al 2018)
}
\value{
The density.
}
\references{
Zhang, Yan, et al. "Estimation of complex effect-size distributions using summary-level statistics from genome-wide association studies across 32 complex traits" Nature Genetics (2018).
}
\author{
Yan Zhang <yandorazhang@gmail.com>
}
\examples{
x_seq <- seq(-0.02,0.02,length.out = 1000); 
y_seq <- apply(matrix(x_seq,ncol=1),1,function(t) dmixssnp(t,est=c(8.899809e-03, 9.476025e-02, 1.458650e-04, 2.227118e-05, 1.567643e-06)));
plot(x_seq, y_seq, type="l"); 
}
