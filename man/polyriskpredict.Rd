\name{polyriskpredict}
\alias{polyriskpredict}
\title{
Polygenic Risk Prediction Function
}
\description{
This function allows to make polygenic risk prediction at given sample size with SNPs included at optimum p-value threshold or genome-wide significance level (5e-8). 
}
\usage{
polyriskpredict(N, Ps, Sig2s, M=1070777, M1, type="optimum", alp.GWAS=5*10^(-8), k.fold=3:5)
}
\arguments{
\item{\code{N}}{Effective sample size, at which the predictive performance is calculated. It must be a scalar.}
\item{\code{Ps}}{a vector of two mixture weights in the effect-size distribution for susceptibility SNPs. The sum should be equal to one. In a case where the one-component normal distribution is assumed for the effect-size distribution, then Ps=c(0.5,0.5).}
\item{\code{Sig2s}}{a vector of two component-specific variances for susceptibility SNPs. In a case where the one-component normal distribution is assumed for the effect-size distribution, then the two elements of the vector should be set to the same value.}
\item{\code{M}}{the number of independent set of susceptibility SNPs with the default of 1070777.}
\item{\code{M1}}{the estimated number of susceptibility SNPs.}
\item{\code{type}}{which threshold should be applied. Either "optimum" (default) or "GWAS" can be chosen.}
\item{\code{alp.GWAS}}{the (scalar) genome-wide significance level, if type=="GWAS". The default value is 5*10^(-8).}
\item{\code{k.fold}}{a vector of multiplicative constants, at which or higher risk then the average population risk, proportions of population and cases are calculated. The default value is set at 3:5.}
}
\details{
This function allows to make polygenic risk prediction at given sample size with SNPs included at optimum p-value threshold or genome-wide significance level (5e-8). (corresponding with Fig. 4 in Zhang et. al 2018)
}
\value{
\item{\code{alpha}}{at which the predictive performance measures are calculated.}
\item{\code{AUC}}{area under the curve.}
\item{\code{PPI}}{proportion of population identifiable to be at k-fold or higher risk than the average population risk.}
\item{\code{PCI}}{proportion of cases identifiable to be at k-fold or higher risk than the average population risk.}
}
\references{
Zhang, Yan, et al. "Estimation of complex effect-size distributions using summary-level statistics from genome-wide association studies across 32 complex traits" Nature Genetics (2018).
}
\author{
Yan Zhang <yandorazhang@gmail.com>
}
\examples{
### Schizophrenia example at a sample size of 10^6
m<-1070777
m1<-19294
ns<-10^6
ps<-c(0.5, 0.5)
sig2s<-c(0.000108801, 0.000108801)

## at an optimal threshold
polyriskpredict(N=ns, Ps=ps, Sig2s=sig2s, M=m, M1=m1, type="optimum", k.fold=3:5)
## at the GWAS significance
polyriskpredict(N=ns, Ps=ps, Sig2s=sig2s, M=m, M1=m1, type="GWAS", alp.GWAS=5*10^(-8), k.fold=3:5)
}