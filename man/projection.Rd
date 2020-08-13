\name{projection}
\alias{projection}
\title{
Future Projection Function
}
\description{
Function to make future projections about the number of expected discoveries and associated heritability explained in future studies based on analytic formulas according to the fitted model. 
}
\usage{
projection(est,v=NULL,n,gwas.significance=5e-8,tol=c(1e-12,1e-15),M=1070777,CI=FALSE,nsim=10000,CI.coverage=0.95,seeds=123)
}
\arguments{
  \item{\code{est}}{parameter estimates by fitting either 2-component model, i.e., (pic, sigmasq, a); or 3-component model, i.e., (pic, p1, sigmasq1, sigmasq2, a).}
  \item{\code{v}}{covariance matrix of parameter estimates by fitting the 2- or 3-component model.}
  \item{\code{n}}{future GWAS sample size.}
  \item{\code{gwas.significance}}{genome-wide significance level, by default it is 5e-8.}
  \item{\code{tol}}{tolerance accuracy vector for intgrate() function. By default, it is c(1e-12,1e-15).}
  \item{\code{M}}{total number of SNPs in the reference panel. By default, it is the total number of common SNPs in Hapmap3 reference panel, which is equal to 1070777. }
\item{\code{CI}}{logic; indicating whether to calculate confidence interval or not; by default, CI=FALSE.}
  \item{\code{nsim}}{the total number of bootstrap samplers in order to calculate CI; by default, it is 1000.}
  \item{\code{CI.coverage}}{coverage level of confidence interval; by default, it is 0.95, i.e., 95\% CI.}
  \item{\code{seeds}}{numeric; random seeds used in simulation; by default, it is 123.}
}
\details{
This function allows to make future projections about the number of expected discoveries and associated heritability explained in future studies based on analytic formulas according to the fitted model.  (corresponding with Fig. 3  in Zhang et. al 2018)
}
\value{
A list of future projections with components
\item{\code{heritability}}{total heritability.}
\item{\code{Numdicoveries}}{predicted number of significant  discoveries  at GWAS significance level "gwas.significance" in  future GWAS study with sample size n. }
\item{\code{pheno.variance}}{predicted phenotypic variance that can be explained by those discovered significant SNPs  at GWAS significance level "gwas.significance"  in  future GWAS study with sample size n.}
\item{\code{GVpercentage}}{predicted percentage  of genetic variance that can be explained  by those discovered significant SNPs  at GWAS significance level "gwas.significance"  in  future GWAS study with sample size n, i.e., GVpercentage = pheno.variance/heritability*100 }
}
\references{
Zhang, Yan, et al. "Estimation of complex effect-size distributions using summary-level statistics from genome-wide association studies across 32 complex traits" Nature Genetics (2018).
}
\author{
Yan Zhang <yandorazhang@gmail.com>
}
\examples{
v = matrix(c(1.327856e-06, -1.131049e-05, 6.912489e-10, -2.901301e-09, -9.388865e-11, -8.568269e-06, 3.380985e-04, -1.543479e-07, 1.771036e-09,8.989668e-10, -1.846216e-09, -1.353542e-07,  1.517545e-10,  1.166563e-11, 1.901686e-14, -3.068113e-09, 9.520417e-09,  5.633801e-12, 9.427492e-12, 1.220433e-13, -9.111521e-11, 1.063848e-09, -1.526281e-13, 1.042123e-13, 1.406410e-14), 5,5)
projection(est=c(8.899809e-03, 9.476025e-02, 1.458650e-04, 2.227118e-05, 1.567643e-06),v, n=253288, CI=TRUE);
}
