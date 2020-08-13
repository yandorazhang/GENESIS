\name{qqplot.plot}
\alias{qqplot.plot}
\title{
Function to Create QQ Plot
}
\description{
This function allows to create QQ plot comparing observed distributions of association statistics against those expected under the fitted models. 
}
\usage{
qqplot.plot(qqplotdata,seq_inx=1,qqplot.axis=10)
}
\arguments{
\item{\code{qqplotdata}}{
A list containing: 
(i) QQdata, a dataframe needed to generate the QQ plot, i.e., a dataframe containing 4 columns, "log10(observed pvalues)", "expected log10(fitted pvalues)", "lower bound of log10(fitted pvalues)","upper bound of log10(fitted pvalues)".
(ii) observedlambda, lambda genomic control (GC) factor for observed summary GWAS statistics.
(iii) meanEXPlambda,  expected lambda GC factor for simulated summary GWAS statistics under the fitted model.
(iv) lowEXPlambda, "qqplotCI.coverage" CI lower bound of lambda GC factor for simulated summary GWAS statistics under the fitted model. 
(v) highEXPlambda, "qqplotCI.coverage" CI upper bound of lambda GC factor for simulated summary GWAS statistics under the fitted model.
}
\item{\code{seq_inx}}{numeric; QQdata dataframe will be thinned every seq_inx. }
\item{\code{qqplot.axis}}{numeric; the x- and y-axis limits are set from 0 to "qqplot.axis" in the QQ plot. By default, it is 10.}
}
\details{
This function allows to create QQ plot. (corresponding with Fig. 1 in Zhang et. al 2018)
}
\value{
QQ plot. 
}
\references{
Zhang, Yan, et al. "Estimation of complex effect-size distributions using summary-level statistics from genome-wide association studies across 32 complex traits" Nature Genetics (2018).
}
\author{
Yan Zhang <yandorazhang@gmail.com>
}
\seealso{
\code{\link{qqplotdata.simu}}
}
\examples{
data(heightGWAS)
qqplotdata <- qqplotdata.simu(heightGWAS, est=c(9.583307e-03,8.562964e-02,1.487684e-04,2.086576e-05,1.498790e-06))
qqplot.plot(qqplotdata)
}
