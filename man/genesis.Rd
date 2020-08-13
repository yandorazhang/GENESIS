\name{genesis}
\alias{genesis}
\title{
GENetic Effect-Size distribution Inference from Summary-level data
}
\description{
Estimation of complex effect-size distributions using summary-level statistics from genome-wide association studies. 
}
\usage{
genesis(summarydata, filter=FALSE,modelcomponents=2, cores=10, LDcutoff=0.1, LDwindow=1, M=1070777,c0=10, BIC.gamma=0.5,print=TRUE, printfreq=10, starting=NA, startingpic=NA, tolerance=NA,qqplot=TRUE, qqplotCI.coverage=0.8, qqplot.name="", qqplot.axis=10, qqplot.nsim=100,summaryGWASdata.save=FALSE, qqplotdata.save=FALSE,herit.liability=FALSE, sample.prevalence=NA, population.prevalence=NA,stratification=TRUE, seeds=123)
}
\arguments{
  \item{\code{summarydata}}{summay-level GWAS data, containing 3 columns: 
  SNP (SNP rsID), 
  Z (GWAS test z-statistic), 
  N (GWAS study sample size which can be different for different SNPs)}
  \item{\code{filter}}{logical; if TRUE, the input summary data will be filtered.}
  \item{\code{modelcomponents}}{2 or 3, indicating fitting a 2-component or 3-component model.}
  \item{\code{cores}}{number of CPU threads in parallel computing.}
  \item{\code{LDcutoff}}{a number from (0.05, 0.1, 0.2); indicating LD score is calculated based on the particular r^2 cutoff. By default, it is 0.1.}
   \item{\code{LDwindow}}{a number from (0.5, 1, 2); indicating LD score is calculated based on the particular window size (MB). By default, it is 1 MB.}
  \item{\code{M}}{total number of SNPs in the reference panel; by default, it is the total number of common SNPs in Hapmap3 reference panel, i.e., 1070777. }
  \item{\code{c0}}{an assumed maximum number of underlying susceptibility SNPs tagged by any individual GWAS marker. By default, c0 is set at 10.}
  \item{\code{BIC.gamma}}{numeric; tuning parameter when calculating BIC, in (0,1). By default, BIC.gamma is set at 0.5.}
  \item{\code{print}}{logical; if TRUE, the EM algorithm iteration details will be output.}
  \item{\code{printfreq}}{a number indicating every "printfreq" iterations of EM algorithms,  results will be output.}
  \item{\code{starting}}{the starting values for the model. For 2-component model, the starting value of (pic, sigmasq, a); for 3-component model, (pic, p1, sigmasq1, sigmasq2, a).}
  \item{\code{startingpic}}{the starting value for pic when starting==NA.}
  \item{\code{tolerance}}{the accuracy of the tolerance. For 2-component model, it is a 6-dim vector with tolerance value for (pic,sigmasq,a,llk,maxEM,steps). For 3-component model, it is a 8-dim vector with tolerance value for (pic,p1,sigmasq1,sigmasq2,a,llk,maxEM,steps).}
  \item{\code{qqplot}}{logical; if TRUE, the QQ plot will be ploted into a pdf file.}
\item{\code{qqplotCI.coverage}}{the coverage rate of confidence band in the QQ plot. By default, it is 0.8.}
  \item{\code{qqplot.name}}{the name of the QQ plot pdf files.}
\item{\code{qqplot.axis}}{numeric; the x- and y-axis limits is set from 0 to "qqplot.axis" in the QQ plot. By default, it is 10.}
  \item{\code{qqplot.nsim}}{the total number of simulations to generate the expected p-values to get the QQ plot. By default, it is 100.}
  \item{\code{summaryGWASdata.save}}{logical; if TRUE, the filtered summary GWAS data after merged with the LD information will be saved as a dataframe.}
  \item{\code{qqplotdata.save}}{logical; if TRUE, the simulated data  to generate the QQ plot will be saved.}
  \item{\code{herit.liability}}{logical; if TRUE, the heritability in log-odds-ratio scale will be transferred to liability scale.}
  \item{\code{sample.prevelance}}{sample prevelance for the disease trait.}
  \item{\code{population.prevelance}}{population prevelance for the disease trait.}
  \item{\code{stratification}}{logical; if TRUE, the population stratification effect is considered.}
\item{\code{seeds}}{numeric; random seeds used in simulation; by default, it is 123.}
}
\details{
If filter==TRUE, the input summary data will be filtered by 3 steps: (i) If sample size varies from SNP to SNP, remove SNPs with an effective sample size less than 0.67 times the 90th percentile of sample size. (ii) Remove SNPs within the major histocompatibility complex (MHR) region; filter SNPs to Hapmap3 SNPs. (iii) Remove SNPs with extremely large effect sizes (chi^2 > 80).

This function allows to estimate effect size distribution by fitting a 2- or 3-component model based on summary-level statistics as well as external linkage disequilibrium (LD) information. Unless users specify their own LD reference panel, the  default LD reference panel provided in this package  is 1070777 common SNPs in Hapmap3. The LD score and Nstar are calculated based on 9 combinations of LD window size and r^2 threshold.  

This function also has an option to draw the QQ plot comparing observed distributions of association statistics against those expected under the fitted models. The QQ plot will allow the users to visually assess the goodness fit of the model.
}
\value{
\item{\code{estimates}}{list containing the parameter estimates of the fitted model.}
\item{\code{summaryGWASdata}}{dataframe containing the filtered summary GWAS data  after merged with the LD information, with 10 columns
("SNP", "Z", "N", "betahat", "varbetahat", "Nstar", "TaggingSNPs", "LD.score.correct", "CHR", "BP").        
}
\item{\code{qqplotdata}}{
A list containing: 
(i) QQdata, a dataframe needed to generate the QQ plot, i.e., a dataframe containing 4 columns, "log10(observed pvalues)", "expected log10(fitted pvalues)", "lower bound of log10(fitted pvalues)","upper bound of log10(fitted pvalues)".
(ii) observedlambda, lambda genomic control (GC) factor for observed summary GWAS statistics.
(iii) meanEXPlambda,  expected lambda GC factor for simulated summary GWAS statistics under the fitted model.
(iv) lowEXPlambda, "qqplotCI.coverage" CI lower bound of lambda GC factor for simulated summary GWAS statistics under the fitted model. 
(v) highEXPlambda, "qqplotCI.coverage" CI upper bound of lambda GC factor for simulated summary GWAS statistics under the fitted model.}
}
\references{
Zhang, Yan, et al. "Estimation of complex effect-size distributions using summary-level statistics from genome-wide association studies across 32 complex traits" Nature Genetics (2018).
}
\author{
Yan Zhang <yandorazhang@gmail.com>
}
\examples{
data(heightGWAS)
fit = genesis(heightGWAS,modelcomponents=2, cores=2, startingpic=0.01)
fit$estimates

fit3 = genesis(heightGWAS,modelcomponents=3, cores=20, startingpic=0.01)
fit3$estimates
}

