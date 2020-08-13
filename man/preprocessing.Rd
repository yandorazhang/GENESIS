\name{preprocessing}
\alias{preprocessing}
\title{
Function to Pre-process Summary-level GWAS Statistics
}
\description{
This function allows to preprocess the summary-level GWAS statistics. 
}
\usage{
preprocessing(summarydata, LDcutoff=0.1,LDwindow=1,filter=FALSE)
}
\arguments{
   \item{\code{summarydata}}{summay-level GWAS data, containing 3 columns: 
  SNP (SNP rsID), 
  Z (GWAS test z-statistic), 
  N (GWAS study sample size which can be different for different SNPs)}
  \item{\code{LDcutoff}}{a number from (0.05, 0.1, 0.2); indicating LD score is calculated based on the particular r^2 cutoff. By default, it is 0.1.}
   \item{\code{LDwindow}}{a number from (0.5, 1, 2); indicating LD score is calculated based on the particular window size (MB). By default, it is 1 MB.}
  \item{\code{filter}}{logical; if TRUE, the input summary data will be filtered.}
}
\details{
If filter==TRUE, the input summary data will be filtered by 3 steps: (i) If sample size varies from SNP to SNP, remove SNPs with an effective sample size less than 0.67 times the 90th percentile of sample size. (ii) Remove SNPs within the major histocompatibility complex (MHR) region; filter SNPs to Hapmap3 SNPs. (iii) Remove SNPs with extremely large effect sizes (chi^2 > 80).
}
\value{
a dataframe containing the filtered summary GWAS data  after merged with the LD information, with 10 columns
("SNP", "Z", "N", "betahat", "varbetahat", "Nstar", "TaggingSNPs", "LD.score.correct", "CHR", "BP").        
}
\author{
Yan Zhang <yandorazhang@gmail.com>
}
\examples{
data(heightGWAS)
df <- preprocessing(heightGWAS,filter=FALSE)
head(df)
}
