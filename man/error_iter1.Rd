\name{error_iter1}
\alias{error_iter1}
\docType{data}
\title{
Random Errors with Particular LD Structure
}
\description{
Random errors following N(0, R), where R is the LD correlation matrix of the 1000 Genome reference panel. 
}
\usage{data("error_iter1")}
\format{
  1070777 errors corresponding to 1070777 SNPs. 
  \describe{
    \item{\code{error}}{a numeric vector}
    \item{\code{SNP}}{a character vector}
  }
}
\details{
For 489 subjects in the 1000 Genome reference dataset, we simulate 489 pseudo null phenotypes independently from N(0,1). Using the genotype data in the 1000 Genome reference panel, calculate association statistics (z-statistics) for the 1070777 SNPs, which are the errors following N(0,R).  (See Methods (Simulation studies) in Zhang et. al 2018)
}
\references{
Zhang, Yan, et al. "Estimation of complex effect-size distributions using summary-level statistics from genome-wide association studies across 32 complex traits" Nature Genetics (2018).
}
\author{
Yan Zhang <yandorazhang@gmail.com>
}
\examples{
data(error_iter1)
}

