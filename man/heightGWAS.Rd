\name{heightGWAS}
\alias{heightGWAS}
\docType{data}
\title{Example GWAS Summary Data for Height (2010)
}
\description{
The dataset contains GWAS summary data of GIANT_HEIGHT_LangoAllen2010_publicrelease_HapMapCeuFreq.txt after quality control, containing 1037172 rows and 3 columns. Each row represent a SNP. 
}
\usage{data("heightGWAS")}
\format{
   A matrix with 1037172 rows and 3 columns (snp, z, n)
    \describe{
    \item{\code{snp}}{SNP rsID}
    \item{\code{z}}{GWAS z-statistic}
    \item{\code{n}}{Sample size associated with the GWAS study of z-statistic}
  }
}
\source{http://portals.broadinstitute.org/collaboration/giant/index.php/GIANT_consortium_data_files#Height_.28download_GZIP.29}
\references{
Lango Allen, H., Estrada, K., Lettre, G., Berndt, S.I., Weedon, M.N., Rivadeneira, F., Willer, C.J., Jackson, A.U., Vedantam, S., Raychaudhuri, S., et al. (2010). Hundreds of variants clustered in genomic loci and biological pathways affect human height. Nature 467, 832-838.
}
\examples{
data("heightGWAS")
}
\keyword{datasets}
