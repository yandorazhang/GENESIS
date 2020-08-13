\name{LDwindow0.5MB_cutoff0.05}
\alias{LDwindow0.5MB_cutoff0.05}
\docType{data}
\title{LD Information for Hapmap 3 SNPs with LD-Square Threshold 0.05, window size 0.5MB}
\description{The LD information is based on 1000 Genome data, i.e., analysis of 489 individuals of European origin from the 1000 GENOME project Phase 3. All the LD-scores and number of tagged SNPs are evaluated based on a reference panel of ~1.07 million common SNPs that were included in the Hapmap3 panel. Nstar is estimated by the number of SNPs in the reference panel which are within 0.5MB distance and have an estimated LD-square coefficient with the GWAS marker above a fixed threshold 0.05. The corresponding LD-score is calculated by summing up the corresponding squared LD coefficients.
}
\usage{data("LDwindow0.5MB_cutoff0.05")}
\format{A dataframe named "dataLD" with 1070777 rows and 6 columns ("SNPname" , "Nstar", "TaggingSNPs","LD.score.correct", "CHR" ,"BP").}
\examples{data("LDwindow0.5MB_cutoff0.05")}
\keyword{datasets}
