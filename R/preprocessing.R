#' Illustration of preprocessing()
#'
#' This function allows to preprocess the summary-level GWAS statistics. 
#' @param summarydata summay-level GWAS data, containing 3 columns: 
#' SNP (SNP rsID), 
#' Z (GWAS test z-statistic), 
#' N (GWAS study sample size which can be different for different SNPs)
#' @param LDcutoff a number from (0.05, 0.1, 0.2); indicating LD score is calculated based on the particular r^2 cutoff. By default, it is 0.1.
#' @param LDwindow a number from (0.5, 1, 2); indicating LD score is calculated based on the particular window size (MB). By default, it is 1 MB.
#' @param filter logical; if TRUE, the input summary data will be filtered.
#' @keywords 
#' @export
#' @examples preprocessing(summarydata, LDcutoff=0.1,LDwindow=1,filter=FALSE)

preprocessing <- function(summarydata, LDcutoff=0.1,LDwindow=1,filter=FALSE){
  #----------------------------------------------------#----------------------------------------------------
  # I. summary GWAS data format check
  #----------------------------------------------------#----------------------------------------------------
  if(ncol(summarydata)!=3){
    stop("The summary GWAS data should have 3 columns with (SNP rsID, Z-statistic, Sample size)!")
  }
  
  #----------------------------------------------------#----------------------------------------------------
  # II. preliminary summary GWAS data filtering
  #----------------------------------------------------#----------------------------------------------------
  colnames(summarydata) <- c("SNP","Z","N")
  
  if(filter==TRUE){
    #a. If sample size varies from SNP to SNP, remove SNPs with an effective sample size less than 0.67 times the 90th percentile of sample size.
    ikeep1 <- which(as.numeric(as.character(summarydata$N))>=0.67*quantile(as.numeric(as.character(summarydata$N)), 0.9))
    summarydata <- summarydata[ikeep1,]
    
    #b. Remove SNPs within the major histocompatibility complex (MHR) region; filter SNPs to Hapmap3 SNPs.
    data(w_hm3.noMHC.snplist)
    ikeep2 <- which(as.character(summarydata$SNP) %in% w_hm3.noMHC.snplist$SNP)
    summarydata <- summarydata[ikeep2,]
    
    #c. Remove SNPs with extremely large effect sizes (chi^2 > 80).
    ikeep3 <- which(as.numeric(as.character(summarydata$Z))^2 <=80)
    summarydata <- summarydata[ikeep3,]
  }
  
  #----------------------------------------------------#----------------------------------------------------
  # III. merge the summary GWAS data with the LD score data
  #----------------------------------------------------#----------------------------------------------------
  data(list=paste0("LDwindow",LDwindow,"MB_cutoff",LDcutoff))
  
  summarydata$SNP <- as.character(summarydata$SNP)
  summarydata$Z <- as.numeric(as.character(summarydata$Z))
  summarydata$N <- as.numeric(as.character(summarydata$N))
  
  # remove NA values.
  inx <- which(is.na(summarydata$Z))
  if(length(inx)>1) summarydata <- summarydata[-inx,]
  inx <- which(is.na(summarydata$N))
  if(length(inx)>1) summarydata <- summarydata[-inx,]
  
  # get the variable needed for our method.
  summarydata$betahat <- summarydata$Z/sqrt(summarydata$N)
  summarydata$varbetahat <- 1/summarydata$N
  df <- merge(summarydata, dataLD,by.x="SNP",by.y="SNPname",sort=F)
  
  return(df)
}


