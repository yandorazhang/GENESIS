Examples
===

## Load summary GWAS data - Height (allen2010)
The summary GWAS dataframe has 3 columns: (SNPname, z-statistic, sample size).  

```{r height}
library(GENESIS)
data(heightGWAS)
# heightGWAS has been pre-processed according to function preprocessing(). 
```


## Fit the 2-component model
#### 1.  Fitting to the model

Note the startingpic value can be specifided at a list of values, and then the one with the largest log-likelihood is selected as the final model. 

```{r 2-component model}
fit2 <- genesis(heightGWAS, filter=F, modelcomponents=2, cores=2, LDcutoff=0.1, LDwindow=1, c0=10, startingpic=0.005)
fit2$estimates

est <- fit2$estimates$`Parameter (pic, sigmasq, a) estimates` # the model parameter estimates
v <- fit2$estimates$`Covariance matrix of parameter estimates` # the covariance matrix of model parameter estimtaes

# the est and v should have below values
est <- c(4.906797e-03, 5.733005e-05, 1.861491e-06)
v <- matrix(c(2.870170e-07, -2.426575e-09, -6.077209e-11, -2.426575e-09,  2.491765e-11,  4.715283e-13, -6.077209e-11, 4.715283e-13,  1.720162e-14),3,3)
```

####  2. Get the density plot for the susceptibility SNPs 
```{r density plot}
x_seq <- seq(-0.02,0.02,length.out = 1000); 
y_seq <- apply(matrix(x_seq,ncol=1),1,function(t) dmixssnp(t,est))
plot(x_seq, y_seq,type="l",xlab="Joint effect size", ylab="Probability Density")
```

#### 3. Make future projections with specified sample size
```{r future projections}
projection(est,v,n=253288, CI=TRUE)
```

#### 4. Calculate number of SNPs falling in an interval
```{r number of SNPs in an interval}
numInterval(0.005,Inf,est)
```

#### 5.  Predict genomic control factor in future GWAS with specified sample size
```{r prediction}
futuregc(est,n=253288,nsim=1)
```


#### 6.  Calculate AUC of polygenic risk prediction model at given sample size
```{r function}
# PRS is calculated with SNPs included at optimum p-value threshold

polyriskpredict(N=253288, Ps=c(0.5,0.5), Sig2s=c(est[3],est[3]), M=1070777, M1=1070777*est[1], type="optimum", k.fold=3:5)

# PRS is calculated with SNPs included at genome-wide significance level
polyriskpredict(N=253288, Ps=c(0.5,0.5), Sig2s=c(est[3],est[3]), M=1070777, M1=1070777*est[1], type="GWAS",alp.GWAS=5e-8, k.fold=3:5)
```


## Fit the 3-component model
####  1. Fitting to the model

Note the startingpic value can be specifided at a list of values, and then the one with the largest log-likelihood is selected as the final model. 

```{r 3-component model}

# starting value of 3-component model comes from 2-component model estimates. 
starting <- rep(0,5)
starting[1] <- est[1]
starting[2] <- 1/9
starting[3] <- est[2]*5
starting[4] <- starting[3]/10
starting[5] <- est[3]

fit3 <- genesis(heightGWAS, filter=F, modelcomponents=3, cores=24, LDcutoff=0.1, LDwindow=1, c0=10,starting=starting)
fit3$estimates

est <- fit3$estimates$`Parameter (pic, p1, sigmasq1, sigmasq2, a) estimates` # the model parameter estimates
v <- fit3$estimates$`Covariance matrix of parameter estimates` # the covariance matrix of model parameter estimtaes

# est and v should have below values
est <-  c(8.899809e-03, 9.476025e-02, 1.458650e-04, 2.227118e-05, 1.567643e-06)
v <- matrix(c(1.327856e-06, -1.131049e-05, 6.912489e-10, -2.901301e-09, -9.388865e-11, -8.568269e-06, 3.380985e-04, -1.543479e-07, 1.771036e-09,8.989668e-10, -1.846216e-09, -1.353542e-07,  1.517545e-10,  1.166563e-11, 1.901686e-14, -3.068113e-09, 9.520417e-09,  5.633801e-12, 9.427492e-12, 1.220433e-13, -9.111521e-11, 1.063848e-09, -1.526281e-13, 1.042123e-13, 1.406410e-14), 5,5) 
```

#### 2. Get the density plot for the susceptibility SNPs 
```{r density plot}
x_seq = seq(-0.02,0.02,length.out = 1000); 
y_seq = apply(matrix(x_seq,ncol=1),1,function(t) dmixssnp(t,est))
plot(x_seq, y_seq,type="l",xlab="Joint effect size", ylab="Probability Density")
```

#### 3.  Make future projections with specified sample size
```{r future projections}
projection(est,v, n=253288, CI=TRUE);
```

#### 4. Calculate number of SNPs falling in an interval
```{r number of SNPs in an interval}
numInterval(0.005,Inf,est)
```

#### 5. Predict genomic control factor in future GWAS with specified sample size
```{r prediction}
futuregc(est,n=253288,nsim=1)
```


#### 6.  Calculate AUC of polygenic risk prediction model at given sample size 
```{r function}
# PRS is calculated with SNPs included at optimum p-value threshold

polyriskpredict(N=253288, Ps=c(est[2],1-est[2]), Sig2s=c(est[3],est[4]), M=1070777, M1=1070777*est[1], type="optimum", k.fold=3:5)


# PRS is calculated with SNPs included at genome-wide significance level

polyriskpredict(N=253288, Ps=c(est[2],1-est[2]), Sig2s=c(est[3],est[4]), M=1070777, M1=1070777*est[1], type="GWAS",alp.GWAS=5e-8, k.fold=3:5)
```