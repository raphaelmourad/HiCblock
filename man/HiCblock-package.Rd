\name{HiCblock-package}
\alias{HiCblock-package}
\alias{HiCblock}
\docType{package}
\title{
\packageTitle{HiCblock}
}
\description{
\packageDescription{HiCblock}
}
\details{
# To install dependencies:\cr
source("https://bioconductor.org/biocLite.R")\cr
biocLite("IRanges","GenomicRanges","GenomeInfoDb","rtracklayer","HiTC")\cr
install.packages(c("methods","MASS","Matrix","glmnet","S4Vectors"))\cr

# To install HiCblock package:\cr
install.packages("HiCblock")\cr

# To use the package, there are two steps:\cr
- To assess the blocking effects of features such as architectural proteins, one should first annotate Hi-C matrix with bias data (GC-content, mappability, fragment length) and with feature data (for instance, ChIP-seq data) with the function HiCblockProcData(). Note that if the Hi-C matrix has already been corrected for biases, then no bias data need to be used. \cr
- Then, one should compute the negative binomial or Poisson lasso regression with the function HiCblockModel() using the output from HiCblockProcData(). \cr
The blocking effect of a protein (or motif) is the associated beta from the regression.

# Choice between Poisson lasso over binomial negative regressions:\cr
Because of Hi-C count overdispersion, we used negative binomial regression as the most appropriate specification of the generalized linear model. However, Poisson regression with lasso shrinkage can also be used. We believe that the choice between both depends mainly on the number of variables to analyze. On the one hand, if there are a few candidate variables (less than 10), it is interesting to estimate beta parameters together with corresponding p-values to assess significance using negative binomial regression. On the other hand, if there are a large number of variables (10 or more), it is more convenient to use Poisson lasso regression in order to select the key variables and to account for correlations among the variables (frequent in ChIP-seq and motif occurrence data).\cr


}
\author{
Raphael Mourad

Maintainer: Raphael Mourad, raphael.mourad@ibcg.biotoul.fr
}
\references{
Raphael Mourad and Olivier Cuvier. TAD-free analysis of architectural proteins and insulators.  Nucleic Acids Research, Volume 46, Issue 5, 16 March 2018, Pages e27.
}

\examples{
# Load data
# The Hi-C matrix is at 20kb resolution (low resolution only for example)
data(dataExample)
genomicFeatureList.GR=dataExample$GenomicFeatureList.GR
annotNames=dataExample$AnnotNames
HTCList=dataExample$HTC
distInter=c(100e3,140e3)
IBP=c("BEAF32","dCTCF","dTFIIIC","GAF","SuHw")

# Annotate Hi-C data with genomic features
HRPD=HiCblockProcData(genomicFeatureList.GR, annotNames, HTCList, distInter,verbose=TRUE)

# Model 1
modelBlock1=as.formula(paste0("Count~logDist+len+GC+map+I(-BEAF32)"))
HRM_Block1=HiCblockModel(HRPD,modelBlock1,"BEAF32",regressionMode="NB")
print(HRM_Block1)

# Output from model 1
# Blocking effect (beta) of BEAF-32 is 0.21
#Call:
#glm.nb(formula = model, data = dataGLM, init.theta = 2.956475677, 
#    link = log)

#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-3.9346  -0.9018  -0.2606   0.4158   5.0515  

#Coefficients:
#             Estimate Std. Error z value Pr(>|z|)    
#(Intercept) -4.420963   1.102393  -4.010 6.06e-05 ***
#logDist     -1.226667   0.073070 -16.788  < 2e-16 ***
#len          0.662864   0.027635  23.987  < 2e-16 ***
#GC          -0.627353   0.087941  -7.134 9.76e-13 ***
#map          1.258398   0.054150  23.239  < 2e-16 ***
#I(-BEAF32)   0.206394   0.007735  26.683  < 2e-16 ***
#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#(Dispersion parameter for Negative Binomial(2.9565) family taken to be 1)

#    Null deviance: 5377.2  on 3375  degrees of freedom
#Residual deviance: 3575.1  on 3370  degrees of freedom
#AIC: 46123

#Number of Fisher Scoring iterations: 1


#              Theta:  2.9565 
#          Std. Err.:  0.0693 

# 2 x log-likelihood:  -46109.2770 



# Model 2
# Blocking effect (beta) of BEAF-32 is 0.18
if(FALSE){
vars2=paste(paste0("I(-",IBP,")"),collapse='+')
modelBlock2=as.formula(paste0("Count~logDist+len+GC+map+",vars2))
HRM_Block2=HiCblockModel(HRPD,modelBlock2,IBP,regressionMode="NB")
print(HRM_Block2)
}

# Output from model 2
#Call:
#glm.nb(formula = model, data = dataGLM, init.theta = 3.138426428, 
#    link = log)

#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-3.7727  -0.8815  -0.2809   0.4282   4.3644  

#Coefficients:
#             Estimate Std. Error z value Pr(>|z|)    
#(Intercept) -4.982498   1.074279  -4.638 3.52e-06 ***
#logDist     -1.236704   0.070946 -17.432  < 2e-16 ***
#len          0.685146   0.027209  25.181  < 2e-16 ***
#GC          -0.560153   0.086482  -6.477 9.35e-11 ***
#map          1.310488   0.053692  24.407  < 2e-16 ***
#I(-BEAF32)   0.176687   0.008283  21.332  < 2e-16 ***
#I(-dCTCF)    0.028798   0.013665   2.107  0.03508 *  
#I(-dTFIIIC)  0.621412   0.077602   8.008 1.17e-15 ***
#I(-GAF)      0.047998   0.007878   6.092 1.11e-09 ***
#I(-SuHw)     0.072596   0.026326   2.758  0.00582 ** 
#---
#Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

#(Dispersion parameter for Negative Binomial(3.1384) family taken to be 1)

#    Null deviance: 5702.9  on 3375  degrees of freedom
#Residual deviance: 3564.8  on 3366  degrees of freedom
#AIC: 45911

#Number of Fisher Scoring iterations: 1


#              Theta:  3.1384 
#          Std. Err.:  0.0738 

# 2 x log-likelihood:  -45889.1960 



# Model 3
# Blocking effect (beta) of BEAF-32 is 0.22
if(FALSE){
HRM_Block3=HiCblockModel(HRPD,NULL,IBP,regressionMode="PoissonLasso")
print(HRM_Block3)
}

# Output from model 3
#        Variable Coefficient
#logDist  logDist    -1.04649
#len          len     0.56508
#GC            GC    -0.43621
#map          map     0.95623
#BEAF32    BEAF32     0.22077
#dCTCF      dCTCF     0.01680
#dTFIIIC  dTFIIIC     0.73750
#GAF          GAF     0.03878
#SuHw        SuHw     0.12617



}
