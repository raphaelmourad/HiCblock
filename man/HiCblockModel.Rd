\name{HiCblockModel}
\alias{HiCblockModel}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The model.
}
\description{
The main function to compute the generalized linear model.
}
\usage{
HiCblockModel(hrpd, model, facBlock, regressionMode = "NB", scale = F, includeBias = T, 
		sampleSize = NULL, distInter = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hrpd}{
Preprocessed data. It includes Hi-C interaction pairs that have been annotated with genomic feature information. To compute this object, you should use the function HiCglmiProcData.
}
  \item{model}{
A formula object (stats R package). See Section HiCglmi-package, for examples. If regressionMode="PoissonLasso", then the variable model should be set to NULL (i.e. model=NULL).
}
  \item{facBlock}{
A vector of character describing the genomic features. For instance, facModel="BEAF32".
}
  \item{regressionMode}{
If regressionMode="Poisson", a Poisson regression. 
If regressionMode="QP", a quasi Poisson regression. 
If regressionMode="NB", a negative binomial regression. 
If regressionMode="PoissonLasso", a Poisson lasso regression. 
}
  \item{scale}{
If scale=TRUE, blocking variables are standardized. 
}
  \item{includeBias}{
If includeBias=True, then GC, mappability and fragment length biases are included in the regression (default mode).
If includeBias=False, biases are not included in the regression. You can used this option if you are using Hi-C matrix that had been previously corrected for biases (not recommended). 
}
  \item{sampleSize}{
Optional. A numerical value for subsampling Hi-C data.
}
  \item{distInter}{
Optional. A vector of two numerical values to set a distance range for Hi-C analysis.
}
}
\value{
If regressionMode is "Poisson", "QP" or "NB", then a summary(glm) object is returned.
}
\references{
Submitted.
}
\author{
Raphael Mourad
}
