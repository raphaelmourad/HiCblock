# HiCblock: TAD-free analysis of architectural proteins and insulators

![alt text](https://github.com/morphos30/PhyloCTCFLooping/blob/master/approach.png)

**Overview**

Because of Hi-C count overdispersion, we used negative binomial regression as the most appropriate specification of the generalized linear model. However, Poisson regression with lasso shrinkage can also be used. We believe that the choice between both depends mainly on the number of variables to analyze. On the one hand, if there are a few candidate variables (less than 10), it is interesting to estimate beta parameters together with corresponding p-values to assess significance using negative binomial regression. On the other hand, if there are a large number of variables (10 or more), it is more convenient to use Poisson lasso regression in order to select the key variables and to account for correlations among the variables (frequent in ChIP-seq and motif occurrence data).

**References**
Raphael Mourad and Olivier Cuvier. TAD-free analysis of architectural proteins and insulators.  Nucleic Acids Research, Volume 46, Issue 5, 16 March 2018, Pages e27.

**Contact**:
raphael.mourad@univ-tlse3.fr
