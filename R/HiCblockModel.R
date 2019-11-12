# Raphael MOURAD
# Cuvier team, LBME lab
# 10/10/2016





# Function to compute generalized linear model with interactions.
HiCblockModel<-function(hrpd,model,facBlock,regressionMode="NB",scale=F,includeBias=T,sampleSize=NULL,distInter=NULL){


# MODEL VARIABLE COMPUTATION ---------------------------------------------------------------
print("Model variable computation")

# Left and right factors and cofactors # CHECKED
annotNames=hrpd$annotNames
if(!is.null(distInter)){ # filter by distance
testDist=hrpd$data[,1]>=distInter[1] & hrpd$data[,1]<=distInter[2]
}else{
testDist=1:nrow(hrpd$data)
}
if(includeBias){ # if biases are included
HiC_left.Bias=as(hrpd$left.Bias[testDist,],"Matrix")
HiC_right.Bias=as(hrpd$right.Bias[testDist,],"Matrix")
}
HiC_data=hrpd$data[testDist,]
HiC_mat.FacBlock=hrpd$betw.Fac[testDist,]

# Sampling
if(!is.null(sampleSize)){
 if(sampleSize<nrow(HiC_data)){ 
  totSize=nrow(HiC_data)
  idxS=sample(1:totSize,sampleSize)
 }else{
  totSize=nrow(HiC_data)
  idxS=1:totSize
 }
}else{
 totSize=nrow(HiC_data)
 idxS=1:totSize
}
HiC_mat.FacBlock=HiC_mat.FacBlock[idxS,]
if(includeBias){
HiC_left.Bias=HiC_left.Bias[idxS,]
HiC_right.Bias=HiC_right.Bias[idxS,]
}
HiC_data=HiC_data[idxS,]
binSize=hrpd$binSize
print(paste0("Data size: ",nrow(HiC_data)," rows"))
rm(hrpd)

# Remove Left and Right # CHECKED!
if(includeBias){
HiC_mat.bias=as(log(HiC_left.Bias*HiC_right.Bias),"dgCMatrix")
}

# Scaling data (usefull when proteins have very different number of peaks, e.g. from 200 to 60000).
if(scale){
 for(i in 2:ncol(HiC_mat.FacBlock)){
  HiC_mat.FacBlock[,i]=scale(HiC_mat.FacBlock[,i])
 }
}

# Distance # CHECKED!
dist=HiC_data[,1]
dist[HiC_data[,1]==0]=binSize/2
HiC_logDist=log(dist)

# HiC data # CHECKED!
HiC_count=HiC_data[,2]

# All data: Matrix format # CHECKED!
HiC_vecbind=as(cbind(HiC_count,HiC_logDist),"dgCMatrix")
if(includeBias){
 HiC_mat.All=cbind(HiC_vecbind,HiC_mat.bias,HiC_mat.FacBlock)
 rm(HiC_vecbind,HiC_mat.bias,HiC_mat.FacBlock)
}else{
 HiC_mat.All=cbind(HiC_vecbind,HiC_mat.FacBlock)
 rm(HiC_vecbind,HiC_mat.FacBlock)
}
colnames(HiC_mat.All)[1:2]=c("Count","logDist")




# REGRESSION ANALYSIS ---------------------------------------------------------------
print("Regression analysis")
print(paste0(regressionMode," regression"))

# GLM # CHECKED!
if(regressionMode!="PoissonLasso" & regressionMode!="RF"){
dataGLM=as.data.frame(as.matrix(HiC_mat.All))
dataGLM=dataGLM[!is.na(dataGLM[,5]),]
dataGLM=dataGLM[dataGLM[,3]!=-Inf,]
dataGLM=dataGLM[dataGLM[,"Blen"]>100,]
if(regressionMode=="Poisson"){
 GLM=suppressWarnings(glm(model,data=dataGLM, family=poisson()))
 #print(dispersiontest(GLM,trafo=1))
}else if(regressionMode=="QP"){
 GLM=suppressWarnings(glm(model,data=dataGLM, family=quasipoisson()))
}else if(regressionMode=="NB"){
 GLM=suppressWarnings(glm.nb(model,data=dataGLM))
}
toReturn=summary(GLM)
}

# Poisson Lasso # CHECKED!
if(regressionMode=="PoissonLasso"){
HiC_mat.All=HiC_mat.All[!is.na(HiC_mat.All[,5]),]
HiC_mat.All=HiC_mat.All[HiC_mat.All[,3]!=-Inf,]
HiC_mat.All=HiC_mat.All[HiC_mat.All[,"Blen"]>100,]
HiC_mat.All=HiC_mat.All[,colnames(HiC_mat.All)!="Blen"]
HiC_mat.All[,6:ncol(HiC_mat.All)]=-HiC_mat.All[,6:ncol(HiC_mat.All)] # -beta

CVLasso=cv.glmnet(HiC_mat.All[,-1],HiC_mat.All[,1],family="poisson",parallel=F)
lambda=CVLasso$lambda.min # CVLasso$lambda.min or CVLasso$lambda.1se
CVLassoError=CVLasso$cvm[which(CVLasso$lambda==lambda)]
devLasso=deviance.glmnet(CVLasso$glmnet.fit)[which(CVLasso$lambda==lambda)]
coefLasso=CVLasso$glmnet.fit$beta[,which(CVLasso$lambda==lambda)]
coefLassoMat=data.frame(Variable=names(coefLasso),Coefficient=round(coefLasso,5))
toReturn=coefLassoMat
}


return(toReturn)
}# End of function HiCblockModel


