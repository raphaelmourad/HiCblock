# Raphael MOURAD
# LBME lab
# 07/10/2016



# Create HiCDataset object used for HiCglmiProcData function from HTClist object from HiTC package.
# Zeroes in the HiC data matrices are filtered out. CHECKED 10/06/2016.
createHiCDataset<-function(HTCL,distInter=NULL,verbose=F){

 if(class(HTCL)!="HTClist"){print("HTCL is not a HTClist object!"); return(0)}
 if(class(distInter)!="numeric" & length(distInter)==2){print("distInter.V is not a vector of 2 numerical values!"); return(0)}

 chr=seqlevels(HTCL)
 binsize=width(x_intervals(HTCL[[1]]))[1]
 left.GR=NULL
 right.GR=NULL
 data=NULL
 

 for(i in 1:length(chr)){

  if(is.null(distInter)){
   distanceInterBin=c(0,length(x_intervals(HTCL[[i]])))
  }else{
   distanceInterBin=c(distInter[1]/binsize,min(distInter[2]/binsize,length(x_intervals(HTCL[[i]]))))
  }

  mati=intdata(HTCL[[i]])
  if(class(mati)=="dgeMatrix" | class(mati)=="dtCMatrix"){
   tabiTemp=as.matrix(summary(as(mati,"dgCMatrix")))
  }else{
   tabiTemp=as.matrix(summary(as(as(mati,"dsCMatrix"),"dgCMatrix")))
  }
  tabi=tabiTemp[(tabiTemp[,2]-tabiTemp[,1])>=distanceInterBin[1] & (tabiTemp[,2]-tabiTemp[,1])<=distanceInterBin[2],]
  mati.GR=x_intervals(HTCL[[i]])
  rm(mati,tabiTemp)

  left.GRi=mati.GR[tabi[,1]]
  left.GRi$idxBin=tabi[,1]
  right.GRi=mati.GR[tabi[,2]]
  right.GRi$idxBin=tabi[,2]
  datai=data.frame(tabi[,2]-tabi[,1],tabi[,3],rep(chr[i],nrow(tabi)))
  datai[,1]=datai[,1]*binsize
  colnames(datai)<-c("dist","count","chr")
  data=rbind(data,datai)
  if(i==1){
   left.GR=left.GRi
   right.GR=right.GRi
  }else{
   left.GR=suppressWarnings(c(left.GR,left.GRi))
   right.GR=suppressWarnings(c(right.GR,right.GRi))
  }
  rm(tabi,datai,left.GRi,right.GRi)

  if(verbose){print(paste(chr[i]," processed",sep=""))}
 }

 HiC_dataset=list(left.GR=left.GR,right.GR=right.GR,data=data)
 return(HiC_dataset)
}


# Test power law distribution before using the model. CHECKED 10/06/2016.
testDistancePowerLawDistrib<-function(HTCL,distInter){

 HiC_data=createHiCDataset(HTCL,distInter)$data

 data_d=HiC_data[HiC_data[,1] >=distInter[1] & HiC_data[,1] <=distInter[2],]
 mean_d=by(data_d[,2], data_d[,1]+5e3, function(x) mean(x))
 mean_HiC=as.vector(mean_d)
 mean_dist=as.numeric(names(mean_d))
 mean_distkb=mean_dist/1e3

 lm_hic_dist=lm(log10(mean_HiC)~log10(mean_distkb))
 coef_hic_dist=summary(lm_hic_dist)$coefficients[2,1]
 R2_hic_dist=summary(lm_hic_dist)$r.squared
 pval_hic_dist=summary(lm_hic_dist)$coefficients[2,4]

 list2return=list(lm=summary(lm_hic_dist),coefficient=coef_hic_dist,r.squared=R2_hic_dist,p.value=pval_hic_dist,mean_HiC=mean_HiC,mean_distkb=mean_distkb)
 return(list2return)
}


# Map a genomic feature to HiC bins. CHECKED 11/10/2016!
annotateHiCBin <- function(HiC_bin.GR, GenomicFeature.GR)
{
 if(!is(HiC_bin.GR, "GenomicRanges")){stop("'HiC_bin.GR' must be a GenomicRanges object")}
 if(!is(GenomicFeature.GR, "GenomicRanges")){stop("'GenomicFeature.GR' must be a GenomicRanges object")}
 #if(any(is.na(seqlengths(GenomicFeature.GR)))){stop("'seqlengths(GenomicFeature.GR)' contains NAs")}

 if(is.null(GenomicFeature.GR$score)){
  cvg <- coverage(GenomicFeature.GR)
 }else{
  cvg <- coverage(GenomicFeature.GR, weight="score")
 }

 chr=seqnames(seqinfo(HiC_bin.GR))
 HiC_bin.list=list()
 for(i in 1:length(chr)){
  HiC_bin.list[[chr[i]]]=HiC_bin.GR[seqnames(HiC_bin.GR)==chr[i]]
 }
 HiC_bin.GRL=GRangesList(HiC_bin.list)

 if(length(chr)==1){
  averageBin=viewMeans(Views(cvg[[chr]],ranges(HiC_bin.GRL[[chr]])))
 }else{
  views_list <- RleViewsList(lapply(names(cvg),function(seqname)Views(cvg[[seqname]], ranges(HiC_bin.GRL[[seqname]]))))
  averageBin=NULL
  for(i in 1:length(views_list)){
   averageBin=c(averageBin,viewMeans(views_list)[[i]])
  }
 }

 return(averageBin)
}







