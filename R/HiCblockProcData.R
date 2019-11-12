# Raphael MOURAD
# LBME lab
# 07/10/2016




# Function to process Hi-C and genomic feature data (such as ChIP-seq peaks). CHECKED 11/06/2016.
HiCblockProcData<-function(genomicFeatureList.GR,annotNames,HTCList,
		      distInter,overlapmode="signal",verbose=F,includeBias=T){


# CHECK INPUT DATA ------------------------------------------------------------ 

if(verbose){print("Data checking")}

if(class(genomicFeatureList.GR)!="list"){print("genomicFeatureList.GR is not a list object!"); return(0)}
for(i in 1:length(genomicFeatureList.GR)){
 if(class(genomicFeatureList.GR[[i]])!="GRanges"){print("i-th object of genomicFeatureList.GR is not a GenomicRanges object!"); return(0)}
}
if(class(annotNames)!="character"){print("annotNames is not a character object!"); return(0)}
if(class(HTCList)!="HTClist"){print("HTCList is not a HiTClist object!"); return(0)}
if(class(distInter)!="numeric" & length(distInter)==2){print("distInter is not a vector of 2 numerical values!"); return(0)}



#  DATA TREATMENT---------------------------------------------------------
if(verbose){print("Data parsing")}

# Create HiC_dataset # CHECKED!
Chr=seqlevels(HTCList)
HiC_dataset=createHiCDataset(HTCList,distInter,verbose)
HiC_binList=lapply(HTCList,x_intervals)
names(HiC_binList)=Chr
binSize=width(HiC_binList[[1]][1])
rm(HTCList,HiC_binList)

# Filter HiC bin couples based on specific filters.GR (e.g. TAD intervals) # CHECKED!
HiC_left.GR=HiC_dataset$left.GR
HiC_right.GR=HiC_dataset$right.GR
HiC_data=HiC_dataset$data
Idx_IC=1:length(HiC_left.GR)
HiC_left.GR$idxPair=Idx_IC
HiC_right.GR$idxPair=Idx_IC
rm(HiC_dataset,Idx_IC) 

# Biases # CHECKED!
if(includeBias){
HiC_left.Bias=mcols(HiC_left.GR)[,1:3]
HiC_right.Bias=mcols(HiC_right.GR)[,1:3]
}else{
HiC_left.Bias=NULL
HiC_right.Bias=NULL
}

# Count overlapping in-between the two bins (blocking variables)
HiC_betw.Fac=NULL
for(i in 1:length(genomicFeatureList.GR)){
 annot_betwi=NULL
 HiC_betwlen.Fac=NULL
 for(j in 1:length(Chr)){
  HiC_left_right.GRi=GRanges(Chr[j],IRanges(end(HiC_left.GR[seqnames(HiC_left.GR)==Chr[j]])+1*binSize,start(HiC_right.GR[seqnames(HiC_right.GR)==Chr[j]])-1*binSize))
  distij=end(HiC_left_right.GRi)-start(HiC_left_right.GRi)-1

  if(overlapmode=="signal"){
   GRi=genomicFeatureList.GR[[i]]
   GRi$score[GRi$score<0]=0
   annot_betwij=annotateHiCBin(HiC_left_right.GRi,GRi)
  }else if(overlapmode=="occurrence"){
   annot_betwij=countOverlaps(HiC_left_right.GRi,genomicFeatureList.GR[[i]])/(distij/1e3)
  }

  annot_betwi=c(annot_betwi,annot_betwij)
  HiC_betwlen.Fac=c(HiC_betwlen.Fac,width(HiC_left_right.GRi))
 }
 if(is.null(HiC_betw.Fac)){
  HiC_betw.Fac=cbind(HiC_betwlen.Fac,annot_betwi)
 }else{
  HiC_betw.Fac=cbind(HiC_betw.Fac,annot_betwi)
 }
 rm(annot_betwi)

 if(verbose){print(paste0(annotNames[i]," annotated"))}
}
colnames(HiC_betw.Fac)=c("Blen",annotNames)


ProcData=list(left.Bias=HiC_left.Bias,right.Bias=HiC_right.Bias,betw.Fac=HiC_betw.Fac,data=HiC_data,annotNames=annotNames,binSize=binSize)
return(ProcData)

}# End of function HiCblockProcData


