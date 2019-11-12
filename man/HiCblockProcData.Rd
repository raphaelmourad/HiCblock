\name{HiCblockProcData}
\alias{HiCblockProcData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Process data for further generalized linear regression.
}
\description{
This function is used to annotate the blocking regions, i.e. in-between bin pairs, with genomic feature information such as ChIP-seq peak intervals or DNA motif occurrence intervals.
}
\usage{
HiCblockProcData(genomicFeatureList.GR, annotNames, HTCList, distInter, 
		 overlapmode="signal", verbose = F, includeBias = T)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{genomicFeatureList.GR}{
A list of GRanges objects. Each GRanges object has been built from coordinate data using readGFBed function (for instance ChIP-seq peak coordinates or DNA motif occurrence coordinates). 
}
  \item{annotNames}{
A character vector defining the name of each genomic feature. 
Names should not comprise any special character such as ":/+-*^,;!?" because an internal R formula object is created inside the function.
}
  \item{HTCList}{
A HTClist object (HiTC R package). This object contains Hi-C intrachromosomal matrices.
}
  \item{distInter}{
Optional. A vector of two numerical values to set a distance range for Hi-C analysis.
}
  \item{overlapmode}{
If ChIP-seq peaks are used, overlapmode should be "signal" to compute the average of coverage of the peak fold-enrichment within the blocking region, divided by the length of the blocking region.
If DNA motif occurrences are used, overlapmode should be "occurrence" to compute the number of occurrences of DNA motifs within the blocking region, divided by the length of the blocking region.
}
  \item{verbose}{
If true, verbose is output while the program is running.
}
  \item{includeBias}{
If True, biases are included.
If False, biases are not included.
}
}
\value{
An object further used for HiCblockModel function.
}
\author{
Raphael Mourad
}
