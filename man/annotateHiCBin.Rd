\name{annotateHiCBin}
\alias{annotateHiCBin}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Annotate blocking regions with ChIP-seq data.
}
\description{
This is an internal function (should not be used). Function to compute the average of the coverage of ChIP-seq peaks within blocking regions. 
}
\usage{
annotateHiCBin(HiC_bin.GR, GenomicFeature.GR)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{HiC_bin.GR}{
A GRanges object containing the Hi-C blocking regions.
}
  \item{GenomicFeature.GR}{
A GRanges object containing the genomic feature intervals with fold-enrichments (score).
}
}
\value{
A vector object is returned.
}
\author{
Raphael Mourad
}
