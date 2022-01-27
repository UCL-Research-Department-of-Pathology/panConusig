\name{getMatrix} 
\alias{getMatrix} 
\title{panConusig::getMatrix} 
\usage{getMatrix(data)}
\arguments{
\item{data}{Matrix. CN profile segments, with columns sample, chr, start
pos, endpos, nMajor, nMinor.}
}
\value{
48xs matrix of CN summaries, where s is the number of samples in data.
} 
\description{
Function that categorises segments into 48 copy number classes, and generates an input matrix for matrix deconvolution. Segments are categorised by LOH status (homozygous deletion, loss of heterozygosity, heterozygous), total copy number (0, 1, 2, 3-4, 5-8, 9+) and segment size (0-100kb, 100kb-1Mb, 1Mb-10Mb, 10Mb-40Mb,>40Mb). 
}
\details{
  
}

