\name{setupWindows} 
\alias{setupWindows} 
\title{panConusig::setupWindows} 
\usage{setupWindows(chr,startpos,endpos,
windowLength = 1000000,doBiomart=FALSE)}
\arguments{
\item{chr}{Vector of chromesome labels for all segments from all copy numebr profiles.}
\item{startpos}{Vector of start positions for all segments from all copy numebr profiles.}
\item{endpos}{Vector of end positions for all segments from all copy numebr profiles.}
\item{windowLength}{Length of genomic windows to be returned, in bases.}
\item{doBiomart}{Boolean. Whether to return a biomart mart object. Used when plotting gene names in plotRecurrence(), if a gene genomicRanges object is not supplied.}
}
\value{
List of: [[windows]]=matrix of genomic windows, [[windowsGR]]=GenomicRanges object of genomic windows, [[splits]]=genomic windows that denote boundaries between chromosomes, [[mids]]=genomic windows that denot midpoints of chromosomes, [[ensembl]]=biomart mart object.
} 
\description{
Split the genome into windows of a given size (windowLength). Used as input when calculating recurrence of signatures across the genome.
}
\details{
  
}

