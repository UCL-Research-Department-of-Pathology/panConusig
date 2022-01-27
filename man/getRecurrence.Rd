\name{getRecurrence} 
\alias{getRecurrence} 
\title{panConusig::getRecurrence} 
\usage{getRecurrence(sample,chrom,start,end,CNS,
                         windowsGR,sigNames,
                         tumType,
                         ncores=8,returnCounts=FALSE)}
\arguments{
\item{sample}{Vector of sample names for each segment from a set of copy number profiles.}
\item{chrom}{Vector of chromosome names for each segment from a set of copy number profiles.}
\item{start}{Vector of start positions for each segment from a set of copy number profiles.}
\item{end}{Vector of end positions for each segment from a set of copy number profiles.}
\item{CNS}{Vector of signature assignments for each segment from a set of copy number profiles. See runAssignment() for method to assign copy number signatures to segments.}
\item{windowsGR}{GenomicRanges object with windows in which to evaluate recurrence of signatures. Can be generated using setupWindows().}
\item{sigNames}{Vector if signature names for which to evaluate recurrence. Should match those in CNS.}
\item{ncores}{Number of cores to use for calculation. Defaults to 8.}
\item{returnCounts}{Boolean. Whether to return counts of recurrence (TRUE), or proportions (FALSE).}
}
\value{
A matrix giving the proportion or count of samples with a signature recurrence in a given genomic window. Matrix will be of wxn where w is the number of genomic windows and n is the number if signatures listed in sigNames.
} 
\description{
Calculates the recurrence of signatures across the genome.
}
\details{
  
}

