\name{generateProfile} 
\alias{generateProfile} 
\title{panConusig::generateProfile} 
\usage{generateProfile(nAltBeforeGD=20,
                           doGD=FALSE,
                           nAltAfterGD=20,
                           probs=NULL,
                           chroms=NULL,
                           index=NULL,
                           doCN=TRUE)}
\arguments{
\item{nAltBeforeGD}{Integer. Number of alterations to simulate before any genome doubling.}
\item{doGD}{Boolean. Whether to simulate a genome doubling event.}
\item{nAltAfterGD}{Integer. Number of alterations to simulate after any genome doubling.}
\item{probs}{Probabilities to simulate any given process. If given, should be a length 9 vector representing the probabilities for: "focal gain","focal loss","chromothripsis","translocation","chromoplexy","BFB","aneuploidy gain","aneuploidy loss","chromoanasynthesis".}
\item{chroms}{List. Set of chromosomes to begin the simulation with. If NULL an unaltered genome is used to start the simulation.}
\item{index}{Integer vector. If given, this restricts alterations to the chromosomes in chroms that match the index.}
\item{doCN}{Boolean. Should the output be converted from a list of derivative chromosomes to a segment file?}
}
\value{
if doCN==TRUE: A segment file with columns sample, chr, start, end. If doCN==FALSE: A list of derivative chromosomes.
} 
\description{
Simulate a copy number profile from any combination of the following processes: "focal gain","focal loss","chromothripsis","translocation","chromoplexy","BFB","aneuploidy gain","aneuploidy loss","chromoanasynthesis", as well as whole genome doubling. 
}
\details{
  
}

