\name{randomRecurrence} 
\alias{randomRecurrence} 
\title{panConusig::randomRecurrence} 
\usage{randomRecurrence(ascat,
                            windowGR,
                            nrep,
                            nsim)}
\arguments{
\item{ascat}{Copy number profiles from which to simulate new datasets. Must include copy number signature assignment to segments. See runAssignment().}
\item{windowGR}{Windows in which to calculate recurrence of randomly generated datasets. See setupWindows(). Should match those used for getRecurrence().}
\item{nrep}{Number of replicates to run (number of random datasets to generate).}
\item{nsim}{Number of samples to simulate per random dataset.}
\item{whichsig}{Which signature to return values for.}
}
\value{
Matrix of recurrence of a given copy number signature (whichsig) across nrep simulations in genomic windows.
} 
\description{
Randomly samples new copy number profiles based on those provided. Segment boundaries for whole chromosomes are drawn from the samples provided in ascat. Signature assignemnts are sampled from all signature assignments in ascat for which the segment length class matched that of the segment being drawn for. Recurrence of the signature is then calculated across all samples from a given replicate.
}
\details{
  
}

