\name{runAssignment} 
\alias{runAssignment} 
\title{panConusig::runAssignment} 
\usage{runAssignment(data, # matrix of copy number profiles with columns sample, chr, startpos, endpos, nMajor, nMinor
                         sigs, # matrix of 48xn where n is the number of signatures
                         exposures, # matrix of n*m where n is the number of signatures, and m is the number of samples
                         outDir = NULL)}
\arguments{
\item{data}{Matrix of copy number profiles with columns sample, chr, startpos, endpos, nMajor, nMinor.}
\item{sigs}{Matrix of 48xn where n is the number of signatures}
\item{exposures}{Matrix of n*m where n is the number of signatures, and m is the number of samples}
\item{outDir}{Directory to save results to.}
}
\value{
Matrix that includes combines data with new columns: the copy number class for each segment in data, and the  copy number signature assignment for each segment in data.
} 
\description{
Takes copy number profiles, signature assignments per sample, and signature definitions, to assign signatures at a segment level.
}
\details{
  
}

