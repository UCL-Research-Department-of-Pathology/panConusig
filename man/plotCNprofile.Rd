\name{plotCNprofile} 
\alias{plotCNprofile} 
\title{panConusig::plotCNprofile} 
\usage{plotCNprofile(chr,start,end,
	majCN,minCN,main=NA,colours=NULL,
	YLIM=NULL,offsetVal=100,LWD=7,CEX=3,
	custom_legend=NULL)}
\arguments{
\item{chr}{Vector of chromosome names for each segment from a single copy number profile.}
\item{start}{Vector of start positions for each segment from a single copy number profile.}
\item{end}{Vector of end positions for each segment from a single copy number profile.}
\item{majCN}{Vector of major copy number for each segment from a single copy number profile.}
\item{minCN}{Vector of minor copy number for each segment from a single copy number profile.}
\item{main}{Title for plot.}
\item{colours}{Colours for plotting segments. If a vector, colours[1] will be used for major copy number, and colours[2] will be used for minor copy number. If a list, colours[[1]] will be used for major segments, colours[[2]] will be used for minor segments. Defaults to c("navyblue","darkorange") if left NULL.}
\item{YLIM}{Y-axis limits for plotting.}
\item{offsetVal}{Major and minor copy number states will be offset from their integer values by max(YLIM)/offsetVal so that they do not overlap.}
\item{LWD}{Line width for plotting segments.}
\item{CEX}{cex value for increasing text size. See par().}
\item{custom_legend}{List of arguments to be called with legend().}
}
\value{
None.
} 
\description{
Produces a plot of an allele specific copy number profile.
}
\details{
  
}

