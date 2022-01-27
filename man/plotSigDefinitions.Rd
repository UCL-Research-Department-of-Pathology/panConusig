\name{plotSigDefinitions} 
\alias{plotSigDefinitions} 
\title{panConusig::plotSigDefinitions} 
\usage{plotSigDefinitions(sigs,sigLength=2,doTitle=TRUE,doTitleBar=FALSE,doText=FALSE,
                               doYlab=TRUE,doLab=FALSE,ncol=1,nPerCol="all",plotBars=c(1,2),
                               myMAR=NULL,barSize=1.15,returnLayout=FALSE,myLayout=NULL,CEX=1.8,
                               YLIM=NULL,mainAdj=0.5,mainLine=0,doBox=FALSE,HDlab="",
                               extraLabel=NULL,extraAdj=0.5,extraAt=15,extraCex=0.5)}
\arguments{
\item{sigs}{Matrix. 48xn matrix of signature definitions, where n is the number of signatures}
\item{sigLength}{Signature plot width to be specified in layout. Integer.}
\item{doTitle}{Whether to plot signature names as a title.}
\item{doTitleBar}{Whether to plot copy number classes as coloured bars.}
\item{doText}{Whether to add copy number class labels as text in coloured bars..}
\item{doYlab}{Whether to plot signature names as y-axis labels.}
\item{doLab}{Whether to plot segment size classes as x-axis labels.}
\item{ncol}{Number of columns to plot signatures in.}
\item{nPerCol}{Number of signatres per column to plot.}
\item{plotBars}{Vector of any combination of c(1,2). If 1 is in plotBars, coloured bars will be plotted above signatures. If 2 is in plotBars, coloured bars will be plotted below signatures.}
\item{myMAR}{Custom margins.}
\item{barSize}{Size for coloured bars.}
\item{returnLayout}{Boolean. Whether to return the layout matrix used for plotting.}
\item{myLayout}{Custom layout for plotting.}
\item{YLIM}{Custom y-axis limits.}
\item{mainAdj}{X-axis adjustment for signature title.}
\item{mainLine}{y-axis line for signature title.}
\item{doBox}{Boolean. Whether to surround each signature definition with a box.}
\item{HDlab}{Label to use for homozygous deletions.}
\item{extraLabel}{Label for signature subtitle.}
\item{extraAdj}{Adjustment for signature subtitle.}
\item{extraAt}{Position for signature subtitle.}
\item{extraCex}{cex for signature subtitle.}
}
\value{
Matrix of layout used if returnLayout==TRUE.
} 
\description{
Generates barplots of signature definitions.
}
\details{
  
}

