\name{plotRecurrence} 
\alias{plotRecurrence} 
\title{panConusig::plotRecurrence} 
\usage{plotRecurrence(props,counts,tumTypes,windows,
arms,doGenes=FALSE,doArms=TRUE,plotGenes=TRUE,outDir="~/",
filename=NULL,whichSigs="CN19",
colours = c("red","blue","green","black","gray","orange","cyan","purple"),
whichGenes="CN19",genesToTest="MDM2",propThreshold=0.2,countThreshold=5,
ensembl=NULL,splits=NA,mids=NA,saveRes=FALSE,geneGR=NULL,windowGR=NULL,
doCombine=TRUE,doSex=FALSE,doLegend=TRUE,
doMain=TRUE,LTY=1:5,YLIM=NULL,LABY=NULL,
offsetBase=0.01,offsetTweak=0.02,xoverlap=100,
yoverlap=0.02,myMAR=NULL,XLIMadj=0,
geneCex=2,legloc="topright",doGenesAbove=TRUE,
mylayout=NULL,ylineoffset=NULL,
pdfWidth=21,pdfHeight=7,CEX=4,doSplit=FALSE,legnames=NULL,propArm=0.75,
geneYtweak=0,extraMain=NULL)}
\arguments{
\item{props}{Matrix of proportions of signature recurrence across genomic bins. See getRecurrence().}
\item{counts}{Matrix of counts of signature recurrence across genomic bins. See getRecurrence().}
\item{tumTypes}{Vector of tumour types for each sample.}
\item{windows}{Genomic windows that were used to generate props and counts. See setupWindows().}
\item{arms}{Chromosome arm information for plotting labels of recurrent signatures at an arm/chromosome level. See loadCytoband().}
\item{doGenes}{Whether to calculate gene names of genes overlapping regions of recurrence.}
\item{doArms}{Whether to plot chromosome arm names of arms overlapping regions of recurrence.}
\item{plotGenes}{Whether to plot gene names of genes overlapping regions of recurrence.}
\item{outDir}{Directory for saving plot as pdf.}
\item{filename}{Filename for saving plot as pdf. If NULL no pdf will be produced.}
\item{whichSigs}{Vector of signature names to plot.}
\item{colours}{Vector of colours for each signature.}
\item{whichGenes}{Vector of signature names for which to plot gene names.}
\item{genesToTest}{Vector of gene names to consider plotting.}
\item{propThreshold}{Proportion threshold above which to plot gene/arm labels.}
\item{countThreshold}{Count threshold above which to plot gene/arm labels.}
\item{ensembl}{Biomart ensembl object to look up gene coordinates. If NULL, geneGR and whichGenes will be used instead.}
\item{splits}{Index of windows that split separate chromosomes. See setupWindows().}
\item{mids}{index of windows that mark the midpoints of chromosomes. See setupWindows().}
\item{saveRes}{Boolean. Whether to save the window and recurrence information as a csv.}
\item{geneGR}{GenomicRanges object denoting genomic positions of genes listed in genesToTest.}
\item{windowGR}{GenomicRanges object denoting genomic windows in which recurrence has been evaluated. See setupWindows().}
\item{doCombine}{Output of getRecurrence is sometimes a list. If so, doCombine=TRUE will convert to a matrix internally.}
\item{doSex}{Boolean. Whether to plot sex chromosomes.}
\item{doLegend}{Boolean. Whether to plot a legend.}
\item{doMain}{Boolean. Whether to plot title.}
\item{LTY}{Vector of line types for each signature to be plotted.}
\item{YLIM}{Custom values for plotting y-axis limits.}
\item{LABY}{Custom label for the y-axis.}
\item{offsetBase}{Y offset for plotting gene labels.}
\item{offsetTweak}{X offset for plotting gene labels.}
\item{xoverlap}{Overlap required for automated moving of overlapping gene labels (x-axis).}
\item{yoverlap}{Overlap required for automated moving of overlapping gene labels (y-axis).}
\item{myMAR}{Custom margins for plotting.}
\item{XLIMadj}{Extension to the right hand xlim value.}
\item{geneCex}{cex value for plotting gene labels.}
\item{legloc}{Location for legend.}
\item{doGenesAbove}{Boolean. Whether to plot gene names on plot (FALSE), or ina separate panel above the plot (TRUE).}
\item{mylayout}{Custom layout for multi-panel figure (When doGenesAbove==TRUE).}
\item{ylineoffset}{Offset for lines connecting gene names to plot.}
\item{pdfWidth}{Width for output pdf.}
\item{pdfHeight}{Height for output pdf.}
\item{CEX}{cex value for plot text.}
\item{doSplit}{Boolean. Whether to split plot into two panels (Only use when plotting two signatures).}
\item{legnames}{Names for legens.}
\item{propArm}{Proportion of chromosome arm that must be above propThreshold to plot arm label.}
\item{geneYtweak}{Offset for plotting gene labels (y-axis).}
\item{extraMain}{Subtitle.}
}
\value{
Pdf if filename is not NULL, csv if saveRes is not FALSE.
} 
\description{
Plot recurrence of copy number signature across the genome, labelling genes of interest or chromosome arms that show highly recurrent signature attribution. 
}
\details{
  
}

