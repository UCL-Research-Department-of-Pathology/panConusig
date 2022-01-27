# =====================================================================================
#			function to plot an allele specific copy number profile
# =====================================================================================
plotCNprofile = function(chr,start,end,majCN,minCN,main=NA,colours=NULL,YLIM=NULL,offsetVal=100,LWD=7,CEX=3,
                         custom_legend=NULL)
{
  layout(matrix(c(1,1,1,1,1,1,1,1,2),nrow=1))
  # setup
  par(mar=c(1,7,3,0))
  offset = max(YLIM)/offsetVal
  lengths = end-start
  cumlengths = cumsum(as.numeric(lengths))
  xlim=c(1,max(cumlengths))
  if(is.null(YLIM)) YLIM=c(0,quantile(majCN,0.95)+1.5)
  plotstarts = c(1,cumlengths[1:(length(cumlengths)-1)])
  plotends = cumlengths
  if(is.null(colours))
  {
    coloursMaj=rep("navyblue",length(start))
    coloursMin=rep("darkorange",length(start))
  } else if(!is.list(colours)) {
    coloursMaj = colours[1]
    coloursMin = colours[2]#apply(col2rgb(colours),MARGIN=2,FUN=function(x) rgb(x[1]/255,x[2]/255,x[3]/255,0.5))
  } else {
    coloursMaj = colours[[1]]
    coloursMin = colours[[2]]
  }
  if(!is.na(main)) main = paste("ASCAT profile: ",main)
  # plot CN
  plot(NA,xlim=xlim,ylim=YLIM,ylab="Copy number",
       bty="n",xaxt="n",xlab=NA,main=main,
       cex=CEX,cex.main=CEX,cex.lab=CEX,cex.axis=CEX,xaxs="i")
  # chromsome labels
chroms = unique(chr)
col = c("white",rgb(0.3,0.3,0.3,0.3))
colLabs = c("white","black")
sapply(chroms,FUN=function(x)
	{
	index = which(chr==x)
	X = c(plotstarts[min(index)],plotends[max(index)])
	X = c(X,rev(X))
	polygon(x=X,y=c(-100,-100,100,100),col=col[((which(chroms==x)%%2)==0)+1],border=NA)
	text(x=mean(X),y=YLIM[2]-(diff(range(YLIM))*0.05),labels=x,col=colLabs[(!(which(chroms==x)%%2)==0)+1],
	     cex=CEX/1.75)
	})
  # add segments
if(!is.list(colours))
  {
  sapply(1:length(start),FUN=function(x) 
  {
    lines(x=c(plotstarts[x],plotends[x]),y=rep(majCN[x],2)-offset,col=coloursMaj,lwd=LWD,lend="butt")
    lines(x=c(plotstarts[x],plotends[x]),y=rep(minCN[x],2)+offset,col=coloursMin,lwd=LWD,lend="butt")
  })
} else {
  sapply(1:length(start),FUN=function(x) 
  {
    lines(x=c(plotstarts[x],plotends[x]),y=rep(majCN[x],2)-offset,col=coloursMaj[x],lwd=LWD,lend="butt")
    lines(x=c(plotstarts[x],plotends[x]),y=rep(minCN[x],2)+offset,col=coloursMin[x],lwd=LWD,lend="butt")
  })
}
  #abline(v=sapply(unique(chr),FUN=function(x) plotstarts[min(which(chr==x))]),
  #       col="gray")
  abline(h=0:100,lty=3,col="gray")
  # plot legend
  par(mar=c(0,0,0,0))
  plot(NA,xlim=0:1,ylim=0:1,axes=FALSE,xlab=NA,ylab=NA)
  if(is.null(custom_legend))
    {
    legend(x=0,y=0.8,legend=c("Major allele","Minor allele"),
         x.intersp=0.7,y.intersp=2,pch=15,col=c("navyblue","darkorange"),
         cex=CEX,bty="n")
  } else {
    do.call(legend,custom_legend)
  }
}

