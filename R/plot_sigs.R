####################################################################
#     functions to plot signature definitions
####################################################################

# function to plot a group of signatures
plotSig = function(sigs, # signatures to plot
                   doText=TRUE, # whether to add text labels to coloured bars (coloured bars = legend)
                   doBar=TRUE, # whether to plot coloured bars on each plot, do not use if plotting legend separately
                   doTitleBar=FALSE,# whether to plot title as bar
                   doTitle=TRUE, # whether to plot title for each sig, deprecated since horizontal
                   doMulti=TRUE, # whether to plot multiple signatures
                   doLab=FALSE, # whether to plot x labels (CN class labels)
                   barSize=1.15, # size of coloured bars, if plotted
                   doYlab=TRUE, # whether to ploot y-axis label (name of signature)
                   myMAR = NULL,
                   YLIM=NULL,
                   CEX=1.8,
                   mainAdj=0.5,mainLine=0,doBox=FALSE,
                   extraLabel=NULL,extraAdj=0.5,extraAt=15,extraCex=0.5
)
{
  print(colnames(sigs))
  print(doLab)
  if(any(is.na(doLab))) doLab[which(is.na(doLab))]=FALSE
  # roeorder and get colours for plotting
  factors = sapply(rownames(sigs),FUN=function(x) strsplit(x,split=":")[[1]])
  sigs = sigs[order(factor(factors[2,],levels=c("homdel","LOH","het")),
                    factor(factors[1,],levels=c("0","1","2","3-4","5-8","9+")),
                    factor(factors[1,],levels=c("0-100kb","100kb-1Mb","1Mb-10Mb","10Mb-40Mb",">40Mb",">1Mb"))),,drop=FALSE]
  colours = vector(length=nrow(sigs))
  delColours = colorRampPalette(c("gray92","gray33"))(5)
  homdelColours = colorRampPalette(c("aliceblue","blue3"))(3)
  neutColours = colorRampPalette(c("mintcream","forestgreen"))(5)
  dupColours = colorRampPalette(c("lavenderblush","purple3"))(5)
  quadColours = colorRampPalette(c("floralwhite","orange3"))(5)
  ampColours = colorRampPalette(c("mistyrose","deeppink4"))(5)
  colours[grep("^1[:]LOH",rownames(sigs))] = delColours
  colours[grep("^0",rownames(sigs))] = homdelColours
  colours[grep("^2",rownames(sigs))] = neutColours
  colours[grep("^3-4",rownames(sigs))] = dupColours
  colours[grep("^5-8",rownames(sigs))] = quadColours
  colours[grep("^9+",rownames(sigs))] = ampColours
  if(doMulti)
  {
    MFROW=c(ncol(sigs),1)
    par(mfrow=MFROW)
  } 
  CNbreaks = c(1.2*c(0,3,(5*(1:9))+3))
  LOHbreaks = CNbreaks[c(1,2,7,11)]
  # loop over signatures
  sapply(1:ncol(sigs),FUN=function(x)
  {
    # set up plotting parameters
    if(doYlab&doLab[x])
       {
         if(is.null(myMAR)) MAR=c(6,2,0,0)
         MGP = c(0,1,0)
    } else if(doYlab) {
      if(is.null(myMAR)) MAR=c(0,2,0,0)
      MGP = c(0,1,0)
    } else if(doLab[x]) {
      if(is.null(myMAR)) MAR=c(6,0,0,0)
      MGP = c(3,1,0)
    } else {
      if(is.null(myMAR)) MAR=c(0,0,0,0)
      MGP = c(3,1,0)  
    }
    if(!is.null(myMAR)) MAR = myMAR
    par(mar=MAR,mgp=MGP)
    # plot signatures
    if(all(is.na(sigs[,x])))
      {
      plot(NA,xlim=0:1,ylim=0:1,axes=FALSE,xlab=NA,ylab=NA)
      return(NA)
      }
    # setup xlim so can plot coloured bar
    if(is.null(YLIM)) YLIM = c(0,max(sigs[,x]))
    if(doBar|doTitleBar|doTitle)
    {
      YLIM[2] = c(YLIM[2]*barSize)
      barPos = (barSize-1)/3
      barPos = barPos*c(1:3)
      barPos = barPos+1
    }
    # setup title
    if(doTitle|doTitleBar)
    {
      titlestring=colnames(sigs)[x]
      if(!is.null(extraLabel))
        {
        extrastring = extraLabel[x]
      } else {
      extrastring = NA  
      }
    } else {
      titlestring=NA 
      extrastring=NA
    }
    # setup labels
    if(doLab[x])
    {
      labnames=rownames(sigs)
      labnames = sapply(labnames,FUN=function(x) strsplit(x,split=":")[[1]][3])
      library(stringr)
      replacements = c("0-100kb","100kb-1Mb",">1Mb","1Mb-10Mb","10Mb-40Mb",">40Mb")
      names(replacements) = c("[(]-0.01,0.1[]]","[(]0.1,1[]]","[(]1,Inf[]]","[(]1,10[]]","[(]10,40[]]","[(]40,Inf[]]")
      labnames = str_replace_all(string=labnames,
                                 pattern=replacements)
    } else {
      labnames=NA
    }
    # setup y label
    if(doYlab)
    {
      #YLAB=colnames(sigs)[x]
      YLAB=NA
    } else {
      YLAB=NA
    }
    # plot sigs
    bp = barplot(sigs[,x],horiz=FALSE,las=1,col=colours,
                 main=NA,names.arg=labnames,
                 las=3,yaxt="n",ylab=YLAB,ylim=YLIM,
                 cex.main=CEX,cex.lab=CEX,cex.names=CEX)
    if(doBox) box()
    if(doTitle)
      {
      title(main=titlestring,adj=mainAdj,line=mainLine,cex.main=CEX*1.5)
      mtext(text=extrastring,adj=extraAdj,line=mainLine,cex=CEX*extraCex,font=1,side=3,at=extraAt)
      }
    if(doYlab)
      {
      maxval = max(sigs[,x])
      maxval = round(maxval/5,2)*5
      yticks = seq(from=0,to=maxval,by=0.05)
      ylabels = yticks
      ylabels[seq(from=2,to=length(yticks),by=2)] = NA
      axis(side=2,cex=CEX,cex.axis=CEX,cex.lab=CEX,
           at=yticks[-which(is.na(ylabels))],
           labels=ylabels[-which(is.na(ylabels))],
           tck=-0.05
           )
      axis(side=2,cex=CEX,cex.axis=CEX,cex.lab=CEX,
           at=yticks,labels=NA,tck=-0.025)
      }
    print(bp)
    #if(!doTitle) abline(v=LOHbreaks[2:3])
    # set up coloured basr parameters
    breaks=max(sigs[,x])
    if(doBar|doTitleBar) breaks = c(breaks*barPos)
    if(doBar)
    {
      # plot LOH info bars
      #LOHcolours = c("blue3","black","gray")
      LOHcolours = c("blue3","white","gray")
      sapply(1:(length(LOHbreaks)-1),FUN=function(x)
      {
        polygon(y=c(breaks[2],breaks[3],breaks[3],breaks[2]),
                x=c(LOHbreaks[x],LOHbreaks[x],LOHbreaks[x+1],LOHbreaks[x+1]),
                col=LOHcolours[x])
      })
      # plot CN info bars
      CNcolours = c("blue3","gray33","forestgreen","purple3","orange3","deeppink4","forestgreen","purple3","orange3","deeppink4")
      sapply(1:(length(CNbreaks)-1),FUN=function(x)
      {
        X = c(breaks[1],breaks[2],breaks[2],breaks[1])
        if(x==1) X = c(breaks[1],breaks[3],breaks[3],breaks[1])
        polygon(y=X,
                x=c(CNbreaks[x],CNbreaks[x],CNbreaks[x+1],CNbreaks[x+1]),
                col=CNcolours[x])
      })
      if(doText)
      {
        # plot LOH bar labels
        LOHlabs = c("Hom","LOH","Het")
        sapply(1:(length(LOHbreaks)-1),FUN=function(i)
        {
          y = mean(c(LOHbreaks[c(i,i+1)]))
          x = breaks[2]+(breaks[3]-breaks[2])/2
          #labcol="white"
          text(x=y,y=x,labels=LOHlabs[i],col=labcol,cex=CEX)#,srt=270
        })
        # plot CN bar labels
        CNlabs = CNcolours = c("del","Del","Neut","Gain","Amp","Amp+","Neut","Gain","Amp","Amp+")
        sapply(1:(length(CNbreaks)-1),FUN=function(i)
        {
          y = mean(c(CNbreaks[c(i,i+1)]))
          x = breaks[1]+(breaks[2]-breaks[1])/2
          text(x=y,y=x,labels=CNlabs[i],col="white",cex=CEX)#,srt=270
        })
      }
    }
    if(doTitleBar)
    {
      # plot title bar
      #polygon(y=c(breaks[1],breaks[3],breaks[3],breaks[1]),
      #        x=LOHbreaks[c(1,1,4,4)],
      #        col="gray",border=NA)
      text(x=mean(LOHbreaks[c(1,4)]),y=mean(breaks[c(1,3)]),labels=title,cex=CEX*1.5,col="black",font=1)
      #lines(x=rep(LOHbreaks[2],2),y=c(0,breaks[1]))
      #lines(x=rep(LOHbreaks[3],2),y=c(0,breaks[1]))
      #abline(h=breaks[1])
    }
  })
}

# function to plot a legend for the signatures
plotLegend = function(sigs,doText=TRUE,doYlab=TRUE,doRev=FALSE,doTitle=FALSE,myMAR=NULL,CEX=1.8,HDlab="")
{
  # roeorder and get colours for plotting
  factors = sapply(rownames(sigs),FUN=function(x) strsplit(x,split=":")[[1]])
  sigs = sigs[order(factor(factors[2,],levels=c("homdel","LOH","het")),
                    factor(factors[1,],levels=c("0","1","2","3-4","5-8","9+")),
                    factor(factors[1,],levels=c("0-100kb","100kb-1Mb","1Mb-10Mb","10Mb-40Mb",">40Mb",">1Mb"))),,drop=FALSE]
  colours = vector(length=nrow(sigs))
  delColours = colorRampPalette(c("gray92","gray33"))(5)
  homdelColours = colorRampPalette(c("aliceblue","blue3"))(3)
  neutColours = colorRampPalette(c("mintcream","forestgreen"))(5)
  dupColours = colorRampPalette(c("lavenderblush","purple3"))(5)
  quadColours = colorRampPalette(c("floralwhite","orange3"))(5)
  ampColours = colorRampPalette(c("mistyrose","deeppink4"))(5)
  colours[grep("^1[:]LOH",rownames(sigs))] = delColours
  colours[grep("^0",rownames(sigs))] = homdelColours
  colours[grep("^2",rownames(sigs))] = neutColours
  colours[grep("^3-4",rownames(sigs))] = dupColours
  colours[grep("^5-8",rownames(sigs))] = quadColours
  colours[grep("^9+",rownames(sigs))] = ampColours
  # set up plotting parameters
  if(doYlab)
  {
    if(is.null(myMAR)) MAR=c(0,2,0,0)
  } else {
    if(is.null(myMAR)) MAR=c(0,0,0,0)  
  }
  if(!is.null(myMAR)) MAR = myMAR
  # set up plotting parameters
  par(mar=MAR)
  CNbreaks = c(1.2*c(0,3,(5*(1:9))+3))
  CNbreaks[1] = 0.1
  LOHbreaks = CNbreaks[c(1,2,7,11)]
  XLIM = range(CNbreaks)
  #XLIM[1] = XLIM[1]+0.01
  plot(NA,ylim=0:1,xlim=XLIM,xlab=NA,ylab=NA,axes=FALSE,bty="n")
  # set up coloured bars parameters
  breaks=c(0,0.5,1)
  if(doRev) breaks=rev(breaks)
  # plot LOH info bars
  #LOHcolours = c("blue3","black","gray")
  LOHcolours = c("blue3","white","gray10")
  sapply(1:(length(LOHbreaks)-1),FUN=function(x)
  {
    polygon(y=c(breaks[2],breaks[3],breaks[3],breaks[2]),
            x=c(LOHbreaks[x],LOHbreaks[x],LOHbreaks[x+1],LOHbreaks[x+1]),
            col=LOHcolours[x])
  })
  # plot CN info bars
  CNcolours = c("blue3","gray33","forestgreen","purple3","orange3","deeppink4","forestgreen","purple3","orange3","deeppink4")
  sapply(1:(length(CNbreaks)-1),FUN=function(x)
  {
    X = c(breaks[1],breaks[2],breaks[2],breaks[1])
    if(x==1) X = c(breaks[1],breaks[3],breaks[3],breaks[1])
    polygon(y=X,
            x=c(CNbreaks[x],CNbreaks[x],CNbreaks[x+1],CNbreaks[x+1]),
            col=CNcolours[x])
  })
  if(doText)
  {
    # plot LOH bar labels
    #LOHlabs = c(ifelse(doRev,"del","Hom"),"LOH","Het")
    LOHlabs = c(HDlab,"LOH","Het")
    sapply(1:(length(LOHbreaks)-1),FUN=function(i)
    {
      y = mean(c(LOHbreaks[c(i,i+1)]))
      x = breaks[2]+(breaks[3]-breaks[2])/2
      #labcol = "white"
      labcol = c("white","black","white")
      text(y=x,x=y,labels=LOHlabs[i],col=labcol[i],cex=CEX)
    })
    # plot CN bar labels
    #CNlabs = CNcolours = c(ifelse(doRev,"Hom","del"),"Del","Neut","Gain","Amp","Amp+","Neut","Gain","Amp","Amp+")
    CNlabs = CNcolours = c("0","1","2","3-4","5-8","9+","2","3-4","5-8","9+")
    sapply(1:(length(CNbreaks)-1),FUN=function(i)
    {
      y = mean(c(CNbreaks[c(i,i+1)]))
      x = breaks[1]+(breaks[2]-breaks[1])/2
      text(y=x,x=y,labels=CNlabs[i],col="white",cex=CEX)
    })
  }
}

# function that plots multiple signatures, with the legends separately
plotSigDefinitions = function(sigs,sigLength=2,doTitle=TRUE,doTitleBar=FALSE,doText=FALSE,
                               doYlab=TRUE,doLab=FALSE,ncol=1,nPerCol="all",plotBars=c(1,2),
                               myMAR=NULL,barSize=1.15,returnLayout=FALSE,myLayout=NULL,CEX=1.8,
                               YLIM=NULL,mainAdj=0.5,mainLine=0,doBox=FALSE,HDlab="",
                               extraLabel=NULL,extraAdj=0.5,extraAt=15,extraCex=0.5)
{
  if(length(doLab)==1) doLab=rep(doLab,ncol(sigs))
  if(length(doText)==1) doText=rep(doText,2)
  if(ncol>1)  # plotting for multiple columns
  {
    if(!is.numeric(nPerCol)) stop("Specify a numeric nPerRow if nRow>1")
    # setup layout of plots
    if(is.null(myLayout))
    {
      if(1%in%plotBars)
      {
        layoutVec = c(1)
      } else {
        layoutVec = c()
      }
      layoutVec = c(layoutVec,rep((1:nPerCol)+(1%in%plotBars),each=sigLength))
      if(2%in%plotBars) layoutVec = c(layoutVec,nPerCol+length(plotBars))
      layoutMat = matrix(layoutVec,ncol=1)
      layoutMatAll = layoutMat
      for(i in 2:ncol)
      {
        #print(i)
        layoutMatAll=cbind(layoutMatAll,
                           layoutMat+((i-1)*max(layoutMat))
        )
      }
      layoutMat = layoutMatAll
      if(returnLayout) return(layoutMatAll)
    } else {
      layoutMat = myLayout
    }
    layout(layoutMat)
    # loop over columns
    for(i in 1:ncol)
    {
      # plot first legend
      if(1%in%plotBars) plotLegend(sigs,doText=doText[1],doYlab=doYlab,doRev=FALSE,myMAR=myMAR,CEX=CEX,HDlab=HDlab)
      # plot signatures
      plotindex = (1:nPerCol)+((i-1)*nPerCol)
      plotindex = plotindex[which(plotindex<=ncol(sigs))]
      plotSig(sigs[,plotindex,drop=FALSE],doMulti=FALSE,
              doTitleBar=doTitleBar,doBar=FALSE,
              doTitle=doTitle,doYlab=doYlab,
              doLab=doLab[plotindex],myMAR=myMAR,
              barSize=barSize,CEX=CEX,YLIM=YLIM,
              mainAdj=mainAdj,mainLine=mainLine,
              doBox=doBox,
              extraLabel=extraLabel[plotindex],extraAdj=extraAdj,
              extraAt=extraAt,extraCex=extraCex)
      # plot secong legend
      if(2%in%plotBars) plotLegend(sigs,doText=doText[2],doYlab=doYlab,doRev=TRUE,myMAR=myMAR,CEX=CEX,HDlab=HDlab)
    }
  } else { # plotting for a single column
    if(is.null(myLayout))
    {
    # setup layout first
    if(1%in%plotBars)
    {
      layoutVec = c(1)
    } else {
      layoutVec = c()
    }
    layoutVec = c(layoutVec,rep((1:ncol(sigs))+(1%in%plotBars),each=sigLength))
    if(2%in%plotBars) layoutVec = c(layoutVec,ncol(sigs)+length(plotBars))
    print(layoutVec)
    layoutMat = matrix(layoutVec,ncol=1)
    layoutMat = cbind(layoutMat,layoutMat,layoutMat)
    } else {
    layoutMat = myLayout  
    }
    layout(layoutMat)
    # plot first legend
    if(1%in%plotBars) plotLegend(sigs,doText=doText[1],doTitle=doTitle,doRev=FALSE,doYlab=doYlab,CEX=CEX,myMAR=c(0,myMAR[2],0,myMAR[4]),HDlab=HDlab)
    # plot signatures
    plotSig(sigs,doMulti=FALSE,doBar=FALSE,doTitle=doTitle,
            doYlab=doYlab,doLab=doLab,CEX=CEX,myMAR=myMAR,
            YLIM=YLIM,mainAdj=mainAdj,mainLine=mainLine,
            doBox=doBox,
            extraLabel=extraLabel,extraAdj=extraAdj,
            extraAt=extraAt,extraCex=extraCex)
    # plot second legend
    if(2%in%plotBars) plotLegend(sigs,doText=doText[2],doTitle=doTitle,doRev=TRUE,doYlab=doYlab,CEX=CEX,myMAR=c(0,myMAR[2],0,myMAR[4]),HDlab=HDlab)
  }
  #layoutMat = rbind(layoutMat,layoutMat,layoutMat)
  
}

