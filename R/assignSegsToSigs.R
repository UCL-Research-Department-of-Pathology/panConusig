###############################################################################
#functions to assign individual segments of a copy number profile to signatures
###############################################################################


# function to plot a CN profile with signature attributions below
plotCNprofile = function(chr,start,end,majCN,minCN,CNS,nSigs=10,main=NA,colours=NULL,YLIM=NULL,offsetVal=100)
{
  par(mfrow=c(2,1))
  # setup
  par(mar=c(1,5,1,1))
  offset = max(majCN)/offsetVal
  lengths = end-start
  cumlengths = cumsum(as.numeric(lengths))
  xlim=c(1,max(cumlengths))
  if(is.null(YLIM)) YLIM=c(0,quantile(majCN,0.95)+1.5)
  plotstarts = c(1,cumlengths[1:(length(cumlengths)-1)])
  plotends = cumlengths
  if(is.null(colours))
  {
    coloursMaj=rep("red",length(start))
    coloursMin=rep("green",length(start))
  } else {
    coloursMaj = colours
    coloursMin = colours#apply(col2rgb(colours),MARGIN=2,FUN=function(x) rgb(x[1]/255,x[2]/255,x[3]/255,0.5))
  }
  # plot CN
  plot(NA,xlim=xlim,ylim=YLIM,ylab="CN",bty="n",xaxt="n",xlab=NA,main=main,cex=2,cex.main=2,cex.lab=2,cex.axis=2)
  sapply(1:length(start),FUN=function(x) 
  {
    lines(x=c(plotstarts[x],plotends[x]),y=rep(majCN[x],2)+offset,col=coloursMaj[x],lwd=7,lend="butt")
    lines(x=c(plotstarts[x],plotends[x]),y=rep(minCN[x],2)-offset,col=coloursMin[x],lwd=7,lend="butt")
  })
  abline(v=sapply(unique(chr),FUN=function(x) plotstarts[min(which(chr==x))]),
         col="gray")
  abline(h=0:100,lty=3,col="gray")
  getMid = function(x)
  {
    index = which(chr==x)
    plotx = min(plotstarts[index])+(max(plotends[index])-min(plotstarts[index]))/2
  }
  text(x=sapply(unique(chr),FUN=getMid),y=YLIM[2]-0.25,labels=unique(chr))
  # plot sig
  par(mar=c(1,5,1,1))
  plot(NA,xlim=xlim,ylim=c(0.5,nSigs+0.5),ylab=NA,yaxt="n",bty="n",xaxt="n",xlab=NA,main=main,cex=2,cex.main=2,cex.lab=2,cex.axis=2)
  sapply(1:length(start),FUN=function(x) 
  {
    thisSig = which(colnames(sigs)==CNS[x])
    polygon(x=c(plotstarts[x],plotends[x],plotends[x],plotstarts[x]),y=c(rep(thisSig-0.5,2),rep(thisSig+0.5,2)),col="black")
  })
  abline(v=sapply(unique(chr),FUN=function(x) plotstarts[min(which(chr==x))]),
         col="gray")
  abline(h=(0:nSigs)+0.5,col="gray")
  flags = xlim
  xlimlength=abs(diff(xlim))
  flags[1]=flags[1]+(0.05*xlimlength)
  flags[2]=flags[2]-(0.05*xlimlength)
  flag = flags[1]
  for(i in 1:nSigs)
  {
    text(x=flag,y=i,labels=colnames(sigs)[i],col="gray",cex=2)
    flag=flags[which(!flags%in%flag)]
  }
}


# assign segs to signatures per chromosome
getSigsPerChrom = function(samp,sigs,exposures,data,doLength=TRUE,outDir=NULL)
{
  print(samp)
  sampIndex = which(data$sample==samp)
  nsegs=length(sampIndex)
  # plot CN profile
  #if(!is.null(outDir)) pdf(paste0(outDir,"/",samp,"-CNprofile.pdf"),width=14)
  #plotCNprofile(data$chr[sampIndex],
  #            data$startpos[sampIndex],
  #            data$endpos[sampIndex],
  #            data$nMajor[sampIndex],
  #            data$nMinor[sampIndex])
  #if(!is.null(outDir)) dev.off()
  # number of segments attributable to each CN class
  nummuts = sapply(1:nrow(sigs),FUN=function(x) sum(sigs[x,]*exposures[,samp]*nsegs))
  sum(nummuts)
  sum(data[,1]==samp)
  # probability of each CN class to come from each signature
  probs = sapply(1:nrow(sigs),FUN=function(x) 
  {
    numerator=(sigs[x,]*exposures[,samp]*nsegs)
    denominator=sum(numerator)
    if(denominator==0) return(rep(0,length(numerator)))
    numerator/denominator
  })
  rownames(probs) = colnames(sigs)
  colnames(probs) = rownames(sigs)
  # probability of each signature in each chromosome
  chromProbs = sapply(unique(paste0(data$chr[sampIndex])),FUN=function(x)
  {
    #print(x)
    index = which(data$sample==samp&data$chr==x)
    # for each segment, the probability of each signature producing that segment type
    # normalised by length of segment so a few tiny segs dont swamp the chrom result
    sepProbs = sapply(index,FUN=function(y)
    {
      out=unlist(probs[,data$CNvar[y]])
      if(doLength) out = out*data$lengths[y]
      out
    })
    # sum across the whole chromosome for each sig
    rowSums(sepProbs)
  })
  #print(chromProbs)
  # plot chrom probs
  #library(ComplexHeatmap)
  #if(!is.null(outDir)) pdf(paste0(outDir,samp,"-chromProbs-heatmap.pdf"))
  #hm = Heatmap(chromProbs,name="CNS\nprobability")
  #draw(hm)
  #if(!is.null(outDir)) dev.off()
  write.csv(chromProbs,file=paste0(outDir,samp,"-chromProbs.csv"),quote=FALSE)
  chromProbs
}

# function to run assignment
# warning can take a long time with many samples
runAssignment = function(data, # matrix of copy number profiles with columns sample, chr, startpos, endpos, nMajor, nMinor
                         sigs, # matrix of 48xn where n is the number of signatures
                         exposures, # matrix of n*m where n is the number of signatures, and m is the number of samples
                         outDir = NULL
                         )
  {
  # get CN class for each segment
  variables = sapply(unique(data$sample),FUN=function(x)
    {
    print(x)
    index = which(data$sample==x)
    copyNumberInfo(chrom=data$chr[index],
                 start=data$startpos[index],
                 end=data$endpos[index],
                 CN=data$nMajor[index]+data$nMinor[index],
                 CNb=data$nMinor[index],
                 doVariables=TRUE,
                 returnSep=TRUE)
    })
  variables = unlist(variables)
  variables = gsub("0[-]1[:]homdel","0:homdel",variables)
  variables = gsub("0[-]1[:]LOH","1:LOH",variables)
  data$CNvar = variables
  # get probabilities of each CNclass in each chromosome
  require(parallel)
  chromProbs = mclapply(colnames(exposures),FUN=function(x) 
    {
    getSigsPerChrom(samp=x,
                  sigs=sigs,
                  exposures=exposures,
                  outDir=outDir,
                  doLength = FALSE)
    },mc.cores=3)#,simplify=FALSE)#,mc.cores=8)
  names(chromProbs) = colnames(exposures)
  # now should be able to assign each segment in a chrom to a signature, based on its chrom level probs
  assignment = vector(length=nrow(data))
  for(samp in colnames(exposures))
    {
    print(samp) 
    for(chrom in unique(data[which(data[,"sample"]==samp),"chr"]))
      {
      print(chrom)
      index = data[,"sample"]==samp&data[,"chr"]==chrom
      nummutsPerSig = chromProbs[[samp]][,chrom]*sum(index)
      typeProbs = sapply(1:length(nummutsPerSig),FUN=function(x) nummutsPerSig[x]*sigs[,x])
      segAssignments = colnames(sigs)[apply(typeProbs,MARGIN=1,FUN=which.max)]
      names(segAssignments) = rownames(sigs)
      #cbind(data[which(index),],segAssignments[data$CNvar[which(index)]])
      assignment[which(index)] = segAssignments[data$CNvar[which(index)]]
      }
    }
  data$CNS = assignment
  write.table(data,file=paste0(outDir,"/seg-with-CNS-assignments.csv"),
              quote=FALSE,row.names=FALSE,sep="\t")
  data
  }

