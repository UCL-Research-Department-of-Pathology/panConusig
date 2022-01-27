# =============================================================================
#	functions to look for signatures recurrence across the genome
# depends on results from assignSegsToSigs.R
# =============================================================================
# setup windows
setupWindows = function(chr,startpos,endpos,windowLength = 1000000,doBiomart=FALSE)
	{
	# get windows
  print("get windows")
	windows = sapply(unique(chr),FUN=function(x)
		{
 		index = which(chr==x)
  		start = min(startpos[index])
  		end = max(endpos[index])
  		starts = seq(from=start,to=end,by=windowLength)
  		ends = starts[2:length(starts)]
  		starts=starts[-length(starts)]
  		cbind(x,starts,ends)
		})
	windows = do.call(rbind,windows)
	# convert to GRanges
	print("window GRs")
	library(GenomicRanges)
	windowsGR = as(paste0(windows[,1],":",windows[,2],"-",windows[,3]),"GRanges")
	# plotting params
	print("plotting params")
	splits = which(diff(as.numeric(factor(windows[,1],levels=c(1:22,"X","Y"))))>0)
	splits = c(0,splits,nrow(windows))
	mids = sapply(1:(length(splits)-1),FUN=function(x) ((splits[x+1]-splits[x])/2)+splits[x])
	# setup biomart
	print("biomart")
	library(biomaRt)
	if(doBiomart)
		{
		ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
		} else {
		ensembl = NULL
		}
	# return
	print("return")
	return(list(windows=windows,
			windowsGR=windowsGR,
			splits=splits,
			mids=mids,
			ensembl=ensembl))
	}


# get recurrence
getRecurrence = function(sample,chrom,start,end,CNS,
                         windowsGR,sigNames,
                         tumType, # named vector of tumour type per sample
                         ncores=8,returnCounts=FALSE)
	{
  # initial setup
  #print("segGR and overlap")
	segGR = as(paste0(chrom,":",start,"-",end),"GRanges")
	overlaps = findOverlaps(windowsGR,segGR)
	# counting
	#print("count matches in windows")
	library(parallel)
	# loop over windows
	counts = mclapply(1:length(windowsGR),FUN=function(x) {
  		#print(x)
 		index = overlaps@to[which(overlaps@from==x)]
 		# loop over signatures
  		counts = sapply(sigNames, FUN=function(y)
  			{
  		    # loop over tumour types
    			sapply(unique(tumType),FUN=function(z) 
      				{
    			    # which segs match signature
      				matchCNindex = CNS[index]==y
      				# which segs match tumour type
      				matchTumTypeIndex = sample[index]%in%names(tumType)[which(tumType==z)]
      				# how many samples match both tumour type and signature
      				matches = length(unique(sample[index[which(matchCNindex&matchTumTypeIndex)]]))
      				# how many samples match tumour type only
      				total = length(unique(sample[index[which(matchTumTypeIndex)]]))
      				if(returnCounts)
					{
					return(matches)
					} else {
					matches/total
					}
    				})
  			})
  		counts
		},mc.cores=ncores)
	counts
	}

# get result into useable format
getRes = function(x,tumType="osteo") {x[tumType,]}

# plot results
plotRecurrence = function(props,counts,tumTypes,windows,arms,doGenes=FALSE,doArms=TRUE,plotGenes=TRUE,outDir="~/",
		filename=NULL,whichSigs="CN19",
		colours = c("red","blue","green","black","gray","orange","cyan","purple"),
		whichGenes="CN19",genesToTest="MDM2",propThreshold=0.2,countThreshold=5,
		ensembl=NULL,splits=NA,mids=NA,saveRes=FALSE,geneGR=NULL,windowGR=NULL,
		doCombine=TRUE,doSex=FALSE,doLegend=TRUE,doMain=TRUE,LTY=1:5,YLIM=NULL,LABY=NULL,
		offsetBase=0.01,offsetTweak=0.02,xoverlap=100,yoverlap=0.02,myMAR=NULL,XLIMadj=0,
		geneCex=2,legloc="topright",doGenesAbove=TRUE,mylayout=NULL,ylineoffset=NULL,
		pdfWidth=21,pdfHeight=7,CEX=4,doSplit=FALSE,legnames=NULL,propArm=0.75,
		geneYtweak=0,extraMain=NULL
		) 
	{
  print("setup")
	if(!is.null(filename)) pdf(paste0(outDir,filename),width=pdfWidth,height=pdfHeight)
  if(doGenesAbove)
    {
    if(is.null(mylayout))
      {
      mylayout = matrix(c(1,1,1),ncol=1)
      if(doGenes)
      {
        mylayout = rbind(2,mylayout)
      if(doSplit) mylayout = rbind(mylayout,c(3,4))
      }
      }
    layout(mylayout)
    }
	if(is.null(myMAR)) myMAR=c(3,8,0,0)
	if(doSplit) myMAR[1] = 0
	par(mar=myMAR,xaxs="i")
	# per tumour type
	for(i in sort(tumTypes))
		{
  		print(i)
		if(doMain) {MAIN=i} else {MAIN=NA}
		# get result into useable format
		print("Combine props")
		if(doCombine)
		  {
		  res = t(sapply(props,FUN=getRes,tumType=i))
		  } else {
		  res = props
		  }
		print("Combine counts")
		if(doCombine)
		  {
		  resC = t(sapply(counts,FUN=getRes,tumType=i))
		  } else {
		  resC = counts
		  }
		print("Check for NAs")
		if(all(is.na(res[,whichSigs])))
		  {
		  print("Skipping: all NA")
		  next
		}
		print("Sort out sex chrom")
		# plot sex chrom?
		labs = unique(windows[,1])
		if(doSex)
			{
		  print("If plot sex")
			index = 1:nrow(res)
			} else {
			print("If not plot sex")
			index = (1:nrow(res))[-which(windows[,1]%in%c("X","Y"))]
			remove = which(labs%in%c("X","Y"))
			labs = labs[-remove]
			mids = mids[-remove]
			splits = splits[-(remove+1)]
			}
		# plot sig recurrences
		print("plot recurrence")
		XLIM = c(0,nrow(res[index,]))
		if(doLegend) XLIM[2] = XLIM[2]+55
		XLIM[2] = XLIM[2]+XLIMadj
		par(mgp=c(3,1,0))
		if(is.null(YLIM))
		  {
		  YLIM = c(0,max(res[index,whichSigs],na.rm=TRUE)+0.05)
		  if(doSplit)
		    {
		    YLIM[1] = -YLIM[2]
		    }
		  }
		par(mgp=c(5,1,0))
		if(doSplit) res[,whichSigs[2]] = -res[,whichSigs[2]]
		# empty plot
		plot(NA,
		        type="l",ylab="Proportion\nsamples",xlab=NA,
		        xaxt="n",bty="n",xlim=XLIM,cex=CEX,
		        cex.axis=CEX,cex.lab=CEX,ylim=YLIM,
		        main=MAIN,yaxt="n")
		# chromsome labels
		if(is.null(LABY)) LABY=max(res[index,whichSigs],na.rm=TRUE)*0.95
		chroms = unique(windows[index,1])
		colchrom = c("white",rgb(0.3,0.3,0.3,0.3))
		colchromtext = c("white","black")
		sapply(chroms,FUN=function(x)
		{
		  chromindex = which(windows[index,1]==x)
		  X = range(chromindex)
		  X = c(X,rev(X))
		  polygon(x=X,y=c(-100,-100,100,100),col=colchrom[((which(chroms==x)%%2)==0)+1],border=NA)
		  text(x=mean(X),y=LABY,labels=x,
		       col=colchromtext[(!(which(chroms==x)%%2)==0)+1],
		       cex=CEX*0.5)
		})
		# add data to plot
		matlines(res[index,whichSigs],
        		lwd=2,cex=CEX,
        		col=colours,lty=LTY)
		if(doSplit) abline(h=0,lty=2,lwd=1,col="gray")
		# y-axis
		axPos = axTicks(side=2)
		axis(side=2,at=axPos,labels=abs(axPos),cex=CEX,cex.lab=CEX,cex.axis=CEX)
		par(mgp=c(3,1,0))
		# gene threshold
		if(doGenes) abline(h=c(propThreshold,-propThreshold),lty="dotted",col="gray")
		# plot legend
		if(doLegend)
		  {
		  if(is.null(legnames)) legnames = whichSigs
		  legend(legloc,lwd=2,lty=LTY,col=colours,legend=legnames,
			bty="n",cex=CEX,ncol=1,x.intersp=0.7,y.intersp=0.7)
		}
		# x-axis label
		par(mgp=c(2,1,0))
		if(!doSplit) title(xlab="Genomic window (1Mb)",cex.axis=CEX,cex.lab=CEX,cex=CEX)
		if(doSplit) res[,whichSigs[2]] = -res[,whichSigs[2]]
		# find genes that overlap peaks
		if(doGenes)
			{
			print("get genes")
		  # loop over each signature we want to find overlaps for
		  if(length(offsetBase)==1) offsetBase = rep(offsetBase,length(whichGenes))
			geneInfo = sapply(whichGenes,FUN=function(x)
				{
				print(x)
				if(all(is.na(res[index,x]))) return(NA)
				if(is.list(genesToTest))
					{
					thisGenes = genesToTest[[x]]
					thisGeneGR = geneGR[[x]]
					thisoffset = offsetBase[which(names(genesToTest)==x)]
					} else {
					thisGenes = genesToTest
					thisGeneGR = geneGR
					thisoffset = offsetBase
					}
			  print(paste0("get genes for ",x))
			 	getGenes(res=res[index,],counts=resC[index,],windows=windows,genesToTest=thisGenes,
					sig=x,propThreshold=propThreshold,countThreshold=countThreshold,
					col=rep(colours,10)[which(colnames(res[,whichSigs,drop=FALSE])==x)],
					ensembl=ensembl,geneGR=thisGeneGR,windowGR=windowGR,
					offsetBase=thisoffset,offsetTweak=offsetTweak,
					xoverlap=xoverlap,yoverlap=yoverlap,geneCex=geneCex,doGenesAbove=doGenesAbove
					)
				},simplify=FALSE)
			GENEINFO <<- geneInfo
			if(doArms)
			  {
			  # get chrom arms
			  print("get arms")
			  ARMIN <<- res[index,]
			armInfo = sapply(whichGenes,FUN=function(x)
			{
			  print(x)
			  if(all(is.na(res[index,x]))) return(NA)
			  print(paste0("get arms for ",x))
			  THISARMIN <<- res[index,x]
			  getArms(res[index,x],propThreshold,propArm=propArm,
			          arms=arms$arms,armsGR=arms$armsGR,overlaps=arms$armsOverlap)
			},simplify=FALSE)
			ARMINFO <<- armInfo
			# remove genes from chrom arm regions
			print("remove genes from arms")
			geneInfo=sapply(whichGenes,FUN=function(x)
			   {
			  print(paste0("remove genes from arms for ",x))
			  toremove = sapply(1:nrow(armInfo[[x]]),FUN=function(y)
			                    {
			                 poscheck = as.numeric(geneInfo[[x]][,1])
			                which(poscheck<=as.numeric(armInfo[[x]][y,3])&poscheck>=as.numeric(armInfo[[x]][y,2]))
			                  })
			  toremove = unique(unlist(toremove))
			  out = geneInfo[[x]]
			  if(length(toremove)>0) out = out[-toremove,,drop=FALSE]
			  out
			  },simplify=FALSE)
			GENEINFO2 <<- geneInfo
			}
			GENEINFOCHECK <<- geneInfo
			print("some checks")
			if(length(geneInfo)==0) next
			if(!doArms&all(sapply(geneInfo,FUN=function(x) all(is.na(x))))) next
			# plot genes
			print("plot genes")
			newMar = myMAR
			newMar[c(1,3)]=0
			par(mar=newMar)
			plot(NA,xlim=XLIM,ylim=c(-0.25,1),axes=FALSE,xlab=NA,ylab=NA)
			if(!is.null(extraMain)) title(main=extraMain,cex=CEX,cex.main=CEX,adj=0.02,line=-15)
			yvals = seq(from=0,to=1,length.out=length(geneInfo)+2)
			yvals = yvals[-c(1,length(yvals))]
			#print(geneInfo)
			whichgenestoplot = 1:length(geneInfo)
			if(doSplit) whichgenestoplot = 1
			for(x in whichgenestoplot)
			  {
			  if(length(geneInfo[[x]])==1) if(is.na(geneInfo[[x]])) next
			  if(!doArms&nrow(geneInfo[[x]])==0) next
			  if(doArms) if(nrow(geneInfo[[x]])+nrow(armInfo[[x]])==0) next
			  X = as.numeric(geneInfo[[x]][,"x"])
			  Y = rep(as.numeric(yvals[x]),length(geneInfo[[x]][,"x"]))
			  if(length(X)>0)
			  {
			  LAB = geneInfo[[x]][,"lab",drop=FALSE]
			  tmp = sapply(1:length(X),FUN=function(j) LAB[which(abs(X[j]-X)<xoverlap)])
			  names(tmp) = LAB
			  if(!all(sapply(tmp,FUN=length)==1))
			    {
			    for(j in names(tmp)) tmp[[j]]=tmp[[j]][-which(tmp[[j]]==j)]
			    # remove duplicates
			    toRemove = NULL
			    for(j in 1:length(tmp))
			      {
			      if(names(tmp)[j]%in%toRemove) next
			      toRemove = c(toRemove,tmp[[j]])
			      }
			    tmp = tmp[-which(names(tmp)%in%unique(toRemove))]
			    # remove empty
			    emptyindex = which(sapply(tmp,FUN=length)==0)
			    if(length(emptyindex)>0) tmp = tmp[-emptyindex]
			    print(tmp)
			    if(length(tmp)==0)
			      {
			      print("Skipping: 0 length"); next
			      }
			    # reorder by genomic position
			    for(ij in 1:length(tmp))
			      {
			      index = which(LAB%in%c(tmp[[ij]],names(tmp)[[ij]]))
			      index = index[order(X[index])]
			      X[index] = X[index]+(offsetTweak*(0:(length(index)-1)))
			      X[index] = X[index]+offsetBase[x]
			      }
			    }
			  }
			  # plot chrom arms
			  if(doArms)
			  {
			    ARMINFOCHECK <<- armInfo
			      thisarms = armInfo[[x]]
			      thiscolour = colours[which(sigsToPlot==names(armInfo)[x])]
			      THISARMS <<- thisarms
			      if(nrow(thisarms)>0)
			      {
			        apply(thisarms,MARGIN=1,FUN=function(i)
			        {
			          lines(x=c(i[2:3][c(1,1,2,2)]),
			                y=c(-0.25,-0.15,-0.15,-0.25),
			                col=thiscolour,lwd=2
			          )
			          lines(x=rep(mean(as.numeric(i[2:3])),2),y=c(-0.15,-0.05),col=thiscolour,lwd=2)
			          text(x=mean(as.numeric(i[2:3])),y=-0.05,labels=gsub(" ","",i[1]),cex=CEX,col=thiscolour,srt=90,adj=c(0,0.5))
			        })
			      }
			    }
			if(plotGenes&length(X)>0)
			{
			  # plot gene text labels
			  text(x=X,
			       y=rep(yvals[x],length(geneInfo[[x]][,"x"]))-geneYtweak,
			       col=geneInfo[[x]][,"col"],
			       labels=geneInfo[[x]][,"lab"],
			       cex=CEX*0.75,srt=90,font=3)
			  # plot lines to gene names
			  if(is.null(ylineoffset))
			    {
			    if(length(yvals)==2) ylineoffset = rev(c(10,5))
			    if(length(yvals)==1) ylineoffset = 5
			    if(length(yvals)>2) ylineoffset = seq(from=5,to=15,length.out=length(yvals))
			    }
			  mapply(FUN=function(x1,x2,y,col) {lines(x=c(x1,x2),y=c(-0.25,y-geneYtweak),col=col)},
			         x1=geneInfo[[x]][,"x"],
			         x2=X,
			         y=yvals[x]-(yvals[x]*nchar(geneInfo[[x]][,"lab"])/ylineoffset[x]),
			         col=geneInfo[[x]][,"col"])
			    }
			  }
			if(doSplit)
			  {
			  plot(NA,xlim=XLIM,ylim=rev(-c(-0.25,1)),xaxt="n",#yaxt="n",
			       axes=FALSE,bty="n",
			       xlab=NA,ylab=NA)
			  whichgenestoplot = 2
			for(x in whichgenestoplot)
			{
			  if(nrow(geneInfo[[x]])==0) next
			  X = as.numeric(geneInfo[[x]][,"x"])
			  Y = rep(as.numeric(yvals[x]),length(geneInfo[[x]][,"x"]))
			  LAB = geneInfo[[x]][,"lab"]
			  tmp = sapply(1:length(X),FUN=function(j) LAB[which(abs(X[j]-X)<xoverlap)])
			  names(tmp) = LAB
			  if(!all(sapply(tmp,FUN=length)==1))
			  {
			    for(j in names(tmp)) tmp[[j]]=tmp[[j]][-which(tmp[[j]]==j)]
			    # remove duplicates
			    toRemove = NULL
			    for(j in 1:length(tmp))
			    {
			      if(names(tmp)[j]%in%toRemove) next
			      toRemove = c(toRemove,tmp[[j]])
			    }
			    tmp = tmp[-which(names(tmp)%in%unique(toRemove))]
			    # remove empty
			    emptyindex = which(sapply(tmp,FUN=length)==0)
			    if(length(emptyindex)>0) tmp = tmp[-emptyindex]
			    print(tmp)
			    if(length(tmp)==0)
			    {
			      print("Skipping: 0 length"); next
			    }
			    # reorder by genomic position
			    for(ij in 1:length(tmp))
			    {
			      index = which(LAB%in%c(tmp[[ij]],names(tmp)[[ij]]))
			      index = index[order(X[index])]
			      X[index] = X[index]+(offsetTweak*(0:(length(index)-1)))
			      X[index] = X[index]+offsetBase[x]
			    }
			  }
			  if(plotGenes)
			  {
			  # plot gene text labels
			  text(x=X,
			       y=-rep(yvals[1],length(geneInfo[[x]][,"x"])),
			       col=geneInfo[[x]][,"col"],
			       labels=geneInfo[[x]][,"lab"],
			       cex=CEX*0.75,srt=90,font=3)
			  # plot lines to gene names
			  if(is.null(ylineoffset))
			  {
			    if(length(yvals)==2) ylineoffset = rev(c(10,5))
			    if(length(yvals)==1) ylineoffset = 5
			    if(length(yvals)>2) ylineoffset = seq(from=5,to=15,length.out=length(yvals))
			  }
			  mapply(FUN=function(x1,x2,y,col) {lines(x=c(x1,x2),y=c(0.25,y),col=col)},
			         x1=geneInfo[[x]][,"x"],
			         x2=X,
			         y=-yvals[x]-(yvals[x]*nchar(geneInfo[[x]][,"lab"])/ylineoffset[x]),
			         col=geneInfo[[x]][,"col"])
			  }
			}
			}
		}
		if(doSplit)
		  {
		  plot(NA,xlim=XLIM,ylim=0:1,xlab=NA,ylab=NA,bty="n",axes=FALSE)
		  text(x=mean(XLIM),y=0.5,labels="Genomic window (1Mb)",cex.axis=CEX,cex.lab=CEX,cex=CEX)
		  }
	# output formatting
	print("out formatting")	
	out = cbind(windows,res)
	colnames(out)[1] = "chrom"
	if(saveRes) write.csv(out,file=paste0(outDir,i,"-recurrence.csv"),row.names=FALSE,quote=FALSE)
	}
	if(!is.null(filename)) dev.off()
}


# get genes that overlap peaks
getGenes = function(res,counts,windows,genesToTest,sig="CN19",
                    propThreshold=0.2,countThreshold=5,colour="red",
                    ensembl=NULL,windowGR=NULL,geneGR=NULL,
                    offsetBase=NULL,offsetTweak=NULL,
                    xoverlap=100,yoverlap=0.02,geneCex=2,doGenesAbove=TRUE)
	{
	if(all(is.na(res[,sig]))) return(NA)
	# which windows have a proportion higher than threshold
  print("get above thresh")
  RES<<-res
  COUNTS<<-counts
  SIG<<-sig
  PROPTHRESH<<-propThreshold
  COUNTTHRESH<<-countThreshold
	index = which(res[,sig]>=propThreshold&counts[,sig]>=countThreshold)
	#which(RES[,SIG]>=PROPTHRESH&COUNTS[,SIG]>=COUNTTHRESH)
	if(length(index)==0) return(NA)
	toTest = windows[index,]
	if(!is.null(ensembl))
		{
		# get genes in those windows using biomart
		values = list(chrom=toTest[,1],start=as.numeric(toTest[,2]),end=as.numeric(toTest[,3]))
		attributes = c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position","end_position")
		filters = c("chromosome_name","start","end")
		genes = mapply(FUN=function(a,b,c) 
			{
  			getBM(attributes=attributes,
        			filters=filters,
        			values=list(a,b,c),
        			mart=ensembl)
			},a=values$chrom,b=values$start,c=values$end,SIMPLIFY=FALSE)
		genes = do.call(rbind,genes)
		# which window genes are in test list of genes
		genesMatch = genes[which(genes[,2]%in%genesToTest),]
		# plot those genes
		genesGR = as(paste0(genesMatch$chromosome_name,":",genesMatch$start_position,"-",genesMatch$end_position),"GRanges")
		overlaps=findOverlaps(genesGR,windowsGR)
		print(genesGR)
		return(sapply(1:length(genesGR),FUN=function(x)
			{
			#abline(v=overlaps@to[which(overlaps@from==x)],col="black",lty=2)
		  out = c(sig=genesMatch[x],x=overlaps@to[which(overlaps@from==x)],
		          y=res[overlaps@to[which(overlaps@from==x)],sig]+0.05,
		          col=colour,labels=genesMatch[x,"hgnc_symbol"])
  			text(x=out["x"],
       				y=out["y"],
       				col=out["colour"],labels=out["labels"],cex=2)
		  return(out)
			}))
		} else if(!is.null(geneGR)) {
		# get genes in those windows using GRanges
		print("get genes in windows")
		overlaps = findOverlaps(geneGR,windowGR)
		overlaps = overlaps[which(overlaps@to%in%index)]
		if(length(overlaps@to)==0) return(NA)
		print(overlaps)
		# plot
		#if(is.null(labelOffsets)) labelOffsets = rep(0.05,unique(overlaps@from))
		X = sapply(unique(overlaps@from),FUN=function(x) mean(overlaps@to[which(overlaps@from==x)]))
		Y = sapply(unique(overlaps@from),FUN=function(x) max(res[overlaps@to[which(overlaps@from==x)],sig]))
		LAB = sapply(unique(overlaps@from),FUN=function(x) genesToTest[x])
		if(doGenesAbove) return(cbind(x=X,y=Y,col=colour,lab=LAB))
		Xsave <<- X; Ysave <<- Y; LABsave <<- LAB
		tmp = sapply(1:length(X),FUN=function(x) LAB[which(abs(X[x]-X)<xoverlap&abs(Y[x]-Y)<yoverlap)],simplify=FALSE)
		names(tmp) = LAB
		print(tmp)
		for(i in names(tmp)) tmp[[i]]=tmp[[i]][-which(tmp[[i]]==i)]
		toRemove = NULL
		for(i in 1:length(tmp))
		  {
		  if(names(tmp)[i]%in%toRemove) next
		  toRemove = c(toRemove,tmp[[i]])
		  }
		tmp = tmp[-which(names(tmp)%in%unique(toRemove))]
		tmp = tmp[-which(sapply(tmp,FUN=length)==0)]
		print(tmp)
		for(i in 1:length(tmp))
		  {
		  index = which(LAB%in%tmp[[i]])
		  Y[index] = Y[index]+(offsetTweak*(1:length(index)))
		  }
		Y = Y + offsetBase
		if(!doGenesAbove) text(x=X,
		     y=Y,
		     col=colour,labels=LAB,cex=geneCex,font=3)
		return(cbind(x=X,y=Y,col=colour,lab=LAB))
		#sapply(unique(overlaps@from),FUN=function(x)
		#	{
		#	text(x=mean(overlaps@to[which(overlaps@from==x)]),
     #  				y=max(res[overlaps@to[which(overlaps@from==x)],sig])+
			#              labelOffsets[which(unique(overlaps@from)==x)],
       #				col=colour,labels=genesToTest[x],cex=2)
		  #print(c(genesToTest[x],
		  #        mean(overlaps@to[which(overlaps@from==x)]),
		  #        max(res[overlaps@to[which(overlaps@from==x)],sig])+0.05,
		  #        colour))
		#	})
		} else {
		return(NA)
		}


	}

# function to randomly draw samples and get recurrence
randomRecurr = function(n,seg)
	{
	samples = unique(seg[,1])
	tokeep = sample(samples,n)
	thisseg = seg[which(seg[,1]%in%tokeep),]
	tumType = rep("tumour",length=length(unique(thisseg$sample)))
	names(tumType) = unique(thisseg$sample)
	counts = getRecurrence(sample=thisseg$sample,chrom=thisseg$chr,
                       start=thisseg$startpos,end=thisseg$endpos,
                       CNS=thisseg$CNS,windowsGR=windows$windowsGR,
                       sigNames=sigNames,
                       tumType=tumType,ncores=4)
	do.call(rbind,counts)
	}


# function to convert monte carlo Ps to Qs
# see Sandve et al, 2011. Sequential Monte Carlo multiple testing
Qval = function(Ps) # p values for each test (each window)
        {
        Ps = Ps[which(!is.na(Ps))]
        n = length(Ps)
        # null proportion estimate
        nullProp = min(1,(2/n)*sum(Ps))
        # ordered p values
        index = 1:length(Ps)
        newOrder = order(Ps)
        ordP = Ps[newOrder]
        index = index[newOrder] 
        # q estimate
        qs = sapply(1:length(ordP),FUN=function(x)
                {
                min(n*nullProp*(ordP[x:n]/(x:n)))
                })
        # reorder qs
        qs = qs[order(index)]
        return(qs)
}



# dealing with plotting arms
setupArms = function(windowGR)
  {
  library(GenomicRanges)
  arms = read.table(paste0(dirs$homeDir,"Dropbox/PostDoc/Rlibs/steeleLib/steeleLib/data/cytoBand.txt"),
                  sep="\t",head=TRUE)
  arms = sapply(unique(arms[,1]),FUN=function(x) sapply(c("p","q"),FUN=function(y)
    {
    index = which(arms[,1]==x&grepl(y,arms[,4]))
    c(paste0(x,y),min(arms[index,2]),max(arms[index,3]))
    }),simplify=FALSE)
  arms = t(do.call(cbind,arms))
  arms[,1] = gsub("chr","",arms[,1])
  #arms[,1] = gsub("chr","",arms[,1])
  armsGR = as(paste0(gsub("p|q","",arms[,1]),":",arms[,2],"-",arms[,3]),"GRanges")
  overlaps = findOverlaps(armsGR,windowGR)
  return(list(arms=arms,armsGR=armsGR,armsOverlap=overlaps))
  }



# function to get arms to plot
getArms = function(thisdata,arms,overlaps,armsGR,threshold=0.25,propArm=0.75)
{
  out = sapply(1:length(armsGR),FUN=function(x)
    {
    armindex = overlaps@to[which(overlaps@from==x)]
    if(length(armindex)==0) return(c(chrom=arms[x,1],indexStart=NA,indexEnd=NA,score=NA,start=NA,end=NA))
    score=mean(thisdata[armindex]>threshold,na.rm=TRUE)
    return(c(chrom=arms[x,1],
             indexStart=range(armindex)[1],
             indexEnd=range(armindex)[2],
             score=score,
             start=arms[x,2],
             end=arms[x,3]))
    },simplify=FALSE)
  out = do.call(rbind,out)
  out = out[which(!is.na(as.numeric(out[,"score"]))),,drop=FALSE]
  out = out[which(out[,"score"]>propArm),,drop=FALSE]
  multiTab = table(sapply(out[,1],FUN=function(x) substr(x,1,nchar(x)-1)))
  if(any(multiTab>1))
    {
    toremove = c()
    toAdd = NULL
    index = which(multiTab>1)
    for(i in index)
      {
      thischrom = names(multiTab)[i]
      thisarms = paste0(thischrom,c("p","q"))
      thisremove = which(out[,1]%in%thisarms)
      toremove = c(toremove,thisremove)
      new = out[thisremove[1],,drop=FALSE]
      new[1,1] = gsub("p|q","",new[1,1])
      new[1,2] = min(as.numeric(out[thisremove,2]))
      new[1,3] = max(as.numeric(out[thisremove,3]))
      new[1,4] = mean(as.numeric(out[thisremove,4]))
      new[1,5] = min(as.numeric(out[thisremove,5]))
      new[1,6] = max(as.numeric(out[thisremove,6]))
      toAdd = rbind(toAdd,new)
      }
    out = out[-toremove,]
    out = rbind(out,toAdd)
    }
  out
}


