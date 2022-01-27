######################################################################
#   function to simulate new CN profiles from processes of:
#   CIN, chromothripsis, genome doubling, amplification
######################################################################

# ====================================================================
#                   MISC FUNCTIONS
# ====================================================================

# read in a file
readFile = function(file,head)
        {
        tmp = gsub("[.]gz","",file)
        ending = rev(strsplit(tmp,split="[.]")[[1]])[1]
        if(ending=="csv") return(read.csv(file,head=head,as.is=TRUE))
        if(ending%in%c("txt","tsv")) return(read.table(file,sep="\t",head=head,as.is=TRUE))
        return(read.table(file,sep="\t",head=head,as.is=TRUE))  
        }


# load cytoBand file
loadCytoBand = function(file=NULL)
        {
        if(is.null(file))
                {
                tmpEnv = new.env()
                data(list="cytoBand", package='panConusig',envir=tmpEnv)
                return(tmpEnv[["cytoBand"]])
                } else {
                return(readFile(file,head=TRUE))
                }
        }


# get chromosome lengths
getChromInfo = function(cytoFile=NULL)
        {
        data = loadCytoBand(cytoFile)
        data[,1] = gsub("chr","",data[,1])
        sapply(c(1:22,"X","Y"),FUN=function(x) 
                c(min(data[which(data[,1]==x),2]),
                max(data[which(data[,1]==x),3])))
        }

# function to setup chromosomes
setupChrom = function(start=0,end=10000,identity="1",allele="A")
  {
  data.frame(identity=identity,allele=allele,start=start,end=end,cumstart=0,cumend=end)
  }

# function to create chromosomes
createChroms = function(chromInfo,sex="male")
	{
	chroms = sapply(colnames(chromInfo),FUN=function(x) 
			{
			setupChrom(start=chromInfo[1,x],end=chromInfo[2,x],identity=x,allele="A")
			},simplify=FALSE)
    	chroms = c(chroms,
			sapply(colnames(chromInfo),FUN=function(x) 
			{
			setupChrom(start=chromInfo[1,x],end=chromInfo[2,x],identity=x,allele="B")
			},simplify=FALSE)
			)
	if(sex=="male")
		{
		chroms = chroms[-c(47:48)]
		} else {
		chroms = chroms[-which(names(chroms)=="Y")]
		}
	chroms
	}


# ====================================================================
#                   BACKEND FUNCTIONS
# ====================================================================
# function to do a breakpoint
doBreakpoint = function(chrom,breakpoint)
	{
	# break
  	index = which(chrom$cumstart<breakpoint&chrom$cumend>breakpoint)
	# chunk before breakpoint
	beforebreakpoint = chrom[1:index,]
	# if section of chromosome is not reversed
	if(beforebreakpoint$end[nrow(beforebreakpoint)]>beforebreakpoint$start[nrow(beforebreakpoint)])
		{
		newstart = beforebreakpoint$start[nrow(beforebreakpoint)]+breakpoint-beforebreakpoint$cumstart[nrow(beforebreakpoint)]
		beforebreakpoint$end[nrow(beforebreakpoint)]=newstart-1
		} else {
		# if it is reversed
		newstart = beforebreakpoint$start[nrow(beforebreakpoint)]-(breakpoint-beforebreakpoint$cumstart[nrow(beforebreakpoint)])
		beforebreakpoint$end[nrow(beforebreakpoint)]=newstart-1
		}
	# chunk after breakpoint
	afterbreakpoint = chrom[index:nrow(chrom),]
	# if section of chromosome isnot reversed
	if(afterbreakpoint$end[1]>afterbreakpoint$start[1])
		{
		newend = afterbreakpoint$end[1]-(afterbreakpoint$cumend[1]-breakpoint)
		afterbreakpoint$start[1]=newend+1
		} else {
		# if it is reversed
		newend = beforebreakpoint$start[nrow(beforebreakpoint)]-(breakpoint-beforebreakpoint$cumstart[nrow(beforebreakpoint)])
		afterbreakpoint$start[1]=newend-1
		}
	return(list(beforebreakpoint=beforebreakpoint,afterbreakpoint=afterbreakpoint))
	}

# function to redo cumulative lengths
redoCumLengths = function(chrom)
	{
  	# redo cumulative lengths
  	lengths = abs(chrom$end-chrom$start)
  	chrom$cumend=cumsum(lengths)
  	chrom$cumstart=c(0,cumsum(lengths)[-length(lengths)])
	chrom
	}


# ====================================================================
#                   ALTERATION FUNCTIONS
# ====================================================================

# function for simple gains (tandem duplications) and losses
# add distribution for break lengths?
focal = function(chrom,breakpoints=NULL,type="gain")
  {
  # get breakpoints
  if(is.null(breakpoints)) #breakpoints = sort(sample(min(chrom$cumstart):max(chrom$cumend),2))
	{
	# gains
	if(type=="gain")
		{
		# gains
		means = c(5.961351, 7.786183)
		sds = sqrt(c(0.4199448, 0.1068539))
		props = c(0.7360366, 0.2639634)
		length=Inf		
		while(length>chrom$cumend[nrow(chrom)])
			{
			# draw type
			whichMix = rbinom(1,1,props[2])
			# do draw
			length = 10^rnorm(1,mean=means[whichMix+1],sd=sds[whichMix+1])
			}
		} else {
		# losses
		means = c(6.188331, 7.588125)
		sds = sqrt(c(0.5686788, 0.1326166))
		props = c(0.6472512, 0.3527488)
		length=Inf		
		while(length>chrom$cumend[nrow(chrom)])
			{
			# draw type
			whichMix = rbinom(1,1,props[2])
			# do draw
			length = 10^rnorm(1,mean=means[whichMix+1],sd=sds[whichMix+1])
			}
		}
	start = runif(1,min=1,max=chrom$cumend[nrow(chrom)]-length)
	breakpoints=c(start,start+length)
	}
  # where in chromsome segments?
  index1 = which(chrom$cumstart<breakpoints[1]&chrom$cumend>breakpoints[1])
  index2 = which(chrom$cumstart<breakpoints[2]&chrom$cumend>breakpoints[2])
  # split section before first breakpoint
  beforebreakpoint = doBreakpoint(chrom,breakpoints[1])$beforebreakpoint
  # split section after second breakpoint
  afterbreakpoint = doBreakpoint(chrom,breakpoints[2])$afterbreakpoint
  # if gain add middle section again (tandem duplication)
  if(type=="gain")
    {
    addition = chrom[index1:index2,]
    changeFactor = ifelse(beforebreakpoint$end[nrow(beforebreakpoint)]<beforebreakpoint$start[nrow(beforebreakpoint)],-1,1)
    addition$start[1] = beforebreakpoint$end[nrow(beforebreakpoint)]+changeFactor
    addition$end[nrow(addition)] = afterbreakpoint$start[1]-changeFactor
    out = rbind(beforebreakpoint,addition,addition,afterbreakpoint)
    }
  # if loss remove middle section
  if(type=="loss")
    {
    out = rbind(beforebreakpoint,afterbreakpoint)
    }
  # redo cumulative lengths
  out = redoCumLengths(out)
  out
}


# function to do translocation
# need to add in possibility that chromosome becomes flipped
# e.g. start of chrom1 can be translocated to start of chrom2
tloc = function(chrom1,chrom2,breakpoints=NULL)
  {
  if(is.null(breakpoints)) breakpoints = c(sample(min(chrom1$cumstart):max(chrom1$cumend),1),
                                           sample(min(chrom2$cumstart):max(chrom2$cumend),1))
  # break on first chromosome
  split1 = doBreakpoint(chrom1,breakpoints[1])
  beforebreakpoint1 = split1$beforebreakpoint
  afterbreakpoint1 = split1$afterbreakpoint
  changeFactor = ifelse(beforebreakpoint1$end[nrow(beforebreakpoint1)]<beforebreakpoint1$start[nrow(beforebreakpoint1)],-1,1)
  beforebreakpoint1$end[nrow(beforebreakpoint1$end)] = beforebreakpoint1$end[nrow(beforebreakpoint1$end)]+changeFactor
  # break on second chromsome
  split2 = doBreakpoint(chrom2,breakpoints[2])
  beforebreakpoint2 = split2$beforebreakpoint
  afterbreakpoint2 = split2$afterbreakpoint
  changeFactor = ifelse(beforebreakpoint2$end[nrow(beforebreakpoint2)]<beforebreakpoint2$start[nrow(beforebreakpoint2)],-1,1)
  beforebreakpoint2$end[nrow(beforebreakpoint2$end)] = beforebreakpoint2$end[nrow(beforebreakpoint2$end)]+changeFactor
  # choose first join
  #choice = c(sample(2,1),sample(2,1))
  #if(choice=c(1,1))
  #    {
  #    tmp = beforebreakpoint2
  #    tmp = tmp[rev(1:nrow(tmp)),]
  #    tmp[,c(3,4)] = tmp[,c(4,3)]
  #    out1 = rbind(beforebreakpoint1,tmp)
  #    }
  #if(choice=c(2,2))
  #  {
  #  tmp = afterbreakpoint1
  #  tmp = tmp[rev(1:nrow(tmp)),]
  #  tmp[,c(3,4)] = tmp[,c(4,3)]
  #  out1 = rbind(tmp,afterbreakpoint2)
  #  }
  #if(choice=c(1,2))
  #  {
  #  out1
  #  }
  # first derivative chromosome
  out1 = rbind(beforebreakpoint1,afterbreakpoint2)
  out1 = redoCumLengths(out1)
  # second derivative chromosome
  out2 = rbind(beforebreakpoint2,afterbreakpoint1)
  out2 = redoCumLengths(out2)
  list(out1,out2)
  }

# function to do chromothripsis
chromothripsis = function(chrom,breakpoints=NULL,nbreaks=20)
  {
  if(is.null(breakpoints))
	{
	#breakpoints = sort(sample(min(chrom$cumstart):max(chrom$cumend),nbreaks))
	lengths=Inf
	flag = 1
	while(sum(lengths)>chrom$cumend[nrow(chrom)])
		{
		lengths = round(10^rnorm(nbreaks,mean=6,sd=0.7))
		flag = flag+1
		if(flag>20)
			{
			lengths = lengths[which(cumsum(lengths)<chrom$cumend[nrow(chrom)])]
			break
			}
		}
	start = runif(n=1,min=1,max=chrom$cumend[nrow(chrom)]-sum(lengths))
	breakpoints = start+cumsum(lengths)
	}
  chrombroken = list(chrom)
  # shatter the chromosome, do sequentially
  for(i in breakpoints)
	{
  n = length(chrombroken)
	split = doBreakpoint(chrombroken[[n]],i)
	beforebreakpoint = split$beforebreakpoint
	afterbreakpoint = split$afterbreakpoint
	beforebreakpoint$cumend[nrow(beforebreakpoint)]=i
	afterbreakpoint$cumstart[1]=i+1
	chrombroken = append(chrombroken,list(beforebreakpoint,afterbreakpoint))
	chrombroken = chrombroken[-n]
	}
  # lose some pieces
  whichLose = which(rbinom(length(chrombroken),1,prob=0.5)==1)
  if(length(whichLose)>0) chrombroken = chrombroken[-whichLose]
  # reorient some pieces
  whichReverse = which(rbinom(length(chrombroken),1,prob=0.5)==1)
  if(length(whichReverse)>0)
	{
  	for(i in whichReverse)
		{
		chrombroken[[i]] = chrombroken[[i]][rev(1:nrow(chrombroken[[i]])),]
		saveStart = chrombroken[[i]]$start
		chrombroken[[i]]$start = chrombroken[[i]]$end
		chrombroken[[i]]$end = saveStart
		}
	}
  # restitch in random order
  combineOrder = sample(length(chrombroken),length(chrombroken),replace=FALSE)
  chrombroken = chrombroken[combineOrder]
  out = do.call(rbind,chrombroken)
  # redo cumulative lengths
  out = redoCumLengths(out)
  out
  }


# function to do aneuploidy
aneuploidy = function(chroms,whichAlter=NULL,type="gain")
  {
  if(is.null(whichAlter)) whichAlter = sample(1:length(chroms),1)
  if(type=="gain") chroms = c(chroms,chroms[whichAlter])
  if(type=="loss") chroms = chroms[-whichAlter]
  chroms
  }

# function to do WGD
WGD = function(chrom)
  {
  chroms = c(chrom,chrom)
  chroms
  }

# function to do BFB
# should we only BFB the longer section (more likely to contain centromere?)
# should we check for centromere?
BFB = function(chrom,breakpoint=NULL)
	{
	# combine
	#neworder = sample(2,2)
	#chrom = rbind(chroms[[neworder[1]]],chroms[[neworder[2]]])
	# redo cum lengths
	#lengths = abs(chrom$end-chrom$start)
	#chrom$cumend=cumsum(lengths)
	#chrom$cumstart=c(0,cumsum(lengths)[-length(lengths)])
	if(is.null(breakpoint)) breakpoint = sample(min(chrom$cumstart):max(chrom$cumend),1)
	BP <<- breakpoint      
	# break
	split = doBreakpoint(chrom,breakpoint)
	beforebreakpoint = split$beforebreakpoint
	afterbreakpoint=split$afterbreakpoint
  # duplicate and reverse
	whichTake = which.max(c(sum(abs(beforebreakpoint$cumend-beforebreakpoint$cumstart)),
	                      sum(abs(afterbreakpoint$cumend-afterbreakpoint$cumstart))))
	#rbinom(1,1,0.5)==1
	# keep largest piece
	if(whichTake==1) added = beforebreakpoint
	if(whichTake!=1) added = afterbreakpoint
	added = added[rev(1:nrow(added)),]
	saveStart = added$start
	added$start=added$end
	added$end = saveStart
	if(whichTake) out = rbind(beforebreakpoint,added)
	if(!whichTake) out = rbind(added,afterbreakpoint)
	# redo cum lengths
	lengths = abs(out$end-out$start)
  	out$cumend=cumsum(lengths)
  	out$cumstart=c(0,cumsum(lengths)[-length(lengths)])
  	out
	}

# function to do chromoplexy
# need to add in possibility that chromosome ends become flipped 
# e.g. can have start of chrom1 and start of chrom2 joined
chromoplexy = function(chroms,breakpoints=NULL)
  {
  if(is.null(breakpoints)) breakpoints = sapply(chroms,FUN=function(x)
          {
          sample(min(x$cumstart):max(x$cumend),1)
          })
  # broken chromosomes
  brokenchroms = mapply(FUN=function(chrom,breakpoint)
        {
        split = doBreakpoint(chrom,breakpoint)
        beforebreakpoint = split$beforebreakpoint
        afterbreakpoint = split$afterbreakpoint
        changeFactor = ifelse(beforebreakpoint$end[nrow(beforebreakpoint)]<beforebreakpoint$start[nrow(beforebreakpoint)],-1,1)
        beforebreakpoint$end[nrow(beforebreakpoint$end)] = beforebreakpoint$end[nrow(beforebreakpoint$end)]+changeFactor
        list(beforebreakpoint,afterbreakpoint)
        },chrom=chroms,breakpoint=breakpoints,SIMPLIFY=FALSE)
  # which broken chroms join to which?
  combs = cbind(sample(length(chroms),length(chroms),replace=FALSE),
                sample(length(chroms),length(chroms),replace=FALSE))
  # do the joining
  sapply(1:nrow(combs),FUN = function(x)
      {
      outchrom = rbind(brokenchroms[[combs[x,1]]][[1]],brokenchroms[[combs[x,2]]][[2]])
      outchrom = redoCumLengths(outchrom)
      outchrom
      },simplify=FALSE)
  }


# ====================================================================
#                   WRAPPER FUNCTIONS
# ====================================================================  

# get CN from chroms
getCN = function(chroms)
  {
  library(GenomicRanges)
  combined = do.call(rbind,chroms)
  #newstart = pmin(combined$start,combined$end)
  #newend = pmin(combined$start,combined$end)
  #combined$start = newstart
  #combined$end = newend
  chromosomes = colnames(chromInfo)
  # breakpoints genomic ranges
  breakpoints = sapply(chromosomes,FUN=function(x)
      {
      index = which(combined$identity==x)
      bp = sort(unique(c(chromInfo[,x],combined$start[index],combined$end[index])))
      as(paste0(x,":",bp[-length(bp)]+1,"-",bp[-1]-1),"GRanges")
      })
  names(breakpoints) = chromosomes
  # A allele genomic ranges
  index = which(combined$allele=="A")
  A = as(paste0(combined$identity[index],":",
		pmin(combined$start[index],combined$end[index]),"-",
		pmax(combined$start[index],combined$end[index])),"GRanges")
  # B allele genomic ranges
  index = which(combined$allele=="B")
  B = as(paste0(combined$identity[index],":",
		pmin(combined$start[index],combined$end[index]),"-",
		pmax(combined$start[index],combined$end[index])),"GRanges")
  # A allele overlaps
  overlapsA = sapply(breakpoints,FUN=function(x) table(findOverlaps(A,x)@to),simplify=FALSE)
  names(overlapsA) = chromosomes
  # B allele overlap
  overlapsB = sapply(breakpoints,FUN=function(x) table(findOverlaps(B,x)@to),simplify=FALSE)
  names(overlapsB) = chromosomes
  # counts
  combined = sapply(chromosomes,FUN=function(x)
	{
	df = data.frame(chrom=x,
		start=breakpoints[[x]]@ranges@start,
		end=breakpoints[[x]]@ranges@start+breakpoints[[x]]@ranges@width,
		CNa=0,
		CNb=0)
	rownames(df) = 1:nrow(df)
	df[as.numeric(names(overlapsA[[x]])),"CNa"] = overlapsA[[x]]
	df[as.numeric(names(overlapsB[[x]])),"CNb"] = overlapsB[[x]]
	df
	},simplify=FALSE)
  do.call(rbind,combined)
  }

# function to do a single alteration
doAlteration = function(alteration,chroms,index=NULL)
  {
  # focal gain
  if(alteration=="focal gain")
  {
    if(is.null(index)) index = sample(length(chroms),1)
    chroms[[index]]=focal(chroms[[index]],type="gain")
  }
  # focal loss
  if(alteration=="focal loss")
  {
    if(is.null(index)) index = sample(length(chroms),1)
    chroms[[index]]=focal(chroms[[index]],type="loss")    
  }
  # chromothripsis
  if(alteration=="chromothripsis")
  {
    #nbreaks = rpois(1,50) # sampling assumption here
    nbreaks = round(10^(rnorm(1,mean=1.3,sd=0.3)))
    if(is.null(index)) index = sample(length(chroms),1)
    chroms[[index]]=chromothripsis(chroms[[index]],nbreaks=nbreaks)     
  }
  # chromoanasynthesis
  if(alteration=="chromoanasynthesis")
  {
    #nbreaks = rpois(1,50) # sampling assumption here
    nbreaks = round(10^(rnorm(1,mean=1.3,sd=0.3)))
    if(is.null(index)) index = sample(length(chroms),1)
    chroms=append(chroms,list(chromothripsis(chroms[[index]],nbreaks=nbreaks)))  
  }
  # translocation
  if(alteration=="translocation")
  {
    if(is.null(index)) index = sample(length(chroms),2,replace=FALSE)
    chroms[index]=tloc(chroms[[index[1]]],chroms[[index[2]]])    
  }
  # chromoplexy
  if(alteration=="chromoplexy")
  {
    nsamp = rpois(1,3)+1 # sampling assumptions here
    if(is.null(index)) index = sample(length(chroms),nsamp,replace=FALSE)
    chroms[index]=chromoplexy(chroms[index])     
  }
  # BFB
  if(alteration=="BFB")
  {
    if(is.null(index)) index = sample(length(chroms),1)
    chroms[[index]]=BFB(chroms[[index]])       
  }
  # aneuploidy gain
  if(alteration=="aneuploidy gain")
  {
    if(is.null(index)) index = sample(length(chroms),1)
    chroms = aneuploidy(chroms,type="gain",whichAlter=index)
  }
  # aneuploidy loss
  if(alteration=="aneuploidy loss")
  {
    if(is.null(index)) index = sample(length(chroms),1)
    chroms = aneuploidy(chroms,type="loss",whichAlter=index)
  }
  chroms
  }

# function to generate a CN profile from various processes
generateProfile = function(nAltBeforeGD=20,
                           doGD=FALSE,
                           nAltAfterGD=20,
                           probs=NULL,
                           chroms=NULL,
                           index=NULL,
                           doCN=TRUE
                           )
  {
  alterationTypes = c("focal gain","focal loss","chromothripsis","translocation","chromoplexy","BFB","aneuploidy gain","aneuploidy loss","chromoanasynthesis")
  if(is.null(probs)) probs = rep(1,length(alterationTypes))
  
  # setup chromosomes
  if(is.null(chroms))
    {
    chroms = createChroms(chromInfo)
    }
  
  # alteration before GD
  if(nAltBeforeGD>0)
  {
  for(i in 1:nAltBeforeGD)
  {
    # choose alteration
    alteration = sample(alterationTypes,1,prob=probs)
    print(alteration)
    chroms = doAlteration(alteration,chroms,index=index)
  }
  }
  
  # genome doubling
  if(doGD) chroms = WGD(chroms)
  
  # alterations after GD
  if(nAltAfterGD>0)
  {
  for(i in 1:nAltAfterGD)
    {
    # choose alteration
    alteration = sample(alterationTypes,1,prob=probs)
    print(alteration)
    chroms = doAlteration(alteration,chroms,index=index)
    }
  }
  # get CN profile from chroms
  if(doCN) return(getCN(chroms))
  chroms
  }


