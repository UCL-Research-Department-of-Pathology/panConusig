# ===========================================================================
#			SIMULATE RANDOM DATASETS TO TEST ENRICHMENT OF RECURRENCE OF SIGNATURES
# ===========================================================================
require(GenomicRanges)
randomRecurrence = function(ascat, # ascat profile from assignSigsToSegs.R
                            windowGR, # GenomicRanges for genomic windows
                            nrep,
                            nsim,
                            whichsig
                            )
  {
  # setup
  lengthsForSim = ascat[,c("CNvar","CNS")]
  lengthsForSim[,1] = sapply(lengthsForSim[,1],FUN=function(x)
                {
                strsplit(x,split=":")[[1]][[3]]
                })
  lengthsForSim[,1][which(lengthsForSim[,1]==">1Mb")] = "1Mb-10Mb"
  allchroms = unique(ascat[,2])
  allsamps = unique(ascat[,1])

  # perform simulations
  simSample = function(whichN)
    {
    chroms = sapply(c(1:22,"X"),FUN=function(chrom)
      {
      thisindex = sample(allsamps,1)
      chromindex = which(ascat[,1]==thisindex&ascat[,2]==chrom)
      thischrom = ascat[chromindex,2:4,drop=FALSE]
      thislengths = lengthsForSim[chromindex,,drop=FALSE]
      newCNS = sapply(thislengths$CNvar,FUN=function(y)
        {
        sample(lengthsForSim[which(lengthsForSim[,1]==y),2],1)
        })
      newCNS = cbind(thischrom,CNS=newCNS)
      newCNS
      },simplify=FALSE)
    out = do.call(rbind,chroms)
    out = cbind(sample=whichN,out)
    out
    }

  simDataset = function(n=1000,whichsig=8)
    {
    newsamps = sapply(1:n,FUN=simSample,simplify=FALSE)
    newsamps = do.call(rbind,newsamps)
    newsampsGR = as(paste0(newsamps$chr,":",newsamps$startpos,"-",newsamps$endpos),"GRanges")
    overlaps = findOverlaps(newsampsGR,windowGR)
    simres = sapply(1:length(windowGR),FUN=function(x)
      {
      #print(x)
      thisseg=newsamps[overlaps@from[which(overlaps@to==x)],]
      if(nrow(thisseg)==0) return(rep(NA,19))
      subres = sapply(unique(thisseg$sample),FUN=function(y)
        {
        index = which(thisseg$sample==y)
        paste0("CN",1:22)%in%thisseg$CNS[index]
        })
    rowMeans(subres)
    })
  simres = t(simres)
  simres[,whichsig] # only return CN8
  }


  require(pbmcapply)
  multiSim = function(nrep=1000,nsim=1000)
    {
    allres = pbmclapply(1:nrep,FUN=function(x) simDataset(nsim),mc.cores=4)
    out = do.call(cbind,allres)
    colnames(out) = paste0("rep",1:nrep)
    out
    }
  
  multiSim(nrep,nsim)
  }


