# function to artificially genome double a set of signatures
# sigs is a 48xn matrix of signature definitions, where n is the number of signatures
getGDsig = function(sigs) 
{
  factors = sapply(rownames(sigs),FUN=function(x) strsplit(x,split=":")[[1]])
  sigs = sigs[order(factor(factors[2,],levels=c("homdel","LOH","het")),
                    factor(factors[1,],levels=c("0","1","2","3-4","5-8","9+")),
                    factor(factors[3,],levels=c("0-100kb","100kb-1Mb","1Mb-10Mb","10Mb-40Mb",">40Mb",">1Mb"))),]
  sapply(1:ncol(sigs),FUN = function(n)
    {
    # create empty new signature
    tmpsig = matrix(0,ncol=1,nrow=nrow(sigs))
    rownames(tmpsig) = rownames(sigs)
    # keep 9+ as is
    tmpsig[grep("9+",rownames(tmpsig)),] = sigs[grep("9+",rownames(tmpsig)),n]
    # keep homdels as are
    tmpsig[grep("homdel",rownames(tmpsig)),] = sigs[grep("homdel",rownames(tmpsig)),n]
    # alter LOH classes
    index = rownames(sigs)[grep("LOH",rownames(sigs))]
    # only alter anything that is not 9+
    indexBefore = index[-grep("9+",index)]
    # only add from anything that is not bottom
    indexAfter = index[-grep("1[:]LOH",index)]
    # do alteration
    tmpsig[indexAfter,] = tmpsig[indexAfter,]+sigs[indexBefore,n]
    # alter het classes
    index = rownames(sigs)[grep("het",rownames(sigs))]
    # only alter anything that is not 9+
    indexBefore = index[-grep("9+",index)]
    indexAfter = index[-grep("2[:]",index)]
    tmpsig[indexAfter,] = tmpsig[indexAfter,]+sigs[indexBefore,n]  
    tmpsig
  })
}

