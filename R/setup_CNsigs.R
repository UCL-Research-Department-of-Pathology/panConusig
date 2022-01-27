###################################################################
# function to generate CN summary matrix from CN profiles
###################################################################

# function to get NMF variables
copyNumberInfo = function(chrom,start,end,CN,CNb=NULL,doVariables=TRUE,returnSep=FALSE)
{ 
  # segment lengths
  print("lengths")
  lengths = (end-start)/1000000
  # allelic imbalance
  print("imbalance")
  CNa = CN-CNb
  # LOH
  LOH=pmin(CNa,CNb)
  # combine
  print("combine")
  if(!doVariables)
  { 
    combined = list(CN=CN,
                    lengths=lengths,
                    LOH=LOH)
    return(combined)
  }
  LOHstatus = ifelse(LOH==0,"LOH","het")
  LOHstatus[which(CN==0)] = "homdel"
  varCN = cut(CN,
              breaks=c(-0.5,0.5,1.5,2.5,4.5,8.5,Inf), # 0,1,2,3-4,5-8,9+
              labels=c("0","1","2","3-4","5-8","9+"))
  varLength = rep(NA,length=length(varCN))
  hdIndex = which(LOHstatus=="homdel")
  if(length(hdIndex)>0)
  { 
    varLength[hdIndex] = paste0(cut(lengths[hdIndex],breaks=c(-0.01,0.1,1,Inf)))
    varLength[-hdIndex] = paste0(cut(lengths[-hdIndex],breaks=c(-0.01,0.1,1,10,40,Inf))) 
  } else {
    varLength = paste0(cut(lengths,breaks=c(-0.01,0.1,1,10,40,Inf)))
  }
  renameVarLengths = c("(-0.01,0.1]"="0-100kb","(0.1,1]"="100kb-1Mb","(1,10]"="1Mb-10Mb","(10,40]"="10Mb-40Mb","(40,Inf]"=">40Mb","(1,Inf]"=">1Mb")
  varLength = renameVarLengths[paste0(varLength)]
  sepVars = paste(varCN,
                  LOHstatus,
                  varLength,
                  sep=":")
  if(returnSep) return(sepVars)
  variables = table(sepVars)
  return(variables)
}

getMatrix = function(data # CN profile segments, with columns sample, chr, startpos, endpos, nMajor, nMinor
                     )
  {
  # get variables
  # variable names may need to change here for different input files
  variables = sapply(unique(data$sample),FUN=function(x)
    {
    index = which(data$sample==x)
    copyNumberInfo(chrom=data$chr[index], # chromosome
               start=data$startpos[index], # start position
               end=data$endpos[index], # end position
               CN=data$nMajor[index]+data$nMinor[index], # total CN
               CNb=data$nMinor[index], # minor CN
               doVariables=TRUE)
    })
  names(variables)=unique(data$sample)
  # create all possible CN classes
  CNclasses = c("1","2","3-4","5-8","9+") # different total CN states
  lengths = c("0-100kb","100kb-1Mb","1Mb-10Mb","10Mb-40Mb",">40Mb") # different lengths
  lohs = c("LOH","het") # loh statuses
  allNames = sapply(lohs,FUN=function(x) sapply(CNclasses,FUN=function(y) sapply(lengths,FUN=function(z) paste0(y,":",x,":",z))))
  # remove het deletions (impossible)
  allNames = allNames[-grep("1[:]het",allNames)]
  # add hom dels
  homdelclasses = "0"
  homdelLengths = c("0-100kb","100kb-1Mb",">1Mb")
  allNames = c(sapply(homdelclasses,FUN=function(x) sapply(homdelLengths,FUN=function(z) paste0(x,":homdel:",z))),allNames)
  # create matrix
  myvars = matrix(0,ncol=length(variables),nrow=length(allNames),dimnames=list(allNames,names(variables)))
  # fill matrix
  for(i in names(variables)) myvars[names(variables[[i]]),i] = variables[[i]]
  myvars
  }

