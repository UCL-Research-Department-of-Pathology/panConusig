#####################################################################
#       functions for random signature assignments
#####################################################################


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

