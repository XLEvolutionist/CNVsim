
#####
# A quick function to examine False positive and True positive rate
#####

# Simon Renny-Byfield

testRates<-function(calls,trues,intervals, geno=0, sample, sizes=1000, repeats=NULL) {
  #####
  # true positive rate
  #####

  # subset the true data by size
  trues<-subset(trues, subset=width(trues) %in% sizes)
  library(data.table)
  # load in a needed function
  source("/Users/simonrenny-byfield/GitHubRepos/cnvwin/scripts/data.frame2GRanges.R")
  # first extract the down CNVs
  down<-subset(calls, subset=sim_CNV %in% geno)
  down<-data.frame2GRanges(as.data.frame(down),keepColumns = TRUE)
  calls<-data.frame2GRanges(as.data.frame(calls),keepColumns = TRUE)
  
  # a subset of ranges that are not CNV
  notCNV<-setdiff(calls,trues)
  
  # see the overlaps among true CNV and calls
  # these are true positves
  match<-subsetByOverlaps(trues,down)
  # now the true postive matches
  TP<-length(match)
  P<-length(trues)
  TPR<-(TP/P)
  
  #####
  # false positive rate
  #####
  
  # see where the overlaps are between calls and "true" deletions
  vecFP<-overlapsAny(down,notCNV)
  # where "down" does not match an simulated deletion -> FALSE POSTIVE
  FP<-length(vecFP[vecFP==TRUE])
  # now we need to know the number of intervals that are NOT CNV
  # how many of our down CNVs overlap normal regions
  vecNeg<-overlapsAny(notCNV,down)
  #the number that do not i.e. FALSE are the true negatives
  TN<-length(vecNeg[vecNeg==FALSE])
  FPR<-(FP/(FP+TN))
  
  ####
  # False Discovery Rate and Precision
  ####
  
  FDR<- (FP/(TP+FP))
  PPV<- (TP/(TP+FP))
  
  return(list(truePos=TPR,falsePos=FPR,FDR=FDR,PPV=PPV))
}#