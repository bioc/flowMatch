
## ========================================================================================
## An S4 class used to store matching of a pair of Templates or  clusteredSamples 
## An object of this class is only created internally and returned 
## Note1: I am not storing cluster distances here because they are easily 
##      obtaied from the diatnce matrix. 
## Note2: I am not storing matching as a matrix since it is really sparse
## ========================================================================================
setClass("ClusterMatch",
         representation(
           match12 = "list", # Matching of clusters (meta-clusters) from sample1 (template1) to sample2 (template2)
           match21 = "list", # Matching of clusters (meta-clusters) from sample2 (template2) to sample1 (template1)
           matching.cost = "numeric", # summation of distances for each pair of matched clusters (meta-clusters)
           unmatch.penalty = "numeric" # the penalty charged for each unmatched cluster (meta-cluster)
         ),
         validity=function(object){
           #cat("***** ClusterMatch: inspector ***** \n")
           # Check the symmetry of the matching in list match12 and match 21
           for(i in 1:length(object@match12))
           {
             mates = object@match12[[i]]
             for(mate in mates)
             {
                 if(!i %in%  object@match21[[mate]])
                 {
                   stop(paste("Matching symmetry error: matching (", i,mate, ") is present in match12 but not in match21!", sep='' ))
                 }
             }
           }
           
           for(i in 1:length(object@match21))
           {
             mates = object@match21[[i]]
             for(mate in mates)
             {
               if(!i %in%  object@match12[[mate]])
               {
                 stop(paste("Matching symmetry error: matching (", mate, i, ") is present in match21 but not in match12!", sep='' ))
               }
             }
           }
           
           return(TRUE)
         } 
)

## ========================================================================================
## Constructor for object of class 'ClusterMatch', mainly used internally
## ========================================================================================
ClusterMatch <- function(match12, match21, matching.cost, unmatch.penalty)
{
  new("ClusterMatch", match12 = match12, match21 = match21, matching.cost = matching.cost, unmatch.penalty=unmatch.penalty)
}

## ========================================================================================
## For internal use only  
## Generate a matching matrix from an object of class 'ClusterMatch'
## Return a matrix mat; 
## mat[i,j] = 1 if ith cluster form sample 1 is matched to jth cluster from sample 2
## ========================================================================================

matching.matrix = function(cluster.match)
{
  num.cluster1 = length(get.match12(cluster.match))
  num.cluster2 = length(get.match21(cluster.match))
  mat = matrix(0, nrow=num.cluster1, ncol=num.cluster2)
  
  match12 = get.match12(cluster.match)
  for(i in 1:length(match12))
  {
    mates = match12[[i]]
    for(mate in mates)
    {
      mat[i,mate] = 1
    }
  }
  rownames(mat) = paste('[1,', 1:num.cluster1, ']', sep='')
  colnames(mat) = paste('[2,', 1:num.cluster2, ']', sep='')
  return (mat)
}


## ========================================================================================
## Summary method object of class 'ClusterMatch'
## ========================================================================================
setMethod("summary", signature(object="ClusterMatch"), function(object){  

    match12 = get.match12(object)
    cat("==========================================================\n")
    cat("clusters/meta-clusters      matched clusters/meta-clusters\n")
    cat("from sample1/template1      sample2/template2 \n")
    cat("==========================================================\n")
    for(i in 1:length(match12))
    {
      cat(paste('           ', i, '                            ', sep=''))
      mates = match12[[i]]
      if(length(mates)!=0)
        cat(mates)
      else
        cat('x')
      cat('\n')
    }
    
    for(i in 1:length(get.match21(object)))
    {
      mates = get.match21(object)[[i]]
      if(length(mates)==0)
        cat(paste('           x                            \n',i ,sep=""))
    }
    cat("==========================================================\n")
    
#     cat("\n=====================Alternative View:====================\n")
#     cat("Matching is shown as a matrix : \n")
#     cat("[i,j] in row/column name denote jth cluster/meta-cluster\n")
#     cat("from the ith sample/template\n")
#     cat("a matched pair of clusters/meta-clusters is denoted with 1\n")
#     cat("==========================================================\n")
#     mat = matching.matrix(object)
#    print(mat)
})



## ========================================================================================
## This function matches a pair of samples or templates using the mixed edge cover
## using already computed dissimilarity matrix
## ========================================================================================

match.clusters.dist <- function(d.matrix, unmatch.penalty=999999)
{
  ## call a c++ function in the module
  mec = computeMEC(d.matrix, unmatch.penalty)
  match = ClusterMatch(match12=mec$match12, match21=mec$match21, matching.cost=mec$matching.cost, unmatch.penalty=unmatch.penalty)
  return(match); 
}



## ========================================================================================
## This function matches a pair of samples or templates using the mixed edge cover
## ========================================================================================

match.clusters <- function(object1, object2, dist.type='Mahalanobis', unmatch.penalty=999999)
{
  if(!is(object1,'ClusteredSample') && !is(object1,'Template'))
  {
    stop('object1 must be an object of class ClusteredSample or Template')
  }
  if(!is(object2,'ClusteredSample') && !is(object2,'Template'))
  {
    stop('object2 must be an object of class ClusteredSample or Template')
  }
  
  dM = dist.matrix(object1, object2, dist.type=dist.type)
  ## call a c++ function in the module
  mec = computeMEC(dM, unmatch.penalty)
  match = ClusterMatch(match12=mec$match12, match21=mec$match21, matching.cost=mec$matching.cost, unmatch.penalty=unmatch.penalty)
  return(match); 
}



## ========================================================================================
## accessor methods for ClusterMatch object
## ========================================================================================

setMethod("get.match12", signature(object="ClusterMatch"), function(object) object@match12)
setMethod("get.match21", signature(object="ClusterMatch"), function(object) object@match21)
setMethod("get.matching.cost", signature(object="ClusterMatch"), function(object) object@matching.cost)
setMethod("get.unmatch.penalty", signature(object="ClusterMatch"), function(object) object@unmatch.penalty)


