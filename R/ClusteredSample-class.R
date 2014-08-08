## ========================================================================================
## An S4 class used to store a clustered sample
## ========================================================================================

setClass("ClusteredSample",
         representation(
           num.clusters = "integer", # Number of clusters 
           labels = "integer", # A vector of integers (from 1:num.clusters) indicating the cluster to which each point is allocated
           clusters = "list", # a list of objects of class 'Cluster'
           dimension = "integer", # dimensionality of the sample (number of markers)
           size = "integer", # number of cells in the sample (in all clusters)
           sample.id = "integer" # index of the sample (relative to other samples of a cohort)
           ),
         validity=function(object){
           #cat("***** ClusteredSample: inspector ***** \n")
           if(max(object@labels) != object@num.clusters)
           {
             stop ("num.clusters is not same as the number of distinct cluster labels \n")
           }
           if(length(object@clusters) != object@num.clusters)
           {
             stop ("num.clusters is not same as the length of list 'clusters' \n")
           }
           
           sample.size = 0
           for(i in 1:length(object@clusters))
           {
             sample.size = sample.size + get.size(object@clusters[[i]])
             if(length(get.center(object@clusters[[i]])) != object@dimension)
              stop("dimension of a cluster does not match the dimension of the sample")
           }
           if(sample.size != object@size)
           {
             stop("Sample size is not equal to the summation of all cluster size!\n")
           }
           
           return(TRUE)
         } 
)


## ========================================================================================
## Constructor method for S4 class 'ClusteredSample'
## Argument: take cluster label, centers, covariance matrices 
##          If cluster centers or covariance matrices are not provided, sample must
##          be passed as argument so that the centers or covariance matrices can be 
##          computed from the data
## return value: An object of class 'ClusteredSample' 
## ========================================================================================

ClusteredSample = function(labels, centers=list(), covs=list(), sample=NULL, sample.id = NA_integer_)
{
  #cat("***** ClusteredSample: constructor ***** \n")
    num.clusters = max(labels)
    sample.size = length(labels)
    clusters=list()
    
    if(length(centers) != length(covs) )
    {
        stop("Number of centers and number of covariance matrices are not same!\n")
    }
    else if(length(centers)==0 || length(covs) ==0)
    {
      if(is.null(sample))
      {
          stop("Please provide the sample in the argument when centers of the cluster are not provided so that the centers can be calculated from data\n")
      }
      if(is(sample,"flowFrame"))
      {
        sample = exprs(sample)
      }
      else if(!is(sample,"matrix") && !is(sample,"data.frame"))
      {
        stop(paste("\"sample\" is not an object of class: flowFrame, matrix or data.frame\n"))
      }
      sample = as.matrix(sample)
    }
    
    
    for(i in 1:num.clusters)
    {
      cluster.idx = which(labels==i)
      
      if(length(centers) == 0)
      {
          selected.rows = matrix(sample[cluster.idx,], nrow=length(cluster.idx))
          cluster.center = colMeans(selected.rows)
      }
      else
        cluster.center = centers[[i]]
      if(length(covs) == 0)
      {
          selected.rows = matrix(sample[cluster.idx,], nrow=length(cluster.idx))
          cluster.cov = cov(selected.rows)
      }
      else
        cluster.cov = covs[[i]]
    
      cluster.cov[is.na(cluster.cov)] = 0
      cluster.size = length(cluster.idx)
      clusters[[i]] = Cluster(size=cluster.size, center=cluster.center, cov=cluster.cov, cluster.id = i, sample.id = sample.id)
    }
  
    dimension = length(get.center(clusters[[1]]))
  new('ClusteredSample', num.clusters=num.clusters, labels=labels,  clusters=clusters, dimension=dimension, size=sample.size, sample.id = sample.id)
}

## ========================================================================================
## Plot method for S4 class 'ClusteredSample'
## ========================================================================================

setMethod("plot", signature(x="flowFrame", y="ClusteredSample"), function(x, y, var.names=NULL, ...){
  plot(exprs(x), y, var.names, ...);
})
  
setMethod("plot", signature(x="ANY", y="ClusteredSample"), function(x, y, var.names=NULL, ...){
  if(is(x,"matrix") || is(x,"data.frame"))
  {
    sample = data.frame(x)
  }
  else
  {
    stop(paste("Object ", as.character(x)," is not of class flowFrame, matrix or data.frame"))
  }

  if(length(var.names)!=0)
  {
    sample = sample[,var.names]
  }
  else if(ncol(sample) == 1) # for proper plot of 1 dimensional data 
  {
    sample = as.numeric(sample[,1])
  }
  
  ## for more than two dimension it calls pairs function
  ## see plot.data.frame for more detail
  plot(sample, col=get.labels(y), pch='.', ...) 
})

#setMethod("show", signature(object="ClusteredSample"), function(object){
#  cat("An Object of class 'ClusteredSample' with the following slots: \n")
#  cat("num.clusters, labels, sizes, centers, cov, df\n")
#})


## ========================================================================================
## summary method for S4 class 'ClusteredSample'
## ========================================================================================

setMethod("summary", signature(object="ClusteredSample"), function(object,...){
  cat("An Object of class 'ClusteredSample'\n")
  if(!is.na(object@sample.id))
  {
    cat(paste("Sample Id: ", object@sample.id, "\n",sep=""))
  }
  cat(sprintf("Number of clusters: %d\n", object@num.clusters))
  for (i in 1:object@num.clusters)
  {
    cluster.size = get.size(object@clusters[[i]])
    cat(sprintf("Number of cells in cluster %d: %d [ %2.1f %% ]\n", i, cluster.size, 100*cluster.size/object@size))
  }
})


## ========================================================================================
## accessor methods for ClusteredSample object
## ========================================================================================

setMethod("get.size", signature(object="ClusteredSample"), function(object) object@size)
setMethod("get.num.clusters", signature(object="ClusteredSample"), function(object) object@num.clusters)
setMethod("get.clusters", signature(object="ClusteredSample"), function(object) object@clusters)
setMethod("get.dimension", signature(object="ClusteredSample"), function(object) object@dimension)
setMethod("get.sample.id", signature(object="ClusteredSample"), function(object) object@sample.id)
setMethod("get.labels", signature(object="ClusteredSample"), function(object) object@labels)


