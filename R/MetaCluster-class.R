## ========================================================================================
## An S4 class used to store statistical parameters for a metacluster
## ========================================================================================
setClass("MetaCluster",
         representation(
           num.clusters = "integer", # Number of clusters in the meta-cluster
           clusters = "list", # a list of objects of class 'Cluster' 
           size = "integer", # The number of points in the metacluster
           center = "numeric", # A vector of metacluster centre
           cov = "matrix" # A matrix storing the covariance matrix of the metacluster.
           #homogeneity.index = "numeric" # Homogeneity of the clusters (for example SCS or cohen's f)
         ),
         validity=function(object){
           #cat("***** MetaCluster: inspector ***** \n")
           if(object@num.clusters != length(object@clusters))
           {
             stop("num.clusters not equal to the number of clusters in the metacluster stored in list clusters\n")
           }
           if(nrow(object@cov) != ncol(object@cov))
           {
             stop ("The covariance matrix is not square !")
           }
           if(nrow(object@cov) != length(object@center))
           {
             stop ("The dimension of covariance matrix is not same as the metacluster center !")
           }
           return(TRUE)
         } 
)

## ========================================================================================
## Constructor for object of class 'MetaCluster', 
## mainly used internally when creating a template
## 
## ========================================================================================
MetaCluster <- function(clusters)
{
  if(!is(clusters, "list"))
  {
    stop("The argument clusters must be a list of object of class Cluster\n")
  }
  if(length(clusters)<1)
  {
    stop("The number of cluster must be at least one\n")
  }
    
  num.clusters = length(clusters)
  dim = length(get.center(clusters[[1]]))
  mc.center = rep(0,dim)
  mc.cov = matrix(0, nrow=dim, ncol=dim)
  mc.size = integer(1)
  
  for(i in 1:num.clusters)
  {
    cluster = clusters[[i]]
    if(!is(cluster, "Cluster"))
    {
      stop("The argument clusters must be a list of object of class Cluster\n")
    }
    mc.size = mc.size + get.size(cluster)
    mc.center = mc.center + get.size(cluster) * get.center(cluster);
    mc.cov = mc.cov + (get.size(cluster) - 1) * get.cov(cluster);
  }
  mc.center = mc.center / mc.size;
  mc.cov = mc.cov / (mc.size - num.clusters - 1);
  
  new("MetaCluster", num.clusters = num.clusters, clusters=clusters, size = mc.size, center = mc.center, cov = mc.cov)
}

## ========================================================================================
## Summary method object of class 'MetaCluster'
## ========================================================================================
setMethod("show", signature(object="MetaCluster"), function(object){
  summary(object)
})

## ========================================================================================
## Summary method object of class 'MetaCluster'
## ========================================================================================
setMethod("summary", signature(object="MetaCluster"), function(object){
  cat("An Object of class 'MetaCluster' : \n")
  cat(sprintf("Number of clusters in this MetaCluster: %d\n", get.num.clusters(object)))
  cat(sprintf("MetaCluster size = %d\n", get.size(object)))
  cat("\nMetaCluster center: \n")
  print(get.center(object))
  cat("\nMetaCluster covariance matrix: \n")
  print(get.cov(object))  
})


## ========================================================================================
## Plot method for S4 class 'Metacluster'
## plot.mc: draw the meta-cluster, not the individual clusters
## by de
## ========================================================================================

setMethod("plot", signature(x="MetaCluster", y="missing"), function(x, y, alpha=.05, plot.mc=FALSE, ...){
  if(!is(x,"MetaCluster"))
  {
    stop(paste("Object ", as.character(x)," is not of class MetaCluster"))
  }
  
  meta.cluster = x
  colors = integer(0)
  clusters = list()
  if(plot.mc)
  {
      clusters = list( meta.cluster)
      colors = 'red'
  }
  else
  {
      clusters = get.clusters(meta.cluster)
      colors = 1:length(get.clusters(meta.cluster))
  }
    

  
  args = match.call()
  if(plot.mc)
  {
    if( !("lwd" %in% names(args)) && !("lty" %in% names(args)) && !("col" %in% names(args))) 
      plot.cluster.contours(clusters, alpha=alpha, col=colors, lwd=3, lty=2, ...)
    else if( !("col" %in% names(args))) 
      plot.cluster.contours(clusters, alpha=alpha, col=colors, ...)
    else
      plot.cluster.contours(clusters,alpha=alpha, ...)
  }  
  else if(!("col" %in% names(args)))
    plot.cluster.contours(clusters, alpha=alpha, col=colors, ...)
  else
    plot.cluster.contours(clusters,alpha=alpha,  ...)
  
})


## ========================================================================================
## accessor methods for MetaCluster object
## ========================================================================================

setMethod("get.size", signature(object="MetaCluster"), function(object) object@size)
setMethod("get.num.clusters", signature(object="MetaCluster"), function(object) object@num.clusters)
setMethod("get.clusters", signature(object="MetaCluster"), function(object) object@clusters)
setMethod("get.center", signature(object="MetaCluster"), function(object) object@center)
setMethod("get.cov", signature(object="MetaCluster"), function(object) object@cov)

