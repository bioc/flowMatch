## ========================================================================================
## An S4 class used to store statistical parameters for a cluster/cell populaiton 
## ========================================================================================
setClass("Cluster",
         representation(
           size = "integer", # The number of points in the cluster
           center = "numeric", # A vector of cluster centre
           cov = "matrix", # A matrix storing the covariance matrix of the cluster.
           cluster.id = "integer", # the index of the cluster (relative to other clusters in same sample)
           sample.id = "integer" # the index of sample in which the cluster belongs to
         ),
         validity=function(object){
           #cat("***** Cluster: inspector ***** \n")
           
           if(nrow(object@cov) != ncol(object@cov))
           {
             stop ("The covariance matrix is not square !")
           }
           if(nrow(object@cov) != length(object@center))
           {
             stop ("The dimension of covariance matrix is not same as the cluster center !")
           }
           return(TRUE)
         } 
)

## ========================================================================================
## Constructor for object of class 'Cluster', mainly used internally
## ========================================================================================
Cluster <- function(size, center, cov, cluster.id=NA_integer_, sample.id=NA_integer_)
{
  new("Cluster", size = size, center = center, cov = cov, cluster.id=cluster.id, sample.id=sample.id)
}


## ========================================================================================
## Summary method object of class 'Cluster'
## ========================================================================================
setMethod("summary", signature(object="Cluster"), function(object){
  cat("An Object of class 'Cluster' : \n")
  cat(paste("Cluster Id: ",object@cluster.id, "  sample.id: ", object@sample.id, "\n",sep="" ))
  cat(sprintf("Cluster size = %d\n", object@size))
  cat("\nCluster center: \n")
  print(object@center)
  cat("\nCluster covariance matrix: \n")
  print(object@cov)

})


## ========================================================================================
## accessor methods for Cluster object
## ========================================================================================

setMethod("get.size", signature(object="Cluster"), function(object) object@size)
setMethod("get.center", signature(object="Cluster"), function(object) object@center)
setMethod("get.cov", signature(object="Cluster"), function(object) object@cov)
setMethod("get.cluster.id", signature(object="Cluster"), function(object) object@cluster.id)
setMethod("get.sample.id", signature(object="Cluster"), function(object) object@sample.id)
setReplaceMethod("sample.id", signature(object="Cluster"),
                 function(object, value) {
                   object@sample.id <- value
                   object } )
