## =========================================================================##
##                      Generic method definitions                          ##
## =========================================================================##


## ---------------------------------------------------------------------------
## Generic function for summary
## ---------------------------------------------------------------------------
setGeneric("summary", function(object,...) standardGeneric("summary"))



## ---------------------------------------------------------------------------
## Generic function for plotting of R objects
## ---------------------------------------------------------------------------
setGeneric("plot",function(x,y,...) standardGeneric("plot"))



## ---------------------------------------------------------------------------
## Generic function for getting size of different objects
## ---------------------------------------------------------------------------
setGeneric("get.size",function(object) standardGeneric("get.size"))
setGeneric("get.dimension",function(object) standardGeneric("get.dimension"))


## ---------------------------------------------------------------------------
## Generic function for accessors of Cluster & MetaCluster classes
## ---------------------------------------------------------------------------
setGeneric("get.center",function(object) standardGeneric("get.center"))
setGeneric("get.cov",function(object) standardGeneric("get.cov"))

## ---------------------------------------------------------------------------
## Generic function for accessors of Cluster & ClusteredSample classes
## ---------------------------------------------------------------------------
setGeneric("get.sample.id",function(object) standardGeneric("get.sample.id"))

## ---------------------------------------------------------------------------
## Generic function for accessors of Cluster class
## ---------------------------------------------------------------------------
setGeneric("get.cluster.id",function(object) standardGeneric("get.cluster.id"))
#setGeneric("set.sample.id",function(object) standardGeneric("set.sample.id"))
setGeneric("sample.id<-", function(object, value) standardGeneric("sample.id<-"))


## ---------------------------------------------------------------------------
## Generic function for accessors of MetaCluster & ClusteredSample classes
## ---------------------------------------------------------------------------
setGeneric("get.num.clusters",function(object) standardGeneric("get.num.clusters"))
setGeneric("get.clusters",function(object) standardGeneric("get.clusters"))

## ---------------------------------------------------------------------------
## Generic function for accessors of ClusteredSample class
## ---------------------------------------------------------------------------
setGeneric("get.labels",function(object) standardGeneric("get.labels"))



## ---------------------------------------------------------------------------
## Generic function for accessors of Template class
## ---------------------------------------------------------------------------
setGeneric("get.num.metaclusters",function(object) standardGeneric("get.num.metaclusters"))
setGeneric("get.template.id",function(object) standardGeneric("get.template.id"))
setGeneric("get.metaClusters",function(object) standardGeneric("get.metaClusters"))
setGeneric("get.tree",function(object) standardGeneric("get.tree"))


## ---------------------------------------------------------------------------
## Generic function for accessors of ClusterMatch class
## ---------------------------------------------------------------------------
setGeneric("get.match12",function(object) standardGeneric("get.match12"))
setGeneric("get.match21",function(object) standardGeneric("get.match21"))
setGeneric("get.matching.cost",function(object) standardGeneric("get.matching.cost"))
setGeneric("get.unmatch.penalty",function(object) standardGeneric("get.unmatch.penalty"))


