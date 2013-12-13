## ========================================================================================
## An S4 class used to store a template
## ========================================================================================

setClass("Template",
         representation(
           num.metaclusters = "integer", # Number of metaclusters 
           metaClusters = "list", # a list of objects of class 'MetaCluster'
           dimension = "integer", # dimensionality of the template (number of markers)
           size = "integer", # number of cells in the template (in all metaclusters)
           tree = "list", # the template tree stored as a hierarchical clustering object  "hclust" not working 
           template.id = "integer" # index of the template (relative to other templates)
         ),
         validity=function(object){
           #cat("***** Template: inspector ***** \n")
           
           if(length(object@metaClusters) != object@num.metaclusters)
           {
             stop ("num.metaclusters is not same as the length of list 'metaClusters' \n")
           }
           
           template.size = integer(1)
           for(i in 1:length(object@metaClusters))
           {
             template.size = template.size + get.size(object@metaClusters[[i]])
             if(length(get.center(object@metaClusters[[i]])) != object@dimension)
               stop("dimension of a metacluster does not match the dimension of the template")
           }
           if(template.size != object@size)
           {
             stop("Template size is not equal to the summation of all metacluster size!\n")
           }
           
           return(TRUE)
         } 
)


## ========================================================================================
## Constructor method for S4 class 'Template'
## is used internally from function create.template
## Argument: take cluster label, centers, covariance matrices 
##          If cluster centers or covariance matrices are not provided, sample must
##          be passed as argument so that the centers or covariance matrices can be 
##          computed from the data
## return value: An object of class 'ClusteredSample' 
## ========================================================================================

Template = function(meta.clusters, tree, template.id = NA_integer_)
{
  #cat("***** Template: constructor ***** \n")
  num.metaclusters = length(meta.clusters)
  if(num.metaclusters==0)
  {
    stop("Number of metaclusters must be at least one in a template \n")
  }
  template.size = integer(1)
  for(i in 1:num.metaclusters)
  {
    mc = meta.clusters[[i]]
    template.size = template.size + get.size(mc)
  }

  dimension = length(get.center(meta.clusters[[1]]))
  new('Template', num.metaclusters=num.metaclusters, metaClusters=meta.clusters,  dimension=dimension, size=template.size, tree=tree, template.id = template.id)
}

## ========================================================================================
## Plot method for S4 class 'Template'
## plot.mc: draw the meta-cluster, not the individual clusters
## color.mc: a vector of length num.mc , if we want to choose color of each mc .. all cluster for ith mc will be drawn by ith entry 
## ========================================================================================

setMethod("plot", signature(x="Template", y="missing"), function(x, y, alpha=.05, plot.mc=FALSE, color.mc=NULL, colorbysample=FALSE, ...){
  if(!is(x,"Template"))
  {
    stop(paste("Object ", as.character(x)," is not of class Template"))
  }
  
  template = x
  colors = integer(0)
  clusters = list()
  for(i in 1:length(get.metaClusters(template)))
  {
    mc = get.metaClusters(template)[[i]]
    if(plot.mc)
    {
      clusters = c(clusters, mc)
      if(!is.null(color.mc) && length(color.mc == length(get.metaClusters(template)) ))
        colors = c(colors, color.mc[i])
      else
       colors = c(colors, i)
    }
    else
    {
      clusters = c(clusters, get.clusters(mc))
      if(colorbysample==TRUE)
      {
        for(j in 1:length(get.clusters(mc)))
          colors = c(colors, get.sample.id(get.clusters(mc)[[j]]))
      }
      else if(!is.null(color.mc) && length(color.mc == length(get.metaClusters(template)) ))
        colors = c(colors, rep(color.mc[i], length(get.clusters(mc))))
      else
        colors = c(colors, rep(i, length(get.clusters(mc))))
    }
    
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
## summary method for S4 class 'Template'
## ========================================================================================

setMethod("show", signature(object="Template"), function(object){
  summary(object)
})



## ========================================================================================
## summary method for S4 class 'Template'
## ========================================================================================

setMethod("summary", signature(object="Template"), function(object,...){
  cat("An Object of class 'Template'\n")
  if(!is.na(get.template.id(object)))
  {
    cat(paste("Template Id: ", get.template.id(object), "\n",sep=""))
  }
  cat(sprintf("Number of metaclusters: %d\n", get.num.metaclusters(object)))
  for (i in 1:get.num.metaclusters(object))
  {
    mc.size = get.size(get.metaClusters(object)[[i]])
    cat(sprintf("Number of cells in metacluster %d: %d [ %2.2f %% ]\n", i, mc.size, 100*mc.size/get.size(object)))
  }
})



## ===========================================================================
## Function to plot and return the template tree
## ===========================================================================
setGeneric("template.tree", 
           function(object, ...) standardGeneric("template.tree"))

setMethod("template.tree", signature=signature(object="Template"), definition=function(object, ...){
  tree = get.tree(object)
  class(tree) = 'hclust';
  dendo = as.dendrogram(tree)
  plot(dendo, ...);
  return(tree)
})


## ========================================================================================
## accessor methods for Template object
## ========================================================================================

setMethod("get.size", signature(object="Template"), function(object) object@size)
setMethod("get.num.metaclusters", signature(object="Template"), function(object) object@num.metaclusters)
setMethod("get.dimension", signature(object="Template"), function(object) object@dimension)
setMethod("get.template.id", signature(object="Template"), function(object) object@template.id)
setMethod("get.metaClusters", signature(object="Template"), function(object) object@metaClusters)
setMethod("get.tree", signature(object="Template"), function(object) object@tree)
