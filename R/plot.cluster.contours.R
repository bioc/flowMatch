# ==========================================================================================
# Draw an ellipse in two dimension for a bivariate normal distribution (for internal use only)
# ==========================================================================================

ellipse <- function(mu, sigma, alpha=.05, npoints=250, ...)         
{
  es <- eigen(sigma)
  e1 <- es$vec%*%diag(sqrt(es$val))
  r1 <- sqrt(qchisq(1-alpha,2))
  theta <- seq(0,2*pi,len=npoints)
  v1 <- cbind(r1*cos(theta),r1*sin(theta))
  pts=t(mu-(e1%*%t(v1)))
  lines(pts, ...)

}

# ==========================================================================================
# Calculate axis limits for a bivariate contour plot (for internal use only)
# input: the list of objects of class 'Cluster', alpha for (1-alpha)*100% contour
# The length of an axis is computed by sqrt of eigen value multiplied by the qhi-square 
# distribution
# this function returns limits (min and max) for each dimension
# ==========================================================================================

limitcalc = function(clust.list, alpha=.05)   
{
  r1 <- sqrt(qchisq(1-alpha,1)) #
  dim = length(get.center(clust.list[[1]]))
  #merge all means and all variance
  means = matrix(0, nrow=0, ncol=dim)
  vars = matrix(0, nrow=0, ncol=dim)
  
  for(i in 1:length(clust.list))
  {
    means = rbind(means, get.center(clust.list[[i]]))
    vars = rbind(vars, diag(get.cov(clust.list[[i]])))
  }
  
  width = r1 * sqrt(apply(vars, 2, max))
  min.pos = apply(means, 2, min) - width
  max.pos = apply(means, 2, max) + width
  diff = max.pos - min.pos
  max.pos = max.pos + diff/10 # a little extra space 
  min.pos = min.pos - diff/10
  minmax = rbind(min.pos, max.pos)
  return(minmax)
}

# ==========================================================================================
# this is an internal function used to plot cluster, meta-cluster or template
# draw Ellipse countour plot matrix for high dimensional clusters
# This will produce a matrix of plot similar to the functions splom or pairs or flowPlot
# Input: 
#   clusters : a list of object of clusss 'Cluster' defined in possibly more than one dimensions
#   axis.labels : a vector of length equal to the dimension of the clsuters
#                 ith entry of the vector denotes the label of ith diemnsion
#                 if not supplied we use FL1, FL2 etc.
#   alpha: (1-alpha)*100% quantile is used for plotting 
#   add: True/False denoting if the classed should be added to an exisitng plot
#   col: the color used to draw the ellipses, if not provided use foreground colors


plot.cluster.contours = function(clusters, axis.labels=NULL, alpha=.05, col=par("fg"), ...)
{
  if(!is(clusters, 'list'))
    stop("The input argument 'clusters' to function plot.cluster.contours must be a list of object 'Clusters'\n")
  len=length(clusters)
  if(len == 0)
  {
    stop("The list of of clusters passed for plotting is  empty\n")
  }
  dim.data = length(get.center(clusters[[1]]))
  if(dim.data < 2)
    stop("The data should be at least two dimensional\n")
  # repeat color if needed
  if(length(col) < length(clusters))
  {
    col = rep(col, length(clusters))
  }
  

  if(is.null(axis.labels))
  {
    if(!is.null(names(get.center(clusters[[1]]))))
    {
      axis.labels = names(get.center(clusters[[1]]))
    }
    else
    {
      axis.labels = paste('FL', 1:dim.data, sep='')
    }
  }

  if(length(axis.labels) != dim.data)
  {
    stop("The number of axis labels are not equal to the dimension of data\n")
  }
  limits=limitcalc(clusters, alpha=alpha)
  
  op <- par(no.readonly = TRUE)
  par(mfrow=c(dim.data-1,dim.data-1));
  par(mar=c(.5,.5,.5,.5)); 
  par(oma=c(4,5,3,3))
   
  for(y in 2: dim.data )
  {
    for(x in 1:(dim.data-1))
    {
      if(x>=y)
      {
        plot.new()
      }
      else 
      {
        if(y==dim.data && x!=1) 
        {
          plot(0,type='n',xlim=limits[,x],ylim=limits[,y], yaxt='n', xlab=axis.labels[x])
          mtext( axis.labels[x], side=1, line=2.5)
          #cex.lab=axis.lab.size, las=0
        }
        else if(x==1 && y!=dim.data) 
        {
          plot(0,type='n',xlim=limits[,x],ylim=limits[,y], xaxt='n', ylab=axis.labels[y])
          mtext(axis.labels[y], side=2, line=2.5)
          
        }
        else if(y==dim.data && x==1)
        {
          plot(0,type='n',xlim=limits[,x],ylim=limits[,y], xlab=axis.labels[x], ylab=axis.labels[y])
          mtext(axis.labels[y], side=2, line=2.5)
          mtext(axis.labels[x], side=1, line=2.5)
        }
        else
        {
          plot(0,type='n',xlim=limits[,x],ylim=limits[,y], xaxt='n', yaxt='n')
        }
        for(k in 1:len)
        {         
          cl=clusters[[k]]
          mu=get.center(cl)[c(x,y)]
          sigma=get.cov(cl)[c(x,y),c(x,y)]      
          ellipse(mu, sigma, alpha, type='l',col=col[k], ... )  
          
        }
      }
    }
  }
  par(op) # restore old parameter
}  

## same function with option of adding ellipses later
## does not work, it seems that the ellipses are shifted left in several subplots
# plot.cluster.contours = function(clusters, axis.labels=NULL, alpha=.05, add=FALSE, col=par("fg"), ...)
# {
#   if(!is(clustSamples, 'list'))
#     stop("The input argument 'clusters' to function plot.cluster.contours must be a list of object 'Clusters'\n")
#   len=length(clusters)
#   if(len == 0)
#   {
#     stop("The list of of clusters passed for plotting is  empty\n")
#   }
#   dim.data = length(clusters[[1]]@center)
#   if(dim.data < 2)
#     stop("The data should be at least two dimensional\n")
#   # repeat color if needed
#   if(length(col) < length(clusters))
#   {
#     col = rep(col, length(clusters))
#   }
#   
#   
#   if(add==TRUE)
#   {
#     dim.plot = par('mfrow')
#     if(dim.plot[1]!=(dim.data-1) && dim.plot[2]!=(dim.data-1))
#       stop("The clusters can not be added because the dimension of plot and the cluster do not match. \n")
#   }
#   else # otherwise create new axis systems
#   {
#     if(is.null(axis.labels))
#     {
#       if(!is.null(names(clusters[[1]]@center)))
#       {
#         axis.labels = names(clusters[[1]]@center)
#       }
#       else
#       {
#         axis.labels = paste('FL', 1:dim.data, sep='')
#       }
#     }
#     
#     if(length(axis.labels) != dim.data)
#     {
#       stop("The number of axis labels are not equal to the dimension of data\n")
#     }
#     limits=limitcalc(clusters, alpha=alpha)
#     
#     
#     par(mfrow=c(dim.data-1,dim.data-1));
#     par(mar=c(.5,.5,.5,.5)); 
#     par(oma=c(4,5,3,3))
#   }
#   
#   for(y in 2: dim.data )
#   {
#     for(x in 1:(dim.data-1))
#     {
#       if(x>=y)
#       {
#         if(add==FALSE)
#           plot.new()
#       }
#       else 
#       {
#         if(add==FALSE) # set up axis labels and others
#         {
#           if(y==dim.data && x!=1) 
#           {
#             plot(0,type='n',xlim=limits[,x],ylim=limits[,y], yaxt='n', xlab=axis.labels[x])
#             mtext( axis.labels[x], side=1, line=2.5)
#             #cex.lab=axis.lab.size, las=0
#           }
#           else if(x==1 && y!=dim.data) 
#           {
#             plot(0,type='n',xlim=limits[,x],ylim=limits[,y], xaxt='n', ylab=axis.labels[y])
#             mtext(axis.labels[y], side=2, line=2.5)
#             
#           }
#           else if(y==dim.data && x==1)
#           {
#             plot(0,type='n',xlim=limits[,x],ylim=limits[,y], xlab=axis.labels[x], ylab=axis.labels[y])
#             mtext(axis.labels[y], side=2, line=2.5)
#             mtext(axis.labels[x], side=1, line=2.5)
#           }
#           else
#           {
#             plot(0,type='n',xlim=limits[,x],ylim=limits[,y], xaxt='n', yaxt='n')
#           }
#         }
#         else
#           par(mfg=c(y-1,x))
#         for(k in 1:len)
#         {         
#           cl=clusters[[k]]
#           mu=cl@center[c(x,y)]
#           print(mu)
#           sigma=cl@cov[c(x,y),c(x,y)]      
#           ellipse(mu, sigma, alpha, type='l',col=col[k], ... )  
#           
#         }
#       }
#     }
#   }
# }  
# 
# 
