#' @title Dendrogram Dataframe Creator
#'
#' @description Transforms the output of a hierarchical algorithm \code{hc} in a dataframe
#'   suitable for ggplot2 visualization
#'
#' @param hc a hierarchy resulting from DIANA
#'
#' @return A hierarchy dataframe suitable for ggplot2 visualization
#'
#' @export
#'
#' @examples
#' mat <- replicate(5, rnorm(5))
#' hc  <- funBI_diana(mat)
#' dendro_df(hc)

dendro_df <- function(hc) {

  position <- hc$position
  copy <- hc$merge
  copyvec <- as.vector(t(copy))
  # vectorize merge order
  or_name <- as.data.frame(copyvec)
  colnames(or_name) <- 'or_name'

  # label ordering and position
  ordering <- cbind(-position,1:length(position))
  colnames(ordering) <- c('or_name','position')

  #create dict to know or_name, new_name and position in the tree
  dict <- merge(or_name,ordering, by='or_name', all.x=TRUE)
  dict <- dict[order(dict$position),]

  dict <- dict[,c('or_name','position')]


  # create merge matrix with position
  posvec <- copyvec
  for(i in 1:length(copyvec[which(copyvec<0)])){
    posvec[which(posvec<0)][i] <- -dict[which(dict$or_name==copyvec[which(copyvec<0)][i]),'position']
  }


  poscopy <- matrix(posvec, ncol=dim(hc$merge)[2], nrow=dim(hc$merge)[1], byrow=TRUE)
  poscopy

  # X-START
  x_start <- dict$position

  positive <- as.data.frame(which(poscopy>0, arr.ind = TRUE))
  positive <- positive[order(positive$row),]
  for( i in 1:dim(positive)[1]){
    poscopy[positive[i,1],positive[i,2]] <- sum((abs(poscopy[poscopy[positive[i,1],positive[i,2]],])))/2
  }

  x_start <- (as.vector(t(poscopy)))

  # vector of x_end
  x_end <- x_start

  # vector of y_start
  y_start <- x_start
  y_start[x_start<0] <- 0
  y_start[x_start>0] <- hc$height[copyvec[which(x_start>0)]]
  y_start

  # vector of y_end
  y_end <- rep(hc$height,each=2)
  y_end

  # vertical segments dataframe
  segments_v <- as.data.frame(cbind(abs(x_start),abs(x_end),y_start,y_end))
  colnames(segments_v) <- c('x_start','x_end','y_start','y_end')

  # horizontal segments
  temp <- segments_v

  x_start <- temp$x_start[seq(1,dim(temp)[1],by=2)]
  x_end <- temp$x_end[seq(2,dim(temp)[1],by=2)]
  y_start <- temp$y_end[seq(2,dim(temp)[1],by=2)]
  y_end <- y_start

  segments_h <- as.data.frame(cbind(abs(x_start),abs(x_end),y_start,y_end))
  colnames(segments_h) <- c('x_start','x_end','y_start','y_end')

  segments_data <- as.data.frame(rbind(segments_h,segments_v))
  return(segments_data)
}


#' @title ggplot2 Hscore Dendrogram Plot
#'
#' @description Creates from the output of a hierarchical algorithm \code{hc} a
#'   dataframe suitable for a ggplot2 visualization using dendro_df. Then it
#'   plots a ggplot2 dendrogram having hscore as height
#'
#' @param hc a hierarchy resulting from DIANA
#'
#' @return ggplot object
#'
#' @export
#'
#' @examples
#'

hscoredendroplot <- function(hc){

  hc_df <- dendro_df(hc)

  p <- ggplot2::ggplot(hc_df) +
    ggplot2::geom_segment(ggplot2::aes(x = x_start, y = y_start, xend = x_end, yend = y_end, colour = "segment"), data = hc_df, col='black') +
    ggplot2::labs(x = "", y="H score") +
    # TO DELETE X AXIS NAMES
    ggplot2::annotate(geom = "text", x = 1:length(hc$position),  y = 0 - min(hc_df$y_end), label = hc$position, size = 2.5, color="#9BC2C7", angle=90) +
    ggplot2::theme(plot.margin = ggplot2::unit(c(1, 1, 1, 1), "lines"),
          legend.position="none",
          # box and guidelines
          #panel.border = element_rect(fill=NA, color="#9BC2C7", size=1.5 ),
          panel.background = ggplot2::element_rect(fill = "white",
                                          colour = "#CDCDCD",
                                          size = 0, linetype = "solid"),
          panel.grid.major = ggplot2::element_blank(),
          panel.grid.minor = ggplot2::element_blank(),
          # no x axis elemet because I fill it by myself
          axis.title.x = ggplot2::element_blank(),
          axis.text.x = ggplot2::element_blank(),
          # axis y
          axis.text.y = ggplot2::element_text(face= 'bold', hjust = 1, size= 9, color='#9BC2C7' ),
          axis.title.y = ggplot2::element_text(face= 'bold', hjust = 1, size= 0.05, color='#9BC2C7' )
    )
  p
}


#' @title Biclusters Map
#'
#' @description Creates a two-dimensional plot where every dot is a bicluster
#'   whose length of the sub-interval is represented by the x-axis, numbers of
#'   functions by the y-axis and Hscore by the dimension of the dot itself, i.e.
#'   small dot refers to low Hscore while big dot is for results having greater
#'   scores.
#'
#' @param mat a matrix
#' @param bic_list a list of biclusters from \code{hc}
#'
#' @return
#'
#' @export
#'
#' @examples
#'
#'

bic_map <- function(mat, bic_list){
  # processing data starting from list of biclusters
  card <- unlist(lapply(bic_list, function(x) length(x[[1]]) ))
  hsc <- unlist(lapply(bic_list, function(x) x[[3]] ))
  L <- unlist(lapply(bic_list, function(x) x[[2]][2] - x[[2]][1] +1 ))

  # data frame created
  mybic <- as.data.frame(cbind(L,card,hsc))

  # x is L, y is cardinality and size is hscore
  p <- ggplot2::ggplot(mybic, ggplot2::aes(x=L, y=card,size = hsc)) +
       ggplot2::geom_point(alpha=0.3,fill="#37AEBA", colour="#37AEBA",pch=21) + ggplot2::ylim(0,dim(mat)[1])
  return(p)
}





