#' @title Hscore
#'
#' @description Calculates Hscore with both alpha and beta different to 0
#'
#' @param mat a matrix
#'
#' @return The hscore of \code{mat}
#'
#' @export
#'
#' @examples
#' mat <- replicate(5, rnorm(5))
#' ccscore(mat)

ccscore<-function(mat)
{
  if(nrow(mat)>1){
    score <- sum((mat - rowMeans(mat) - matrix(colMeans(mat),nrow=nrow(mat),ncol=ncol(mat),byrow=TRUE) + mean(mat)) ^2)/(nrow(mat)*ncol(mat))
  }else{
    score <- 0
  }
  return(score)
}


#' @title Hscore alpha = 0
#'
#' @description Calculates Hscore with alpha equal to 0
#'
#' @param mat a matrix
#'
#' @return The hscore of \code{mat}
#'
#' @export
#'
#' @examples
#' mat <- replicate(5, rnorm(5))
#' alphazeroccscore(mat)

alphazeroccscore <- function(mat)
{
  if(nrow(mat)>1){
    score <- sum((mat - matrix(colMeans(mat),nrow=nrow(mat),ncol=ncol(mat),byrow=TRUE)) ^2)/(nrow(mat)*ncol(mat))
  }else{
    score <- 0
  }
  return(score)
}

#' @title Hscore beta = 0
#'
#' @description Calculates Hscore with beta equal to 0
#'
#' @param mat a matrix
#'
#' @return The hscore of \code{mat}
#'
#' @export
#'
#' @examples
#' mat <- replicate(5, rnorm(5))
#' betazeroccscore(mat)

betazeroccscore <- function(mat)
{
  score <- sum((mat - rowMeans(mat)  + mean(mat)) ^2)/(nrow(mat)*ncol(mat))
  score
}

#' @title Two functions Hscore
#'
#' @description Calculates the Hscore of the matrix composed by function \code{x}
#'   and function \code{y} as rows
#'
#' @param x a function
#' @param y a function
#'
#' @return The hscore of a matrix having \code{x} and \code{y} as rows
#'
#' @examples
#' x <- runif(5, min=0, max=100)
#' y <- runif(5, min=0, max=100)
#' d_ccscore(x,y)

d_ccscore <- function(x, y)
{
  mat <- as.matrix(rbind(x,y))
  ccscore(mat)
}




