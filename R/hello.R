#' @title Hscore
#'
#' @description Calculate Hscore with both alpha and beta both different to 0
#'
#' @param mat a matrix
#'
#' @return The hscore of \code{mat}
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


