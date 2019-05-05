#' @title Delta-Threshold Cutter
#'
#' @description Cuts a Hscore dendrogram according to the threshold \code{delta}
#'   finding biclusters whose hscore is lesser than \code{delta}.
#'
#' @param biclist a list of candidate biclusters, resulting from
#'   \code{funBI_harvesting}
#'
#' @return A reduced version of \code{biclist} containing all the candidates
#'   bicluster which are not nested into other candidates biclusters.
#'
#' @export
#'
#' @examples


funBI_tasting <- function(biclist){
  delete <- c()

  for(i in 1:(length(biclist) - 1)){
    # if one is already deleted it is not useful to use it for comparison
    for(j in ((i + 1):length(biclist))[!(i + 1):length(biclist) %in% delete]){

      # verify if the domain of i is inside the domain of j
      if(biclist[[i]][[2]][1]>=biclist[[j]][[2]][1] & biclist[[i]][[2]][2]<=biclist[[j]][[2]][2]){

        # calculate the intersection
        A <- intersect(biclist[[i]][[1]], biclist[[j]][[1]])
        # if there are some element in common
        if(length(A)!=0){
          # if the interesection is equal to the considered element
          if(all(A %in% biclist[[i]][[1]]) && (length(biclist[[i]][[1]])==length(A))){
            delete <- c(delete,i)
            break
          }
        }
      }
    }
  }
  return(delete)
}
