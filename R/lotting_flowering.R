#' @title Lotting
#'
#' @description Divides the domain in sub-intervals of lenght greater than
#'   \code{int}. The smallest sub-intervals have length equal to \code{int} while the
#'   others are created by shifting by \code{int} along the continuous dimension and by
#'   enlarging them by multiples of \code{int} in order to cover the whole T
#'
#' @param mat a matrix
#' @param int the minimum length of sub-intervarls
#'
#' @return A list containing all the possible starts and lengths of the sub-intervals
#'
#' @export
#'
#' @examples
#' mat <- replicate(5, rnorm(5))
#' int <- 2
#' funBI_lotting(mat, int)

funBI_lotting <- function(mat, int){
  positions <-  seq(1, by=int-1,ncol(mat))
  elem <- unique(c(int * (1:floor(ncol(mat) / int)), ncol(mat)))
  return(list(positions,elem))
}

#' @title Diana Tree
#'
#' @description Performs hscore based DIANA (DIvisive ANAlysis) algorithm
#'
#' @param mat a matrix
#'
#' @return A DIANA hierarchy, i.e. a DIANA list with hscore height
#'
#' @export
#'
#' @examples
#' mat <- replicate(20, rnorm(15))
#' funBI_diana(mat)

funBI_diana <- function(mat){
  mydist <- proxy::dist(mat,d_ccscore)
  hc <- cluster::diana(mydist, diss=TRUE, keep.diss = TRUE)
  addinfo <- score_merging(mat,hc)
  hc$height <- addinfo[[1]]
  hc$position <- addinfo[[2]]
  return(hc)
}

#' @title Score Merger
#'
#' @description Gets the height of every single DIANA merging step
#'
#' @param mat a matrix
#' @param hc a hierarchy resulting from DIANA
#'
#' @return A list containing the hscores of every single DIANA merging steps and
#'   the resulting groups of those steps
#'
#' @export
#'
#' @examples


score_merging <- function(mat, hc){

  # create list of all biclusters coming from the merging procedure
  mergelist <- list()
  for(i in 1:dim(hc$merge)[1])
  {
    if(sum(hc$merge[i, ] < 0) == 2){ #both negative
      mergelist[[i]] <- abs(hc$merge[i, ][which(hc$merge[i, ]<0)])

    }else if(sum(hc$merge[i, ] < 0) == 1){ # one negative
      mergelist[[i]] <- c(abs(hc$merge[i, ][which(hc$merge[i, ] < 0)]),
                           mergelist[[hc$merge[i, ][which(hc$merge[i, ] > 0)] ]] )
    }else{ #both positive
      mergelist[[i]] <- c(mergelist[[hc$merge[i, ][1]]], mergelist[[hc$merge[i, ][2]]])
    }
  }

  mergehscore <- c()
  for(i in 1:length(mergelist)){
    mergehscore[i] <- ccscore(as.matrix(mat[mergelist[[i]], ]))
  }

  return(list(mergehscore,mergelist[[i]]))
}


#' @title Lotting and Flowering
#'
#' @description Performs the Lotting and Flowering steps by dividing the whole
#'   domain in sub-intervals of length greater than \code{int} and by applying
#'   DIANA in every sub-interval
#'
#' @param mat a matrix
#' @param int the minimum length of sub-intervarls
#'
#' @return A list of DIANA hierarchy, one for each sub-intervals, with the
#'   reference of starting and ending of the corresponding intervals
#'
#' @export
#'
#' @examples
#' mat <- replicate(5, rnorm(5))
#' int <- 2
#' funBI_lot_and_flower(mat,int)

funBI_lot_and_flower <- function(mat,int){

  list_hc <- list()
  k <- 1

  # Lotting
  lotting <- funBI_lotting(mat, int)
  positions <- lotting[[1]]
  elem <- lotting[[2]]

  # Flower
  for (i in 1: length(elem)){
    for (start in positions){
      if(start + elem[i]-1 <= ncol(mat)){ # check to be inside the original columns
        list_hc[[k]] <- funBI_diana(mat[, start:(start+elem[i]-1)])
        list_hc[[k]]$elements <- c(start, start+elem[i]-1) # to know which column
        k <- k + 1
      }
    }
  }
  return(list_hc)
}

