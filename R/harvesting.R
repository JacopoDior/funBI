#' @title Delta-Threshold Cutter
#'
#' @description Cuts a Hscore dendrogram according to the threshold \code{delta}
#'   finding biclusters whose hscore is lesser than \code{delta}.
#'
#' @param mat a matrix
#' @param hc a hierarchy resulting from DIANA built from \code{mat}
#' @param delta a scalar value indicating the delta-threshold
#'
#' @return A list containing all the candidate biclusters whose Hscore is lesser
#'   than \code{delta}. Every single element of the list is a candidate
#'   bicluster and the elements belonging to it are reprented using the
#'   row index/label of the matrix \code{mat}.
#'
#' @export
#'
#' @examples


delta_cutter <- function(mat, hc, delta){

  # create list of all biclusters coming from the merging procedure
  mergelist <- list()
  for( i in 1:dim(hc$merge)[1])
  {
    if(sum(hc$merge[i,]<0)==2){ #both negative
      mergelist[[i]] <- abs(hc$merge[i,][which(hc$merge[i,]<0)])

    }else if(sum(hc$merge[i,]<0)==1){ # one negative
      mergelist[[i]] <- c( abs(hc$merge[ i,][which(hc$merge[i,]<0)]),
                           mergelist[[ hc$merge[i,][which(hc$merge[i,]>0)] ]] )
    }else{ #both positive
      mergelist[[i]] <- c( mergelist[[ hc$merge[i,][1] ]], mergelist[[ hc$merge[i,][2] ]])
    }
  }

  # delete dangerous behaviours
  new_height <- deletedanger(hc)[,3]

  # find delta-biclusters
  cluster <- hcut(new_height, mergelist, delta, mat)
  return(cluster)
}


#' @title Delete Inversions
#'
#' @description Identifies and solves branches Inversions
#'
#' @param hc a hierarchy resulting from DIANA
#'
#' @return A vector with the new hscore-based height to avoid having branches
#'   inversions in the dendogram
#'
#' @examples
#'

deletedanger <- function(hc){
  # create df composed by merge e hscore
  temp <- cbind(hc$merge,hc$height)

  danger <- c()
  for(i in dim(temp)[1]:1){
    check <- which(temp[i,c(1,2)]>0)
    if(length(check)!=0){
      if(any(temp[i,3] < temp[temp[i,check],3])){
        danger <- c(danger,i)
        case1 <- temp[i,1]
        case2 <- temp[i,2]
        if(temp[case1,3]>temp[i,3])
          temp[case1,3] <- temp[i,3]*(1-0.000001)

        if(temp[case2,3]>temp[i,3])
          temp[case2,3] <- temp[i,3]*(1-0.000001)
      }
    }
  }
  return(temp)
}

#' @title Cutter
#'
#' @description Cuts a Hscore dendrogram according to the threshold \code{delta}
#'   finding biclusters whose hscore is lesser than \code{delta}.
#'
#' @param new_height a vector with the hscore-based height (modified or not) of
#'   the dendogram built from \code{delta}
#' @param mergelist a list of all the possible biclusters one can get from \code{hc}
#' @param delta a scalar value indicating the delta-threshold
#' @param mat a matrix
#'
#' @return A list containing all the candidate biclusters whose Hscore is lesser
#'   than \code{delta}. Every single element of the list is a candidate
#'   bicluster and the elements belonging to it are reprented using the
#'   row index/label of the matrix \code{mat}.
#'
#' @examples
#'

hcut <- function(new_height, mergelist, delta,mat){
  temp <- c()
  cluster <- rep(0,dim(mat)[1])
  check <- sort(which(new_height < delta), decreasing = TRUE)

  if(length(check)!=0){
    j <- 1
    cluster[mergelist[[check[1]]]] <- j
    temp <- mergelist[[check[1]]]
    for( i in check){
      if(all(mergelist[[i]] %in% temp)==FALSE){
        temp <- c(temp,mergelist[[i]])
        j <- j+1
        cluster[mergelist[[i]]] <- j
      }
    }
    cluster[which(cluster==0)] <- (j+1):(j+length(cluster[which(cluster==0)]))
  }else{
    cluster <- 1:dim(mat)[1]
  }
  return(cluster)
}


#' @title funBI Harvesting
#'
#' @description Perform Harvesting to a list of DIANA hierarcy obetained using
#'   the \code{funBI_lot_and_flower} function. The available harvesting
#'   strategies are three: a \code{delta}-threshold strategy, a dynamic strategy
#'   with a threshold obtained as a percentage of the total hscore of a
#'   considered hierarchical dendrogram, and an automatic strategy using
#'   \code{gap} statistics.
#'
#' @param list_hc a list of hierarchy resulting from DIANA built from \code{mat}
#' @param delta a scalar value indicating the delta-threshold in order to
#'   perform a delta-bicluster. In this case \code{perc} and \code{gap} should
#'   be missing.
#' @param perc a scalar value indicated a percentage of the total score in order
#'   to perform a delta-bicluster cut  using as delta-threshold the
#'   percentage \code{perc} of the total hscore of a particular dendrogram. In
#'   this case \code{delta} and \code{gap} should be missing.
#' @param gap TRUE/FALSE indicating if gap statistics should be used as
#'   Harvesting strategy.In this case \code{delta} and \code{perc} should be
#'   missing.
#' @param mat a matrix
#' @param Kmax the maximum number of clusters to consider when performing gap
#'   statistics (\code{gap} = TRUE), must be at least two.
#' @param B integer value, number of Monte Carlo (“bootstrap”) samples to
#'   consider when performing gap statistics (\code{gap} = TRUE).
#'
#' @return \code{biclist}, a list containing all the candidate biclusters
#'   according to the used Harvesting strategy. Every element of the list is a
#'   bicluster. Every bicluster is characterized by the following components:
#'   \enumerate{ \item{a vector containing the elements belonging to the
#'   bicluster reprented using the row index/label of the matrix \code{mat}}
#'   \item{a vector containing the boundaries of the sub-interval
#'   considered represented usign the colum index/label of the matrix
#'   \code{mat}} \item{the Hscore of the bicluster} \item{The index
#'   of \code{list_hc}, representing which dendogram the bicluster was cut} }
#'
#' @export
#'
#' @examples
#'



funBI_harvesting <- function(list_hc, delta, perc, gap, mat, Kmax, B) {

  if(missing(delta)==FALSE){
    # creo lista all biclusters
    biclist <- list()
    k <- 1
    n <- length(list_hc)
    pb <- progress::progress_bar$new(total = n)

    for (i in 1:n){
      pb$tick()
      bic <- delta_cutter(mat, list_hc[[i]], delta)

      for(j in 1:max(bic)){
        biclist[[k]] <- list(which(bic==j),
                             list_hc[[i]]$elements,
                             ccscore(as.matrix(mat[which(bic==j),list_hc[[i]]$elements[1]:list_hc[[i]]$elements[2]])),
                             i)
        k <- k + 1
      }
      Sys.sleep(1 / length(n))
    }
  }

  if(missing(perc)==FALSE){
    # creo lista all biclusters
    biclist <- list()
    k <- 1
    n <- length(list_hc)
    pb <- progress::progress_bar$new(total = n)
    for (i in 1:n){
      pb$tick()
      totalscore <- ccscore(as.matrix(mat[ ,list_hc[[i]]$elements[1]:list_hc[[i]]$elements[2] ]))
      delta <- totalscore*perc
      bic <- delta_cutter(mat, list_hc[[i]], delta)

      for(j in 1:max(bic)){
        biclist[[k]] <- list(which(bic==j),
                             list_hc[[i]]$elements,
                             ccscore(as.matrix(mat[which(bic==j),list_hc[[i]]$elements[1]:list_hc[[i]]$elements[2]])),
                             i)
        k <- k + 1
      }
      Sys.sleep(1 / length(n))
    }
  }

  if(missing(gap)==FALSE){
    biclist <- list()
    k <- 1
    n <- length(list_hc)
    if(missing(Kmax))
      Kmax <- 20

    if(missing(B))
      B <- 100

    pb <- progress::progress_bar$new(total = n)

    for (i in 1:n){
      pb$tick()
      gapstat <- gap_clust(x=mat, tree=list_hc[[i]], K.max = Kmax, B=B, d.power=1)
      bestk <- with(gapstat,maxSE(Tab[,"gap"],Tab[,"SE.sim"]))
      bic <- cutree(list_hc[[i]], k = bestk)

      for(j in 1:max(bic)){
        biclist[[k]] <- list(which(bic==j),
                             list_hc[[i]]$elements,
                             ccscore(as.matrix(mat[which(bic==j),list_hc[[i]]$elements[1]:list_hc[[i]]$elements[2]])),
                             i)
        k <- k + 1
      }
      Sys.sleep(1 / length(n))
    }
  }
  return(biclist)
}

#' @title Gap Statistics
#'
#' @description calculates a goodness of clustering measure, the \code{gap} statistic.
#'
#' @param mat a matrix
#' @param hc a hierarchy resulting from DIANA
#' @param K.max the maximum number of clusters to consider, must be at least
#'   two.
#' @param B integer value, number of Monte Carlo (“bootstrap”) samples
#'
#' @return An object of S3 class "clusGap", basically a list with components:
#'   \describe{ \item{Tab}{a matrix with K.max rows and 4 columns, named "logW",
#'   "E.logW", "gap", and "SE.sim"} \item{call}{the function call}
#'   \item{spaceH0}{the \code{spaceH0} argument} \item{n}{number of observation}
#'  \item{B}{input \code{B}}}
#'
#' @examples
#'



gap_clus <- function(mat, hc, K.max, B, d.power = 1, spaceH0 = c("scaledPCA", "original")){
  stopifnot(length(dim(mat)) == 2, K.max >= 2, (n <- nrow(mat)) >= 1, ncol(mat) >= 1)
  if (B != (B. <- as.integer(B)) || (B <- B.) <= 0) {
    stop("'B' has to be a positive integer")}
  cl. <- match.call()

  if (is.data.frame(mat)){
    mat <- as.matrimat(mat)}

  ii <- seq_len(n)

  W.k <- function(mat, hc, kk) {
    clus <- cuhc(hc, k = kk)
    0.5 * sum(vapply(split(ii, clus), function(I) {
      mats <- mat[I, , drop = FALSE]
      sum(promaty::dist(mats,newccscore)^d.power/nrow(mats))
    }, 0))
  }

  logW <- numeric(K.max)
  E.logW <- logW
  SE.sim <- logW

  for (k in 1:K.max){ logW[k] <- log(W.k(mat, hc, k)) }

  spaceH0 <- match.arg(spaceH0)
  mats <- scale(mat, center = TRUE, scale = FALSE)
  m.mat <- rep(attr(mats, "scaled:center"), each = n)

  switch(spaceH0, scaledPCA = {
    V.smat <- svd(mats, nu = 0)$v
    mats <- mats %*% V.smat
  }, original = {
  }, stop("invalid 'spaceH0':", spaceH0))

  rng.mat1 <- apply(mats, 2L, range)
  logWks <- matrimat(0, B, K.max)

  for (b in 1:B) {
    z1 <- apply(rng.mat1, 2, function(M, nn) runif(nn, min = M[1],
                                                 mamat = M[2]), nn = n)
    z <- switch(spaceH0, scaledPCA = tcrossprod(z1, V.smat),
                original = z1) + m.mat
    hc.z <- diana(z)
    for (k in 1:K.max) {
      logWks[b, k] <- log(W.k(mat = z, hc = hc.z, kk = k))
    }
  }

  E.logW <- colMeans(logWks)
  SE.sim <- sqrt((1 + 1/B) * apply(logWks, 2, var))
  #return(list(Tab = cbind(logW, E.logW,  gap = E.logW - logW, SE.sim)))
  structure(class = "clusGap", list(Tab = cbind(logW, E.logW,
                                                gap = E.logW - logW, SE.sim), call = cl., spaceH0 = spaceH0,
                                    n = n, B = B))
}
