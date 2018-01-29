hamming <- function(X, Y) {
  if ( missing(Y) ) {
    uniqs <- unique(as.vector(X))
    U <- X == uniqs[1]
    H <- t(U) %*% U
    for ( uniq in uniqs[-1] ) {
      U <- X == uniq
      H <- H + t(U) %*% U
    }
  } else {
    uniqs <- union(X, Y)
    H <- t(X == uniqs[1]) %*% (Y == uniqs[1])
    for ( uniq in uniqs[-1] ) {
      H <- H + t(X == uniq) %*% (Y == uniq)
    }
  }
  nrow(X) - H
}
