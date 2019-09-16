clean.data <- function(x) {

  if (!(grepl("matrix", class(x), ignore.case = TRUE))) {

    x <- Matrix::Matrix(as.matrix(x))

    message("Converting x to matrix.")

    if (!is.numeric(x)) {

      warning("Make sure x is numeric.")

    }

  }

  np <- dim(x)

  size <- as.numeric(np[1])*as.numeric(np[2])

  if(size > 2^31-1){

    inds <- split(1:np[2], ceiling(1:np[2]/1000))

    for(i in 1:length(inds)){

      x[, inds[[i]]][x[, inds[[i]]] < 0.001] <- 0

    }

  } else {

    x[x < 0.001] <- 0

  }

  if (is.null(np) | (np[2] <= 1))

    stop("x should be a matrix with 2 or more columns")

  if (min(Matrix::colSums(x)) == 0) {

    nzerocells <- sum(Matrix::colSums(x) == 0)

    x <- x[, Matrix::colSums(x) != 0]

    message("Removing ", nzerocells, " cell(s) with zero expression.")

  }

  if (is.null(rownames(x))) {

    rownames(x) <- 1:np[1]

  }

  x

}





calc.size.factor <- function(x) {

  sf <- Matrix::colSums(x)/mean(Matrix::colSums(x))

  scale.sf <- 1

  list(unname(sf), scale.sf)

}



#update.output <- function(f, ind, start, stop, out, x, sf, scale.sf, mu,

#                          nworkers, estimates.only) {

#  n <- stop-start+1

#  ind1 <- ind[start:stop]

#  results <- f(x[ind1, , drop = FALSE], sf, scale.sf,

#               mu[ind1, , drop = FALSE], nworkers, estimates.only)

#  out$estimate[ind1, ] <- results$est

#  if (!estimates.only) {

#    out$se[ind1, ] <- results$se

#  }

#  out$info[[2]][ind1] <- results[[8]]

#  return(out)

#}



calc.lambda <- function(X_count, percent){

  X <- log_normalization(X_count, percent, preprocess.only = FALSE)

  lambda <- sd(X)

}




shrinkage <- function(X, lambda, penalize_diagonal){

  zeros = matrix(0, nrow(X), ncol(X))

  Bbar <- pmax(zeros, X-lambda) - pmax(zeros, -X-lambda)

  if (!penalize_diagonal){

    diag(Bbar) <- 0

  }

  Bbar <- pmax(Bbar, zeros)

  return(Bbar)

}


