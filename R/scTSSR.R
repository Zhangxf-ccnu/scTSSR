#' use scTSSR to impute dropout values in scRNA-seq data
#'
#' @param X_count An expression count matrix. The rows correspond to genes and
#' the columns correspond to cells. Can be sparse.
#'
#' @param lambda1 Tuning parameter to facilitate feature selection and regularization.
#'
#' @param lambda2 Tuning parameter to penalize the diagonal elements of the parameter to
#' eliminate the trivial solution of representing an expression level as a linear combination
#' of itself.
#'
#' @param initA The initionlization of A. The elements of A represent the similarities
#' between genes.
#'
#' @param initB The initionlization of B. The elements of B represent the similarities
#' between cells.
#'
#' @param percent The expression count matrix is preprocessed by filtering out the genes
#' expressed in at most percent*\eqn{100\%} of the cells.
#'
#' @param ncores Number of cores to use. Default is 1.
#'
#' @param MAX_ITER Maximum iteration of the external circulation of scTSSR.
#'
#' @param ABSTOL Absolute tolerance of the external circulation.
#'
#' @param learning_rate A hyper-parameter that controls the speed of adjusting the weights of the network
#' with respect to the loss gradient.
#'
#' @param epochs The number of the entire training set going through the entire network.
#'
#' @param batch_size The number of examples that are fed to the algorithm at a time.
#'
#' @param run_batch Whether to use batch or to set the number of all the samples as the value of the batch size.
#' Default is TRUE.
#'
#' @param verbose Whether to output the value of metrics at the end of each epoch. Default is TRUE.
#'
#'
#'
#' @return If `estimates.only = TRUE', then a matrix of scTSSR estimates.
#'
#' If `estimates.only = FALSE', a list with the following components
#'
#' \item{\code{estimate}}{Recovered (normalized) expression.}
#'
#' \item{\code{se}}{Standard error of estimates.}
#'
#' \item{\code{info}}{Information about dataset.}
#'
#'
#' The \code{info} element is a list with the following components:
#'
#' \item{\code{size.factor}}{Size factor used for normalization.}
#'
#' \item{\code{pred.time}}{Time taken to generate predictions.}
#'
#' \item{\code{posterior.time}}{Time taken to compute the posterior distribution.}
#'
#' \item{\code{total.time}}{Total time for scTSSR estimation.}
#'
#'
#' @export
#'
#'
#' @import keras
#'
#' @import tensorflow
#'
#' @import SAVER
#'
#' @author Ke Jin, \email{kej13@mails.ccnu.edu.cn}
#'
#' @examples
#'
#' data("baron")
#'
#' baron_imputation_result = scTSSR(baron$count.samp)
#'
scTSSR <- function(X_count, lambda1 = NULL, lambda2 = 1e10, initA = NULL, initB = NULL,

                   percent = 0.05, ncores = 1, MAX_ITER = 4, ABSTOL = 1e-3, learning_rate = 0.0001,

                   epochs = 100, batch_size = 128, run_batch = TRUE, verbose = TRUE, estimates.only = FALSE){

  # calculating the tunning parameter lambda

  if (is.null(lambda1)) {

    message("Calculating the penalty parameter lambda ...")

    lambda1 <- calc.lambda(X_count, percent)

    message("Done!")

  }


  # preprocessing and log-normalization

  message("Starting preprocessing and log-normalization ...")

  X <- log_normalization(X_count, percent, preprocess.only = FALSE)

  X_count <- log_normalization(X_count, percent, preprocess.only = TRUE)

  message("Done!")


  X_count <- clean.data(X_count)

  ngenes <- nrow(X_count)

  ncells <- ncol(X_count)

  gene.names <- rownames(X_count)

  cell.names <- colnames(X_count)


  # assign size factor

  sf.out <- calc.size.factor(X_count)

  sf <- sf.out[[1]]

  scale.sf <- sf.out[[2]]



  result <- list()

  result$estimate <- matrix(0, ngenes, ncells, dimnames = list(gene.names, cell.names))

  result$A <- matrix(0, ngenes, ngenes)

  result$B <- matrix(0, ncells, ncells)
  
  result$obj.prior <- c()  #############
  
  

  if (!estimates.only) {

    result$se <- matrix(0, ngenes, ncells, dimnames = list(gene.names, cell.names))

  } else {

    result$se <- NA

  }

  result$info <- c(list(0), list(0), list(0), list(0))

  names(result$info) <- c("size.factor", "pred.time", "posterior.time", "total.time")

  result$info$size.factor <- scale.sf*sf



  # initialize A and B

  if (is.null(initA)) {

    initA <- matrix(0, ngenes, ngenes)

  }

  if (is.null(initB)) {

    initB <- matrix(0, ncells, ncells)

  }

  A <- initA

  B <- initB


  message("Imputation starts ...")

  message("Calculating the prior mean for each gene in each cell ...")

  message("iter", ' objective')

  pred.st <- Sys.time()

  k <- 1

  history.objval <- c()

  zeros = matrix(0, ngenes, ncells)

  while (k <= MAX_ITER){

    Z <- pmax((diag(ngenes)+A)%*%X, zeros)

    Y <- X - A%*%X

    B_old <- B

    if (!run_batch){

      batch_size = nrow(Z)

    } else{

      batch_size = batch_size

    }

    # using keras with SGD algorithm to calculate B or A

    B <- keras_lasso_regression(Y, Z, epochs = epochs, batch_size = batch_size, lambda1 = lambda1, lambda2 = lambda2,

                                learning_rate = learning_rate, verbose = verbose)

    if (k > 1){

      B <- (B_old + B)/2

    } else {

      B <- B

    }

    Z <- t(pmax(X%*%(diag(ncells)+B), zeros))

    Y <- X - X%*%B

    A_old <- A

    if (!run_batch){

      batch_size = nrow(Z)

    } else{

      batch_size = batch_size

    }


    A <- t(keras_lasso_regression(t(Y), Z, epochs = epochs, batch_size = batch_size, lambda1 = lambda1, lambda2 = lambda2,

                                  learning_rate = learning_rate, verbose = verbose))


    if (k > 1){

      A <- (A_old + A)/2

    } else {

      A <- A

    }

    history.objval[k] <- objective(X, lambda1, A, B)

    message(k, '    ', round(history.objval[k], 3))

    if (k > 1 && abs(history.objval[k] - history.objval[k-1]) < ABSTOL) break

    k <- k + 1

  }

  Xhat <- pmax(A%*%X + X%*%B + A%*%X%*%B, zeros)
  
  obj.prior <- history.objval #########
  

  pred.time <- Sys.time() - pred.st

  message("Calculating the posterior means with ", ncores, " worker(s)")

  posterior.st <- Sys.time()

  out <- SAVER::saver(X_count, ncores = ncores, mu =  exp(as.matrix(Xhat)) )

  posterior.time <- Sys.time() - posterior.st

  total.time <- Sys.time() - pred.st

  result$estimate <- out$estimate

  result$se <- out$se

  result$info$pred.time <- pred.time

  result$info$posterior.time <- posterior.time

  result$info$total.time <- total.time

  result$A <- A

  result$B <- B
  result$obj_prior <- obj.prior    ###########
  

  message("Done!")

  message("Finish time: ", Sys.time())

  message("Total time: ", format(result$info$total.time))

  if (!estimates.only) {

    class(result) <- "scTSSR"

    result

  } else {

    result$estimate

  }

}
