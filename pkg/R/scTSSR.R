#' use scTSSR to impute dropout values in scRNA-seq data
#'
#'
#' @param X_count An expression count matrix. The rows correspond to genes and
#' the columns correspond to cells. Can be sparse.
#'
#' @param lambda Tuning parameter to facilitate feature selection and regularization.
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
#' @param gamma The step size.
#'
#' @param beta The line search parameter.
#'
#' @param max_iter Maximum iteration of the accelerated proximal gradient descent algorithm.
#'
#' @param abstol Absolute tolerance of the accelerated proximal gradient descent algorithm.
#'
#' @param penalize_diagonal Whether penalize the diagonal elements of
#' the regression coefficient matrix. Default is FALSE.
#'
#' @param estimates.only Only return the imputed data. Default is FALSE.
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
#' @import SAVER
#'
#' @author Ke Jin, \email{kej13@mails.ccnu.edu.cn}
#'
#' @examples
#'
#' data("baron")
#' baron_imputation_result = scTSSR(baron$count.samp)
#'


scTSSR <- function(X_count, lambda = NULL, initA = NULL, initB = NULL,

                   percent = 0.1, ncores = 1, MAX_ITER = 4, ABSTOL = 1e-3,

                   gamma = 1, beta = 0.5, max_iter = 100, abstol = 1e-4,

                   penalize_diagonal = FALSE, estimates.only = FALSE){

  # calculating the tunning parameter lambda

  if (is.null(lambda)) {

    message("Calculating the penalty parameter lambda ...")

    lambda <- calc.lambda(X_count, percent)

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

    # using accelerated proximal gradient descent to calculate B or A

    B <- proximal_gradient(Y, Z, lambda, B_old, gamma = gamma, beta = beta,

                           MAX_ITER = max_iter, ABSTOL = abstol, penalize_diagonal = penalize_diagonal)

    if (k > 1){

      B <- (B_old + B)/2

    } else {

      B <- B

    }

    Z <- t(pmax(X%*%(diag(ncells)+B), zeros))

    A_old <- A

    A <- t(proximal_gradient(t(Y), Z, lambda, t(A_old), gamma = gamma, beta = beta,

                             MAX_ITER = max_iter, ABSTOL = abstol, penalize_diagonal = penalize_diagonal))

    if (k > 1){

      A <- (A_old + A)/2

    } else {

      A <- A

    }

    history.objval[k] <- objective(X, lambda, A, B)

    message(k, '    ', round(history.objval[k], 3))

    if (k > 1 && abs(history.objval[k] - history.objval[k-1]) < ABSTOL) break

    k <- k + 1

  }

  Xhat <- pmax(A%*%X + X%*%B + A%*%X%*%B, zeros)

  pred.time <- Sys.time() - pred.st

  message("Calculating the posterior means with ", ncores, " worker(s)")

  posterior.st <- Sys.time()

  out <- SAVER::saver(X_count, ncores = ncores, mu =  exp(as.matrix(Xhat)))

  posterior.time <- Sys.time() - posterior.st

  total.time <- Sys.time() - pred.st

  result$estimate <- out$estimate

  result$se <- out$se

  result$info$pred.time <- pred.time

  result$info$posterior.time <- posterior.time

  result$info$total.time <- total.time

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


