# Preprocession and log-normalization

log_normalization = function(x, percent, preprocess.only = FALSE){

  x <- as.matrix(x)

  if (preprocess.only){

    n <- dim(x)[2]

    gene.exprs.count <- rowSums(x != 0)

    x <- x[gene.exprs.count > n * percent, ]

    return(x)

  } else {

    n <- dim(x)[2]

    gene.exprs.count <- rowSums(x != 0)

    x <- x[gene.exprs.count > n * percent, ]

    sf <- colSums(x)/median(colSums(x))

    return(log(sweep(x, 2, sf, '/')+1))

  }

}


