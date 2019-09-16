
objective <- function(X, lambda, A, B){

  objval <- (0.5*norm(X - A%*%X - X%*%B - A%*%X%*%B, 'F')^2 + lambda*sum(abs(A)) + lambda*sum(abs(B)))

  return(objval)
}





pg_objective <- function(X, Z, lambda, B){

  objval <- (0.5*norm(X - Z%*%B, 'F')^2 + lambda*sum(abs(B)))

  return(objval)

}




f <- function(X, Z, B){

  cost <- 0.5*norm(X - Z%*%B, 'F')^2

  return(cost)

}


