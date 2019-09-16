proximal_gradient <- function(X, Z, lambda, pre_B, gamma, beta,

                              MAX_ITER, ABSTOL, penalize_diagonal){

  ngenes <- nrow(Z)

  ncells <- ncol(Z)

  k <- 1

  history.objval <- c()

  # warm start

  B <- pre_B

  Bprev <- B

  ZtX <- t(Z)%*%X

  if(ngenes < ncells){

    while (k <= MAX_ITER){

      Q <- B + (k/(k+3))*(B-Bprev)

      grad_Q <- t(Z)%*%(Z%*%Q) - ZtX

      while (1){

        Bbar <- shrinkage(Q - gamma*grad_Q, lambda*gamma, penalize_diagonal)

        loss_B <- Bbar - Q

        if (f(X, Z, Bbar) <= f(X, Z, Q) + sum(grad_Q*loss_B) + (1/(2*gamma))*norm(loss_B, 'F')^2) break

        gamma <- beta*gamma

      }

      Bprev <- B

      B <- Bbar

      history.objval[k] <- pg_objective(X, Z, lambda, B)

      if (k > 1 && abs(history.objval[k] - history.objval[k-1]) < ABSTOL) break

      k <- k + 1

    }

    return(B)

  } else {

    ZtZ <- t(Z)%*%Z

    while (k <= MAX_ITER){

      Q <- B + (k/(k+3))*(B-Bprev)

      grad_Q <- ZtZ%*%Q - ZtX

      while (1){

        Bbar <- shrinkage(Q - gamma*grad_Q, lambda*gamma, penalize_diagonal)

        loss_B <- Bbar - Q

        if (f(X, Z, Bbar) <= f(X, Z, Q) + sum(grad_Q*loss_B) + (1/(2*gamma))*norm(loss_B, 'F')^2) break

        gamma <- beta*gamma

      }

      Bprev <- B

      B <- Bbar

      history.objval[k] <- pg_objective(X, Z, lambda, B)

      if (k > 1 && abs(history.objval[k] - history.objval[k-1]) < ABSTOL) break

      k <- k + 1

    }

    return(B)

  }

}
