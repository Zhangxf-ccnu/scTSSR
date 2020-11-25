keras_lasso_regression <- function(Y, X,  epochs = 100, batch_size = 128, lambda1 = 1.0, lambda2 = 1e10,
                                   learning_rate = 0.0001, lasso_threshold = 0, verbose = TRUE){

  p <- ncol(X)
  n <- nrow(X)

  X <- array_reshape(X, dim = c(nrow(X), ncol(X)))
  Y <- array_reshape(Y, dim = c(nrow(X), ncol(X)))


  network <- keras::keras_model_sequential() %>%
    layer_dense(units = p, input_shape = p, activation = "linear", kernel_regularizer = function(weight_matrix)
      regularizer_define(weight_matrix = weight_matrix ,lambda1 = lambda1,  lambda2 = lambda2),
      use_bias = FALSE, kernel_constraint = constraint_nonneg(), kernel_initializer = "zeros")


  network %>% keras::compile(loss = loss_define, optimizer = optimizer_rmsprop(lr = learning_rate), metrics = list("mean_squared_error"))
  history <- network %>% keras::fit(x = X, y = Y, epochs = epochs, batch_size = batch_size, validation_split = 0,
                                    verbose = verbose)


  weight <- keras::get_weights(network)[[1]]

  weight_final <- weight

  weight_final[weight_final <= lasso_threshold] <- 0

  return(weight_final)

}

