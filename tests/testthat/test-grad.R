compare_anlytical_and_numerical_grad <- function(
  model_name, n_obs = 32, n_pred = 4, n_test = 10, data_seed = 1918, loc_seed = 615
) {
  n_obs <- 32; n_pred <- 4
  data <- simulate_data(n_obs, n_pred, model_name, seed = data_seed)
  design <- data$design; outcome <- data$outcome
  model <- new_regression_model(design, outcome, model_name)
  loglik_func <- function (coef) { 
    model$reg_coef = coef
    return(calc_loglik(model))
  }
  grad_func <- function (coef) {
    if (model_name == "linear") {
      model$reg_coef = coef
      calc_grad(model, list(noise_var=1))
    } else {
      model$reg_coef = coef
      calc_grad(model)
    }
  }
  set.seed(loc_seed)
  grads_are_close <- TRUE
  for (i in 1:n_test) {
    if (!grads_are_close) break
    regcoef <- rnorm(n_pred)
    analytical_grad <- grad_func(regcoef)
    numerical_grad <- approx_grad_via_finite_diff(loglik_func, regcoef)
    grads_are_close <- are_all_close(
      analytical_grad, numerical_grad, abs_tol = Inf, rel_tol = 1e-3
    )
  }
  return(grads_are_close)
}

test_that("linear model's analytical gradient is close to numerical one", {
  expect_true(
    compare_anlytical_and_numerical_grad("linear")
  )
})

test_that("logit model's analytical gradient is close to numerical one", {
  expect_true(
    compare_anlytical_and_numerical_grad("logit")
  )
})

test_that("poisson model's analytical gradient is close to numerical one", {
  skip_if_not(Sys.getenv("MANUAL_TEST") != "")
  expect_true(
    compare_anlytical_and_numerical_grad("poisson")
  )
})