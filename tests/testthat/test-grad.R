compare_anlytical_and_numerical_grad <- function(
  model_name, n_obs = 32, n_pred = 4, n_test = 10, data_seed = 1918, loc_seed = 615
) {
  n_obs <- 32; n_pred <- 4
  data <- simulate_data(n_obs, n_pred, model_name, seed = data_seed)
  design <- data$design; outcome <- data$outcome
  # The `do.call` trick below might seem like a clever solution and a 
  # simpler-to-implement alternative to S3 methods. However, it tends to make 
  # the code less readable and harder to maintain. For example, you wouldn't 
  # know the exact function being called by `sprintf("calc_%s_loglik", model)` 
  # without reading other parts of code. It also obscures the usage of those 
  # functions within codebase; e.g. `git grep "calc_logit_grad"` would fail to 
  # detect its usage here. The trick is used here only for didactic/illustrative 
  # purpose --- I advise against its use in general.
  loglik_func <- function (coef) { 
    do.call(
      sprintf("calc_%s_loglik", model_name), 
      list(coef, design, outcome)
    )
  }
  grad_func <- function (coef) {
    do.call(
      sprintf("calc_%s_grad", model_name), 
      list(coef, design, outcome)
    )
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
  expect_true(
    compare_anlytical_and_numerical_grad("poisson")
  )
})