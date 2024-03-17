get_default_and_bfgs_optimizer_outputs <- function(
  model_name, n_obs = 32, n_pred = 4, data_seed = 1918
) {
  data <- simulate_data(n_obs, n_pred, model_name, seed = data_seed)
  design <- data$design; outcome <- data$outcome
  via_default_out <- hiper_glm(design, outcome, model_name)
  via_bfgs_out <- hiper_glm(
    design, outcome, model_name, option = list(mle_solver = 'BFGS')
  )
  return(list(via_default = via_default_out, via_bfgs = via_bfgs_out))
}

test_that("linalg and optim least-sq coincide", {
  out <- get_default_and_bfgs_optimizer_outputs("linear")
  expect_true(are_all_close(
    coef(out$via_default), coef(out$via_bfgs), abs_tol = 1e-2, rel_tol = 1e-2
  ))
})

test_that("newton and bfgs outputs coincide on logit model", {
  out <- get_default_and_bfgs_optimizer_outputs("logit")
  expect_true(are_all_close(
    coef(out$via_default), coef(out$via_bfgs), abs_tol = 1e-2, rel_tol = 1e-2
  ))
})

test_that("newton and bfgs outputs coincide on poisson model", {
  out <- get_default_and_bfgs_optimizer_outputs("poisson")
  expect_true(are_all_close(
    coef(out$via_default), coef(out$via_bfgs), abs_tol = 1e-2, rel_tol = 1e-2
  ))
})

test_that("vanilla/weighted least-sq Newton updates coincide", {
  n_obs <- 32; n_pred <- 4
  data <- simulate_data(n_obs, n_pred, model_name = 'logit', seed = 1918)
  design <- data$design; outcome <- data$outcome
  set.seed(615)
  init_coef <- rnorm(n_pred)
  wls_updated_coef <- 
    take_one_newton_step(init_coef, design, outcome, solver = "weighted-leqst-sq")
  ne_updated_coef <- 
    take_one_newton_step(init_coef, design, outcome, solver = "normal-eq")
  expect_true(are_all_close(wls_updated_coef, ne_updated_coef))
})

test_that("direct/via-QR inversion of the Gram matrix coincide", {
  set.seed(1918)
  n_row <- 32; n_col <- 4
  X <- matrix(rnorm(n_row * n_col), nrow = n_row, ncol = n_col)
  direct_inverse <- solve(t(X) %*% X)
  R <- qr_wrapper(X)$R
  qr_inverse <- invert_gram_mat_from_qr(R)
  expect_true(are_all_close(
    as.vector(direct_inverse), as.vector(qr_inverse)
  ))
})

test_that("least square via built-in and Eigen QR coincide", {
  set.seed(1918)
  n_row <- 32; n_col <- 4
  X <- matrix(rnorm(n_row * n_col), nrow = n_row, ncol = n_col)
  y <- rnorm(n_row)
  eigen_qr <- solve_least_sq_via_qr_cpp_eig(X, y)
  eigen_sol <- eigen_qr$solution
  eigen_gram_mat <- invert_gram_mat_from_qr(eigen_qr$R)
  lapack_qr <- solve_least_sq_via_qr(X, y)
  lapack_sol <- lapack_qr$solution
  lapack_gram_mat <- invert_gram_mat_from_qr(lapack_qr$R)
  expect_true(are_all_close(eigen_sol, lapack_sol))
  expect_true(are_all_close(
    as.vector(eigen_gram_mat), as.vector(lapack_gram_mat)
  ))
})