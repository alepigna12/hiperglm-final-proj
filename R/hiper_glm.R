#' @export
hiper_glm <- function(design, outcome, model_name = "linear", option = list()) {
  supported_model <- c("linear", "logit")
  if (!(model_name %in% supported_model)) {
    stop(sprintf("The model %s is not supported.", model_name))
  }
  model <- new_regression_model(design, outcome, model_name)
  hglm_out <- find_mle(model, option)
  class(hglm_out) <- "hglm"
  return(hglm_out)
}

new_regression_model <- function(design, outcome, model_name) {
  model <- list(design = design, outcome = outcome, name = model_name)
  class(model) <- paste(model_name, "model", sep = "_")
  return(model)
}

find_mle <- function(model, option) {
  if (is.null(option$mle_solver)) {
    if (model$name == 'linear') {
      result <- solve_via_least_sq(model$design, model$outcome)
    } else {
      result <- solve_via_newton(
        model$design, model$outcome, option$n_max_iter, option$rel_tol, option$abs_tol
      )
    }
  } else {
    result <- solve_via_optim(model$design, model$outcome, model$name, option$mle_solver)
  }
  return(result)
}

solve_via_least_sq <- function(design, outcome) {
  ls_result <- solve_least_sq_via_qr(design, outcome)
  mle_coef <- ls_result$solution
  noise_var <- mean((outcome - design %*% mle_coef)^2)
  n_obs <- nrow(design); n_pred <- ncol(design)
  noise_var <- noise_var / (1 - n_pred / n_obs) 
    # Use the same nearly-unbiased estimator as in `stats::lm`
  cov_est <- noise_var * invert_gram_mat_from_qr(ls_result$R)
  return(list(coef = mle_coef, cov_est = cov_est))
}

solve_via_newton <- function(design, outcome, n_max_iter, rel_tol, abs_tol) {
  if (is.null(n_max_iter)) { n_max_iter <- 25L }
  if (is.null(rel_tol)) { rel_tol <- 1e-6 }
  if (is.null(abs_tol)) { abs_tol <- 1e-6 }
  coef_est <- rep(0, ncol(design))
  n_iter <- 0L
  max_iter_reached <- FALSE
  converged <- FALSE
  curr_loglik <- calc_logit_loglik(coef_est, design, outcome)
  while (!(converged || max_iter_reached)) {
    prev_loglik <- curr_loglik
    coef_est <- take_one_newton_step(coef_est, design, outcome)
    curr_loglik <- calc_logit_loglik(coef_est, design, outcome)
    converged <- (
      2 * abs(curr_loglik - prev_loglik) < (abs_tol + rel_tol * abs(curr_loglik))
    )
    n_iter <- n_iter + 1L
    max_iter_reached <- (n_iter == n_max_iter)
  }
  if (max_iter_reached && !converged) {
    warning("Newton's method did not converge. The estimates may be meaningless.")
  }
  cov_est <- - calc_logit_hessian_inverse(coef_est, design, outcome)
  return(list(
    coef = coef_est, cov = cov_est, 
    converged = converged, n_iter = n_iter
  ))
}

take_one_newton_step <- function(
    coef_est, design, outcome, solver = "weighted-leqst-sq"
) {
  if (solver == "weighted-leqst-sq") {
    loglink_grad <- 
      calc_logit_loglink_deriv(coef_est, design, outcome, order = 1)
    weight <- calc_logit_loglink_deriv(coef_est, design, outcome, order = 2)
    if (any(weight == 0)) {
      stop("Exact 0 or 1 found in predicted probability while solving for MLE.")
        # TODO: pursue alternative path forward in this case. Maybe just fall 
        # back on a Newton step with explicit computation of weighted Hessian. 
    }
    ls_target_vec <- loglink_grad / weight
    coef_update <- solve_least_sq_via_qr(design, ls_target_vec, weight)$solution
  } else {
    grad <- calc_logit_grad(coef_est, design, outcome)
    hess <- calc_logit_hessian(coef_est, design, outcome)
    coef_update <- - solve(hess, grad)
  }
  coef_est <- coef_est + coef_update
  return(coef_est = coef_est)
}

solve_via_optim <- function(design, outcome, model_name, method) {
  init_coef <- rep(0, ncol(design))
  if (model_name == 'linear') {
    obj_fn <- function (coef) {
      calc_linear_loglik(coef, design, outcome) 
    }
    obj_grad <- function (coef) {
      calc_linear_grad(coef, design, outcome)
    }
  } else {
    obj_fn <- function (coef) {
      calc_logit_loglik(coef, design, outcome) 
    }
    obj_grad <- function (coef) {
      calc_logit_grad(coef, design, outcome)
    }
  }
  optim_result <- stats::optim(
    init_coef, obj_fn, obj_grad, method = method,
    control = list(fnscale = -1) # Maximize the function
  )
  optim_converged <- (optim_result$convergence == 0L)
  if (!optim_converged) {
    warning("Optimization did not converge. The estimates may be meaningless.")
  }
  return(list(coef = optim_result$par))
}
