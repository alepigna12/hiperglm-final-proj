calc_loglik <- function(model, ...) UseMethod("calc_loglik")

calc_loglik.default <- function(model){
  stop("Model type not implemented yet")
}

calc_loglik.linear_model <- function(model, noise_var=1){
  design = model$design
  reg_coef = model$reg_coef
  outcome = model$outcome
  predicted_val <- design %*% reg_coef
  loglik <- - 0.5 * sum((outcome - predicted_val)^2) / noise_var
  return(loglik)
}

calc_loglik.logit_model <- function(model){
  design = model$design
  reg_coef = model$reg_coef
  outcome = model$outcome
  if (is.list(outcome)) {
    n_success <- outcome$n_success
    n_trial <- outcome$n_trial
  } else {
    n_success <- outcome
    n_trial <- rep(1, length(n_success)) # Assume binary outcome
  }
  logit_prob <- design %*% reg_coef
  loglik <- sum(n_success * logit_prob - n_trial * log(1 + exp(logit_prob)))
  # TODO: improve numerical stability for logit_prob >> 1
  return(loglik)
}

calc_loglik.poisson_model <- function(model){
  design = model$design
  reg_coef = model$reg_coef
  outcome = model$outcome
  nu <- design %*% reg_coef
  loglik <- sum(outcome * nu - exp(nu))
  return(loglik)
}

calc_grad <- function(model, ...){
  design <- model$design
  if (model$name == "linear") {
    args <- list(...)
    noise_var <- if ("noise_var" %in% names(args)) args$noise_var else 1
    reg_coef = model$reg_coef
    outcome = model$outcome
    predicted_val <- design %*% reg_coef
    grad <- t(design) %*% (outcome - predicted_val) / noise_var
    grad <- as.vector(grad)
    return(grad)
  } else {
    loglink_grad <- calc_loglink_deriv(model, order = 1)
    grad <- t(design) %*% loglink_grad
    grad <- as.vector(grad)
    return(grad)
  }
}

calc_loglink_deriv <- function(model, ...) UseMethod("calc_loglink_deriv")

calc_loglink_deriv.logit_model <- function(model, order) {
  design = model$design
  reg_coef = model$reg_coef
  outcome = model$outcome
  if (is.list(outcome)) {
    n_success <- outcome$n_success
    n_trial <- outcome$n_trial
  } else {
    n_success <- outcome
    n_trial <- rep(1, length(n_success)) # Assume binary outcome
  }
  logit_prob <- as.vector(design %*% reg_coef)
  predicted_prob <- 1 / (1 + exp(-logit_prob))
  if (order == 1) {
    deriv <- n_success - n_trial * predicted_prob
  } else if (order == 2) {
    deriv <- n_trial * predicted_prob * (1 - predicted_prob)
  } else {
    stop("3rd+ order derivative calculations are not supported")
  }
  deriv <- as.vector(deriv)
  return(deriv)
}

calc_loglink_deriv.poisson_model <- function(model, order) {
  design = model$design
  reg_coef = model$reg_coef
  outcome = model$outcome
  lambda <- as.vector(exp(design %*% reg_coef))
  if (order == 1) {
    deriv <- outcome - lambda
  } else if (order == 2) {
    deriv <- lambda
  } else {
    stop("3rd+ order derivative calculations are not supported")
  }
  deriv <- as.vector(deriv)
  return(deriv)
}

calc_logit_hessian <- function(model) {
  weight <- calc_loglink_deriv(model, order = 2)
  design <- model$design
  hess <- - t(design) %*% (outer(weight, rep(1, ncol(design))) * design)
  return(hess)
}

calc_logit_hessian_inverse <- function(model) {
  weight <- calc_loglink_deriv(model, order = 2)
  design = model$design
  sqrt_weighted_design <- outer(sqrt(weight), rep(1, ncol(design))) * design
  R <- qr_wrapper(sqrt_weighted_design)$R
  inverse <- - invert_gram_mat_from_qr(R)
  return(inverse)
}
