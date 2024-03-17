test_that("`hglm` linear model est coincides with `lm` one", {
  data <- simulate_data(32, 4, intercept = 1, seed = 1918)
  design <- data$design; outcome <- data$outcome
  lm_out <- stats::lm(outcome ~ design + 0)
  hglm_out <- hiper_glm(design, outcome)
  expect_true(are_all_close(coef(hglm_out), coef(lm_out)))
  expect_true(are_all_close(  
    as.vector(vcov(hglm_out)), as.vector(vcov(lm_out))
  ))
})

test_that("`hglm` logit model est coincides with `glm` one for binary outcome", {
  data <- simulate_data(32, 4, model_name = 'logit', seed = 1918)
  design <- data$design; outcome <- data$outcome
  glm_out <- stats::glm(outcome ~ design + 0, family = binomial('logit'))
  hglm_out <- hiper_glm(design, outcome, model_name = 'logit')
  expect_true(are_all_close(coef(hglm_out), coef(glm_out)))
  expect_true(are_all_close(
    as.vector(vcov(hglm_out)), as.vector(vcov(glm_out)), 
    abs_tol = Inf, rel_tol = 1e-3
  ))
})

test_that("`hglm` logit model est coincides with `glm` one for binomial outcome", {
  set.seed(615)
  n_obs <- 32; n_pred <- 4
  n_trial <- 1L + rpois(n_obs, lambda = 1)
  data <- simulate_data(
    n_obs, n_pred, model_name = 'logit', seed = 1918, 
    option = list(n_trial = n_trial)
  )
  design <- data$design; outcome <- data$outcome
  n_success <- outcome$n_success; n_trial <- outcome$n_trial
  glm_out <- stats::glm(
    cbind(n_success, n_trial - n_success) ~ design + 0, 
    family = binomial('logit')
  )
  hglm_out <- hiper_glm(design, outcome, model_name = 'logit')
  expect_true(are_all_close(coef(hglm_out), coef(glm_out)))
  expect_true(are_all_close(
    as.vector(vcov(hglm_out)), as.vector(vcov(glm_out)), 
    abs_tol = Inf, rel_tol = 1e-3
  ))
})

test_that("`hglm` poisson model est coincides with `glm` one", {
  data <- simulate_data(32, 4, model_name = 'poisson', seed = 1918)
  design <- data$design; outcome <- data$outcome
  glm_out <- stats::glm(outcome ~ design + 0, family = poisson)
  hglm_out <- hiper_glm(design, outcome, model_name = 'poisson')
  expect_true(are_all_close(coef(hglm_out), coef(glm_out)))
  expect_true(are_all_close(
    as.vector(vcov(hglm_out)), as.vector(vcov(glm_out)), 
    abs_tol = Inf, rel_tol = 1e-3
  ))
})