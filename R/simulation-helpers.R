model_nullnde <- function(n, k, ...) {
  # exactly the block2 model
  m <- model_mediator_block2(n = n, k = k, ...)
  dim_w <- ncol(m$W)
  o <- model_outcome_t(m, beta_w = rep(0, dim_w))
  o$nde <- o$beta_w["trt"]
  o$nie <- drop(o$mediator$Theta %*% o$beta_x)["trt"]
  o$model_name <- "block_nullnde"
  o
}

model_nullnie <- function(n, k, ...) {
  # exactly the block2 model
  m <- model_mediator_block2(n = n, k = k, ...)
  dim_x <- ncol(m$X)
  o <- model_outcome_t(m, beta_x = rep(0, dim_x))
  o$nde <- o$beta_w["trt"]
  o$nie <- drop(o$mediator$Theta %*% o$beta_x)["trt"]
  o$model_name <- "block_nullnie"
  o
}

embed_adjacency_matrix <- function(graph, rank) {
  A <- as_adjacency_matrix(graph)
  US(A, rank)
}

logspace <- function(min, max, num, base = 10) {
  round(base^seq(log(min, base = base), log(max, base = base), length.out = num))
}

subspace_loss <- function(u, v) {
  # see [1] Vu and Lei 2013 section 2.3
  # and [2] Rohe, Chatterjee, Yu 2011 Annals of Statistics page 1908

  u <- normalize.cols(u, tol = 1e-10)
  v <- normalize.cols(v, tol = 1e-10)

  s <- svd(crossprod(u, v))
  ncol(u) - sum(s$d^2)
}

# rotate Y to align it to X
mypro <- function(X, Y) {
  XY <- crossprod(X, Y)
  s <- svd(XY)
  rotation <- s$v %*% t(s$u)
  Yrot <- Y %*% rotation
  diff <- X - Yrot

  colwise_loss <- Matrix::colSums(diff^2)
  two_infty_loss <- max(sqrt(Matrix::rowSums(diff^2)))

  list(
    rotation = rotation,
    Yrot = Yrot,
    colwise_loss = colwise_loss,
    two_infty_loss = two_infty_loss,
    sum_squared_error = sum(colwise_loss)
  )
}

spectral_loss <- function(model, xhat) {

  xhat_pop <- ASE(model$mediator$A_model)
  sin_theta_loss <- subspace_loss(xhat, xhat_pop)

  pro <- mypro(xhat, xhat_pop)

  list(
    procrustes_rotation = pro$rotation,
    model = model,
    sin_theta_loss = sin_theta_loss,
    colwise_loss = pro$colwise_loss,
    two_infty_loss = pro$two_infty_loss,
    model_name = model$model_name
  )
}

u_loss_helper <- function(rotated_model, params, chunk) {

  enframe(rotated_model$colwise_loss, name = "column", value = "loss") |>
    mutate(
      frob_loss = sum(rotated_model$colwise_loss),
      two_infty_loss = rotated_model$two_infty_loss,
      sin_theta_loss = rotated_model$sin_theta_loss,
      n = params$n,
      rank = params$rank,
      model_name = rotated_model$model_name,
      chunk = chunk
    )
}

fit_models <- function(node_data, US) {

  # NOTE: do not pass graph to nodelm_robust() to separate construct the
  # m_fit and the o_fit here because US(A, k) used in that context did
  # not given the same singular vectors (in particular, I somehow ended up
  # with some sort of sign flipping). this results in a different Q for the
  # m_fit and the o_fit and also possibly explains issues with finding
  # a rotation using X and Xhat

  m_fit <- lm_robust(US ~ . - name - y - 1, data = node_data)
  m_fit$terms <- NULL

  o_fit <- lm_robust(y ~ . - name - 1 + US, data = node_data)
  o_fit$terms <- NULL

  list(
    m_fit = m_fit,
    o_fit = o_fit
  )
}

# rotated model is not actually rotated, just includes information about
# a rotation between population US and sample US
mediator_loss <- function(estimates, rotated_model) {

  model <- rotated_model$model$mediator
  prediction_error <- mean(estimates$r.squared)

  coef_estimates <- coef(estimates)

  coef_model <- coef(model)
  rownames(coef_model) <- rownames(coef_estimates)

  procrustes_coefs <- mypro(coef_estimates, coef_model)

  # procrustes learned on coefficients
  # coef_rot <- procrustes_coefs$Yrot

  # procrustes learned on ASE
  coef_rot <- coef_model %*% rotated_model$procrustes_rotation

  rownames(coef_rot) <- rownames(coef_estimates)
  colnames(coef_rot) <- 1:model$k

  tidy_estimates <- tidy(estimates, conf.int = TRUE)

  coverage_table <- coef_rot |>
    as_tibble(rownames = "term") |>
    pivot_longer(
      -term,
      names_to = "outcome",
      values_to = "aligned_true_coef"
    ) |>
    left_join(tidy_estimates, by = c("term", "outcome")) |>
    mutate(
      aligned_coef_covered = conf.low <= aligned_true_coef & aligned_true_coef <= conf.high
    )

  list(
    frob_coef_loss = procrustes_coefs$sum_squared_error,
    coverage_table = coverage_table,
    prediction_error = prediction_error,
    model_name = rotated_model$model_name
  )
}

mediator_loss_helper <- function(loss, params, chunk) {
  loss$coverage_table |>
    mutate(
      total_coef_loss = loss$frob_coef_loss,
      prediction_error = loss$prediction_error,
      n = params$n,
      rank = params$rank,
      model_name = loss$model_name,
      chunk = chunk
    )
}

clean_term_names <- function(term_names) {
  str_replace(term_names, "\\(A, [0-9]+\\)", "")
}

outcome_loss <- function(estimates, rotated) {

  prediction_error <- estimates$r.squared

  model <- rotated$model

  coef_estimates <- coef(estimates)
  names(coef_estimates) <- clean_term_names(names(coef_estimates))

  ase_coef_names <- names(coef_estimates)[
    str_detect(names(coef_estimates), "US[0-9]+")
  ]

  coef_model <- coef(model)

  # important that these are row vectors here even though we typically
  # think of them as column vectors

  ase_coef_model <- coef_model[ase_coef_names]
  # ase_coef_model <- matrix(coef_model[ase_coef_names], nrow = 1)
  # ase_coef_estimates <- matrix(coef_estimates[ase_coef_names], nrow = 1)

  # ase_coef_pro <- mypro(ase_coef_estimates, ase_coef_model)
  # ase_coef_pro$Yrot

  # print(ase_coef_model)
  # print(dim(rotated$procrustes_rotation))

  # procrustes_ase <- rotated$procrustes_ase
  coef_rot <-  t(rotated$procrustes_rotation) %*% ase_coef_model

  # ase_coef_estimates
  # coef_rot
  coef_model_aligned <- coef_model
  coef_model_aligned[ase_coef_names] <- drop(coef_rot) # as.numeric(ase_coef_pro$Yrot)

  coef_model_df <- enframe(
    coef_model_aligned,
    name = "term",
    value = "aligned_true_coef"
  )

  coverage_table <- tidy(estimates, conf.int = TRUE) |>
    mutate(
      term = clean_term_names(term)
    ) |>
    left_join(coef_model_df, by = "term") |>
    mutate(
      aligned_coef_covered = conf.low <= aligned_true_coef & aligned_true_coef <= conf.high
    )

  list(
    frob_coef_loss = sum((coef_estimates - coef_model)^2),
    coverage_table = coverage_table,
    prediction_error = prediction_error,
    model_name = model$model_name,
    kappa = kappa(estimates$vcov)
  )
}

outcome_loss_helper <- function(loss, params, chunk) {
  loss$coverage_table |>
    mutate(
      total_coef_loss = loss$frob_coef_loss,
      prediction_error = loss$prediction_error,
      kappa = loss$kappa,
      n = params$n,
      rank = params$rank,
      model_name = loss$model_name,
      chunk = chunk
    )
}

causal_loss <- function(estimates, model) {

  o_fit <- estimates$o_fit
  m_fit <- estimates$m_fit

  # may need modification if you use simulation models
  # other than the uninformative/block2 classes

  num_coefs <- nrow(coef(m_fit))

  nde_table <- tidy(o_fit, conf.int = TRUE)[1:num_coefs, ] |>
    select(term, estimate, conf.low, conf.high) |>
    mutate(
      covered = conf.low <= model$nde & model$nde <= conf.high
    )

  betaw_hat <- coef(o_fit)[1:num_coefs]
  betax_hat <- coef(o_fit)[-c(1:num_coefs)]
  theta_hat <- coef(m_fit)

  nie_hat <- drop(theta_hat %*% betax_hat)

  nie_table <- enframe(nie_hat, name = "term", value = "estimate") |>
    select(term, estimate)

  sigmabetax_hat <- vcov(o_fit)[-c(1:num_coefs), -c(1:num_coefs)]
  sigmatheta_hat <- vcov(m_fit)

  # need to re-arrange sigmatheta_hat from enormous square into something
  # more tensor-y / considering each covariate one at a time

  coef_names <- names(betaw_hat)

  # everything following is under the assumption that Theta_tc = 0
  # for convenience since it's a pain to handle Theta_tc != 0



  for (i in seq_along(coef_names)) {

    nm <- coef_names[i]

    indices <- which(
      str_detect(
        colnames(sigmatheta_hat),
        coll(nm)
      )
    )

    thetat_hat <- theta_hat[i, ]
    sigmathetat_hat <- sigmatheta_hat[indices, indices, drop = FALSE]


    # delta method. in the misspecification simulations, this fails for some
    # reason. currently just ignore these failures.

    nie_var <- tryCatch({
      t(betax_hat) %*% sigmathetat_hat %*% betax_hat +
        t(thetat_hat) %*% sigmabetax_hat %*% thetat_hat
    },
    error = function(e) {
      0
    })

    nie_table[i, "conf.low"] <- nie_hat[i] - 1.96 * sqrt(nie_var)
    nie_table[i, "conf.high"] <- nie_hat[i] + 1.96 * sqrt(nie_var)
  }

  nie_table <- nie_table |>
    mutate(
      covered = conf.low <= model$nie & model$nie <= conf.high
    )

  nde_hat <- nde_table |>
    filter(term == "trt") |>
    pull(estimate)

  nie_hat <- nie_table |>
    filter(term == "trt") |>
    pull(estimate)

  nde_covered <- nde_table |>
    filter(term == "trt") |>
    pull(covered)

  nie_covered <- nie_table |>
    filter(term == "trt") |>
    pull(covered)

  list(
    nde = (nde_hat - model$nde)^2,
    nie = (nie_hat - model$nie)^2,
    nde_covered = nde_covered,
    nie_covered = nie_covered,
    model_name = model$model_name
  )
}

causal_loss_helper <- function(loss, params, chunk) {
  as_tibble(loss) |>
    mutate(
      n = params$n,
      rank = params$rank,
      chunk = chunk
    )
}

