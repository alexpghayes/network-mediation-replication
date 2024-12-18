# eventual plot cleanup plan

# plots i want:
#
#  - side by side average coefficient error and individual absolute coefficient error
#  - coefficient convergence summaries:
#    - frobenius loss
#    - prediction error [0, 1] scale
#    - coverage
#  - Xhat to X convergence
#    - two to infinity loss (x model), frobenius loss, sin-theta loss
#    - column-wise frobenius loss of X vs Xhat (rank x model)

plot_spectral_loss <- function(loss, file_type = "png", width = NA, height = NA) {
  model <- loss$model_name[1] |> stringr::str_remove("_normal")
  model <- model <- ifelse(model == "block", "informative", model)

  sliced <- loss |>
    group_by(chunk, rep_within_chunk, n, rank) |>
    select(two_infty_loss, frob_loss, sin_theta_loss, chunk, rep_within_chunk, n, rank) |>
    slice(1) |>
    ungroup()

  num_sims_per_paramset <- nrow(distinct(sliced, chunk, rep_within_chunk))

  plot1 <- sliced |>
    pivot_longer(
      cols = c("two_infty_loss", "frob_loss", "sin_theta_loss"),
      names_to = "loss_type",
      values_to = "loss"
    ) |>
    mutate(
      loss_type_nice = case_when(
        loss_type == "two_infty_loss" ~ "Two-to-infinity norm of Xhat - X Q",
        loss_type == "frob_loss" ~ "Frobenius norm of Xhat - X Q",
        loss_type == "sin_theta_loss" ~ "Sin(Theta) distance between Xhat and X"
      )
    ) |>
    ggplot(aes(n, loss, color = as.factor(rank))) +
    geom_point() +
    scale_x_log10() +
    scale_color_viridis_d() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    facet_grid(
      rows = vars(loss_type_nice),
      cols = vars(rank),
      scales = "free",
      labeller = labeller(rank = label_both, loss_type_nice = label_value)
    ) +
    labs(
      x = "log(num nodes)",
      title = glue("Spectral distances between X and Xhat ({model})")
    ) +
    theme_minimal(base_size = 16) +
    theme(legend.position = "none")

  plot1

  if (!dir.exists(here("figures/simulations/ase"))) {
    dir.create(here("figures/simulations/ase"))
  }

  path1 <- here(
    glue("figures/simulations/ase/{model}_loss.{file_type}")
  )

  ggsave(
    path1,
    plot = plot1,
    width = width,
    height = height
  )

  plot2 <- loss |>
    ggplot(aes(n, loss, color = as.factor(rank))) +
    geom_point() +
    # geom_smooth() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    # scale_y_log10() +
    scale_color_viridis_d() +
    facet_grid(
      rows = vars(column),
      cols = vars(rank),
      labeller = "label_both",
      scales = "free_y"
    ) +
    labs(
      y = "Loss",
      x = "Number of nodes",
      title = glue("Sum square error between X and Xhat, by column ({model})")
    ) +
    theme_minimal(base_size = 16) +
    theme(legend.position = "none")

  plot2

  path2 <- here(
    glue("figures/simulations/ase/{model}_pc_loss.{file_type}")
  )

  ggsave(
    path2,
    plot = plot2,
    width = width,
    height = height
  )

  c(path1, path2)
}

plot_mediator_glance <- function(loss, file_type = "png", width = NA, height = NA) {
  model <- loss$model_name[1] |> stringr::str_remove("_normal")
  model <- ifelse(model == "block", "informative", model)

  num_sims_per_paramset <- nrow(distinct(loss, chunk, rep_within_chunk))

  # coverage is summarized across reps, but total_coef_loss and prediction error are not

  coverage <- loss |>
    group_by(n, rank) |>
    summarize(
      coverage = mean(aligned_coef_covered),
      .groups = "drop"
    )

  plot <- loss |>
    mutate_at(vars(aligned_coef_covered), as.numeric) |>
    group_by(chunk, rep_within_chunk, n, rank) |>
    select(total_coef_loss, prediction_error, chunk, rep_within_chunk, n, rank) |>
    slice(1) |>
    ungroup() |>
    left_join(coverage, by = c("n", "rank")) |>
    pivot_longer(
      cols = c("total_coef_loss", "prediction_error", "coverage"),
      names_to = "loss_type",
      values_to = "loss"
    ) |>
    mutate(
      loss_type_nice = case_when(
        loss_type == "total_coef_loss" ~ "Coefficient frobenius loss\n(after procrustes alignment)",
        loss_type == "prediction_error" ~ "Average variance explained (R^2)\nacross all outcomes",
        loss_type == "coverage" ~ "Confidence interval coverage\n(after procrustes alignment)"
      )
    ) |>
    ggplot(aes(n, loss, color = as.factor(rank))) +
    geom_point() +
    scale_color_viridis_d() +
    scale_x_log10() +
    facet_grid(
      rows = vars(loss_type_nice),
      cols = vars(rank),
      scales = "free",
      labeller = labeller(rank = label_both, loss_type_nice = label_value)
    ) +
    labs(
      y = "",
      x = "Number of nodes (log scale)",
      title = glue("Convergence of mediator regression ({model})")
    ) +
    theme_minimal(base_size = 16) +
    theme(legend.position = "none")

  plot

  if (!dir.exists(here("figures/simulations/mediator/"))) {
    dir.create(here("figures/simulations/mediator/"))
  }

  path <- here(
    glue("figures/simulations/mediator/{model}_glance.{file_type}")
  )

  ggsave(
    path,
    plot = plot,
    width = width,
    height = height
  )

  path
}

plot_mediator_elementwise_loss <- function(loss, file_type = "png", width = NA, height = NA) {
  model <- loss$model_name[1] |> stringr::str_remove("_normal")
  model <- ifelse(model == "block", "informative", model)

  num_sims_per_paramset <- nrow(distinct(loss, chunk, rep_within_chunk))

  plot <- loss |>
    mutate(
      error = estimate - aligned_true_coef
    ) |>
    group_by(n, rank, term, outcome) |>
    summarize(
      avg_l2_error = mean(error^2)
    ) |>
    ungroup() |>
    ggplot(aes(n, avg_l2_error, color = outcome)) +
    # geom_smooth(se = FALSE) +
    geom_point() +
    geom_line() +
    scale_x_log10(labels = scales::label_log(digits = 2)) +
    scale_y_log10(labels = scales::label_log(digits = 2)) +
    scale_color_brewer(palette = "Dark2") +
    facet_grid(
      cols = vars(rank),
      rows = vars(term),
      scales = "free_y",
      labeller = labeller(.rows = label_value, .cols = label_both)
    ) +
    labs(
      y = "Mean squared error\n(log scale, after procrustes alignment)",
      x = "Number of nodes (log scale)",
      color = "Dim"
    ) +
    theme_minimal(base_size = 16)

  if (!dir.exists(here("figures/simulations/mediator/"))) {
    dir.create(here("figures/simulations/mediator/"))
  }

  path <- here(
    glue("figures/simulations/mediator/{model}_theta_loss.{file_type}")
  )

  ggsave(
    path,
    plot = plot,
    width = width,
    height = 8
  )

  path
}


plot_outcome_glance <- function(loss, file_type = "png", width = NA, height = NA) {
  model <- loss$model_name[1] |> stringr::str_remove("_normal")
  model <- ifelse(model == "block", "informative", model)

  num_sims_per_paramset <- nrow(distinct(loss, chunk, rep_within_chunk))

  # coverage is summarized across reps, but total_coef_loss and prediction error are not

  coverage <- loss |>
    group_by(n, rank) |>
    summarize(
      coverage = mean(aligned_coef_covered),
      .groups = "drop"
    )

  plot <- loss |>
    mutate_at(vars(aligned_coef_covered), as.numeric) |>
    group_by(rep_within_chunk, chunk, n, rank) |>
    select(total_coef_loss, prediction_error, kappa, rep_within_chunk, chunk, n, rank) |>
    slice(1) |>
    ungroup() |>
    left_join(coverage, by = c("n", "rank")) |>
    pivot_longer(
      cols = c("total_coef_loss", "prediction_error", "coverage", "kappa"),
      names_to = "loss_type",
      values_to = "loss"
    ) |>
    mutate(
      loss_type_nice = case_when(
        loss_type == "total_coef_loss" ~ "Coefficient frobenius loss\n(after procrustes alignment)",
        loss_type == "prediction_error" ~ "Variance explained",
        loss_type == "coverage" ~ "Confidence interval coverage\n(after procrustes alignment)",
        loss_type == "kappa" ~ "Condition number of vcov(betahat)"
      )
    ) |>
    ggplot(aes(n, loss, color = as.factor(rank))) +
    geom_point() +
    scale_x_log10() +
    scale_y_log10() +
    scale_color_brewer(palette = "Dark2") +
    facet_grid(
      rows = vars(loss_type_nice),
      cols = vars(rank),
      scales = "free",
      labeller = labeller(rank = label_both, loss_type_nice = label_value)
    ) +
    labs(
      y = "",
      x = "Number of nodes (log scale)",
      title = glue("Convergence of outcome regression ({model})")
    ) +
    theme_minimal(base_size = 16) +
    theme(legend.position = "none")

  plot

  if (!dir.exists(here("figures/simulations/outcome/"))) {
    dir.create(here("figures/simulations/outcome/"))
  }

  path <- here(
    glue("figures/simulations/outcome/{model}_glance.{file_type}")
  )

  ggsave(
    path,
    plot = plot,
    width = width,
    height = 8
  )

  path
}

plot_outcome_elementwise_loss <- function(loss, file_type = "png", width = NA, height = NA) {
  model <- loss$model_name[1] |> stringr::str_remove("_normal")
  model <- ifelse(model == "block", "informative", model)

  num_sims_per_paramset <- nrow(distinct(loss, chunk, rep_within_chunk))

  plot <- loss |>
    mutate(
      error = estimate - aligned_true_coef,
    ) |>
    group_by(n, rank, term, outcome) |>
    summarize(
      avg_l2_error = mean(error^2)
    ) |>
    mutate(
      partition = ifelse(stringr::str_detect(term, "US"), "beta[x]", "beta[w]"),
    ) |>
    group_by(n, rank, outcome, partition) |>
    mutate(
      index = as.factor(row_number())
    ) |>
    ungroup() |>
    ggplot(aes(n, avg_l2_error, color = index)) +
    # geom_smooth(se = FALSE) +
    geom_point() +
    geom_line() +
    scale_x_log10(labels = scales::label_log(digits = 2)) +
    scale_y_log10(labels = scales::label_log(digits = 2)) +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    facet_grid(
      cols = vars(rank),
      rows = vars(partition),
      scales = "free_y",
      labeller = labeller(.rows = label_parsed, .cols = label_both)
    ) +
    labs(
      y = "Mean squared error\n(log scale, after procrustes alignment)",
      x = "Number of nodes (log scale)"
    ) +
    theme_minimal(base_size = 16)

  if (!dir.exists(here("figures/simulations/outcome/"))) {
    dir.create(here("figures/simulations/outcome/"))
  }

  path <- here(
    glue("figures/simulations/outcome/{model}_beta_loss_average.{file_type}")
  )

  ggsave(
    path,
    plot = plot,
    width = width,
    height = height
  )

  plot2 <- loss |>
    mutate(
      error = (estimate - aligned_true_coef)^2,
      partition = ifelse(stringr::str_detect(term, "US"), "beta[x]", "beta[w]")
    ) |>
    ggplot(aes(n, error, color = term)) +
    # geom_smooth(se = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point() +
    scale_x_log10(labels = scales::label_log(digits = 2)) +
    scale_y_log10(labels = scales::label_log(digits = 2)) +
    scale_color_brewer(palette = "Dark2") +
    facet_grid(
      cols = vars(rank),
      rows = vars(partition),
      scales = "free_y",
      labeller = labeller(.rows = label_parsed, .cols = label_both)
    ) +
    labs(
      y = "Squared error\n(log scale, after procrustes alignment)",
      x = "Number of nodes (log scale)",
      title = glue("Outcome coefficient error"),
      subtitle = glue("{model} model")
    ) +
    theme_minimal(base_size = 16)

  path2 <- here(
    glue("figures/simulations/outcome/{model}_beta_loss.{file_type}")
  )

  ggsave(
    path2,
    plot = plot2,
    width = width,
    height = height
  )

  plot3 <- loss |>
    mutate(
      error = estimate - aligned_true_coef
    ) |>
    group_by(n, rank, term, outcome) |>
    summarize(
      bias = mean(abs(error))
    ) |>
    mutate(
      partition = ifelse(stringr::str_detect(term, "US"), "beta[x]", "beta[w]")
    ) |>
    group_by(n, rank, outcome, partition) |>
    mutate(
      index = as.factor(row_number())
    ) |>
    ungroup() |>
    ggplot(aes(n, bias, color = index)) +
    geom_line() +
    # geom_smooth(se = FALSE) +
    geom_point() +
    scale_x_log10(labels = scales::label_log(digits = 2)) +
    scale_y_log10(labels = scales::label_log(digits = 2)) +
    scale_color_brewer(palette = "Dark2", guide = FALSE) +
    facet_grid(
      cols = vars(rank),
      rows = vars(partition),
      scales = "free_y",
      labeller = labeller(.rows = label_parsed, .cols = label_both)
    ) +
    labs(
      y = "Coefficient-wise absolute bias\n(log scale, after procrustes alignment)",
      x = "Number of nodes (log scale)"
    ) +
    theme_minimal(base_size = 16)

  path3 <- here(
    glue("figures/simulations/outcome/{model}_beta_bias.{file_type}")
  )

  ggsave(
    path3,
    plot = plot3,
    width = width,
    height = height
  )

  c(path, path2, path3)
}

plot_causal_loss <- function(loss, file_type = "png", width = 9, height = 9 * 9 / 16) {
  plot <- loss |>
    filter(model_name %in% c("block_normal", "uninformative_normal")) |>
    mutate(
      model = stringr::str_remove(model_name, "_normal"),
      # slightly fragile
      model = ifelse(model == "block", "informative", "uninformative")
    ) |>
    pivot_longer(
      cols = c("nie", "nde"),
      names_to = "estimand",
      values_to = "loss"
    ) |>
    group_by(n, rank, estimand, model) |>
    summarize(
      mse = mean(loss),
      sd = sd(loss),
      dist = dist_truncated(dist_normal(mse, sd = sd), lower = 0)
    ) |>
    ungroup() |>
    mutate(
      estimand = toupper(estimand)
    ) |>
    ggplot(aes(x = n, y = mse, ymin = mse - sd, ymax = mse + sd, color = estimand, ydist = dist)) +
    # geom_pointrange() +
    # geom_ribbon() +
    # stat_pointinterval(
    #   .width = c(0.5, 0.80, 0.95)
    # ) +
    geom_point() +
    geom_line() +
    scale_x_log10(labels = scales::label_log(digits = 2)) +
    scale_y_log10(labels = scales::label_log(digits = 2)) +
    scale_color_brewer(palette = "Dark2") +
    facet_grid(
      cols = vars(rank),
      rows = vars(model),
      scales = "free_y",
      labeller = labeller(
        .rows = label_value,
        .cols = label_both
      )
    ) +
    labs(
      y = "Mean squared error (log scale)",
      x = "Number of nodes (log scale)"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      legend.title = element_blank()
    )

  plot

  if (!dir.exists(here("figures/simulations/causal/"))) {
    dir.create(here("figures/simulations/causal/"))
  }

  path <- here(
    glue("figures/simulations/causal/loss_average.{file_type}")
  )

  ggsave(
    path,
    plot = plot,
    width = width,
    height = height
  )

  plot2 <- loss |>
    filter(model_name %in% c("block_normal", "uninformative_normal")) |>
    mutate(
      model = stringr::str_remove(model_name, "_normal"),
      # slightly fragile
      model = ifelse(model == "block", "informative", "uninformative")
    ) |>
    pivot_longer(
      cols = c("nie_covered", "nde_covered"),
      names_to = "estimand",
      values_to = "covered"
    ) |>
    mutate(
      estimand = stringr::str_remove(estimand, "_covered"),
      estimand = toupper(estimand)
    ) |>
    group_by(n, rank, estimand, model) |>
    summarize(
      coverage = mean(covered),
      conf.low = mean(covered) - 1.96 * sd(covered) / n(),
      conf.high = mean(covered) + 1.96 * sd(covered) / n()
    ) |>
    ungroup() |>
    ggplot() +
    aes(
      x = n,
      ymin = conf.low,
      y = coverage,
      ymax = conf.high,
      color = estimand,
      fill = estimand
    ) +
    geom_point() +
    # geom_ribbon(alpha = 0.5) +
    geom_line() +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    expand_limits(x = 0, y = 0.7) +
    scale_x_log10() +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    facet_grid(
      cols = vars(rank),
      rows = vars(model),
      scales = "free_y",
      labeller = labeller(
        .rows = label_value,
        .cols = label_both
      )
    ) +
    labs(
      y = "Coverage rate",
      x = "Number of nodes (log scale)"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      legend.title = element_blank()
    )

  path2 <- here(
    glue("figures/simulations/causal/coverage.{file_type}")
  )

  ggsave(
    path2,
    plot = plot2,
    width = width,
    height = height
  )

  c(path, path2)
}


plot_null_causal_loss <- function(loss, file_type = "png", width = 9, height = 9 * 9 / 16) {
  plot <- loss |>
    filter(model_name %in% c("block_nullnde", "block_nullnie")) |>
    mutate(
      # slightly fragile
      model = ifelse(model_name == "block_nullnde", "null nde", "null nie")
    ) |>
    pivot_longer(
      cols = c("nie", "nde"),
      names_to = "estimand",
      values_to = "loss"
    ) |>
    group_by(n, rank, estimand, model) |>
    summarize(
      mse = mean(loss)
    ) |>
    ungroup() |>
    mutate(
      estimand = toupper(estimand)
    ) |>
    ggplot(aes(n, mse, color = estimand)) +
    geom_point() +
    geom_line() +
    scale_x_log10(labels = scales::label_log(digits = 2)) +
    scale_y_log10(labels = scales::label_log(digits = 2)) +
    scale_color_brewer(palette = "Dark2") +
    facet_grid(
      cols = vars(rank),
      rows = vars(model),
      scales = "free_y",
      labeller = labeller(
        .rows = label_value,
        .cols = label_both
      )
    ) +
    labs(
      y = "Mean squared error (log scale)",
      x = "Number of nodes (log scale)"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      legend.title = element_blank()
    )

  if (!dir.exists(here("figures/simulations/causal/"))) {
    dir.create(here("figures/simulations/causal/"))
  }

  path <- here(
    glue("figures/simulations/causal/loss_average_null.{file_type}")
  )

  ggsave(
    path,
    plot = plot,
    width = width,
    height = height
  )

  plot2 <- loss |>
    filter(model_name %in% c("block_nullnde", "block_nullnie")) |>
    mutate(
      # slightly fragile
      model = ifelse(model_name == "block_nullnde", "null nde", "null nie")
    ) |>
    pivot_longer(
      cols = c("nie_covered", "nde_covered"),
      names_to = "estimand",
      values_to = "covered"
    ) |>
    mutate(
      estimand = stringr::str_remove(estimand, "_covered"),
      estimand = toupper(estimand)
    ) |>
    group_by(n, rank, estimand, model) |>
    summarize(
      coverage = mean(covered),
      conf.low = mean(covered) - 1.96 * sd(covered) / n(),
      conf.high = mean(covered) + 1.96 * sd(covered) / n()
    ) |>
    ungroup() |>
    ggplot() +
    aes(
      x = n,
      ymin = conf.low,
      y = coverage,
      ymax = conf.high,
      color = estimand,
      fill = estimand
    ) +
    geom_point() +
    # geom_ribbon(alpha = 0.5) +
    geom_line() +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    expand_limits(x = 0, y = 0.7) +
    scale_x_log10() +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    facet_grid(
      cols = vars(rank),
      rows = vars(model),
      scales = "free_y",
      labeller = labeller(
        .rows = label_value,
        .cols = label_both
      )
    ) +
    labs(
      y = "Coverage rate",
      x = "Number of nodes (log scale)"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      legend.title = element_blank()
    )

  path2 <- here(
    glue("figures/simulations/causal/coverage_null.{file_type}")
  )

  ggsave(
    path2,
    plot = plot2,
    width = width,
    height = height
  )

  c(path, path2)
}