misspecification_loss <- function(est, pop) {
  nde_loss_table <- est |>
    dplyr::filter(term == "trt", estimand == "nde") |>
    dplyr::mutate(
      truth = pop$nde
    )

  nie_loss_table <- est |>
    dplyr::filter(term == "trt", estimand == "nie") |>
    dplyr::mutate(
      truth = pop$nie
    )

  bind_rows(nde_loss_table, nie_loss_table) |>
    dplyr::mutate(
      true_rank = pop$mediator$k,
      n = pop$mediator$n,
      squared_error = (estimate - truth)^2,
      covered = conf.low <= truth & truth <= conf.high,
      model_name = pop$model_name
    )
}

plot_misspecification_loss <- function(losses, num_sims_per_paramset, file_type = "png", ...) {
  # https://twitter.com/rdataberlin/status/1329045024357703681
  str_rank_natural <- function(x) {
    match(x, stringr::str_sort(x, numeric = TRUE))
  }

  summarized_losses <- losses |>
    mutate(
      model = stringr::str_remove(model_name, "_normal"),
      model = ifelse(model == "block", "informative", model),
      modeln = glue("{model}\n(n = {n})")
    ) |>
    arrange(str_rank_natural(modeln)) |>
    mutate(
      modeln = forcats::fct_inorder(modeln)
    ) |>
    group_by(rank, n, modeln, estimand, true_rank) |>
    summarize(
      mse = mean(squared_error),
      coverage = mean(covered)
    ) |>
    ungroup() |>
    mutate(
      `True rank` = true_rank,
      estimand = toupper(estimand)
    ) |>
    filter(rank <= 2 * true_rank)

  plot <- summarized_losses |>
    ggplot(aes(rank, mse, color = estimand)) +
    # geom_smooth(se = FALSE) +
    geom_vline(aes(xintercept = true_rank), linetype = "dashed") +
    geom_point() +
    geom_line() +
    # scale_x_log10() +
    # scale_y_log10() +
    scale_color_brewer(palette = "Dark2") +
    facet_grid(
      cols = vars(`True rank`),
      rows = vars(modeln),
      scales = "free",
      labeller = labeller(
        .rows = label_value,
        .cols = label_both
      )
    ) +
    labs(
      y = "Mean squared error",
      x = "Embedding dimension"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      # legend.position = "bottom",
      legend.title = element_blank()
    )

  if (!dir.exists(here("figures/simulations/misspecification/"))) {
    dir.create(here("figures/simulations/misspecification/"))
  }

  path <- here(
    glue("figures/simulations/misspecification/mean_squared_error.{file_type}")
  )

  ggsave(
    path,
    plot = plot,
    width = 7.5,
    height = 8,
    ...
  )

  plot2 <- summarized_losses |>
    ggplot(aes(rank, coverage, color = estimand)) +
    # geom_smooth(se = FALSE) +
    geom_hline(yintercept = 0.95, linetype = "dashed") +
    geom_vline(aes(xintercept = true_rank), linetype = "dashed") +
    geom_point() +
    geom_line() +
    # scale_x_log10() +
    # scale_y_log10() +
    scale_color_brewer(palette = "Dark2") +
    facet_grid(
      cols = vars(`True rank`),
      rows = vars(modeln),
      scales = "free",
      labeller = labeller(
        .rows = label_value,
        .cols = label_both
      )
    ) +
    expand_limits(y = c(0, 1)) +
    labs(
      y = "Coverage rate",
      x = "Embedding dimension"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      # legend.position = "bottom",
      legend.title = element_blank()
    )

  path2 <- here(
    glue("figures/simulations/misspecification/coverage.{file_type}")
  )

  ggsave(
    path2,
    plot = plot2,
    width = 7.5,
    height = 8,
    ...
  )


  plot3 <- losses |>
    mutate(
      model = stringr::str_remove(model_name, "_normal"),
      model = ifelse(model == "block", "informative", model),
      modeln = glue("{model}\n(n = {n})")
    ) |>
    arrange(str_rank_natural(modeln)) |>
    mutate(
      `True rank` = true_rank,
      estimand = toupper(estimand)
    ) |>
    filter(rank <= 2 * true_rank) |>
    ggplot(aes(rank, estimate - truth, color = estimand, group = paste0(rep_within_chunk, truth))) +
    geom_line(alpha = 0.1) +
    geom_vline(aes(xintercept = true_rank), linetype = "dashed") +
    scale_color_brewer(palette = "Dark2") +
    facet_grid(
      cols = vars(`True rank`),
      rows = vars(modeln),
      scales = "free",
      labeller = labeller(
        .rows = label_value,
        .cols = label_both
      )
    ) +
    labs(
      y = "Bias = Estimate - Estimand",
      x = "Embedding dimension"
    ) +
    guides(color = guide_legend(override.aes = list(alpha = 1))) +
    theme_minimal(base_size = 16) +
    theme(
      legend.title = element_blank()
    )

  path3 <- here(
    glue("figures/simulations/misspecification/bias-trajectories.{file_type}")
  )

  ggsave(
    path3,
    plot = plot3,
    width = 7.5,
    height = 8,
    ...
  )

  plot4 <- losses |>
    mutate(
      model = stringr::str_remove(model_name, "_normal"),
      model = ifelse(model == "block", "informative", model),
      modeln = glue("{model}\n(n = {n})")
    ) |>
    arrange(str_rank_natural(modeln)) |>
    mutate(
      `True rank` = true_rank,
      estimand = toupper(estimand)
    ) |>
    filter(rank <= 2 * true_rank) |>
    ggplot(
      aes(
        x = rank,
        y = c(NA, diff(estimate - truth)),
        color = estimand,
        group = paste0(rep_within_chunk, truth)
      )
    ) +
    geom_line(alpha = 0.1) +
    geom_vline(aes(xintercept = true_rank), linetype = "dashed") +
    scale_color_brewer(palette = "Dark2") +
    facet_grid(
      cols = vars(`True rank`),
      rows = vars(modeln),
      scales = "free",
      labeller = labeller(
        .rows = label_value,
        .cols = label_both
      )
    ) +
    labs(
      y = "Change in bias as a function of embedding embedding dimension",
      x = "Embedding dimension"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      legend.title = element_blank()
    )

  path4 <- here(
    glue("figures/simulations/misspecification/bias-trajectory-slopes.{file_type}")
  )

  ggsave(
    path4,
    plot = plot4,
    width = 7.5,
    height = 8,
    ...
  )

  c(path, path2, path3, path4)
}
