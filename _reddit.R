library(targets)
library(here)

data(reddit, package = "netmediate")

tar_option_set(
  packages = c(
    "glue",
    "broom",
    "tidygraph",
    "netmediate",
    "invertiforms",
    "here",
    "ggplot2",
    "dplyr",
    "igraph",
    "stringr",
    "ggraph",
    "scales"
  ),
  imports = c("netmediate", "fastRG") # changes to these packages will invalidate pipeline
)

options(clustermq.scheduler = "multiprocess")

subreddit_network <- function(subreddit) {
  reddit |>
    activate(nodes) |>
    filter(subreddit == !!subreddit | node_type == "word")
}

as_csparse <- function(ig) {
  ig |>
    as_incidence_matrix(sparse = TRUE) |>
    as("CsparseMatrix")
}


plot_eigcv <- function(eigcv, subreddit, laplacian) {
  matrix <- if (laplacian) {
    "L"
  } else {
    "A"
  }

  plot <- gdim:::plot.eigcv(eigcv) +
    labs(
      y = "Z-score for cross-validated eigenvalue",
      title = glue("Rank estimates for {subreddit} based on {matrix}")
    )

  path <-
    here(
      "figures",
      "reddit",
      "rank-estimates",
      glue("{subreddit}-{matrix}.png")
    )

  ggsave(
    path,
    plot = plot,
    dpi = 600,
    width = 8,
    height = 8 * 9 / 16
  )

  path
}

plot_localization <- function(A, subreddit, rank, num_tau = 20) {
  tau <- 10^seq(0, 4, length.out = num_tau)

  print("here")

  laplacians <- tibble(tau = tau) %>%
    mutate(
      scaler = purrr::map(tau, ~ RegularizedLaplacian(A, .x, .x)),
      L_tau = purrr::map(scaler, ~ transform(.x, A))
    )

  print("here 2")

  decomposed <- laplacians %>%
    mutate(
      svd = purrr::map(L_tau, RSpectra::svds, k = rank)
    )

  cumulative_participation <- function(U) {
    sum(rowSums(U^2)^2)
  }

  print("here 3")

  localization <- decomposed %>%
    mutate(
      u = purrr::map_dbl(svd, ~ cumulative_participation(.x$u)),
      v = purrr::map_dbl(svd, ~ cumulative_participation(.x$v))
    )

  print("here 4")

  # these are the same now but might differ for rectangular A
  mean_rs <- mean(rowSums(A))
  mean_cs <- mean(colSums(A))

  plot <- localization %>%
    select(-scaler, -L_tau, -svd) %>%
    tidyr::gather(subspace, localization, u, v) %>%
    ggplot() +
    aes(tau, localization, color = subspace, group = subspace) +
    geom_line() +
    geom_vline(xintercept = mean_rs, linetype = "dashed") +
    geom_vline(xintercept = mean_cs, linetype = "dashed", color = "purple") +
    scale_color_viridis_d(begin = 0.15, end = 0.85) +
    scale_x_log10(
      breaks = trans_breaks("log10", function(x) 10^x),
      labels = trans_format("log10", math_format(10^.x))
    ) +
    labs(
      title = glue("Cumulative localization of rank {rank} decomposition of L_tau for r/{subreddit}")
    )

  path <-
    here(
      "figures",
      "reddit",
      "localization",
      glue("{subreddit}.png")
    )

  ggsave(
    path,
    plot = plot,
    dpi = 600,
    width = 8,
    height = 8 * 9 / 16
  )

  path
}

fit_curve <- function(A, subnetwork, max_rank, tau_row, tau_col) {
  iform <- RegularizedLaplacian(A, tau_row, tau_col)
  L <- transform(iform, A)

  s_max <- RSpectra::svds(L, max_rank, max_rank)
  X_max <- s_max$u %*% diag(sqrt(s_max$d))

  # and now we plug them into the product-of-coefs estimator

  sensitivity_curve_custom(subnetwork, score ~ flair, X_max)
}


plot_curve <- function(curve,
                       subreddit,
                       plot_file_type = "png") {
  plot <- plot(curve) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_minimal(base_size = 16) +
    labs(
      title = glue("Estimated effect of male flair on comment score in r/{subreddit}"),
      subtitle = "Using embeddings from the regularized graph Laplacian"
    )

  path <-
    here(
      "figures",
      "reddit",
      "sensitivity-curves",
      glue("{subreddit}.{plot_file_type}")
    )

  ggsave(
    path,
    plot = plot,
    dpi = 600,
    width = 8,
    height = 8 * 9 / 16
  )

  path
}

add_subreddit_column <- function(curve, subreddit) {
  curve |>
    mutate(subreddit = subreddit)
}



plot_curves_for_manuscript <- function(all, plot_file_type = "png") {
  plot <- all |>
    filter(rank <= 50) |>
    dplyr::mutate(
      estimand = dplyr::recode(
        estimand,
        nde = "NDE",
        nie = "NIE"
      )
    ) |>
    ggplot() +
    aes(
      x = rank,
      y = estimate,
      ymin = conf.low,
      ymax = conf.high,
      color = estimand,
      fill = estimand
    ) +
    geom_ribbon(alpha = 0.3) +
    geom_line() +
    geom_point() +
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    labs(
      x = "Number of estimated latent topics",
      y = "Causal effect of flair\n on comment score"
    ) +
    facet_grid(
      rows = vars(subreddit),
      scales = "free_y"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      legend.title = element_blank()
    )

  path <-
    here(
      "figures",
      "reddit",
      glue("effects.{plot_file_type}")
    )

  ggsave(
    path,
    plot = plot,
    dpi = 600,
    width = 8,
    height = 8 * 9 / 16
  )

  path
}


list(
  tar_target(
    plot_file_type,
    c("png", "pdf")
  ),
  tar_target(
    subreddit,
    c("OkCupid", "keto", "childfree")
  ),
  tar_target(
    subnetwork,
    subreddit_network(subreddit),
    pattern = map(subreddit),
    iteration = "list"
  ),
  tar_target(A,
    as_csparse(subnetwork),
    pattern = map(subnetwork)
  ),
  tar_target(
    A_eigcv,
    gdim::eigcv(A, k_max = 75, laplacian = FALSE),
    pattern = map(A)
  ),
  tar_target(
    A_eigcv_plots,
    plot_eigcv(A_eigcv, subreddit, FALSE),
    pattern = map(A_eigcv, subreddit),
    format = "file"
  ),
  tar_target(
    L_eigcv,
    gdim::eigcv(A, k_max = 150, laplacian = TRUE),
    pattern = map(A)
  ),
  tar_target(
    L_eigcv_plots,
    plot_eigcv(L_eigcv, subreddit, TRUE),
    pattern = map(L_eigcv, subreddit),
    format = "file"
  ),

  # tar_target(
  #   localization_plots,
  #   plot_localization(A, subreddit, rank = 75),
  #   pattern = map(A, subreddit),
  #   format = "file"
  # ),

  tar_target(
    fa,
    vsp::vsp(A, rank = 50, degree_normalize = TRUE, tau_row = 100, tau_col = 100),
    pattern = map(A)
  ),
  tar_target(
    curve,
    fit_curve(A, subnetwork, max_rank = 150, tau_row = 200, tau_col = 200),
    pattern = map(A, subnetwork),
    iteration = "list"
  ),
  tar_target(
    curve_plots,
    plot_curve(curve, subreddit, plot_file_type),
    pattern = cross(map(curve, subreddit), plot_file_type),
    iteration = "list",
    format = "file"
  ),
  tar_target(
    all_subreddit_curves,
    purrr::map2_dfr(curve, subreddit, add_subreddit_column)
  ),
  tar_target(
    reddit_plot_manuscript,
    plot_curves_for_manuscript(all_subreddit_curves, plot_file_type),
    pattern = map(plot_file_type),
    format = "file"
  )
)
