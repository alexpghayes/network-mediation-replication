clean_network <- function(graph) {
  graph |>
    activate(edges) |>
    filter(friendship != "Structurally missing") |>
    activate(nodes) |>
    mutate(
      in_degree = centrality_degree(mode = "in"),
      out_degree = centrality_degree(mode = "out")
    ) |>
    filter(in_degree > 0 | out_degree > 0) |>
    mutate(
      tobacco_dimaria = as.numeric(tobacco_int > 1),
      alcohol_dimaria = as.numeric(alcohol_int > 2),
      cannabis_dimaria = as.numeric(cannabis_int > 2)
    )
}


plot_network <- function(graph, time, file_type) {
  layout <- graph |>
    create_layout(layout = "stress")

  plot <- layout |>
    ggraph() +
    geom_edge_fan(
      arrow = arrow(length = unit(1, "mm")),
      end_cap = circle(2.5, "mm")
    ) +
    geom_node_point(aes(size = in_degree)) +
    labs(
      title = "Friendships in a secondary school in Glasgow",
      subtitle = glue("Teenage Friends and Lifestyle Study, 1995 (wave {time})"),
      caption = "Each node represents one student",
      size = "Popularity (in-degree)"
    ) +
    theme_graph()

  path <-
    here(
      "figures",
      "glasgow",
      "network",
      glue("time{time}.{file_type}")
    )

  ggsave(
    path,
    plot = plot,
    device = if (file_type == "pdf") cairo_pdf else NULL,
    dpi = 600,
    width = 8,
    height = 8 * 9 / 16,
    create.dir = TRUE
  )

  path
}

plot_measure <- function(graph, time, measure, file_type) {
  layout <- graph |>
    create_layout(layout = "stress")

  plot <- layout |>
    ggraph() +
    geom_edge_fan(
      arrow = arrow(length = unit(1, "mm")),
      end_cap = circle(2.5, "mm")
    ) +
    geom_node_point(aes_string(size = "in_degree", color = measure)) +
    labs(
      title = "Friendships in a secondary school in Glasgow",
      subtitle = glue("Teenage Friends and Lifestyle Study, 1995 (wave {time})"),
      caption = "Each node represents one student",
      size = "Popularity (in-degree)"
    ) +
    theme_graph()

  path <-
    here(
      "figures",
      "glasgow",
      "network",
      glue("time{time}-{measure}.{file_type}")
    )

  ggsave(
    path,
    plot = plot,
    dpi = 600,
    device = if (file_type == "pdf") cairo_pdf else NULL,
    width = 8,
    height = 8 * 9 / 16,
    create.dir = TRUE
  )

  path
}

eigcv_rank_estimate <- function(tbl_graph, ...) {
  A <- igraph::as_adj(tbl_graph)

  gdim::eigcv(A, ...)
}

plot_eigcv <- function(eigcv, time, laplacian) {
  plot <- gdim:::plot.eigcv(eigcv) +
    labs(y = "Z-score for cross-validated eigenvalue")

  matrix <- if (laplacian) {
    "L"
  } else {
    "A"
  }

  path <-
    here(
      "figures",
      "glasgow",
      "rank-estimates",
      glue("eigcv-time{time}-{matrix}.png")
    )

  ggsave(
    path,
    plot = plot,
    dpi = 600,
    width = 8,
    height = 8 * 9 / 16,
    create.dir = TRUE
  )

  path
}

plot_diagnostics <-
  function(tbl_graph, time, rank, laplacian) {
    # TODO
    library(vsp)

    fa1 <- vsp(tbl_graph, rank = 9, degree_normalize = FALSE)
    fa1

    plot_mixing_matrix(fa1)
    screeplot(fa1)
  }

fit_curve <-
  function(tbl_graph,
           max_rank,
           laplacian = c(FALSE, TRUE),
           coembedding = c("U", "V", "symmetrized"),
           outcome,
           treatment) {
    coembedding <- rlang::arg_match(coembedding)

    # first, we construct embeddings

    if (coembedding == "symmetrized") {
      tbl_graph <- tbl_graph |>
        igraph::as.undirected() |>
        as_tbl_graph()
    }

    A <- igraph::as_adj(tbl_graph)

    if (laplacian) {
      iform <- RegularizedLaplacian(A)
      A <- transform(iform, A)
    }

    s_max <- RSpectra::svds(A, max_rank, max_rank)

    if (coembedding == "V") {
      X_max <- s_max$v %*% diag(sqrt(s_max$d))
    } else {
      X_max <- s_max$u %*% diag(sqrt(s_max$d))
    }

    formula <- as.formula(glue("{outcome} ~ {treatment}"))

    sensitivity_curve_custom(tbl_graph, formula, X_max)
  }

plot_curve <- function(curve,
                       time,
                       laplacian = c(FALSE, TRUE),
                       coembedding = c("U", "V", "symmetrized"),
                       outcome,
                       treatment,
                       plot_file_type = "png") {
  clean_trt <- stringr::str_remove(treatment, "_fct")
  clean_out <- stringr::str_remove(outcome, "_dimaria")

  if (coembedding == "U") {
    clean_embedding <- "left spectral embedding of "
  } else if (coembedding == "V") {
    clean_embedding <- "right spectral embedding of"
  } else {
    clean_embedding <- "spectral embedding of symmetrized"
  }

  graph <-
    if (laplacian) {
      "regularized graph Laplacian"
    } else {
      "adjacency matrix"
    }

  plot <- plot(curve) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    theme_classic() +
    labs(
      title = glue(
        "Estimated effects of {clean_trt} on {clean_out} as function of latent dimension at time {time}"
      ),
      subtitle = glue("Using {clean_embedding} {graph}")
    )

  matrix <- if (laplacian) {
    "L"
  } else {
    "A"
  }

  path <-
    here(
      "figures",
      "glasgow",
      "sensitivity-curves",
      glue(
        "time{time}-{clean_trt}-{clean_out}-{matrix}-{coembedding}.{plot_file_type}"
      )
    )

  ggsave(
    path,
    plot = plot,
    dpi = 600,
    width = 8,
    height = 8 * 9 / 16,
    create.dir = TRUE
  )

  path
}

plot_curve_manuscript <- function(graph, plot_file_type) {
  curve_v <- sensitivity_curve(
    graph,
    tobacco_dimaria ~ sex_fct + age + leisure_church,
    coembedding = "V",
    max_rank = 25,
    ranks_to_consider = 24
  )

  plot_v <- curve_v |>
    dplyr::filter(str_detect(term, "sex_fct")) |>
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
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_ribbon(alpha = 0.3) +
    geom_line() +
    geom_point() +
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    labs(
      color = "Natural Effect",
      fill = "Natural Effect",
      x = "Latent dimension of network",
      y = "Causal effect of sex\non probability of smoking"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      legend.title = element_blank()
    )

  path_v <-
    here(
      "figures",
      "glasgow",
      glue("effects_v.{plot_file_type}")
    )

  ggsave(
    path_v,
    plot = plot_v,
    dpi = 600,
    width = 8,
    height = 8 * 9 / 16,
    create.dir = TRUE
  )

  curve_u <- sensitivity_curve(
    graph,
    tobacco_dimaria ~ sex_fct + age + leisure_church,
    coembedding = "U",
    max_rank = 25,
    ranks_to_consider = 24
  )

  plot_u <- curve_u |>
    dplyr::filter(str_detect(term, "sex_fct")) |>
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
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_ribbon(alpha = 0.3) +
    geom_line() +
    geom_point() +
    scale_fill_brewer(palette = "Dark2") +
    scale_color_brewer(palette = "Dark2") +
    labs(
      color = "Natural Effect",
      fill = "Natural Effect",
      x = "Latent dimension of network",
      y = "Causal effect of sex\non probability of smoking"
    ) +
    theme_minimal(base_size = 16) +
    theme(
      legend.title = element_blank()
    )

  path_u <-
    here(
      "figures",
      "glasgow",
      glue("effects_u.{plot_file_type}")
    )

  ggsave(
    path_u,
    plot = plot_u,
    dpi = 600,
    width = 8,
    height = 8 * 9 / 16,
    create.dir = TRUE
  )

  c(path_u, path_v)
}

plot_tobacco_manuscript <- function(graph, plot_file_type) {
  set.seed(28)

  tobacco_plot <- ggraph(graph) +
    geom_edge_fan(
      arrow = arrow(length = unit(1, "mm")),
      end_cap = circle(2.5, "mm"), alpha = 0.2
    ) +
    geom_node_point(aes(size = in_degree, color = tobacco_fct)) +
    scale_color_brewer(type = "seq", palette = "Set1", direction = -1) +
    scale_size_continuous(guide = "none") +
    labs(
      color = "Tobacco use",
    ) +
    theme_void(base_size = 24)

  path <- here("figures", "glasgow", glue("tobacco.{plot_file_type}"))

  ggsave(
    path,
    plot = tobacco_plot,
    dpi = 600,
    width = 9,
    height = 6,
    create.dir = TRUE
  )

  path
}

plot_sex_manuscript <- function(graph, plot_file_type) {
  set.seed(28)

  sex_plot <- ggraph(graph) +
    geom_edge_fan(
      arrow = arrow(length = unit(1, "mm")),
      end_cap = circle(2.5, "mm"), alpha = 0.2
    ) +
    geom_node_point(aes(size = in_degree, color = sex_fct)) +
    scale_color_brewer(type = "seq", palette = "Set2", direction = -1) +
    scale_size_continuous(guide = "none") +
    labs(
      color = "Recorded sex",
    ) +
    theme_void(base_size = 24)

  path <- here("figures", "glasgow", glue("sex.{plot_file_type}"))

  ggsave(
    path,
    plot = sex_plot,
    dpi = 600,
    width = 9,
    height = 6,
    create.dir = TRUE
  )

  path
}

make_sex_positivity_plot <- function(clean, time_index) {
  fa <- vsp(clean, rank = 12, degree_normalize = FALSE)

  plot <- clean |>
    bind_svd_v(fa) |>
    as_tibble() |>
    mutate(Sex = sex_fct) |>
    select(Sex, matches("v[0-9]+")) |>
    ggpairs(
      mapping = aes(color = Sex, fill = Sex, alpha = 0.5),
      columns = colnames(fa$v),
      columnLabels = paste0("widehat(F)[", 1:fa$rank, "]"),
      labeller = "label_parsed",
      legend = c(2, 1),
      upper = "blank",
      axisLabels = "none"
    ) +
    scale_color_brewer(type = "qual") +
    scale_fill_brewer(type = "qual") +
    scale_alpha_continuous(guide = "none") +
    theme_minimal()

  path <- here("figures", "glasgow", glue("positivity-time-{time_index}.pdf"))

  ggsave(
    path,
    plot = plot,
    dpi = 300,
    width = 10,
    height = 10,
    create.dir = TRUE
  )

  path
}
