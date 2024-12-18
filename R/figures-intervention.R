create_intervention_plots <- function() {
  k <- 5
  n <- 100
  B <- matrix(0.0, nrow = k, ncol = k)
  diag(B) <- 0.8

  Z <- matrix()

  z_helper <- function(i, n = 20, k = 5) {
    alpha <- rep(2, k)
    alpha[i] <- 40
    Z <- t(sample_dirichlet(n, alpha))
    Z
  }

  Z <- do.call("rbind", map(1:5, z_helper))

  trt <- rep(1, n)
  Z_post <- Z
  Z_post[, 2] <- Z[, 2] - trt * runif(n, 0.01, 0.05)
  Z_post[, 3] <- Z[, 3] + trt * runif(n, 0.1, 0.3)

  rownames(Z) <- 1:nrow(Z)
  colnames(Z) <- 1:ncol(Z)

  rownames(Z_post) <- 1:nrow(Z_post)
  colnames(Z_post) <- 1:ncol(Z_post)

  Z_post_long <- Z_post |>
    as_tibble(rownames = "row") |>
    gather(col, value, -row) |>
    mutate_all(as.numeric) |>
    mutate(
      type = "Intervention"
    )

  Z_long <- Z |>
    as_tibble(rownames = "row") |>
    gather(col, value, -row) |>
    mutate_all(as.numeric) |>
    mutate(
      type = "No intervention"
    )

  Z_full <- bind_rows(Z_long, Z_post_long) |>
    mutate(
      type = as.factor(type),
      type = fct_relevel(type, "No intervention", "Intervention"),
      col = glue("Z[{col}]")
    )

  latent_plot <- ggplot(Z_full, aes(x = col, y = row, fill = value)) +
    geom_raster() +
    scale_y_reverse() +
    scale_x_discrete(
      labels = label_parse()
    ) +
    scale_fill_gradient(
      low = "white",
      high = "black",
      guide = "none"
    ) +
    facet_grid(
      cols = vars(type)
    ) +
    theme_minimal(
      base_size = 10
    ) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
    ) +
    labs(
      title = "Community membership",
      y = "Node",
      x = "Block"
    )
  
  latent_plot

  path6 <- here("figures", "canonical-intervention", "latent-space.png")

  ggsave(
    filename = path6,
    plot = latent_plot,
    width = 5,
    height = 2.75,
    dpi = 800
  )

  pre <- undirected_factor_model(
    X = Z,
    S = B,
    expected_density = 0.2,
    poisson_edges = FALSE,
    allow_self_loops = FALSE
  )

  post <- pre
  post$X <- pmax(Z_post, 0) # lil hack

  plot1 <- plot_expectation(pre) +
    scale_fill_gradient(low = "white", high = "black", limits = c(0, 1)) +
    theme(
      legend.position = "none"
    )

  path1 <- here("figures", "canonical-intervention", "expected-a-pre-trt.pdf")

  ggsave(
    filename = path1,
    plot = plot1
  )

  plot2 <- plot_expectation(post) +
    scale_fill_gradient(low = "white", high = "black", limits = c(0, 1)) +
    theme(
      legend.position = "none"
    )

  path2 <- here("figures", "canonical-intervention", "expected-a-post-trt.pdf")

  ggsave(
    filename = path2,
    plot = plot2
  )

  A <- sample_sparse(pre)
  A_post <- sample_sparse(post)

  path4 <- here("figures", "canonical-intervention", "a-treated.pdf")

  ggsave(
    filename = path4,
    plot = plot_sparse_matrix(A_post) +
      theme(
        legend.position = "none"
      )
  )

  path5 <- here("figures", "canonical-intervention", "a-untreated.pdf")

  ggsave(
    filename = path5,
    plot = plot_sparse_matrix(A) +
      theme(
        legend.position = "none"
      )
  )

  diff <- as.matrix(expectation(post) - expectation(pre))
  rownames(diff) <- 1:n
  colnames(diff) <- 1:n

  plot3 <- plot_dense_matrix(diff) +
    scale_fill_gradient2(
      name = "Diff"
    )

  path3 <- here("figures", "canonical-intervention", "expected-a-pre-post-diff.pdf")

  ggsave(
    path3,
    plot3
  )

  c(path1, path2, path3, path4, path5, path6)
}
