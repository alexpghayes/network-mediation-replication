# status_colors <- c(
#   "exposure" = "firebrick",
#   "outcome" = "steelblue",
#   "observed" = "black",
#   "latent" = "grey"
# )

status_colors <- c(
  "exposure" = "black",
  "outcome" = "black",
  "observed" = "black",
  "latent" = "black"
)

make_hm_mediating_dag <- function() {
  coords <- tibble(
    name = c("Ci", "Ti", "Xi", "Yi"),
    x = c(1, 0, 2, 1),
    y = c(2, 1, 1, 0)
  )
  
  #  example from the dagitty package
  dag <- dagify(
    Yi ~ Ti + Ci + Xi,
    Xi ~ Ti + Ci,
    Ti ~ Ci,
    coords = coords,
    labels = c(
      Ci = "Age & sex",
      Ti = "Meditation",
      Xi = "Cognitive state",
      Yi = "Anxiety"
    ),
    exposure = "Ti",
    outcome = "Yi"
  )
  
  
  dag |>
    tidy_dagitty() |>
    node_status() |>
    mutate(
      status = replace_na(as.character(status), "observed")
    ) |>
    ggplot(
      aes(
        x = x,
        y = y,
        xend = xend,
        yend = yend
      ) 
    ) +
    geom_dag_label(aes(color = status, label = label), fill = "black", color = "white") +
    scale_color_manual(values = status_colors, guide = "none") +
    geom_dag_edges_fan(
      # arrow = grid::arrow(length = grid::unit(7, "pt"), type = "closed"),
      aes(
        start_cap = label_rect(label, padding = margin(2.5, 2.5, 2.5, 2.5, "mm")),
        end_cap = label_rect(label, padding = margin(10, 10, 10, 10, "mm"))
      )
    ) +
    # geom_dag_text(aes(label = label), size = 5) +
    theme_dag(base_size = 22, base_family = "Fira Sans")
  
  path <- here::here("figures", "dags", "mediating.png")
  
  ggsave(
    path,
    height = 3.5,
    width = 3.5,
    dpi = 500
  )
  
  path
}

make_mediation_dag <- function() {
  coords <- tibble(
    name = c("Ci", "Ti", "Xi", "Yi"),
    x = c(1, 0, 2, 1),
    y = c(2, 1, 1, 0)
  )

  #  example from the dagitty package
  dag <- dagify(
    Yi ~ Ti + Ci + Xi,
    Xi ~ Ti + Ci,
    Ti ~ Ci,
    coords = coords,
    labels = c(
      Ci = "C[i %.% phantom(j)]",
      Ti = "T[i]",
      Xi = "X[i %.% phantom(j)]",
      Yi = "Y[i]"
    ),
    exposure = "Ti",
    latent = "Xi",
    outcome = "Yi"
  )

  dag |>
    tidy_dagitty() |>
    node_status() |>
    mutate(
      status = replace_na(as.character(status), "observed")
    ) |>
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_point(aes(color = status)) +
    scale_color_manual(values = status_colors, guide = "none") +
    geom_dag_edges_diagonal(
      arrow = grid::arrow(length = grid::unit(7, "pt"), type = "closed")
    ) +
    geom_dag_text(aes(label = label), parse = TRUE, size = 5) +
    theme_dag(base_size = 22, base_family = "Computer Modern")

  path <- here::here("figures", "dags", "mediating.png")

  ggsave(
    path,
    height = 3.5,
    width = 3.5,
    dpi = 500
  )

  path
}

make_bipartite_mediation_figure <- function() {
  coords <- tibble(
    name = c("Ci", "Ti", "Xi", "Yi", "Aij", "Vj"),
    x = c(1, 0, 2, 1, 3, 4),
    y = c(2, 1, 1, 0, 0, 1)
  )

  #  example from the dagitty package
  dag <- dagify(
    Yi ~ Ti + Ci + Xi,
    Aij ~ Xi + Vj,
    Xi ~ Ti + Ci,
    Ti ~ Ci,
    coords = coords,
    latent = c("Xi", "Vj"),
    labels = c(
      Ci = "C[i %.% phantom(j)]",
      Ti = "T[i]",
      Xi = "X[i %.% phantom(j)]",
      Yi = "Y[i]",
      Aij = "A[ij]",
      Vj = "L[j %.% phantom(j)]"
    ),
    exposure = "Ti",
    outcome = "Yi"
  )

  dag |>
    tidy_dagitty() |>
    node_status() |>
    mutate(
      status = replace_na(as.character(status), "observed")
    ) |>
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_point(aes(color = status)) +
    scale_color_manual(values = status_colors, guide = "none") +
    geom_dag_edges_diagonal(
      arrow = grid::arrow(length = grid::unit(7, "pt"), type = "closed")
    ) +
    geom_dag_text(aes(label = label), parse = TRUE, size = 5) +
    theme_dag(base_size = 16, base_family = "Fira Sans")

  path <- here::here("figures", "dags", "bipartite-mediation.png")

  ggsave(
    path,
    height = 3.5,
    width = 5,
    dpi = 300
  )

  path
}


make_homophily_mediating_figure <- function() {
  
  coords <- tibble(
    name = c("Ci", "Ti", "Xi", "Yi", "Aij", "Xj", "Cj", "Yj", "Tj"),
    x = c(1, 0, 2, 1, 3, 4, 5, 5, 6),
    y = c(2, 1, 1, 0, 0, 1, 2, 0, 1)
  )
  
  #  example from the dagitty package
  dag <- dagify(
    Yi ~ Ti + Ci + Xi,
    Aij ~ Xi + Xj,
    Xi ~ Ti + Ci,
    Ti ~ Ci,
    Yj ~ Tj + Cj + Xj,
    Xj ~ Tj + Cj,
    Tj ~ Cj,
    coords = coords,
    latent = c("Xi", "Xj"),
    labels = c(
      Ci = "C[i %.% phantom(j)]",
      Ti = "T[i]",
      Xi = "X[i %.% phantom(j)]",
      Yi = "Y[i]",
      Aij = "A[ij]",
      Xj = "X[j %.% phantom(j)]",
      Cj = "C[j %.% phantom(j)]",
      Yj = "Y[j]",
      Tj = "T[j]"
    ),
    exposure = "Ti",
    outcome = "Yi"
  )
  
  dag |>
    tidy_dagitty() |>
    node_status() |>
    mutate(
      status = replace_na(as.character(status), "observed")
    ) |>
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_point(aes(color = status)) +
    scale_color_manual(values = status_colors, guide = "none") +
    geom_dag_edges_diagonal(
      arrow = grid::arrow(length = grid::unit(7, "pt"), type = "closed")
    ) +
    geom_dag_text(aes(label = label), parse = TRUE, size = 5) +
    theme_dag(base_size = 22, base_family = "Computer Modern")
  
  path <- here::here("figures", "dags", "homophily-mediating.png")
  
  ggsave(
    path,
    height = 3.5,
    width = 3.5 * 16/9,
    dpi = 500
  )
  
  path
}

make_confounding_homophily_figure <- function() {
  
  coords <- tibble(
    name = c("Ci", "Ti", "Xi", "Yi", "Aij", "Xj", "Cj", "Yj", "Tj"),
    x = c(1, 0, 2, 1, 3, 4, 5, 5, 6),
    y = c(2, 1, 1, 0, 0, 1, 2, 0, 1)
  )
  
  #  example from the dagitty package
  dag <- dagify(
    Yi ~ Ti + Ci + Xi,
    Aij ~ Xi + Xj,
    Xi ~ Ci,
    Ti ~ Xi + Ci,
    Yj ~ Tj + Cj + Xj,
    Xj ~ Cj,
    Tj ~ Xj + Cj,
    coords = coords,
    latent = c("Xi", "Xj"),
    labels = c(
      Ci = "C[i %.% phantom(j)]",
      Ti = "T[i]",
      Xi = "X[i %.% phantom(j)]",
      Yi = "Y[i]",
      Aij = "A[ij]",
      Xj = "X[j %.% phantom(j)]",
      Cj = "C[j %.% phantom(j)]",
      Yj = "Y[j]",
      Tj = "T[j]"
    ),
    exposure = "Ti",
    outcome = "Yi"
  )
  
  dag |>
    tidy_dagitty() |>
    node_status() |>
    mutate(
      status = replace_na(as.character(status), "observed")
    ) |>
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_point(aes(color = status)) +
    scale_color_manual(values = status_colors, guide = "none") +
    geom_dag_edges_diagonal(
      arrow = grid::arrow(length = grid::unit(7, "pt"), type = "closed")
    ) +
    geom_dag_text(aes(label = label), parse = TRUE, size = 5) +
    theme_dag(base_size = 22, base_family = "Computer Modern")
  
  path <- here::here("figures", "dags", "confounding-homophily.png")
  
  ggsave(
    path,
    height = 3.5,
    width = 3.5 * 16/9,
    dpi = 500
  )
  
  path
}

make_confounding_homophily_interference_figure <- function() {
  
  coords <- tibble(
    name = c("Ci", "Ti", "Xi", "Yi", "Aij", "Xj", "Cj", "Yj", "Tj"),
    x = c(1, 0, 2, 1, 3, 4, 5, 5, 6),
    y = c(2, 1, 1, 0, 0, 1, 2, 0, 1)
  )
  
  #  example from the dagitty package
  dag <- dagify(
    Yi ~ Tj + Aij + Ti + Ci + Xi,
    Aij ~ Xi + Xj,
    Xi ~ Ci,
    Ti ~ Xi + Ci,
    Yj ~ Ti + Aij + Tj + Cj + Xj,
    Xj ~ Cj,
    Tj ~ Xj + Cj,
    coords = coords,
    latent = c("Xi", "Xj"),
    labels = c(
      Ci = "C[i %.% phantom(j)]",
      Ti = "T[i]",
      Xi = "X[i %.% phantom(j)]",
      Yi = "Y[i]",
      Aij = "A[ij]",
      Xj = "X[j %.% phantom(j)]",
      Cj = "C[j %.% phantom(j)]",
      Yj = "Y[j]",
      Tj = "T[j]"
    ),
    exposure = "Ti",
    outcome = "Yi"
  )
  
  dag |>
    tidy_dagitty() |>
    node_status() |>
    mutate(
      status = replace_na(as.character(status), "observed")
    ) |>
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_point(aes(color = status)) +
    scale_color_manual(values = status_colors, guide = "none") +
    geom_dag_edges_diagonal(
      arrow = grid::arrow(length = grid::unit(7, "pt"), type = "closed")
    ) +
    geom_dag_text(aes(label = label), parse = TRUE, size = 5) +
    theme_dag(base_size = 22, base_family = "Computer Modern")
  
  path <- here::here("figures", "dags", "confounding-homophily-interference.png")
  
  ggsave(
    path,
    height = 3.5,
    width = 3.5 * 16/9,
    dpi = 500
  )
  
  path
}
