library(tidygraph)
library(tidyverse)
library(here)
library(ggraph)
library(ggdag)

make_dag_figures <- function() {
  c(
    make_mediating_figure(),
    make_confounding_figure(),
    make_interference_figure(),
    make_contagion_figure(),
    make_full_mediating_figure(),
    make_full_confounding_figure()
  )
}

make_mediating_figure <- function() {

  coords <- tibble(
    name = c("Ci", "Ti", "Xi", "Yi", "Aij", "Xj"),
    x = c(1, 0, 2, 1, 3, 4),
    y = c(2, 1, 1, 0, 0, 1)
  )

  #  example from the dagitty package
  dag <- dagify(
    Yi ~ Ti + Ci + Xi,
    Aij ~ Xi + Xj,
    Xi ~ Ti + Ci,
    Ti ~ Ci,
    coords = coords,
    latent = c("Xi", "Xj"),
    labels = c(
      Ci = "C[i %.% phantom(j)]",
      Ti = "T[i]",
      Xi = "X[i %.% phantom(j)]",
      Yi = "Y[i]",
      Aij = "A[ij]",
      Xj = "X[j %.% phantom(j)]"
    ),
    exposure = "Ti",
    outcome = "Yi"
  )

  dag |>
    tidy_dagitty() |>
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_point() +
    geom_dag_edges_diagonal() +
    geom_dag_text(aes(label = label), parse = TRUE, size = 5) +
    theme_dag(base_size = 22, base_family = "Computer Modern")

  path <- here::here("figures", "dags", "mediating.png")

  ggsave(
    path,
    height = 3.5,
    width = 3.5 * 16/9,
    dpi = 500
  )

  path
}


make_confounding_figure <- function() {

  coords <- tibble(
    name = c("Ci", "Ti", "Xi", "Yi", "Aij", "Xj"),
    x = c(1, 0, 2, 1, 3, 4),
    y = c(2, 1, 1, 0, 0, 1)
  )

  #  example from the dagitty package
  dag <- dagify(
    Yi ~ Ti + Ci + Xi,
    Aij ~ Xi + Xj,
    Xi ~ Ci,
    Ti ~ Ci + Xi,
    coords = coords,
    latent = c("Xi", "Xj"),
    labels = c(
      Ci = "C[i %.% phantom(j)]",
      Ti = "T[i]",
      Xi = "X[i %.% phantom(j)]",
      Yi = "Y[i]",
      Aij = "A[ij]",
      Xj = "X[j %.% phantom(j)]"
    ),
    exposure = "Ti",
    outcome = "Yi"
  )

  dag |>
    tidy_dagitty() |>
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_point() +
    geom_dag_edges_diagonal() +
    geom_dag_text(aes(label = label), parse = TRUE, size = 5) +
    theme_dag(base_size = 22, base_family = "Computer Modern")

  path <- here::here("figures", "dags", "confounding.png")

  ggsave(
    path,
    height = 3.5,
    width = 3.5 * 16/9,
    dpi = 500
  )

  path
}

make_full_mediating_figure <- function() {

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
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_point() +
    geom_dag_edges_diagonal() +
    geom_dag_text(aes(label = label), parse = TRUE, size = 5) +
    theme_dag(base_size = 22, base_family = "Computer Modern")

  path <- here::here("figures", "dags", "full_mediating.png")

  ggsave(
    path,
    height = 3.5,
    width = 3.5 * 16/9,
    dpi = 500
  )

  path
}

make_full_confounding_figure <- function() {

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
    Ti ~ Ci + Xi,
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
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_point() +
    geom_dag_edges_diagonal() +
    geom_dag_text(aes(label = label), parse = TRUE, size = 5) +
    theme_dag(base_size = 22, base_family = "Computer Modern")

  path <- here::here("figures", "dags", "full_confounding.png")

  ggsave(
    path,
    height = 3.5,
    width = 3.5 * 16/9,
    dpi = 500
  )

  path
}

make_interference_figure <- function() {

  coords <- tibble(
    name = c("Ci", "Ti", "Xi", "Yi", "Aij", "Xj", "Cj", "Yj", "Tj"),
    x = c(1, 0, 2, 1, 3, 4, 5, 5, 6),
    y = c(2, 1, 1, 0, 0, 1, 2, 0, 1)
  )

  #  example from the dagitty package
  dag <- dagify(
    Yi ~ Tj + Aij + Ti + Ci + Xi,
    Aij ~ Xi + Xj,
    Xi ~ Ti + Ci,
    Ti ~ Ci,
    Yj ~ Ti + Aij + Tj + Cj + Xj,
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
    exposure = "Tj",
    outcome = "Yi"
  )

  dag |>
    tidy_dagitty() |>
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_point() +
    geom_dag_edges_diagonal() +
    geom_dag_text(aes(label = label), parse = TRUE, size = 5) +
    theme_dag(base_size = 22, base_family = "Computer Modern")

  path <- here::here("figures", "dags", "interference.png")

  ggsave(
    path,
    height = 3.5,
    width = 3.5 * 16/9,
    dpi = 500
  )

  path
}


make_contagion_figure <- function() {

  coords <- tibble(
    name = c("Ci", "Ti", "Xi", "Yi", "Aij", "Xj", "Cj", "Yj", "Tj"),
    x = c(1, 0, 2, 1, 3, 4, 5, 5, 6),
    y = c(2, 1, 1, 0, 0.5, 1, 2, 0, 1)
  )

  #  example from the dagitty package
  dag <- dagify(
    Yi ~ Yj + Aij + Ti + Ci + Xi,
    Aij ~ Xi + Xj,
    Xi ~ Ti + Ci,
    Ti ~ Ci,
    Yj ~ Yi + Aij + Tj + Cj + Xj,
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
    exposure = "Tj",
    outcome = "Yi"
  )

  dag |>
    tidy_dagitty() |>
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_point() +
    geom_dag_edges_diagonal() +
    geom_dag_text(aes(label = label), parse = TRUE, size = 5) +
    theme_dag(base_size = 22, base_family = "Computer Modern")

  path <- here::here("figures", "dags", "contagion.png")

  ggsave(
    path,
    height = 3.5,
    width = 3.5 * 16/9,
    dpi = 500
  )

  path
}

