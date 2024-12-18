library(targets)
library(crew)
library(here)

data(glasgow, package = "latentnetmediate")

source(here("R", "glasgow.R"))

tar_option_set(
  controller = crew_controller_local(workers = 8),
  packages = c(
    "broom",
    "dplyr",
    "GGally",
    "ggraph",
    "ggplot2",
    "glue",
    "here",
    "igraph",
    "invertiforms",
    "latentnetmediate",
    "stringr",
    "tidygraph",
    "vsp"
  ),
  imports = c("latentnetmediate", "fastRG") # changes to these packages will invalidate pipeline
)

list(
  tar_target(
    plot_file_type,
    c("png", "pdf")
  ),
  tar_target(
    treatment,
    c("money", "sex_fct")
  ),
  tar_target(
    outcome_fct,
    c(
      "tobacco_fct",
      "alcohol_fct",
      "cannabis_fct"
    )
  ),
  tar_target(
    outcome_dimaria,
    c(
      "tobacco_dimaria",
      "alcohol_dimaria",
      "cannabis_dimaria"
    )
  ),
  tar_target(
    outcome_int,
    c(
      "tobacco_int",
      "alcohol_int",
      "cannabis_int"
    )
  ),
  tar_target(laplacian, c(FALSE, TRUE)),
  tar_target(coembedding, c("U", "V", "symmetrized")),
  tar_target(
    measures,
    c(treatment, outcome_fct)
  ),
  tar_target(clean,
    purrr::map(glasgow, clean_network),
    iteration = "list"
  ),
  tar_target(
    sex_plot_manuscript,
    plot_sex_manuscript(clean, plot_file_type),
    pattern = cross(slice(clean, index = 1), plot_file_type),
    format = "file"
  ),
  tar_target(
    tobacco_plot_manuscript,
    plot_tobacco_manuscript(clean, plot_file_type),
    pattern = cross(slice(clean, index = 1), plot_file_type),
    format = "file"
  ),
  tar_target(
    network_plots,
    purrr::map2_chr(
      clean,
      1:3,
      \(graph, time) plot_network(graph, time, plot_file_type)
    ),
    pattern = map(plot_file_type),
    format = "file"
  ),
  tar_target(
    network_measure_plots,
    purrr::map2_chr(
      clean,
      1:3,
      \(graph, time) plot_measure(
        graph,
        time,
        measures,
        plot_file_type
      )
    ),
    pattern = cross(plot_file_type, measures),
    format = "file"
  ),
  tar_target(
    eigcv_A,
    purrr::map(
      clean,
      eigcv_rank_estimate,
      k_max = 50,
      laplacian = FALSE
    )
  ),
  tar_target(
    eigcv_L,
    purrr::map(
      clean,
      eigcv_rank_estimate,
      k_max = 50,
      laplacian = TRUE
    )
  ),
  tar_target(
    positivity_plots,
    purrr::map2_chr(
      clean,
      1:3,
      make_sex_positivity_plot
    ),
    format = "file"
  ),
  tar_target(
    eigcv_A_plot,
    purrr::map2_chr(eigcv_A, 1:3, plot_eigcv, laplacian = FALSE),
    format = "file"
  ),
  tar_target(
    eigcv_L_plot,
    purrr::map2_chr(eigcv_L, 1:3, plot_eigcv, laplacian = TRUE),
    format = "file"
  ),
  tar_target(
    curve,
    purrr::map(
      clean,
      fit_curve,
      max_rank = 25,
      laplacian = laplacian,
      coembedding = coembedding,
      outcome = outcome_dimaria,
      treatment = treatment
    ),
    pattern = cross(laplacian, coembedding, outcome_dimaria, treatment)
  ),
  tar_target(
    curve_plots,
    purrr::map2_chr(
      curve,
      1:3,
      plot_curve,
      laplacian = laplacian,
      coembedding = coembedding,
      outcome = outcome_dimaria,
      treatment = treatment
    ),
    pattern = map(
      curve,
      cross(laplacian, coembedding, outcome_dimaria, treatment)
    ),
    format = "file"
  ),
  tar_target(
    manuscript_plot,
    plot_curve_manuscript(clean, plot_file_type),
    pattern = cross(slice(clean, index = 1), plot_file_type),
    format = "file"
  )
)
