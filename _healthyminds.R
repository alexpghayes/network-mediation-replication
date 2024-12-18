library(targets)
library(tarchetypes)
library(crew)
library(here)

tar_option_set(
  controller = crew_controller_local(workers = 8),
  packages = c(
    "broom",
    "distributional",
    "estimatr",
    "forcats",
    "gdim",
    "GGally",
    "ggdag",
    "ggdist",
    "ggraph",
    "glue",
    "here",
    "hms",
    "igraph",
    "lubridate",
    "latentnetmediate",
    "marginaleffects",
    "Matrix",
    "methods",
    "scales",
    "tidygraph",
    "tidyverse",
    "vsp"
  ),
  imports = "latentnetmediate" # changes to these packages will invalidate pipeline
)

tar_source()

get_sparse_survey_responses <- function(df) {
  A <- df |>
    select(
      contains("nih_lonliness"),
      contains("mlq_10"),
      contains("ffmq_8"),
      contains("dds_10")
    ) |>
    as.matrix() |>
    as("sparseMatrix")
  rownames(A) <- df$id
  A
}

get_node_data <- function(df) {
  df |>
    select(
      -contains("nih_lonliness"),
      -contains("mlq_10"),
      -contains("ffmq_8"),
      -contains("dds_10")
    ) |>
    mutate(
      name = as.character(id)
    ) |>
    select(-id)
}

make_tbl_graph <- function(A, node_data) {
  graph_from_biadjacency_matrix(A, weighted = TRUE) |>
    as_tbl_graph() |>
    left_join(
      node_data,
      by = "name"
    )
}

nest_by_event <- function(data, rank = 5, max_rank = 18) {
  data |>
    na.omit() |>
    nest(
      data = -c(redcap_event_name, event_name, event_name_fct)
    ) |>
    mutate(
      A = map(data, get_sparse_survey_responses),
      node_data = map(data, get_node_data),
      tbl_graph = map2(A, node_data, make_tbl_graph),
      fa = map(A, vsp, rank = rank, degree_normalize = FALSE),
      ecv = map(A, eigcv, k_max = 10, laplacian = FALSE, regularize = FALSE),
      meddep = map(tbl_graph, netmediate, depression ~ intervention + age + I(age^2) + sex, rank = rank),
      medanx = map(tbl_graph, netmediate, anxiety ~ intervention + age + I(age^2) + sex, rank = rank),
      curvedep = map(tbl_graph, sensitivity_curve, depression ~ intervention + age + sex, max_rank = max_rank, ranks_to_consider = max_rank - 1),
      curveanx = map(tbl_graph, sensitivity_curve, anxiety ~ intervention + age + sex, max_rank = max_rank, ranks_to_consider = max_rank - 1)
    )
}

list(
  tar_target(
    redcap_data_path,
    # here("data", "healthyminds", "df_camped_alex2.csv"),
    here("data", "healthyminds", "deidentified.csv"),
    format = "file"
  ),
  tar_target(
    data,
    clean_cahmped_data(redcap_data_path)
  ),
  tar_target(
    ate_figure,
    make_ate_figure(data),
    format = "file"
  ),
  tar_target(
    nested_by_event,
    nest_by_event(data)
  ),
  tar_group_by(
    grouped_by_event,
    nested_by_event,
    event_name
  ),
  tar_target(
    eigcv_plots,
    make_eigcv_plots(grouped_by_event$ecv[[1]], grouped_by_event$redcap_event_name),
    pattern = map(grouped_by_event)
  ),
  tar_target(
    mediation_dag,
    make_mediation_dag(),
    format = "file"
  ),
  tar_target(
    week4_scatter_figure,
    make_week4_scatter_figure(data),
    format = "file"
  ),
  tar_target(
    week4_responses_figure,
    make_week4_responses_figure(nested_by_event$A[[5]]),
    format = "file"
  ),
  tar_target(
    bipartite_mediation_figure,
    make_bipartite_mediation_figure(),
    format = "file"
  ),
  tar_target(
    xhat_figure,
    make_xhat_figure(nested_by_event$fa[[5]]),
    format = "file"
  ),
  tar_target(
    vhat_figure,
    make_vhat_figure(nested_by_event$fa[[5]]),
    format = "file"
  ),
  tar_target(
    yhat_figure,
    make_yhat_figure(nested_by_event$fa[[5]]),
    format = "file"
  ),
  tar_target(
    rank_estimate_figures,
    make_rank_estimate_figures(nested_by_event$A[[5]], nested_by_event$ecv[[5]]),
    format = "file"
  ),
  tar_target(
    regression_figures,
    make_regression_figures(nested_by_event),
    format = "file"
  ),
  tar_target(
    regression_figures_depression,
    make_regression_figures_depression(nested_by_event),
    format = "file"
  ),
  tar_target(
    sensitivity_figure,
    make_sensitivity_figure(nested_by_event$curveanx[[5]]),
    format = "file"
  ),
  tar_target(
    trajectory_figure,
    make_mediation_trajectory_figure(nested_by_event),
    format = "file"
  ),
  tar_target(
    trajectory_figure_depression,
    make_mediation_trajectory_figure_depression(nested_by_event),
    format = "file"
  ),
  tar_target(
    positivity_figure,
    make_latent_positivity_plot(nested_by_event),
    format = "file"
  )
)
