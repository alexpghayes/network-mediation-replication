library(targets)
library(tarchetypes)
library(here)
library(crew)
library(tibble)

tar_config_set(
  seconds_reporter = 0.5,
  seconds_meta_append = 15,
  reporter_make = "timestamp_positives",
  script = "_simulations.R",
  store = "_simulations",
  garbage_collection = TRUE
)

tar_option_set(
  packages = c(
    "broom",
    "dplyr",
    "distributional",
    "estimatr",
    "fastRG",
    "forcats",
    "ggdist",
    "ggplot2",
    "glue",
    "here",
    "igraph",
    "latentnetmediate",
    "Matrix",
    "purrr",
    "scales",
    "stats",
    "stringr",
    "tibble",
    "tidygraph",
    "tidyr",
    "vegan",
    "wordspace"
  ),
  imports = c(
    "fastRG",
    "latentnetmediate"
  ),

  # https://books.ropensci.org/targets/performance.html
  controller = crew_controller_local(workers = 12),
  format = "qs",
  memory = "transient",
  garbage_collection = TRUE,
  storage = "worker",
  retrieval = "worker"
)

source(here("R", "figures-intervention.R"))
source(here("R", "simulation-helpers.R"))
source(here("R", "simulation-plotting.R"))

models <- c(
  "model_uninformative",
  "model_block",
  "model_nullnde",
  "model_nullnie"
)

static_branch_targets <- tar_map(
  unlist = FALSE,
  values = tibble(model = rlang::syms(models)),
  tar_target(
    population,
    purrr::map(1:chunk_size, ~ model(n = params$n, k = params$rank)),
    iteration = "list",
    pattern = cross(chunk_indices, map(params))
  ),
  tar_target(
    tbl_graph,
    purrr::map(population, sample_tidygraph),
    iteration = "list",
    pattern = map(population)
  ),
  tar_target(
    xhat,
    purrr::map2(tbl_graph, population, ~ embed_adjacency_matrix(.x, .y$mediator$k)),
    iteration = "list",
    pattern = map(tbl_graph, population)
  ),
  tar_target(
    node_data,
    purrr::map(tbl_graph, tidygraph::as_tibble),
    iteration = "list",
    pattern = map(tbl_graph)
  ),
  tar_target(
    rotated,
    purrr::map2(population, xhat, spectral_loss),
    iteration = "list",
    pattern = map(population, xhat)
  ),
  tar_target(

    # combined within chunk
    combined_u_losses_chunk,
    purrr::map_dfr(rotated, u_loss_helper, params = params, chunk = chunk_indices, .id = "rep_within_chunk"),

    # without chunking here, we would use
    #
    #   pattern = map(rotated, params)
    #
    # but now we have to use cross to replicate the params target for each
    # chunk
    pattern = map(rotated, cross(chunk_indices, map(params)))
  ),


  # there is some redundancy in the following. we could derive xhat from tbl_graph,
  # (i.e. the list of graphs), but this process seems to be slightly non-
  # deterministic, so we fix xhat once early one to avoid having slightly
  # different svds in different places (so far this has shown up as sign
  # flip issues mostly)

  tar_target(
    estimates,
    purrr::map2(node_data, xhat, fit_models),
    iteration = "list",
    pattern = map(node_data, xhat)
  ),
  tar_target(
    mediator_losses,
    purrr::map2(rotated, estimates, ~ mediator_loss(estimates = .y$m_fit, rotated = .x)),
    iteration = "list",
    pattern = map(rotated, estimates)
  ),
  tar_target(
    outcome_losses,
    purrr::map2(rotated, estimates, ~ outcome_loss(estimates = .y$o_fit, rotated = .x)),
    iteration = "list",
    pattern = map(rotated, estimates)
  ),
  tar_target(
    causal_losses,
    purrr::map2(estimates, population, ~ causal_loss(estimates = .x, model = .y)),
    iteration = "list",
    pattern = map(estimates, population)
  ),

  # more aggregation of losses within each chunk

  tar_target(
    combined_mediator_losses,
    purrr::map_dfr(
      mediator_losses,
      mediator_loss_helper,
      params = params,
      chunk = chunk_indices,
      .id = "rep_within_chunk"
    ),
    pattern = map(mediator_losses, cross(chunk_indices, map(params)))
  ),
  tar_target(
    combined_outcome_losses,
    purrr::map_dfr(
      outcome_losses,
      outcome_loss_helper,
      params = params,
      chunk = chunk_indices,
      .id = "rep_within_chunk"
    ),
    pattern = map(outcome_losses, cross(chunk_indices, map(params)))
  ),
  tar_target(
    combined_causal_losses,
    purrr::map_dfr(
      causal_losses,
      causal_loss_helper,
      params = params,
      chunk = chunk_indices,
      .id = "rep_within_chunk"
    ),
    pattern = map(causal_losses, cross(chunk_indices, map(params)))
  )
)

list(
  tar_target(
    canonical_intervention_figures,
    create_intervention_plots(),
    deployment = "main",
    format = "file"
  ),
  tar_target(
    plot_file_type,
    c("png", "pdf"),
    deployment = "main"
  ),
  tar_target(chunk_size, 10, deployment = "main"),
  tar_target(num_chunks, 1, deployment = "main"),  # TURN UP TO 30 TO MATCH PAPER

  # 100  182  331  603 1099 2000
  tar_target(n, c(100, 182, 331, 603, 1099), deployment = "main"),
  tar_target(rank, c(2, 5, 7), deployment = "main"),
  tar_target(chunk_indices, 1:num_chunks, deployment = "main"),
  tar_target(
    params,
    tibble::tibble(
      n = n,
      rank = rank
    ),
    pattern = cross(n, rank),
    deployment = "main"
  ),
  static_branch_targets,
  tar_combine(
    combined_u_losses,
    static_branch_targets[[6]],
    command = bind_rows(!!!.x)
  ),
  tar_group_by(
    grouped_u_losses,
    combined_u_losses,
    model_name
  ),
  tar_target(
    u_loss_plot,
    plot_spectral_loss(grouped_u_losses, file_type = plot_file_type),
    pattern = cross(grouped_u_losses, plot_file_type),
    format = "file"
  ),
  tar_combine(
    all_mediator_losses,
    static_branch_targets[[11]],
    command = bind_rows(!!!.x)
  ),
  tar_combine(
    all_outcome_losses,
    static_branch_targets[[12]],
    command = bind_rows(!!!.x)
  ),
  tar_combine(
    all_causal_losses,
    static_branch_targets[[13]],
    command = bind_rows(!!!.x)
  ),
  tar_group_by(
    grouped_mediator_losses,
    all_mediator_losses,
    model_name
  ),
  tar_group_by(
    grouped_outcome_losses,
    all_outcome_losses,
    model_name
  ),
  tar_target(
    mediator_glance_plot,
    plot_mediator_glance(grouped_mediator_losses, file_type = plot_file_type),
    pattern = cross(grouped_mediator_losses, plot_file_type),
    format = "file"
  ),
  tar_target(
    average_elementwise_mediator_loss_plot,
    plot_mediator_elementwise_loss(
      grouped_mediator_losses,
      file_type = plot_file_type
    ),
    pattern = cross(grouped_mediator_losses, plot_file_type),
    format = "file"
  ),
  tar_target(
    outcome_glance_plot,
    plot_outcome_glance(grouped_outcome_losses, file_type = plot_file_type),
    pattern = cross(grouped_outcome_losses, plot_file_type),
    format = "file"
  ),
  tar_target(
    average_elementwise_outcome_loss_plot,
    plot_outcome_elementwise_loss(grouped_outcome_losses, file_type = plot_file_type),
    pattern = cross(grouped_outcome_losses, plot_file_type),
    format = "file"
  ),
  tar_target(
    causal_loss_plot,
    plot_causal_loss(all_causal_losses, file_type = plot_file_type),
    pattern = map(plot_file_type),
    format = "file"
  ),
  tar_target(
    null_causal_loss_plot,
    plot_null_causal_loss(all_causal_losses, file_type = plot_file_type),
    pattern = map(plot_file_type),
    format = "file"
  )
)
