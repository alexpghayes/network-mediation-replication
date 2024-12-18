library(targets)
library(tarchetypes)
library(here)
library(tibble)
library(crew)

tar_option_set(
  controller = crew_controller_local(workers = 8),
  packages = c("glue", "broom", "tidygraph", "latentnetmediate", "here", "ggplot2", "dplyr", "vegan", "estimatr", "tidyr", "wordspace"),
  imports = c("latentnetmediate", "fastRG")  # changes to these packages will invalidate pipeline,
)

source(here("R", "simulation-helpers.R"))
source(here("R", "simulation-plotting.R"))
source(here("R", "misspecification-helpers.R"))

models <- c(
  "model_uninformative",
  "model_block"
)

static_branch_targets <- tar_map(

  unlist = FALSE,

  values = tibble(model = rlang::syms(models)),

  tar_target(
    population,
    purrr::map(1:chunk_size, ~model(n = params$n, k = params$rank)),
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
    estimates,
    purrr::map(
      tbl_graph,
      latentnetmediate::sensitivity_curve,
      formula = y ~ . - name + 0,
      max_rank = 30,
      ranks_to_consider = 29
    ),
    iteration = "list",
    pattern = map(tbl_graph)
  ),

  # collect losses within each chunk into a data frame
  # note that this does not track which chunk each loss comes from,
  # but we don't actually need that information
  tar_target(
    losses_by_chunk,
    purrr::map2_dfr(estimates, population, misspecification_loss, .id = "rep_within_chunk"),
    pattern = map(estimates, population)
  )
)

list(

  tar_target(
    plot_file_type,
    c("png", "pdf")
  ),

  tar_target(chunk_size, 10),
  tar_target(num_chunks, 10),

  tar_target(n, c(500, 1000)),
  tar_target(rank, c(5, 10, 15)),

  tar_target(chunk_indices, 1:num_chunks),

  tar_target(
    params,
    tibble::tibble(
      n = n,
      rank = rank
    ),
    pattern = cross(n, rank)
  ),

  static_branch_targets,

  tar_combine(
    losses,
    static_branch_targets[[4]],
    command = bind_rows(!!!.x)
  ),

  tar_target(
    misspecification_loss_plot,
    plot_misspecification_loss(
      losses,
      num_sims_per_paramset = chunk_size * num_chunks,
      file_type = plot_file_type,
      dpi = 500
    ),
    pattern = map(plot_file_type),
    format = "file"
  )
)
