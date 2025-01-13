library(targets)
library(here)
library(crew)

data(addhealth, package = "latentnetmediate")

tar_option_set(
  controller = crew_controller_local(workers = 4),
  packages = c(
    "glue",
    "broom",
    "tidygraph",
    "latentnetmediate",
    "here",
    "ggplot2",
    "dplyr",
    "igraph",
    "stringr",
    "ggraph"
  ),
  imports = c("latentnetmediate", "fastRG") # changes to these packages will invalidate pipeline
)

##### AddHealth overcontrol bias example ---------------------------------------

# use our cleaned data in the netmediate and roughly match the Li-Le pre-processing
#
#   - symmetrize graph, subset to largest weakly connected component
#   - impute missing Y to zero, mean impute grade, mean impute sex (dropped NAs)
#   - lump race categories < 0.05 percent (didn't do this)

create_leli_graph <- function() {
  addhealth[[36]] |>
    mutate(
      out_degree = centrality_degree(weights = weight, mode = "out"),
      out_friends = centrality_degree(mode = "out"),
      level = out_degree / out_friends
    ) |>
    as.undirected(mode = "collapse") |>
    as_tbl_graph() |>
    mutate(
      component = group_components(type = "weak"),
      race = relevel(race, "white")
    ) |>
    filter(component == 1, !is.na(level), !is.na(sex), !is.na(race), !is.na(grade)) |>
    select(-component, -out_degree, -out_friends) |>
    activate(edges) |>
    select(-weight) |>
    activate(nodes)
}

compute_rank_curve <- function(graph) {
  sensitivity_curve(graph,
    level ~ sex + race + grade + 0,
    max_rank = 250,
    ranks_to_consider = 25
  )
}

plot_rank_curve <- function(rank_curve, file_type, ...) {
  rank_curve |>
    filter(str_detect(term, "race")) |>
    mutate(term = str_remove(term, "race")) |>
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
      x = "Latent dimension of network",
      y = "Mental health relative to white students"
    ) +
    facet_wrap(
      vars(term),
      nrow = 2
    ) +
    theme_minimal(base_size = 16) +
    theme(
      legend.title = element_blank()
    )

  path <- here("figures", "addhealth", glue("effects.{file_type}"))

  ggsave(path,
    dpi = 600,
    width = 8,
    height = 8 * 9 / 16,
    create.dir = TRUE
  )

  path
}

list(
  tar_target(
    plot_file_type,
    c("png", "pdf")
  ),
  tar_target(
    leli_graph,
    create_leli_graph()
  ),
  tar_target(
    leli_rank_curve,
    compute_rank_curve(leli_graph)
  ),
  tar_target(
    rank_curve_plot,
    plot_rank_curve(
      leli_rank_curve,
      file_type = plot_file_type,
      height = 6,
      width = 6 * 16 / 9,
      dpi = 500
    ),
    pattern = map(plot_file_type),
    format = "file"
  )
)
