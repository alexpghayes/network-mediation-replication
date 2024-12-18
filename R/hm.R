clean_cahmped_data <- function(raw_path) {
  raw_path |>
    read_csv() |>
    select(
      id,
      group, # treatment assignment
      # contains("timestamp"),
      redcap_event_name,
      contains("nih_lonliness_q0"), # mediator indicators
      contains("mlq_10_q"), # mediator indicators
      contains("ffmq_8"), # 8 item-level mediators,
      contains("dds_10_q"), # mediator
      promis_bank_v10_depression_tscore, # outcome 1
      promis_bank_v10_anxiety_tscore, # outcome 2,
      dem_sex_assigned_birth,
      pre_age
    ) |>
    mutate(
      intervention = case_match(
        group,
        1 ~ "Control",
        2 ~ "Meditation"
      ),
      depression = promis_bank_v10_depression_tscore,
      anxiety = promis_bank_v10_anxiety_tscore,
      event_name = case_match(
        redcap_event_name,
        "pre_arm_1" ~ "Week 0",
        "w1_arm_1" ~ "Week 1",
        "w2_arm_1" ~ "Week 2",
        "w3_arm_1" ~ "Week 3",
        "post_arm_1" ~ "Week 4",
        "3mo_arm_1" ~ "Week 16"
      ),
      event_name_fct = fct_inorder(event_name),
      age = pre_age,
      sex = dem_sex_assigned_birth
    ) |>
    group_by(id) |>
    fill(
      age, sex, intervention,
      .direction = "downup"
    ) |>
    ungroup() |>
    filter(!is.na(event_name)) |>
    select(
      id, contains("event_name"), intervention, age, sex, anxiety, depression, everything(),
      -group, -contains("bank_v10"), -pre_age, -dem_sex_assigned_birth
    )
}

make_ate_figure <- function(data) {
  plot <- data |>
    lm(anxiety ~ intervention * event_name_fct, data = _) |>
    augment(interval = "confidence") |>
    ggplot() +
    aes(
      x = event_name_fct,
      ymin = .lower,
      y = anxiety,
      ymax = .upper,
      color = intervention,
      fill = intervention,
      group = intervention
    ) +
    geom_jitter(
      position = position_dodge2(width = 0.3),
      # width = 0.25,
      # alpha = 0.5
    ) +
    geom_ribbon(
      alpha = 0.75
    ) +
    # annotate(
    #   geom = "text",
    #   x = 2,
    #   y = 30.2,
    #   label = "High anxiety",
    #   color = "#444444"
    # ) +
    # annotate(
    #   geom = "text",
    #   x = 2.1,
    #   y = 83,
    #   label = "Low anxiety",
    #   color = "#444444"
    # ) +
    labs(
      y = "Anxiety level (PROMIS)",
      fill = NULL,
      color = NULL,
      title = NULL
    ) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_minimal(
      base_size = 12,
      base_family = ""
    ) +
    theme(
      axis.title.x = element_blank()
    )

  path <- here("figures", "healthyminds", "ate.pdf")

  ggsave(
    path,
    plot = plot,
    height = 2.75,
    width = 5.75,
    dpi = 300
  )

  path
}


make_week4_scatter_figure <- function(data) {
  week4 <- data |>
    filter(event_name == "Week 4")

  plot <- week4 |>
    ggplot() +
    aes(
      x = intervention,
      y = anxiety,
      color = intervention
    ) +
    geom_jitter(
      position = position_dodge2(width = 0.3),
    ) +
    labs(
      title = "Week 4",
      y = "Anxiety level (PROMIS)"
    ) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_minimal(
      base_size = 12,
      base_family = ""
    ) +
    theme(
      axis.title.x = element_blank(),
      legend.position = "none"
    )

  plot

  path <- here("figures", "healthyminds", "week4_scatter.png")

  ggsave(
    path,
    plot = plot,
    height = 2.5,
    width = 2.25,
    dpi = 300
  )

  path
}


make_week4_responses_figure <- function(A) {
  plot <- A |>
    as.matrix() |>
    as_tibble(rownames = "name") |>
    pivot_longer(
      -name,
      names_to = "item",
      values_to = "response"
    ) |>
    mutate(
      survey = case_match(
        substr(item, 1, 3),
        "nih" ~ "Loneliness",
        "mlq" ~ "Purpose",
        "dds" ~ "Defusion",
        "ffm" ~ "Awareness"
      )
    ) |>
    ggplot(aes(item, name, fill = as.factor(response))) +
    geom_raster() +
    scale_fill_viridis_d() +
    labs(
      fill = "Response",
      x = glue("Survey question (m = {ncol(A)})"),
      y = glue("Participant (n = {nrow(A)})")
    ) +
    facet_grid(
      cols = vars(survey),
      scales = "free_x"
    ) +
    theme_minimal(
      base_size = 11,
      base_family = ""
    ) +
    theme(
      axis.text = element_blank(),
      axis.ticks = element_blank()
    )

  path <- here("figures", "healthyminds", "week4-responses.png")

  ggsave(
    path,
    plot = plot,
    height = 2.75,
    width = 5.75,
    dpi = 600
  )

  path
}


make_xhat_figure <- function(fa) {
  plot <- fa |>
    get_svd_u() |>
    set_names(nm = c("participant", paste0("Xhat", 1:fa$rank))) |>
    ggplot(aes(Xhat1 * sqrt(fa$d[1]), Xhat2 * sqrt(fa$d[2]))) +
    geom_point() +
    theme_minimal(
      base_size = 11,
      base_family = ""
    ) +
    labs(
      x = "Xhat1",
      y = "Xhat2"
    ) +
    theme(
      axis.text = element_blank()
    )

  path <- here("figures", "healthyminds", "xhat.png")

  ggsave(
    path,
    plot = plot,
    height = 3.5,
    width = 4 * 16 / 9,
    dpi = 300
  )

  path
}


make_yhat_figure <- function(fa) {
  plot <- fa |>
    get_varimax_y() |>
    set_names(nm = c("item", paste0("Lhat", 1:fa$rank))) |>
    pivot_longer(
      contains("Lhat")
    ) |>
    mutate(
      survey = case_match(
        substr(item, 1, 3),
        "nih" ~ "Loneliness",
        "mlq" ~ "Purpose",
        "dds" ~ "Defusion",
        "ffm" ~ "Awareness"
      ),
      question_num = substr(item, nchar(item) - 1, nchar(item)),
      item_name = paste(survey, question_num),
      symbol = case_match(
        name,
        "Lhat1" ~ "widehat(F)[1]",
        "Lhat2" ~ "widehat(F)[2]",
        "Lhat3" ~ "widehat(F)[3]",
        "Lhat4" ~ "widehat(F)[4]",
        "Lhat5" ~ "widehat(F)[5]"
      ),
      hint = case_match(
        name,
        "Lhat1" ~ "Purposeful",
        "Lhat2" ~ "Purposeless",
        "Lhat3" ~ "Distraction",
        "Lhat4" ~ "Defusion",
        "Lhat5" ~ "Loneliness"
      ),
      name_hinted = glue("atop({symbol}, {hint})")
    ) |>
    ggplot(aes(name_hinted, item_name, fill = value)) +
    geom_tile() +
    scale_fill_gradient2() +
    scale_x_discrete(labels = label_parse()) +
    labs(
      fill = "Loading"
    ) +
    theme_classic(
      base_size = 10.5,
      base_family = ""
    ) +
    theme(
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )

  path <- here("figures", "healthyminds", "yhat.pdf")

  ggsave(
    path,
    plot = plot,
    height = 4.5,
    width = 5.5,
    dpi = 600
  )

  path
}

make_vhat_figure <- function(fa) {
  plot <- fa |>
    get_svd_v() |>
    set_names(nm = c("item", paste0("Vhat", 1:fa$rank))) |>
    pivot_longer(
      contains("Vhat")
    ) |>
    ggplot(aes(name, item, fill = value)) +
    geom_raster() +
    scale_fill_gradient2() +
    labs(
      fill = "Loading"
    ) +
    theme_classic(
      base_size = 9.5,
      base_family = ""
    ) +
    theme(
      axis.title = element_blank(),
      axis.ticks = element_blank()
    )

  path <- here("figures", "healthyminds", "vhat.png")

  ggsave(
    path,
    plot = plot,
    height = 2.9,
    width = 5.35,
    dpi = 300
  )

  path
}

make_eigcv_plots <- function(ecv, event) {
  plot <- plot(ecv) +
    labs(
      title = glue("Z-scores for cross-validated eigs for {event}")
    )

  path <- here("figures", "healthyminds", glue("ecv-{event}.png"))

  ggsave(
    path,
    plot = plot,
    height = 3.5,
    width = 3.5,
    dpi = 500
  )

  path
}

make_rank_estimate_figures <- function(A, ecv) {
  pca <- prcomp(A)

  vars <- pca$sdev^2
  prop <- cumsum(vars / sum(vars))

  plot <- ecv$summary |>
    mutate(
      cum_var_prop = prop[1:10]
    ) |>
    filter(k <= 6) |>
    pivot_longer(
      cols = c(z, cum_var_prop)
    ) |>
    mutate(
      name = if_else(name == "z", "Z-statistic", "Variance\nexplained"),
      name = as.factor(name),
      name = relevel(name, "Z-statistic")
    ) |>
    ggplot() +
    aes(k, value, color = name) +
    geom_line() +
    geom_point() +
    scale_x_continuous(
      breaks = scales::pretty_breaks()
    ) +
    scale_color_brewer(palette = "Dark2") +
    labs(
      caption = "Cross-validated eigenvalue method selects d = 2",
      x = "Rank"
    ) +
    facet_grid(
      rows = vars(name),
      scales = "free"
    ) +
    theme_minimal(
      base_size = 12,
      base_family = ""
    ) +
    theme(
      axis.title.y = element_blank(),
      legend.position = "none"
    )


  path <- here("figures", "healthyminds", "rank-determination.png")

  ggsave(
    path,
    plot = plot,
    height = 2.75,
    width = 5.25,
    dpi = 400
  )

  path
}

make_regression_figures <- function(nested_by_event) {
  A <- nested_by_event$A[[5]]
  fa <- nested_by_event$fa[[5]]

  Xhat <- fa$Z %*% fa$B |>
    as.matrix() |>
    as.data.frame() |>
    set_names(paste0("Xhat", 1:5)) |>
    mutate(
      name = rownames(A)
    )

  node_data <- nested_by_event$node_data[[5]]
  merged <- left_join(node_data, Xhat, by = "name")

  o_fit <- merged |>
    select(anxiety, intervention, age, sex, Xhat1:Xhat5) |>
    lm(
      anxiety ~ .,
      data = _
    )

  m_fit <- lm(
    cbind(Xhat1, Xhat2, Xhat3, Xhat4, Xhat5) ~ intervention + age + sex,
    data = merged
  )

  tidy_m_fit <- tidy(m_fit, conf.int = TRUE) |>
    mutate(
      term = str_replace(term, "interventionMeditation", "Meditation"),
      regression = "Mediator"
    ) |>
    filter(term == "Meditation")

  tidy_o_fit <- tidy(o_fit, conf.int = TRUE) |>
    mutate(
      term = str_replace(term, "interventionMeditation", "Meditation"),
      response = "Anxiety",
      regression = "Outcome"
    ) |>
    filter(term %in% c("Xhat1", "Xhat2", "Xhat3", "Xhat4", "Xhat5", "Meditation"))

  tidy_product <- tidy_o_fit |>
    left_join(tidy_m_fit, by = c("term" = "response")) |>
    mutate(
      estimate = if_else(term == "Meditation", estimate.x, estimate.x * estimate.y),
      std.error = if_else(term == "Meditation", std.error.x, sqrt(std.error.x^2 * estimate.y^2 + std.error.y^2 * estimate.x^2)),
      regression = "Product"
    ) |>
    select(term, estimate, std.error, regression)

  plot <- tidy_m_fit |>
    mutate(
      term = response
    ) |>
    bind_rows(tidy_o_fit, tidy_product) |>
    mutate(
      is_product = as.factor(if_else(regression != "Product", "Coefficients", "Product")),
      is_product = fct_relevel(is_product, "Coefficients", "Product"),
      symbol = case_match(
        term,
        "Xhat1" ~ "widehat(X)[1]",
        "Xhat2" ~ "widehat(X)[2]",
        "Xhat3" ~ "widehat(X)[3]",
        "Xhat4" ~ "widehat(X)[4]",
        "Xhat5" ~ "widehat(X)[5]",
        "Meditation" ~ "phantom(M)"
      ),
      hint = case_match(
        term,
        "Xhat1" ~ "Purposeful",
        "Xhat2" ~ "Purposeless",
        "Xhat3" ~ "Distraction",
        "Xhat4" ~ "Defusion",
        "Xhat5" ~ "Loneliness",
        "Meditation" ~ "Meditation"
      ),
      term_hinted = glue("{symbol} ~~ {hint}")
    ) |>
    ggplot() +
    aes(
      x = term_hinted,
      ydist = dist_normal(mean = estimate, sd = std.error),
      color = regression
    ) +
    stat_pointinterval(
      position = position_dodge(width = 1),
      .width = c(0.5, 0.80, 0.95)
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "darkgrey"
    )  +
    coord_flip() +
    scale_color_manual(
      values = c(
        "Mediator" = "firebrick",
        "Outcome" = "steelblue",
        "Product" = "darkgreen"
      ),
      name = ""
    ) +
    scale_x_discrete(labels = label_parse()) +
    facet_grid(
      cols = vars(is_product),
      scales = "free_x"
    ) +
    theme_minimal(
      base_size = 11,
      base_family = ""
    ) +
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "top"
    )

  path <- here("figures", "healthyminds", "coefficients.pdf")

  ggsave(
    path,
    plot = plot,
    height = 3,
    width = 5.75,
    dpi = 300
  )

  path
}

make_regression_figures_depression <- function(nested_by_event) {
  A <- nested_by_event$A[[5]]
  fa <- nested_by_event$fa[[5]]

  Xhat <- fa$Z %*% fa$B |>
    as.matrix() |>
    as.data.frame() |>
    set_names(paste0("Xhat", 1:5)) |>
    mutate(
      name = rownames(A)
    )

  node_data <- nested_by_event$node_data[[5]]
  merged <- left_join(node_data, Xhat, by = "name")

  o_fit <- merged |>
    select(depression, intervention, age, sex, Xhat1:Xhat5) |>
    lm(
      depression ~ .,
      data = _
    )

  m_fit <- lm(
    cbind(Xhat1, Xhat2, Xhat3, Xhat4, Xhat5) ~ intervention + age + sex,
    data = merged
  )

  tidy_m_fit <- tidy(m_fit, conf.int = TRUE) |>
    mutate(
      term = str_replace(term, "interventionMeditation", "Meditation"),
      regression = "Mediator"
    ) |>
    filter(term == "Meditation")

  tidy_o_fit <- tidy(o_fit, conf.int = TRUE) |>
    mutate(
      term = str_replace(term, "interventionMeditation", "Meditation"),
      response = "Depression",
      regression = "Outcome"
    ) |>
    filter(term %in% c("Xhat1", "Xhat2", "Xhat3", "Xhat4", "Xhat5", "Meditation"))

  tidy_product <- tidy_o_fit |>
    left_join(tidy_m_fit, by = c("term" = "response")) |>
    mutate(
      estimate = if_else(term == "Meditation", estimate.x, estimate.x * estimate.y),
      std.error = if_else(term == "Meditation", std.error.x, sqrt(std.error.x^2 * estimate.y^2 + std.error.y^2 * estimate.x^2)),
      regression = "Product"
    ) |>
    select(term, estimate, std.error, regression)

  plot <- tidy_m_fit |>
    mutate(
      term = response
    ) |>
    bind_rows(tidy_o_fit, tidy_product) |>
    mutate(
      is_product = as.factor(if_else(regression != "Product", "Coefficients", "Product")),
      is_product = fct_relevel(is_product, "Coefficients", "Product"),
      symbol = case_match(
        term,
        "Xhat1" ~ "widehat(X)[1]",
        "Xhat2" ~ "widehat(X)[2]",
        "Xhat3" ~ "widehat(X)[3]",
        "Xhat4" ~ "widehat(X)[4]",
        "Xhat5" ~ "widehat(X)[5]",
        "Meditation" ~ "phantom(M)"
      ),
      hint = case_match(
        term,
        "Xhat1" ~ "Purposeful",
        "Xhat2" ~ "Purposeless",
        "Xhat3" ~ "Distraction",
        "Xhat4" ~ "Defusion",
        "Xhat5" ~ "Loneliness",
        "Meditation" ~ "Meditation"
      ),
      term_hinted = glue("{symbol} ~~ {hint}")
    ) |>
    ggplot() +
    aes(
      x = term_hinted,
      ydist = dist_normal(mean = estimate, sd = std.error),
      color = regression
    ) +
    stat_pointinterval(
      position = position_dodge(width = 1),
      .width = c(0.5, 0.80, 0.95)
    ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "darkgrey"
    )  +
    coord_flip() +
    scale_color_manual(
      values = c(
        "Mediator" = "firebrick",
        "Outcome" = "steelblue",
        "Product" = "darkgreen"
      ),
      name = ""
    ) +
    scale_x_discrete(labels = label_parse()) +
    facet_grid(
      cols = vars(is_product),
      scales = "free_x"
    ) +
    theme_minimal(
      base_size = 11,
      base_family = ""
    ) +
    theme(
      axis.title.y = element_blank(),
      axis.title.x = element_blank(),
      legend.position = "top"
    )

  path <- here("figures", "healthyminds", "coefficients-depression.pdf")

  ggsave(
    path,
    plot = plot,
    height = 3,
    width = 5.75,
    dpi = 300
  )

  path
}

make_sensitivity_figure <- function(curve) {
  plot <- curve |>
    mutate(
      term = str_replace(term, "interventionMeditation", "Meditation"),
      estimand = toupper(estimand)
    ) |>
    filter(term %in% c("Meditation")) |>
    ggplot() +
    aes(
      x = rank,
      ymin = conf.low,
      y = estimate,
      ymax = conf.high,
      color = estimand,
      fill = estimand
    ) +
    geom_point() +
    geom_line() +
    geom_ribbon(alpha = 0.3) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
    labs(
      x = "Number of latent factors",
      y = "Estimate"
    ) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_minimal(
      base_family = "",
      base_size = 11
    ) +
    theme(
      legend.title = element_blank()
    )

  path <- here("figures", "healthyminds", "curve.pdf")

  ggsave(
    path,
    plot = plot,
    height = 2.25,
    width = 5,
    dpi = 300
  )

  path
}

make_mediation_trajectory_figure <- function(nested_by_event) {
  plot <- nested_by_event |>
    mutate(
      effects = map(medanx, pluck, "effects")
    ) |>
    select(
      event_name_fct, effects
    ) |>
    unnest(c(effects)) |>
    mutate(
      estimand = toupper(estimand)
    ) |>
    filter(term == "interventionMeditation") |>
    ggplot() +
    aes(
      x = event_name_fct,
      ymin = conf.low,
      y = estimate,
      ymax = conf.high,
      color = estimand,
      group = estimand,
      fill = estimand
    ) +
    geom_ribbon(alpha = 0.3) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(
      y = "Effect"
    ) +
    theme_minimal() +
    theme(
      base_size = 11,
      base_family = "",
      axis.title.x = element_blank(),
      legend.title = element_blank()
    )

  path <- here("figures", "healthyminds", "mediation-trajectory.pdf")

  ggsave(
    path,
    plot = plot,
    height = 2.75,
    width = 5.25,
    dpi = 300
  )

  path
}


make_mediation_trajectory_figure_depression <- function(nested_by_event) {
  plot <- nested_by_event |>
    mutate(
      effects = map(meddep, pluck, "effects")
    ) |>
    select(
      event_name_fct, effects
    ) |>
    unnest(c(effects)) |>
    mutate(
      estimand = toupper(estimand)
    ) |>
    filter(term == "interventionMeditation") |>
    ggplot() +
    aes(
      x = event_name_fct,
      ymin = conf.low,
      y = estimate,
      ymax = conf.high,
      color = estimand,
      group = estimand,
      fill = estimand
    ) +
    geom_ribbon(alpha = 0.3) +
    geom_line() +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey") +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(
      y = "Effect"
    ) +
    theme_minimal() +
    theme(
      base_size = 11,
      base_family = "",
      axis.title.x = element_blank(),
      legend.title = element_blank()
    )

  path <- here("figures", "healthyminds", "mediation-trajectory-depression.pdf")

  ggsave(
    path,
    plot = plot,
    height = 2.75,
    width = 5.25,
    dpi = 300
  )

  path
}



make_latent_positivity_plot <- function(nested_by_event) {
  fa <- nested_by_event$fa[[5]]
  graph <- nested_by_event$tbl_graph[[5]]

  U_df <- get_svd_u(fa)

  plot <- graph |>
    filter(!type) |>
    as_tibble() |>
    left_join(U_df, by = c("name" = "id")) |>
    mutate(Condition = intervention) |>
    select(Condition, matches("u[0-9]+")) |>
    ggpairs(
      mapping = aes(color = Condition, fill = Condition, alpha = 0.5),
      columns = colnames(fa$u),
      columnLabels = paste0("widehat(X)[", 1:fa$rank, "]"),
      labeller = "label_parsed",
      legend = c(2, 1),
      upper = "blank",
      axisLabels = "none"
    ) +
    scale_color_brewer(type = "qual") +
    scale_fill_brewer(type = "qual") +
    scale_alpha_continuous(guide = "none") +
    theme_minimal()

  path <- here("figures", "healthyminds", "positivity.pdf")

  ggsave(
    path,
    plot = plot,
    dpi = 300,
    width = 10,
    height = 10
  )

  path
}
