# library(targets)
# library(netmediate)

# graph <- tar_read(clean)[[1]]

# node_data <- tidygraph::as_tibble(graph)

# rank <- 5
# A <- igraph::as_adjacency_matrix(graph, sparse = TRUE)
# X <- US(A, rank = rank)  # use scaled left singular vectors

# outcome_model <- stats::lm(
#   tobacco_int ~ sex_fct + age + leisure_church + X,
#   data = node_data
# )

# mediator_model <- stats::lm(
#   X ~ sex_fct * (age + leisure_church),
#   data = node_data
# )

# trt_chr <- "sex_fct"

# nde_table <- broom::tidy(outcome_model, conf.int = TRUE) |>
#   dplyr::filter(stringr::str_detect(term, trt_chr)) |>
#   dplyr::mutate(estimand = "nde") |>
#   dplyr::select(term, estimand, estimate, conf.low, conf.high)

# beta_hat <- stats::coef(outcome_model)
# dim_betaw <- length(beta_hat) - rank

# betaw_hat <- beta_hat[1:dim_betaw]
# betax_hat <- beta_hat[-c(1:dim_betaw)]

# Theta_hat <- stats::coef(mediator_model)

# intercept_index <- 1
# t_ind <- stringr::str_detect(rownames(Theta_hat), trt_chr)
# colon_ind <- stringr::str_detect(rownames(Theta_hat), ":")

# t_index <- which(xor(t_ind, colon_ind))
# tc_indices <- which(colon_ind)

# c_indices <- setdiff(
#   1:nrow(Theta_hat),
#   c(intercept_index, t_index, tc_indices)
# )

# thetat_hat <- Theta_hat[t_index, ]
# Thetatc_hat <- Theta_hat[tc_indices, ]



# nie_hat <- drop(thetat_hat %*% betax_hat + mu_c %*% Thetatc_hat %*% betax_hat)

# nie_table <- tibble::enframe(nie_hat, name = "term", value = "estimate") |>
#   dplyr::mutate(estimand = "nie") |>
#   dplyr::select(term, estimand, estimate)

# sigmabetax_hat <- stats::vcov(outcome_model)[-c(1:num_coefs), -c(1:num_coefs)]
# sigmatheta_hat <- stats::vcov(mediator_model)

# # need to re-arrange sigmatheta_hat from enormous square into something
# # more tensor-y / considering each covariate one at a time

# coef_names <- names(betaw_hat)

# # everything following is under the assumption that Theta_tc = 0
# # for convenience since it's a pain to handle Theta_tc != 0

# for (i in seq_along(coef_names)) {

#   nm <- coef_names[i]

#   indices <- which(
#     stringr::str_detect(
#       colnames(sigmatheta_hat),
#       stringr::coll(nm)
#     )
#   )

#   thetat_hat <- theta_hat[i, ]
#   sigmathetat_hat <- sigmatheta_hat[indices, indices, drop = FALSE]

#   # delta method
#   nie_var <- t(betax_hat) %*% sigmathetat_hat %*% betax_hat +
#     t(thetat_hat) %*% sigmabetax_hat %*% thetat_hat

#   nie_table[i, "conf.low"] <- nie_hat[i] - 1.96 * sqrt(nie_var)
#   nie_table[i, "conf.high"] <- nie_hat[i] + 1.96 * sqrt(nie_var)
# }

# effects <- dplyr::bind_rows(nde_table, nie_table) |>
#   dplyr::filter(!stringr::str_detect(term, "Intercept")) |>
#   dplyr::mutate(
#     term = stringr::str_replace(term, "W", "")
#   )

# object <- list(
#   formula = formula,
#   rank = rank,
#   outcome_model = outcome_model,
#   mediator_model = mediator_model,
#   effects = effects
# )

