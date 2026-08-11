
library(dplyr)
library(ggdist)
library(tidyr)
library(tibble)

sim_table <- list()

for (sim in 1:100) {

    data_source <- sprintf("./data/sim_cc_%d.RData", sim)              # contains sd_x1/sd_x2
  data_res    <- sprintf("./results/result_sim_cc_%03d.rds", sim)    # contains beta/gamma

  if (!file.exists(data_res)) {
    warning(sprintf("Missing result file: %s (skipping sim=%d)", data_res, sim))
    next
  }
  if (!file.exists(data_source)) {
    warning(sprintf("Missing data file: %s (skipping sim=%d)", data_source, sim))
    next
  }

  # load scaling SDs
  load(data_source)   # expects sd_x1, sd_x2 saved in the .RData
  results <- readRDS(data_res)

  # ---- posterior draws selection ----
  idx <- seq(10001, 60000, by = 5)
  sampT <- length(idx)

  # ---- build transformed parameters ----
  # theta rows correspond to:
  # 1 = AgeDifference (beta[2])
  # 2 = SpatialDistance (beta[3])
  # 3 = CombinedAge / combined age effect (gamma[1])
  theta <- matrix(NA_real_, nrow = 3, ncol = sampT)
  rownames(theta) <- c("AgeDifference", "SpatialDistance", "CombinedAge")

  theta["AgeDifference", ]   <- results$beta[2, idx] 
  theta["SpatialDistance", ] <- results$beta[3, idx] 
  theta["CombinedAge", ]          <- results$gamma[1, idx]

  # ---- convert to OR draws ----
  draws_or <- exp(theta)

  # ---- summarize (median + 95% interval) ----
  cov_fin <- as.data.frame(t(draws_or)) %>%
    pivot_longer(cols = everything(), names_to = "key", values_to = "value")

  sum_single <- cov_fin %>%
    group_by(key) %>%
    ggdist::median_hdci(value, .width = 0.95) %>%
    rename(
      value = value,
      .lower = .lower,
      .upper = .upper
    ) %>%
    mutate(source = "SingleRun")

  # ---- truth + bias/coverage ----
  truth <- tibble(
    key = c("AgeDifference", "CombinedAge", "SpatialDistance"),
    true_value = c(0.752, 0.754, 0.516) # the reference value from the benchmark dataset
  )

  comparison <- sum_single %>%
    left_join(truth, by = "key") %>%
    mutate(
      bias = value - true_value,
      coverage = true_value >= .lower & true_value <= .upper,
      ci_width = .upper - .lower
    ) %>%
    select(key, source, value, true_value, bias, coverage, ci_width)

  # add simulation id
  comparison$simulation <- sim
  sim_table[[sim]] <- comparison
}

final_bias <- bind_rows(sim_table)
write.csv2(final_bias, "final_bias_update_new.csv", row.names = FALSE)