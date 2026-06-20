library(tidyverse)

mu <- 4.37e-9
gen <- 2

pops <- c("L01","L02","L03","L04","L05","L06",
          "L07","L08","L09","L10","L11","L12")

main_dir <- "/home/xiongh/2026/Fish/fitered/new_1/msmc2/results_singlepop"
boot_dir <- "/home/xiongh/2026/Fish/fitered/new_1/msmc2/results_singlepop_bootstrap"
out_dir  <- "/home/xiongh/2026/Fish/fitered/new_1/msmc2/Ne_summary"

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

read_msmc2 <- function(file, pop, type = "main", bootstrap = NA) {
  x <- read.table(file, header = TRUE)

  x %>%
    mutate(
      pop = pop,
      type = type,
      bootstrap = bootstrap,
      time_index = as.integer(time_index),
      time_years = left_time_boundary / mu * gen,
      Ne = 1 / (2 * mu * lambda)
    ) %>%
    filter(time_index >= 2)   # 去掉前两个 time segment
}

all_main <- list()
all_boot <- list()
Ne4000_main <- list()
Ne4000_boot <- list()

for (pop in pops) {

  message("Processing ", pop)

  main_file <- file.path(main_dir, paste0(pop, ".final.txt"))

  if (!file.exists(main_file)) {
    warning("Missing main MSMC2 result: ", main_file)
    next
  }

  main <- read_msmc2(main_file, pop, "main")

  all_main[[pop]] <- main

  ne4000 <- approx(
    x = main$time_years,
    y = main$Ne,
    xout = 4000,
    rule = 2
  )$y

  Ne4000_main[[pop]] <- tibble(
    pop = pop,
    Ne_4000_main = ne4000
  )

  boot_files <- file.path(
    boot_dir,
    paste0(pop, ".bootstrap_", 1:20, ".final.txt")
  )

  boot_existing <- boot_files[file.exists(boot_files)]

  if (length(boot_existing) > 0) {

    boot <- map_dfr(seq_along(boot_existing), function(i) {
      read_msmc2(
        file = boot_existing[i],
        pop = pop,
        type = "bootstrap",
        bootstrap = i
      )
    })

    all_boot[[pop]] <- boot

    ne4000_b <- boot %>%
      group_by(pop, bootstrap) %>%
      summarise(
        Ne_4000 = approx(time_years, Ne, xout = 4000, rule = 2)$y,
        .groups = "drop"
      )

    Ne4000_boot[[pop]] <- ne4000_b
  }
}

main_df <- bind_rows(all_main)
boot_df <- bind_rows(all_boot)

Ne4000_main_df <- bind_rows(Ne4000_main)
Ne4000_boot_df <- bind_rows(Ne4000_boot)

Ne4000_ci <- Ne4000_boot_df %>%
  group_by(pop) %>%
  summarise(
    Ne_4000_boot_median = median(Ne_4000, na.rm = TRUE),
    Ne_4000_boot_lower95 = quantile(Ne_4000, 0.025, na.rm = TRUE),
    Ne_4000_boot_upper95 = quantile(Ne_4000, 0.975, na.rm = TRUE),
    n_boot = n(),
    .groups = "drop"
  )

Ne4000_summary <- Ne4000_main_df %>%
  left_join(Ne4000_ci, by = "pop")

write.csv(main_df, file.path(out_dir, "MSMC2_main_all_pops_scaled.csv"), row.names = FALSE)
write.csv(boot_df, file.path(out_dir, "MSMC2_bootstrap_all_pops_scaled.csv"), row.names = FALSE)
write.csv(Ne4000_summary, file.path(out_dir, "MSMC2_Ne4000_summary.csv"), row.names = FALSE)

pdf(file.path(out_dir, "MSMC2_Ne_curves_all_pops.pdf"), width = 8, height = 6)

ggplot() +
  geom_step(
    data = boot_df,
    aes(x = time_years, y = Ne, group = interaction(pop, bootstrap)),
    alpha = 0.15
  ) +
  geom_step(
    data = main_df,
    aes(x = time_years, y = Ne),
    linewidth = 0.7
  ) +
  facet_wrap(~ pop, scales = "free_y") +
  scale_x_log10() +
  theme_bw() +
  labs(
    x = "Years before present",
    y = expression(paste("Effective population size ", italic(N[e]))),
    title = "MSMC2 demographic history"
  )

dev.off()

print(Ne4000_summary)
