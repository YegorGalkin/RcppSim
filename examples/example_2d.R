library(pacman)
p_load(tidyverse)
library(MathBioSim)

sim <- initialize_simulator(
  area_length_x = 15,
  area_length_y = 15,
  dd = 0.01,
  seed = 1235,
  initial_population_x = c(7.5),
  initial_population_y = c(7.5),
  death_r = 5,
  death_y = dnorm(seq(0, 5, length.out = 1001), sd = 1),
  birth_ircdf_y = qnorm(seq(0.5, 1 - 1e-6, length.out = 101), sd = 0.2),
  realtime_limit = 10,
  ndim = 2
)

epochs <- 100
sim_results <- run_simulation(sim, epochs, calculate.pcf = TRUE)

dir.create("plots", showWarnings = FALSE)

ggplot(sim_results$pattern, aes(x = x, y = y)) +
  geom_bin2d(bins = 30) +
  scale_fill_viridis_c() +
  labs(x = "X coordinate", y = "Y coordinate", fill = "Count") +
  coord_fixed(ratio = 1) +
  theme_minimal() +
  theme(
    plot.background = element_rect(fill = "white", color = NA),
  )
ggsave("plots/pattern_density_2d.png") # 2d гистограмма численности популяции


make_animation(
  sim_results = sim_results,
  simulator = sim,
  epochs = epochs,
  output_path = "plots/population_1d.mp4"
)
