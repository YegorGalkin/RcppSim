library(pacman)
p_load(tidyverse)
library(MathBioSim)

sim <- initialize_simulator(
  area_length_x = 100,
  dd = 0.01,
  initial_population_x = c(50),
  death_r = 5,
  death_y = dnorm(seq(0, 5, length.out = 1001), sd = 1),
  birth_ircdf_y = qnorm(seq(0.5, 1 - 1e-6, length.out = 101), sd = 0.2),
  realtime_limit = 10,
  ndim = 1
)

epochs <- 100
sim_results <- run_simulation(sim, epochs, calculate.pcf = TRUE)

dir.create("plots", showWarnings = FALSE)

ggplot(sim_results$pattern, aes(x = x)) +
  geom_histogram(bins = 100, breaks = seq(0, 100, length.out = 101))
ggsave("plots/pattern_histogram.png") # гистограмма численности популяции по x

ggplot(sim_results$population, aes(x = time, y = pop)) +
  geom_line()
ggsave("plots/population_dynamics.png") # динамика численности популяции

ggplot(sim_results$pcf, aes(x = r, y = pcf)) +
  geom_line()
ggsave("plots/pair_correlation.png") # функция парного корреляционного радиуса

ggplot(sim_results$K, aes(x = r, y = K_iso)) +
  geom_line()
ggsave("plots/k_function.png") # функция K


make_animation(
  sim_results = sim_results,
  simulator = sim,
  epochs = epochs,
  output_path = "plots/population_1d.mp4"
)
