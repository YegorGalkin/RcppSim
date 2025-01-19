library(pacman)
p_load(tidyverse, spatstat, future, promises, listenv,
       ggplot2, gganimate, av)
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


animation <- enframe(sim_results$frames$x) %>%
  unnest(value) %>%
  rename(x_coordinate = value, time_point = name) %>%
  ggplot(aes(x = x_coordinate)) +
  geom_histogram(bins = 100, breaks = seq(0, 100, length.out = 101)) +
  labs(x = "X axis", y = "Count") +
  transition_time(time_point) +
  theme_minimal() +
  coord_cartesian(
    xlim = c(0, 100),
    ylim = c(0, 170)
  ) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

anim_save(
  "plots/population_1d.mp4",
  animation = animation,
  renderer = av_renderer(),
  width = 800,
  height = 800,
  fps = epochs / sim$realtime_limit,
  res = 150
)

system2("xdg-open", "plots/population_1d.mp4")
