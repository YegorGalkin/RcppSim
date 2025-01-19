library(pacman)
p_load(tidyverse, spatstat, future, promises, listenv,
       ggplot2, gganimate, ggridges, av)

x_grid_birth <- x_grid_death <- seq(0, 5, length.out = 1001)

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

animation <- enframe(sim_results$frames$x) %>%
  unnest(value) %>%
  rename(x_coordinate = value) %>%
  bind_cols(
    enframe(sim_results$frames$y) %>%
      unnest(value) %>%
      select(value)
  ) %>%
  rename(y_coordinate = value, time_point = name) %>%

  ggplot(aes(x = x_coordinate, y = y_coordinate)) +
  geom_point(size = 0.5, alpha = 0.5) +
  labs(x = "X axis", y = "Y axis") +
  transition_time(time_point) +
  theme_minimal() +
  coord_fixed(ratio = 1) +
  coord_cartesian(
    xlim = c(0, 15),
    ylim = c(0, 15)
  ) +
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14)
  )

anim_save(
  "plots/population_2d.mp4",
  animation = animation,
  renderer = av_renderer(),
  width = 800,
  height = 800,
  fps = epochs / sim$realtime_limit,
  res = 150
)

# # 2D density plot
# ggplot(sim_results$pattern, aes(x = x, y = y)) +
#   geom_bin2d(bins = 30) +
#   scale_fill_viridis_c() +
#   labs(x = "X coordinate", y = "Y coordinate", fill = "Count") +
#   coord_fixed(ratio = 1) +
#   theme_minimal()
# ggsave("plots/pattern_density_2d.png")

system2("xdg-open", "plots/population_2d.mp4")
