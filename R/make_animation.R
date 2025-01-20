#' Create and save animation from simulation results
#'
#' @param sim_results Results from run_simulation
#' @param simulator Simulator object for parameters
#' @param epochs Number of epochs
#' @param output_path Path to save animation
#' @param width Plot width
#' @param height Plot height
#' @param res Plot resolution
#'
#' @export
#'
make_animation <- function(sim_results, simulator, epochs,
                           output_path = "plots/population.mp4",
                           width = 800, height = 800, res = 150) {

  library(ggplot2)
  library(gganimate)
  library(tibble)
  library(dplyr)
  library(tidyr)
  library(av)

  value <- NULL
  name <- NULL
  x_coordinate <- NULL
  time_point <- NULL
  y_coordinate <- NULL

  if (class(simulator)[1] == "Rcpp_poisson_1d") {
    animation <- enframe(sim_results$frames$x) %>%
      unnest(value) %>%
      rename(x_coordinate = value, time_point = name) %>%
      ggplot(aes(x = x_coordinate)) +
      geom_histogram(bins = 100, breaks = seq(0, 100, length.out = 101)) +
      labs(x = "X axis", y = "Count") +
      transition_time(time_point) +
      theme_minimal() +
      coord_cartesian(
        xlim = c(0, simulator$area_length_x),
        ylim = c(0, 170)
      ) +
      theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)
      )

  } else if (class(simulator)[1] == "Rcpp_poisson_2d") {
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
        xlim = c(0, simulator$area_length_x),
        ylim = c(0, simulator$area_length_y)
      ) +
      theme(
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14)
      )

  } else if (class(simulator)[1] == "Rcpp_poisson_3d") {
    # TODO: implement 3d animation
  }

  anim_save(
    output_path,
    animation = animation,
    renderer = av_renderer(),
    width = width,
    height = height,
    fps = epochs / simulator$realtime_limit,
    res = res
  )

  system2("xdg-open", output_path)
}
