#' Simplified simulator runner with animation support
#'
#' @param simulator 1-species simulator object (1d, 2d or 3d)
#' @param epochs amount of population-dependent time units to run
#' @param calculate.pcf if TRUE, adds pcf estimate to output
#' @param pcf_grid grid for pcf calculation
#'
#' @return list with results and animation data
#' @export
#'
run_simulation <- function(simulator, epochs, calculate.pcf = FALSE, pcf_grid) {
  require(dplyr)
  require(glue)

  if (dir.exists("plots")) {
    unlink("plots", recursive = TRUE) # delete old plots
  }
  dir.create("plots", showWarnings = FALSE) # create new plots directory

  time <- numeric(epochs + 1)
  pop <- numeric(epochs + 1)
  frames <- list(x = vector("list", epochs + 1),
                 y = vector("list", epochs + 1),
                 z = vector("list", epochs + 1))


  time[1] <- simulator$time
  pop[1] <- simulator$total_population

  frames$x[[1]] <- simulator$get_all_x_coordinates()
  if (class(simulator)[1] == "Rcpp_poisson_2d") {
    frames$y[[1]] <- simulator$get_all_y_coordinates()
  } else if (class(simulator)[1] == "Rcpp_poisson_3d") {
    frames$y[[1]] <- simulator$get_all_y_coordinates()
    frames$z[[1]] <- simulator$get_all_z_coordinates()
  }

  for (j in 2:(epochs + 1)) {
    simulator$run_events(simulator$total_population)

    time[j] <- simulator$time
    pop[j] <- simulator$total_population

    frames$x[[j]] <- simulator$get_all_x_coordinates()
    if (class(simulator)[1] == "Rcpp_poisson_2d") {
      frames$y[[j]] <- simulator$get_all_y_coordinates()
    } else if (class(simulator)[1] == "Rcpp_poisson_3d") {
      frames$y[[j]] <- simulator$get_all_y_coordinates()
      frames$z[[j]] <- simulator$get_all_z_coordinates()
    }
  }

  if (simulator$total_population == 0) {
    return(
      list(
        "realtime_limit_reached" = simulator$realtime_limit_reached,
        "population" = data.frame(time = time, pop = pop) %>% distinct()
      )
    )
  }

  if (class(simulator)[1] == "Rcpp_poisson_1d") {
    pattern <- data.frame(
      x = simulator$get_all_x_coordinates()
    )
  } else if (class(simulator)[1] == "Rcpp_poisson_2d") {
    pattern <- data.frame(
      x = simulator$get_all_x_coordinates(),
      y = simulator$get_all_y_coordinates()
    )
  } else if (class(simulator)[1] == "Rcpp_poisson_3d") {
    pattern <- data.frame(
      x = simulator$get_all_x_coordinates(),
      y = simulator$get_all_y_coordinates(),
      z = simulator$get_all_z_coordinates()
    )
  }

  result <- list(
    "realtime_limit_reached" = simulator$realtime_limit_reached,
    "population" = data.frame(time = time, pop = pop) %>% distinct(),
    "pattern" = pattern,
    "frames" = frames
  )

  if (calculate.pcf) {
    require(spatstat)
    if (class(simulator)[1] == "Rcpp_poisson_1d") {
      points <- spatstat.geom::ppp(
        simulator$get_all_x_coordinates(),
        rep(0, simulator$total_population),
        c(0, simulator$area_length_x),
        c(-simulator$area_length_x / 2, simulator$area_length_x / 2)
      )

      if (missing(pcf_grid)) {
        K_estimate <- spatstat.explore::Kest(points, correction = "Ripley")
        pcf_grid <- K_estimate$r
      } else {
        K_estimate <- spatstat.explore::Kest(points, r = pcf_grid, correction = "Ripley") # nolint
      }

      pcf_estimate <- data.frame(Kest = K_estimate$iso / 2, x = pcf_grid) %>%
        mutate(pfc = (.data$Kest - lag(.data$Kest)) / (pcf_grid - lag(pcf_grid)) / simulator$area_length_x) %>% # nolint
        pull(.data$pfc)

    } else if (class(simulator)[1] == "Rcpp_poisson_2d") {

      points <- spatstat.geom::ppp(
        simulator$get_all_x_coordinates(),
        simulator$get_all_y_coordinates(),
        c(0, simulator$area_length_x),
        c(0, simulator$area_length_y)
      )

      if (missing(pcf_grid)) {
        K_estimate <- spatstat.explore::Kest(points, correction = "Ripley")
        pcf_grid <- K_estimate$r
      } else {
        K_estimate <- spatstat.explore::Kest(points, r = pcf_grid, correction = "Ripley") # nolint
      }

      pcf_estimate <- spatstat.explore::pcf(K_estimate, method = "b")$pcf

    } else if (class(simulator)[1] == "Rcpp_poisson_3d") {

      points <- spatstat.geom::pp3(
        simulator$get_all_x_coordinates(),
        simulator$get_all_y_coordinates(),
        simulator$get_all_z_coordinates(),
        c(0, simulator$area_length_x),
        c(0, simulator$area_length_y),
        c(0, simulator$area_length_z)
      )

      if (missing(pcf_grid)) {
        K_estimate <- spatstat.explore::Kest(points, correction = "Ripley")
        pcf_grid <- K_estimate$r
      } else {
        K_estimate <- spatstat.explore::Kest(points, r = pcf_grid, correction = "Ripley") # nolint
      }

      pcf_estimate <- spatstat.explore::pcf(K_estimate, method = "b")$pcf

    } else {
      stop(glue(
        "simulator is not supported, object should be poisson_1d, 2d or 3d",
        cl = class(simulator)[1]
      ))
    }

    result[["pcf"]] <- data.frame(r = pcf_grid, pcf = pcf_estimate)
    result[["K"]] <- data.frame(r = pcf_grid, K_iso = K_estimate$iso, K_poisson = K_estimate$theo) # nolint
  }

  return(result)
}
