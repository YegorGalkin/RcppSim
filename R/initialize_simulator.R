#' Creates space-time birth-death poisson point process simulation with given parameters.
#' 
#' @description
#' Length of simulation area should be large enough for birth and death kernels.
#' Requires birth rate to be greater than death rate.
#' Death kernel should be specified as a pair of max death interaction radius
#' and values of death kernel on uniform grid.
#' Birth kernel should be specified as a inverse cumulative distribution function 
#' (function of quantiles) of radius-distribution. In 1d case for normal distributions,
#' quantiles of normal distribution could be used. In 2d and 3d case, Rayleigh and Maxwell 
#' distribution quantiles are required.
#'
#' @param area_length_x Length of simulated area by x axis
#' @param area_length_y Length of simulated area by y axis
#' @param area_length_z Length of simulated area by z axis
#' @param cell_count_x Internal parameter for simulation speed-up
#' @param cell_count_y Internal parameter for simulation speed-up
#' @param cell_count_z Internal parameter for simulation speed-up
#' @param periodic Specifies whether to use periodic or killing boundary condition
#' @param b Birth rate
#' @param d Death rate
#' @param dd Competitive death rate
#' @param seed Simulation RNG seed
#' @param initial_population_x Coordinates of initial population, x axis
#' @param initial_population_y Coordinates of initial population, y axis
#' @param initial_population_z Coordinates of initial population, z axis
#' @param death_r Max radius of death interaction
#' @param death_y Values of death kernel on uniform grid from 0 to death_r
#' @param birth_ircdf_y Inverse radial cumulative distribution function used in birth simulation. 
#' Used to get random variable that corresponds to displacement distance on birth
#' @param realtime_limit Limit on simulation lifetime in seconds
#' @param ndim Dimension count, only 1, 2 and 3 supported
#'
#' @return Simulator object with methods for running
#' @export
#'
#' @examples
#' sim<-initialize_simulator(area_length_x = 100, dd=0.01,
#'                           initial_population_x = c(10),
#'                           death_r = 5,
#'                           death_y = dnorm(seq(0,5,length.out = 1001), sd = 1),
#'                           birth_ircdf_y = qnorm(seq(0.5,1-1e-6,length.out = 101), sd = 0.2),
#'                           realtime_limit = 60)
#' # Runs 1d simulator for a million events or 60 seconds maximum                          
#' sim$run_events(1e6)
#' # Shows population in the end of simulation
#' sim$total_population
initialize_simulator <- 
  function(area_length_x,area_length_y,area_length_z,
           cell_count_x=100,cell_count_y=100,cell_count_z=100,
           periodic=TRUE,
           b=1,d=0,dd,
           seed=1234L,
           initial_population_x,initial_population_y,initial_population_z,
           death_r,death_y,
           birth_ircdf_y,
           realtime_limit=1e6,
           ndim=1){
    
    stopifnot(ndim %in% c(1, 2, 3))
    stopifnot(b>=d)
    stopifnot(all(death_y>=0))
    stopifnot(all(birth_ircdf_y>=0))
    stopifnot(all(birth_ircdf_y<area_length_x))
    stopifnot(death_r<area_length_x)
    
    sim_params <-
      list("area_length_x"=area_length_x, 
           "cell_count_x"=cell_count_x,  
           
           "periodic"=periodic,
           
           "b"=b,    
           "d"=d,    
           "dd"=dd, 
           
           "seed"=seed,  
           "initial_population_x"=initial_population_x,
           
           "death_r"=death_r,
           "death_y"=death_y,
           
           "birth_ircdf_y"=birth_ircdf_y,
           
           "realtime_limit"=realtime_limit
      )
    
    if(ndim == 1){
      return(new(poisson_1d,sim_params))  
    }
    sim_params[['area_length_y']]<-area_length_y
    sim_params[['cell_count_y']]<-cell_count_y
    sim_params[['initial_population_y']]<-initial_population_y
    
    if(ndim == 2){
      return(new(poisson_2d,sim_params))  
    }
    sim_params[['area_length_z']]<-area_length_z
    sim_params[['cell_count_z']]<-cell_count_z
    sim_params[['initial_population_z']]<-initial_population_z
    if(ndim == 3){
      return(new(poisson_3d,sim_params))  
    } 
    
}
