library(pacman)
p_load(tidyverse)
library(MathBioSim)

x_grid_birth <- x_grid_death <- seq(0, 3, length.out = 1001)

sim_params <-
  list("area_length_x"=10,
       "area_length_y"=10, 
       "area_length_z"=10, 
       
       "cell_count_x"=10,
       "cell_count_y"=10,  
       "cell_count_z"=10,  
       "periodic"=TRUE,
       
       "b"=1,    
       "d"=0.1,    
       "dd"=0.03, 
       
       "seed"=123,  
       "initial_population_x"=c(0.1),
       "initial_population_y"=c(0.5),
       "initial_population_z"=c(0.9),
       
       "population_limit" = 1e5L,
       
       "death_kernel_r"=3,
       "death_kernel_y"=dnorm(x_grid_death, sd = 1),
       
       "birth_kernel_r"=3,
       "birth_kernel_y"=dnorm(x_grid_birth, sd = 1), 
       
       "spline_precision" = 1e-9 
  )

sim <- new(poisson_3d,sim_params)
