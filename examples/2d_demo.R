library(pacman)
p_load(tidyverse,spatstat,future,promises,listenv,
       ggplot2,gganimate,ggridges,av)
library(MathBioSim)

x_grid_birth <- x_grid_death <- seq(0,10,length.out = 1001)

sim_params <-
  list("area_length_x"=5, 
       "area_length_y"=5, 
       "cell_count_x"=50,
       "cell_count_y"=50, 
       "periodic"=TRUE,
       
       "b"=1,    
       "d"=0,    
       "dd"=0.01, 
       
       "seed"=1235,  
       "initial_population_x"=c(0.1),
       "initial_population_y"=c(0.1),
       
       "death_kernel_r"=5,
       "death_kernel_y"=dnorm(x_grid_death, sd = 1),
       
       "birth_kernel_r"=5,
       "birth_kernel_y"=dnorm(x_grid_birth, sd = 1)*x_grid_birth, 
       
       "spline_precision" = 1e-9  
  )

sim<-new(poisson_2d,sim_params)
