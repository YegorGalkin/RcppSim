library(pacman)
p_load(tidyverse,spatstat,future,promises,listenv,
       ggplot2,gganimate,ggridges,av)
library(MathBioSim)


x_grid_birth <- x_grid_death <- seq(0,10,length.out = 1001)

sim_params <-
  list("area_length_x"=1000, 
       "cell_count_x"=100,  
       "periodic"=FALSE,
       
       "b"=1,    
       "d"=0,    
       "dd"=0.01, 
       
       "seed"=123,  
       "initial_population"=c(500),
       
       "death_kernel_r"=5,
       "death_kernel_y"=dnorm(x_grid_death, sd = 1),
       
       "birth_kernel_r"=5,
       "birth_kernel_y"=dnorm(x_grid_birth, sd = 1), 
       
       "spline_precision" = 1e-9  
  )

sim<-new(poisson_1d,sim_params)

n_events = 10000
population_list=list()
population_list[[1]]=500
for(j in 2:10000){
  sim$run_events(1)
  population_list[[j]]=sim$get_all_coordinates()
}

animation <- 
  enframe(population_list)%>%
  unnest(value)%>%
  filter(name<100)%>%
  ggplot(aes(x=value,y='Population')) +
  geom_density_ridges(
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7
  ) + 
  transition_time(name)

animate(animation,renderer=av_renderer("population_start.mp4"))

animation <- 
  enframe(population_list)%>%
  unnest(value)%>%
  filter(name%%100==0)%>%
  ggplot(aes(x=value,y='Population')) +
  geom_density_ridges(
    jittered_points = TRUE,
    position = position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 3, point_alpha = 1, alpha = 0.7
  ) + 
  transition_time(name)

animate(animation, fps = 20, duration = 15,renderer=av_renderer("population_samples.mp4"))
