library(pacman)
p_load(tidyverse,spatstat,future,promises,listenv,
       ggplot2,gganimate,ggridges,av)
library(MathBioSim)

x_grid_birth <- x_grid_death <- seq(0,10,length.out = 1001)

sim_params <-
  list("area_length_x"=10, 
       "cell_count_x"=100,  
       "periodic"=TRUE,
       
       "b"=1,    
       "d"=0,    
       "dd"=0.01, 
       
       "seed"=1235,  
       "initial_population_x"=c(0.1),
       "population_limit" = 100,
       
       "death_kernel_r"=5,
       "death_kernel_y"=dnorm(x_grid_death, sd = 1),
       
       "birth_kernel_r"=5,
       "birth_kernel_y"=dnorm(x_grid_birth, sd = 1), 
       
       "spline_precision" = 1e-9  
  )

sim<-new(poisson_1d,sim_params)

n_times = 100
population_list=list()
population_list[[1]]=0.1

for(j in 2:n_times){
  sim$run_for(1)
  population_list[[j]]=sim$get_all_coordinates()
}

data <- enframe(population_list)%>%
  unnest(value)%>%
  rename(time_point=name,coordinate=value)

animation <- data%>%
  bind_rows(data%>%
              filter(coordinate<2)%>%
              mutate(coordinate=coordinate+10))%>%
  bind_rows(data%>%
              filter(coordinate>8)%>%
              mutate(coordinate=coordinate-10))%>%
  ggplot(aes(x=coordinate,y='')) +
  scale_x_continuous(limits = c(-2, 12), expand = c(0, 0))+
  geom_vline(xintercept = 0)+
  geom_vline(xintercept = 10)+
  geom_density_ridges(
    jittered_points = TRUE,
    panel_scaling=FALSE,
    position = position_points_jitter(width = 0.05, height = 0),
    point_shape = '|', point_size = 3, point_alpha = 0.1, alpha = 0.7
  )+
  labs(x='Area',y='Population density')+
  transition_time(time_point)

animate(animation,renderer=av_renderer("population_periodic.mp4"),fps=5)
