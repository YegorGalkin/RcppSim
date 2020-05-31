library(pacman)
p_load(tidyverse,spatstat,future,promises,listenv,
       ggplot2,gganimate,ggridges,av)
library(MathBioSim)

x_grid_birth <- x_grid_death <- seq(0,5,length.out = 1001)

sim_params <-
  list("area_length_x"=15, 
       "area_length_y"=15, 
       "cell_count_x"=100,
       "cell_count_y"=100, 
       "periodic"=TRUE,
       
       "b"=1,    
       "d"=0,    
       "dd"=0.01, 
       
       "seed"=1235,  
       "initial_population_x"=c(0.05),
       "initial_population_y"=c(0.05),
       
       "death_kernel_r"=5,
       "death_kernel_y"=dnorm(x_grid_death, sd = 1),
       
       "birth_kernel_r"=5,
       "birth_kernel_y"=dnorm(x_grid_birth, sd = 1), 
       
       "spline_precision" = 1e-9  
  )

sim<-new(poisson_2d,sim_params)

n_times = 100
population_list_x=list()
population_list_x[[1]]=0.05

population_list_y=list()
population_list_y[[1]]=0.05
for(j in 2:n_times){
  sim$run_for(0.33)
  population_list_x[[j]]=sim$get_all_x_coordinates()
  population_list_y[[j]]=sim$get_all_y_coordinates()
}

data <- enframe(population_list_x)%>%
  unnest(value)%>%
  rename(x_coordinate=value)%>%
  bind_cols(enframe(population_list_y)%>%unnest(value)%>%select(value))%>%
  rename(y_coordinate=value,time_point=name)

animation <- data%>%
  ggplot(aes(x=x_coordinate,y=y_coordinate)) +
  geom_point(size=0.5,alpha=0.5)+
  labs(x='X axis',y='Y axis')+
  transition_time(time_point)

animate(animation,renderer=av_renderer("population_periodic_2d.mp4"),fps=5)
