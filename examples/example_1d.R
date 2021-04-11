library(pacman)
p_load(tidyverse)
library(MathBioSim)

sim<-initialize_simulator(area_length_x = 100,
                          cell_count_x = 100,
                          periodic = TRUE,
                          b=1,d=0,dd=0.01,
                          seed=1234,
                          initial_population_x = c(10),
                          death_r = 5,
                          death_y = dnorm(seq(0,5,length.out = 1001)),
                          birth_ircdf_y = qnorm(seq(0.5,1-1e-6,length.out = 101),sd = 0.2),
                          realtime_limit = 60)
sim$run_events(1000000)

ggplot(data.frame(x=sim$get_all_x_coordinates()),aes(x=x))+
  geom_histogram(bins=100,breaks=seq(0,100,length.out = 101))
  
