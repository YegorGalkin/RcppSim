library(pacman)
p_load(tidyverse)
library(MathBioSim)

sim<-initialize_simulator(area_length_x = 100, dd=0.01,
                          initial_population_x = c(10),
                          death_r = 5,
                          death_y = dnorm(seq(0,5,length.out = 1001), sd = 1),
                          birth_ircdf_y = qnorm(seq(0.5,1-1e-6,length.out = 101), sd = 0.2),
                          realtime_limit = 60)

sim_results <- run_simulation(sim, 1000, calculate.pcf = TRUE)

ggplot(sim_results$pattern,aes(x=x))+
  geom_histogram(bins=100,breaks=seq(0,100,length.out = 101))
  
ggplot(sim_results$population,aes(x=time,y=pop))+
  geom_line()

ggplot(sim_results$pcf,aes(x=r,y=pcf))+
  geom_line()

ggplot(sim_results$K,aes(x=r,y=K_iso))+
  geom_line()
