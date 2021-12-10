library(pacman)
p_load(tidyverse,future,listenv)
library(MathBioSim)

initial_population = list()

# Set up multiprocess parallel execution environment with availableCores() parallel processes
plan(multiprocess(workers = availableCores()))

# Collects results from different processes
sim_results = listenv()

simulation <- function(area_length_x = 100, area_length_y = 100, area_length_z = 100, 
                cell_count_x = 100, cell_count_y = 100, cell_count_z = 100, periodic = TRUE, 
                b = 1, d = 0, dd = 0, seed = 1234L, 
                initial_population_x = c(10), initial_population_y = c(10), initial_population_z = c(10),
                death_r = 1, death_y, birth_ircdf_y, realtime_limit = 1e6, ndim = 1, epochs = 1000, 
                calculate.pcf = FALSE, pcf_grid) {
  # Stores expression as a future, future gets send to different process
  # Futures eat errors so carefully debug code
  # Different starting populations for example
  # Increase area length and runtime for better PCF results
  sim_results[[i]]%<-%{
    sim<-initialize_simulator(area_length_x = area_length_x, area_length_y = area_length_y, 
                                area_length_z = area_length_z, cell_count_x = cell_count_x,
                                cell_count_y = cell_count_y, cell_count_z = cell_count_z, periodic = periodic,
                                b = b, d = d, dd = dd, seed = seed, 
                                initial_population_x = initial_population_x, initial_population_y = initial_population_y,
                                initial_population_z = initial_population_z, death_r = death_r,
                                death_y = death_y, birth_ircdf_y = birth_ircdf_y, realtime_limit = realtime_limit, ndim = ndim)
    # Last expression is a return value of a future
    run_simulation(sim, epochs = epochs, calculate.pcf = calculate.pcf, pcf_grid = pcf_grid)
  }
}

N_l = 1
N_r = 10
step = 1

for (i in seq(N_l, N_r, step)) {
  simulation(area_length_x = 10, dd = 0.01, initial_population_x = seq(0, 10, length.out = 100*i), death_r = 5, 
      death_y = dnorm(seq(0, 5, length.out = 1001), sd = 1), 
      birth_ircdf_y = qnorm(seq(0.5, 1-1e-6, length.out = 101), sd = 0.2), 
      realtime_limit = 60, calculate.pcf = TRUE, 
      pcf_grid = seq(0,5,length.out = 1001))
}


# Waits for all tasks to finish, converts results to list
sim_results <- sim_results%>%as.list()

sim_results%>%
  purrr::map(~.x[['population']])%>% # Pulls population from every result
  bind_rows(.id='ID')%>%             # Concatenates them and adds ID row
  ggplot(aes(x=time,y=pop,color=ID))+
  geom_line()+
  xlim(0,10)                         # Shows different starting positions


sim_results%>%
  purrr::map(~.x[['pcf']])%>%
  bind_rows(.id='ID')%>% 
  group_by(r)%>%
  summarise(pcf=mean(pcf))%>% # Average across 10 runs
  ggplot(aes(x=r,y=pcf))+
  geom_line()


# Same for K function
sim_results%>%
  purrr::map(~.x[['K']])%>%
  bind_rows(.id='ID')%>% 
  group_by(r)%>%
  summarise(K_iso=mean(K_iso))%>%
  ggplot(aes(x=r,y=K_iso))+
  geom_line()
