library(pacman)
p_load(tidyverse,future,listenv)
library(MathBioSim)

initial_population = list()

# Set up multisession parallel execution environment with availableCores() parallel processes
plan(multisession(workers = availableCores()))

# Collects results from different processes
sim_results = listenv()

simulation <- function(area_length_x = 100, 
                cell_count_x = 100, periodic = TRUE, 
                b = 1, d = 0, dd = 0, seed = 1234L, 
                initial_population_x = c(10),
                death_r = 1, death_y, birth_ircdf_y, realtime_limit = 1e6, ndim = 1, epochs = 1000, 
                calculate.pcf = TRUE, pcf_grid) {
  # Stores expression as a future, future gets send to different process
  # Futures eat errors so carefully debug code
  # Different starting populations for example
  # Increase area length and runtime for better PCF results
  sim_results[[i]]%<-%{
    sim<-initialize_simulator(area_length_x = area_length_x, 
                                cell_count_x = cell_count_x,
                                periodic = periodic,
                                b = b, d = d, dd = dd, seed = seed, 
                                initial_population_x = initial_population_x,
                                death_r = death_r,
                                death_y = death_y, birth_ircdf_y = birth_ircdf_y, realtime_limit = realtime_limit, ndim = ndim)
    
    res<-run_simulation(sim, epochs = epochs, calculate.pcf = calculate.pcf, pcf_grid = pcf_grid)
    res[['initial_parameters']]<-data.frame(b=b,d=d,dd=dd)
    # Last expression is a return value of a future
    res
  }
}


i = 0L

b_l = 1.0
b_r = 1.1
b_step = 0.1

d_l = 0
d_r = 0.1
d_step = 0.1

dd_l = 0.001
dd_r = 0.0011
dd_step = 0.0001


start_time <- Sys.time()
for (b_value in seq(b_l, b_r, b_step)) {
  for (d_value in seq(d_l, d_r, d_step)) {
    for (dd_value in seq(dd_l, dd_r, dd_step)) {
      i = i + 1L
      simulation(area_length_x = 10,
          b = b_value,
          d = d_value,
          dd = dd_value,
          initial_population_x = seq(0, 10, length.out = 100), death_r = 5, 
          death_y = dnorm(seq(0, 5, length.out = 1001), sd = 1), 
          birth_ircdf_y = qnorm(seq(0.5, 1-1e-6, length.out = 101), sd = 0.2), 
          realtime_limit = 60, calculate.pcf = TRUE, 
          pcf_grid = seq(0,5,length.out = 1001))
    }
  }
}
finish_time <- Sys.time()
time_taken <- finish_time - start_time
time_taken
count_of_simulations = i
mean_time_per_simulation <- time_taken / count_of_simulations
mean_time_per_simulation

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
  #summarise(pcf=mean(pcf))%>% # Average across 10 runs
  ggplot(aes(x=r,y=pcf,color=ID))+
  geom_line()

# Same for K function
sim_results%>%
  purrr::map(~.x[['K']])%>%
  bind_rows(.id='ID')%>% 
  group_by(r)%>%
  #summarise(K_iso=mean(K_iso))%>%
  ggplot(aes(x=r,y=K_iso,color=ID))+
  geom_line()


work_dir=file.path(getwd(),'sim')
for (i in seq(1, count_of_simulations)) {
  # Exporting to csv
  write_csv(do.call(rbind.data.frame, sim_results[[i]]['pcf']),
            file.path(work_dir,"pcf",paste0(i,'.csv')))
  write_csv(do.call(rbind.data.frame, sim_results[[i]]['initial_parameters']),
            file.path(work_dir,"initial_parameters",paste0(i,'.csv')))
}
