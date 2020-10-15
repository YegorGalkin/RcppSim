# Assuming everything is installed
library(pacman)
p_load(tidyverse, future, promises, listenv, fOptions)
library(MathBioSim)

# Prepare directory for results
result_dir = './14_10_2020_sims_2d/'
dir.create(result_dir, showWarnings = FALSE)
dir.create(paste0(result_dir,"pop/"), showWarnings = FALSE)

# Multiprocessing at simulation run level - single process for single parameter set
plan(multiprocess)

n_dim=2

# This code block prepares parameters for simulations
# Check params_all

grid_points <- expand.grid(comp_strength = seq(0.1,10,by=0.1),
                           comp_width = 10^seq(-1,1,by=0.1))

params_all <- grid_points%>%
  mutate(b=1,d=0,sm=1,sw=comp_width)%>%
  mutate(dd=comp_strength*(sw*sw*2*pi)^(n_dim/2))%>%
  mutate(kernel_trim = 3, spline_precision = 1e-9, spline_nodes=1001L)%>%
  mutate(death_kernel_r=kernel_trim*sw,birth_kernel_r=kernel_trim*sm)%>%
  mutate(id=row_number(),seed=1234L)%>%
  mutate(initial_population = 1e3L, population_limit = 1e5L)%>%
  mutate(area_length_x=10*pmax(sm,sw),periodic=TRUE,cell_count_x=10L)%>%
  mutate(n_samples = 100L)

mv_normal_radius<-function(r,sd,n_dim){
  1/(2*pi*sd*sd)^(n_dim/2)*exp(-r^2/(2*sd^2))
}

# Prepares environment for multiprocess results

all_runs = listenv()

# Checks if any simulations are already done - if simulator is killed, it retrieves 
# already calculated data

params_done <- list.files(paste0(result_dir,"pop/"))%>%str_remove('.csv')%>%as.numeric()

for (i in setdiff(params_all$id,params_done)%>%.[sample(length(.))]) {
  params=params_all[i,]
  all_runs[[i]]%<-%
  {
    x_grid_death = seq(0,params$kernel_trim * params$sw, length.out = params$spline_nodes)
    x_grid_birth = seq(0,params$kernel_trim * params$sm, length.out = params$spline_nodes)
    
    initial_points = fOptions::runif.halton(params$initial_population,n_dim) * params$area_length_x
    
    sim_params <-
      list("area_length_x"=params$area_length_x, 
           "area_length_y"=params$area_length_x, 
           
           "cell_count_x"=params$cell_count_x,  
           "cell_count_y"=params$cell_count_x, 
           
           "periodic"=params$periodic,
           
           "b"=params$b,    
           "d"=params$d,    
           "dd"=params$dd, 
           
           "seed"=params$seed,  
           "initial_population_x"=initial_points[,1],
           "initial_population_y"=initial_points[,2],
           
           "population_limit" = params$population_limit,
           
           "death_kernel_r"=params$death_kernel_r,
           "death_kernel_y"=mv_normal_radius(r=x_grid_death, sd=params$sw, n_dim=n_dim),
           
           "birth_kernel_r"=params$birth_kernel_r,
           "birth_kernel_y"=dnorm(x_grid_birth, sd = params$sm), 
           
           "spline_precision" = params$spline_precision 
      )
    
    sim<-new(poisson_2d,sim_params)
    time<-numeric(params$n_samples)
    pop<-numeric(params$n_samples)
    pop_cap_reached<-numeric(params$n_samples)
    
    for(j in 1:params$n_samples){
      sim$run_events(sim$total_population)
      pop[j]=sim$total_population
      time[j]=sim$time
      cap_reached=sim$pop_cap_reached
    }
    
    pops<-data.frame(id=i,time=time,pop=pop,cap_reached=cap_reached)
    
    write_csv(pops,paste0(result_dir,"pop/",i,".csv"))
  }
}

# Wait for completion
all_runs%>%as.list()
# Prepare summary results table and plots

write_csv(params_all,paste0(result_dir,"params.csv"))

list.files(paste0(result_dir,"pop/"),full.names = TRUE)%>%
  map_dfr(read_csv)%>%
  write_csv(paste0(result_dir,"pop.csv"))

summary_results <-
  read_csv(paste0(result_dir,"pop.csv"))%>%
  group_by(id)%>%
  arrange(id,time)%>%
  slice_tail(prop=0.8)%>%
  summarise(average_pop=mean(pop),cap_reached=any(cap_reached))%>%
  left_join(params_all%>%select(id,area_length_x,comp_strength,comp_width))%>%
  mutate(N = average_pop / area_length_x^n_dim)%>%
  select(id,comp_strength,comp_width,N,everything()) 

summary_results%>%write_csv(paste0(result_dir,"summary.csv"))

ggplot(summary_results, 
       aes(x=log10(comp_width), 
           y=comp_strength, 
           z = log10(N),
           fill=log10(N))) +
  geom_tile() +
  geom_contour() + 
  scale_fill_distiller(palette = "Spectral")
