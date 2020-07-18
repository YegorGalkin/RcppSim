library(pacman)
p_load(tidyverse, future, promises, listenv)
library(MathBioSim)

result_dir = './7_18_2020_sims/'
dir.create(result_dir, showWarnings = FALSE)
dir.create(paste0(result_dir,"pop/"), showWarnings = FALSE)

plan(multiprocess)

grid_points <- expand.grid(comp_strength = seq(0.1,10,by=0.1),
                           comp_width = 10^seq(-2,1,by=0.1))

params_all <- grid_points%>%
  mutate(b=1,d=0,sm=1,sw=comp_width)%>%
  mutate(dd=comp_strength*sw*sqrt(2*pi))%>%
  mutate(kernel_trim = 5, spline_precision = 1e-9, spline_nodes=1001L)%>%
  mutate(death_kernel_r=kernel_trim*sw,birth_kernel_r=kernel_trim*sm)%>%
  mutate(id=row_number(),seed=1234L)%>%
  mutate(initial_population = 1e4L, population_limit = 1e5L)%>%
  mutate(area_length_x=100*pmax(sm,sw),periodic=TRUE,cell_count_x=100L)%>%
  mutate(n_samples = 100L)



all_runs = listenv()

for (i in params_all$id) {
  params=params_all[i,]
  all_runs[[i]]%<-%
  {
    x_grid_death = seq(0,params$kernel_trim * params$sw, length.out = params$spline_nodes)
    x_grid_birth = seq(0,params$kernel_trim * params$sm, length.out = params$spline_nodes)
    
    sim_params <-
      list("area_length_x"=params$area_length_x, 
           "cell_count_x"=params$cell_count_x,  
           "periodic"=params$periodic,
           
           "b"=params$b,    
           "d"=params$d,    
           "dd"=params$dd, 
           
           "seed"=params$seed,  
           "initial_population_x"=seq(0.01,params$area_length_x-0.01, length.out = params$initial_population),
           "population_limit" = params$population_limit,
           
           "death_kernel_r"=params$death_kernel_r,
           "death_kernel_y"=dnorm(x_grid_death, sd = params$sw),
           
           "birth_kernel_r"=params$birth_kernel_r,
           "birth_kernel_y"=dnorm(x_grid_birth, sd = params$sm), 
           
           "spline_precision" = params$spline_precision 
      )
    
    sim<-new(poisson_1d,sim_params)
    time<-numeric(params$n_samples)
    pop<-numeric(params$n_samples)
    
    
    for(j in 1:params$n_samples){
      sim$run_events(sim$total_population)
      pop[j]=sim$total_population
      time[j]=sim$time
    }
    
    pops<-data.frame(id=i,time=time,pop=pop)
    
    write_csv(pops,paste0(result_dir,"pop/",i,".csv"))
    sim$pop_cap_reached
  }
}

pop_cap_reached <- all_runs%>%as.list()%>%unlist()

write_csv(params_all,paste0(result_dir,"params.csv"))

list.files(paste0(result_dir,"pop/"),full.names = TRUE)%>%
  map_dfr(read_csv)%>%
  write_csv(paste0(result_dir,"pop.csv"))

read_csv(paste0(result_dir,"pop.csv"))%>%
  group_by(id)%>%
  arrange(id,time)%>%
  slice_tail(prop=0.9)%>%
  summarise(average_pop=mean(pop))%>%
  left_join(params_all%>%select(id,area_length_x,comp_strength,comp_width))%>%
  mutate(cap_reached = pop_cap_reached, N = average_pop / area_length_x)%>%
  select(id,comp_strength,comp_width,N,cap_reached,everything())%>%
  write_csv(paste0(result_dir,"summary.csv"))
