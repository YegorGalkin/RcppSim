library(pacman)
p_load(tidyverse, future, promises, listenv, fOptions)
library(MathBioSim)

params_all <- data.frame(sm = 0.1,
                         sw = 10^seq(-1,1,length.out = 30),
                         b=1,d=0,dd=0.1)%>%
  mutate(kernel_trim = 5, spline_precision = 1e-9, spline_nodes=1001L)%>%
  mutate(death_kernel_r=kernel_trim*sw,birth_kernel_r=kernel_trim*sm)%>%
  mutate(id=row_number(),seed=1234L)%>%
  mutate(initial_population = 100, population_limit = 30e3L)%>%
  mutate(area_length_x=50*pmax(sm,sw),periodic=TRUE,cell_count_x=100L)%>%
  mutate(n_samples = 2000L)

plan(multiprocess)
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
    pop_cap_reached<-numeric(params$n_samples)
    
    for(j in 1:params$n_samples){
      sim$run_events(sim$total_population)
      pop[j]=sim$total_population
      time[j]=sim$time
      cap_reached=sim$pop_cap_reached
    }
    
    data.frame(id=i,time=time,pop=pop,cap_reached=cap_reached)
  }
}

# Wait for completion
results<-all_runs%>%as.list()%>%bind_rows()

ggplot(results%>%dplyr::filter(id %in% c(1,5,10,15,20,25,30)))+
  geom_line(aes(x=time,y=log2(pop),color=as.factor(id)))

ggplot(results)+
  geom_line(aes(x=time,y=log2(pop),color=as.factor(id)))

results%>%
  arrange(id,time)%>%
  group_by(id)%>%
  slice_tail(prop=0.9)%>%
  summarise(mean_pop=mean(pop))%>%
  full_join(params_all)%>%
  mutate(N=mean_pop/area_length_x)%>%
  select(id,mean_pop,N,everything())%>%
  write_csv('D:/Rstuff/mathbio/MBSruns/2_march_2021_1d.csv')
