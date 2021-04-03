library(pacman)
p_load(tidyverse, future, promises, listenv, fOptions)
library(MathBioSim)

normal_radius<-function(r,sd,n_dim){
  1/(2*pi*sd*sd)^(n_dim/2)*exp(-r^2/(2*sd^2))
}

params_all <- data.frame(sm = 0.1,
                         sw = 10^seq(-1,1,length.out = 30),
                         b=1,d=0,dd=0.1)%>%
  mutate(kernel_trim = 5, spline_precision = 1e-9, spline_nodes=1001L)%>%
  mutate(death_kernel_r=kernel_trim*sw,birth_kernel_r=kernel_trim*sm)%>%
  mutate(id=row_number(),seed=1234L)%>%
  mutate(initial_population = 100, population_limit = 30e3L)%>%
  mutate(area_length_x=10*pmax(sm,sw),periodic=TRUE,cell_count_x=100L)%>%
  mutate(n_samples = 200L)


plan(multiprocess(workers = availableCores()-1))
all_runs = listenv()
n_dim=2

for (i in params_all$id) {
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
           "death_kernel_y"=normal_radius(r=x_grid_death, sd=params$sw, n_dim=n_dim),
           
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
    
    data.frame(id=i,time=time,pop=pop,cap_reached=cap_reached)
  }
}

# Wait for completion
results<-all_runs%>%as.list()%>%bind_rows()

ggplot(results)+
  geom_line(aes(x=time,y=log2(pop),color=as.factor(id)))

results%>%
  arrange(id,time)%>%
  group_by(id)%>%
  dplyr::filter(max(time)>50)%>%
  slice_tail(prop=0.5)%>%
  summarise(mean_pop=mean(pop))%>%
  full_join(params_all)%>%
  mutate(N=mean_pop/area_length_x^n_dim)%>%
  select(id,mean_pop,N,everything())%>%
  write_csv('D:/Rstuff/mathbio/MBSruns/2_march_2021_2d.csv')
