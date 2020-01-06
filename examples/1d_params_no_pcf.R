library(pacman)
p_load(tidyverse,spatstat,future,promises,listenv)
library(MathBioSim)

result_dir = 'D:/Rstuff/SimulationRuns/9_27_2019_pops/'
dir.create(result_dir, showWarnings = FALSE)
dir.create(paste0(result_dir,"pop/"), showWarnings = FALSE)

plan(multiprocess)

params_file = read_csv("D:/Rstuff/params.csv")

n_samples = 1000
initial_population = 10000
time_limit = 3600

params_all=data.frame(id=params_file$id,
                      sm=params_file$sm,
                      sw=params_file$sw,
                      b=1,d=0,dd=0.01,
                      samples=n_samples,
                      start_pop=initial_population,
                      seed=1234, 
                      area=pmax(params_file$sm, params_file$sw) * 100)

all_runs = listenv()

for (i in 1:nrow(params_all)) {
  params=params_all[i,]
  all_runs[[i]]%<-%
  {
    start_time = Sys.time()
    params$death_kernel_r = 10 * params$sw
    params$death_kernel_nodes = 1001
    x_grid_death = seq(0,params$death_kernel_r, 
                       length.out = params$death_kernel_nodes)
    
    
    params$birth_kernel_r = 10 * params$sm
    params$birth_kernel_nodes = 1001
    x_grid_birth = seq(0,params$birth_kernel_r, 
                       length.out = params$birth_kernel_nodes)
    
    params$area_length_x = params$area
    params$init_density = params$start_pop / params$area_length_x
    
    sim_params <-
      list("area_length_x"=params$area_length_x, 
           "cell_count_x"=100,  
           
           "b"=params$b,    
           "d"=params$d,    
           "dd"=params$dd, 
           
           "seed"=params$seed,  
           "init_density"=params$init_density,
           
           "death_kernel_r"=params$death_kernel_r,
           "death_kernel_y"=dnorm(x_grid_death, sd = params$sw),
           
           "birth_kernel_r"=params$birth_kernel_r,
           "birth_kernel_y"=dnorm(x_grid_birth, sd = params$sm), 
           
           "spline_precision" = 1e-9  
      )
    
    pcf_grid = seq(0,max(c(params$sw,params$sm))*10,length.out = 1001)
    
    sim<-new(poisson_1d,sim_params)
    time<-numeric(n_samples)
    pop<-numeric(n_samples)
    calculated_limit = n_samples
    
    for(j in 1:n_samples){
      sim$run_events(sim$total_population)
      pop[j]=sim$total_population
      time[j]=sim$time
      if (Sys.time()-start_time>time_limit){
        calculated_limit = j
        break
      }
    }
    
    pops<-data.frame(id=i,time=time,pop=pop)
    
    write_csv(pops,paste0(result_dir,"pop/",i,".csv"))
  }%stdout%TRUE
}

all_runs%>%as.list()

write_csv(params_all,paste0(result_dir,"params.csv"))

list.files(paste0(result_dir,"pop/"),full.names = TRUE)%>%
  map_dfr(read_csv)%>%
  write_csv(paste0(result_dir,"pop.csv"))


