#' Runs multiple simulations
#'
#' @description
#' Experimental
#' @param area_length_x 
#' @param b 
#' @param d 
#' @param dd 
#' @param sw 
#' @param sm 
#' @param work_dir 
#' @param calculate.pcf 
#' @param pcf.r 
#' @param pcf.nodes 
#' @param ndim 
#' @param periodic 
#' @param seed 
#' @param initial_population 
#' @param realtime_limit 
#' @param birth_tail_cutoff 
#' @param kernel_nodes 
#' @param cell_count_x 
#' @param n_cores 
#' @param epochs 
#' @param death_radius_cutoff 
#'
#' @return
#' @export
#'
#' @examples
batch_run_simulations <- 
  function(area_length_x,
           b=1,d=0,dd,sw,sm=1,
           work_dir=file.path(getwd(),Sys.time()),
           calculate.pcf=FALSE,
           pcf.r = 5*pmax(sm,sw),
           pcf.nodes = 101,
           ndim=1,
           periodic=TRUE,
           seed=12345,
           initial_population=1000,
           realtime_limit=3600,
           epochs=1000,
           birth_tail_cutoff=1e-6,
           death_radius_cutoff = 5*sw,
           kernel_nodes=101,
           cell_count_x=100,
           n_cores=availableCores()
           ){
    require(future)
    require(fOptions)
    require(VGAM)
    
    dir.create(work_dir, showWarnings = FALSE)
    dir.create(file.path(work_dir,"pop"), showWarnings = FALSE)
    dir.create(file.path(work_dir,"pcf"), showWarnings = FALSE)
    dir.create(file.path(work_dir,"K"), showWarnings = FALSE)
    dir.create(file.path(work_dir,"pattern"), showWarnings = FALSE)
    
    params_all <- data.frame(
      b=b,d=d,dd=dd,sw=sw,sm=sm,ndim=ndim,
      pcf.r=pcf.r,pcf.nodes=pcf.nodes,
      periodic=periodic,area_length_x=area_length_x,
      seed=seed,initial_population=initial_population,
      realtime_limit=realtime_limit,epochs=epochs,
      birth_tail_cutoff=birth_tail_cutoff,death_radius_cutoff=death_radius_cutoff,
      kernel_nodes=kernel_nodes,cell_count_x=cell_count_x
      )%>%
      dplyr::mutate(id=row_number())%>%
      dplyr::select(id,everything())
    
    plan(multisession(workers = n_cores))
    
    all_runs = listenv()
    
    for (i in params_all$id) {
      
      #Create asynchronous tasks and store them in all_runs list environment
      all_runs[[i]]%<-%
        {
          initial_points <- fOptions::runif.halton(params$initial_population[i],
                                                   params$n_dim[i]) * params$area_length_x[i]
          initial_points <-cbind(initial_points,initial_points,initial_points)
          
          birth_x_grid <- seq(0.5, 1-params$birth_tail_cutoff[i], length.out = params_all$kernel_nodes[i])
          if(params$n_dim[i]==1){
            birth_ircdf <- qnorm(birth_x_grid, sd = params_all$sm[i])
          }else if(params$n_dim[i]==2){
            birth_ircdf <- VGAM::qrayleigh(birth_x_grid,scale = params_all$sm[i])
          }else{
            birth_ircdf <- VGAM::qmaxwell(birth_x_grid,rate = params_all$sm[i])
          }
          
          sim<-initialize_simulator(
            area_length_x = params_all$area_length_x[i],
            area_length_y =params_all$area_length_x[i],
            area_length_z =params_all$area_length_x[i],
            cell_count_x=params_all$cell_count_x[i],
            cell_count_y=params_all$cell_count_x[i],
            cell_count_z=params_all$cell_count_x[i],
            periodic = params_all$periodic[i],
            b=params_all$b[i],
            d=params_all$d[i],
            dd=params_all$dd[i],
            seed=params_all$seed[i],
            initial_population_x=initial_points[,1],
            initial_population_y=initial_points[,2],
            initial_population_z=initial_points[,3],
            death_r=params_all$death_radius_cutoff[i],
            death_y=dnorm(seq(0,params_all$death_radius_cutoff[i],
                              length.out = params_all$kernel_nodes[i] ), 
                          sd = params_all$sw[i])^params_all$ndim[i],
            birth_ircdf_y=birth_ircdf,
            realtime_limit = params_all$realtime_limit[i],
            ndim = params_all$ndim[i]
          )
          
          results<-run_simulation(sim,
                                  epochs=params_all$epochs[i],
                                  calculate.pcf=TRUE,
                                  pcf_grid=seq(0,params_all$pcf.r[i],
                                               length.out = params_all$pcf.nodes[i])
                                  )
          write_csv(results$population,
                    file.path(work_dir,"pop",paste0(params_all$id[i],'.csv')))
          write_csv(results$pcf,
                    file.path(work_dir,"pcf",paste0(params_all$id[i],'.csv')))
          write_csv(results$K,
                    file.path(work_dir,"K",paste0(params_all$id[i],'.csv')))
          write_csv(results$pattern,
                    file.path(work_dir,"pattern",paste0(params_all$id[i],'.csv')))
          
          results$realtime_limit_reached
        }
    }
    
    # Wait for completion
    all_runs%>%as.list()
  }
