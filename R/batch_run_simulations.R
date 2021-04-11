#' Runs multiple simulations
#'
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
#'
#' @return
#'
#' @examples
batch_run_simulations <- 
  function(area_length_x,
           b=1,d=0,dd,sw,sm=1,
           work_dir=getwd(),
           calculate.pcf=FALSE,
           pcf.r = 5*pmax(sm,sw),
           pcf.nodes = 101,
           ndim=1,
           periodic=TRUE,
           seed=12345,
           initial_population=1000,
           realtime_limit=3600,
           birth_tail_cutoff=1e-6,
           kernel_nodes=101,
           cell_count_x=100,
           n_cores=availableCores()
           ){
    
  }
