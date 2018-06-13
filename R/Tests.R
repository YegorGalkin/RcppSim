if(0){

require(devtools)
require(Rcpp)
load_all()

params<-list("area_length_x"=1,    
             "cell_count_x"=100,  
             "b"=0.3,    
             "d"=0,    
             "dd"=0.01,   
             "seed"=1234,  
             "init_density"=100,
             "death_kernel_x"=(1:100)/20000,
             "death_kernel_y"=dnorm((1:100)/20000,sd=0.001), 
             "birth_kernel_x"=(1:100)/10000,
             "birth_kernel_y"=dnorm((1:100)/10000,sd=0.002), 
             "spline_precision" = 1-1e-6
             )
sim<-new(poisson_1d,params)

values<-c()
for(i in 1:1000){
  values<-c(values,sim$birth_reverse_cdf_spline_at(((1:1000)/100000)[[i]]))
}


qnorm((500:1000)/1000,sd=0.002)
}



