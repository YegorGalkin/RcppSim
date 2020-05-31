context("Testing spline building")

test_that("Death spline is correctly built", {
  
  params<-list("area_length_x"=1,    
               "cell_count_x"=100,  
               "periodic"=FALSE, 
               
               "b"=0.3,    
               "d"=0,    
               "dd"=0.01, 
               
               "seed"=1234,  
               "initial_population_x"=runif(100,min=0,max=1),
               
               "death_kernel_r"=1/200,
               "death_kernel_y"=dnorm((0:100)/20000,sd=0.001),
               
               "birth_kernel_r"=1/100,
               "birth_kernel_y"=dnorm((0:100)/10000,sd=0.002), 
               
               "spline_precision" = 1e-9  #without spline trimming
  )
  
  sim<-new(poisson_1d,params)
  
  death_kernel_x_test=(0:1000)/200000
  death_kernel_y_test=dnorm(death_kernel_x_test,sd=0.001)
  
  total_diff=0
  st=death_kernel_x_test[2]
  for(i in 1:1001)
  {
    total_diff=total_diff+abs(death_kernel_y_test[[i]]-sim$death_spline_at(death_kernel_x_test[[i]]))*st
  }
  expect_true(total_diff<1e-6) #0.0001% ish integral precision for uniform net 
  
  
})


test_that("Birth reverse cdf spline is correctly built", {
  require(Rcpp)
  require(devtools)
  params<-list("area_length_x"=1,    
               "cell_count_x"=100,  
               "periodic"=FALSE,
               
               "b"=0.3,    
               "d"=0,    
               "dd"=0.01, 
               
               "seed"=1234,  
               "initial_population_x"=runif(100,min=0,max=1),
               
               "death_kernel_r"=1/200,
               "death_kernel_y"=dnorm((0:100)/20000,sd=0.001),
               
               "birth_kernel_r"=1/100,
               "birth_kernel_y"=dnorm((0:100)/10000,sd=0.002), 
               
               "spline_precision" = 1e-9 #without spline trimming
  )
  
  sim<-new(poisson_1d,params)
  
  birth_kernel_x_quantile_test=(0:1000)/1000
  
  birth_kernel_y_quantile_test=qnorm(birth_kernel_x_quantile_test/2+0.5,sd=0.002)
  
  total_diff=0
  st=birth_kernel_x_quantile_test[2]
  
  for(i in 1:999)
  {
    total_diff=total_diff+abs(birth_kernel_y_quantile_test[[i]]-
                              sim$birth_reverse_cdf_spline_at(birth_kernel_x_quantile_test[[i]])
                              )*st
  }
  expect_true(total_diff<1e-4) #0.01% ish integral precision for uniform net 
  
  p.values<-sapply(1:100,function(i){
    random_numbers=runif(1000)
    simulated_results=sapply(random_numbers,function(x){sim$birth_reverse_cdf_spline_at(x)})*
      (rbinom(1000,1,0.5)*2-1)
    suppressWarnings(ks.test(simulated_results,"pnorm",0,0.002)$p.value)
  })

  expect_true(mean(p.values)>0.4) #High p-values all around
  
  expect_silent(for(i in 1:100000) sim$make_event())
})


test_that("Spline trimming works, try birth only",{

  params<-list("area_length_x"=1,    
               "cell_count_x"=100,  
               "periodic"=FALSE,  
               
               "b"=0.3,    
               "d"=0,    
               "dd"=0.01, 
               
               "seed"=1234,  
               "initial_population_x"=runif(100,min=0,max=1),
               
               "death_kernel_r"=1/20,
               "death_kernel_y"=dnorm((0:1000)/20000,sd=0.001),
               
               "birth_kernel_r"=1/10,
               "birth_kernel_y"=dnorm((0:1000)/10000,sd=0.002), 
               
               "spline_precision" = 1e-9  
  )
  
  sim<-new(poisson_1d,params)
  
  birth_kernel_x_quantile_test=(0:1000)/1000
  
  birth_kernel_y_quantile_test=qnorm(birth_kernel_x_quantile_test/2+0.5,sd=0.002)
  
  total_diff=0
  st=birth_kernel_x_quantile_test[2]
  
  for(i in 1:999)
  {
    total_diff=total_diff+abs(birth_kernel_y_quantile_test[[i]]-
                                sim$birth_reverse_cdf_spline_at(birth_kernel_x_quantile_test[[i]])
    )*st
  }
  expect_true(total_diff<1e-4) #0.01% ish integral precision for uniform net on reverse CDF
  
  
  p.values<-sapply(1:100,function(i){
    random_numbers=runif(1000)
    simulated_results=sapply(random_numbers,function(x){sim$birth_reverse_cdf_spline_at(x)})*
      (rbinom(1000,1,0.5)*2-1)
    suppressWarnings(ks.test(simulated_results,"pnorm",0,0.002)$p.value)
  })
  
  expect_true(mean(p.values)>0.1) #High p-values all around
  
  expect_silent(for(i in 1:100000) sim$make_event())
})

