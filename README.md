# Poisson point process simulation of spatial population dynamics

This repository contains source code (R/Rcpp) required to run spatial population dynamics simualation.

# Installation
### Windows
  - Install R ([R Installer](https://cran.r-project.org/bin/windows/base/)) - latest version used is 3.6.2
  - Install Rtools ([Rtools Installer](https://cran.r-project.org/bin/windows/Rtools/))
  - Install required R packages from CRAN (devtools, Rcpp, BH)
  - Install this source package using devtools
```R
devtools::install_github("YegorGalkin/RcppSim")
```
### Linux \ Ubuntu
TODO: like windows but without Rtools
### MacOS
TODO: Like Linux

# Simulator usage (1d)
1) Prepare named list of simulation options. Parameters:
    - area_length_x - area size
    - cell_count_x - used to split simulation interval into subintervals for filtering negliable interactions.
    - b - poisson birth rate
    - d - poisson death rate
    - dd - competitive death rate
    - seed - RNG seed
    - init_density - per unit starting density of speciments
    - death_kernel_r - max death kernel interaction radius
    - death_kernel_y - death kernel values on uniform grid on [0, death_kernel_r] interval
    - birth_kernel_r - max birth kernel interaction radius
    - birth_kernel_y - birth kernel values on uniform grid on [0, birth_kernel_r] interval
    - spline_precision - tries to keep spline approximation precise up to this value (for numerical integration and cdf inversal)
2) Run simulation for specified time or number of events
    ```R
    sim$run_events(1e9)
    sim$run_time(100)
    ```
3) Extract and visualize results (using packages spatstat and tidyverse)
    ```R
    sim$time
    sim$total_population
    sim$get_all_coordinates()
    
    points <- ppp(sim$get_all_coordinates(),
                  rep(0,sim$total_population),
                  c(0,sim$area_length_x),
                  c(-sim$area_length_x/2,sim$area_length_x/2)
                  )
                  
    pcf_grid = seq(0,1,length.out = 1001)    
    
    K_estimate <- Kest(points,r=pcf_grid,correction="Ripley")
    
    pcf_estimate <- data.frame(Kest=K_estimate$iso/2,x=pcf_grid)%>%
                    mutate(pfc=(Kest-lag(Kest))/(pcf_grid-lag(pcf_grid))/sim$area_length_x)%>%
                    pull(pfc)
                    
    ggplot(data=data.frame(x=pcf_grid,y=pcf_estimate),aes(x=x,y=y))+
    geom_point()
    ```
4) Check examples in this repository (tests, examples)
