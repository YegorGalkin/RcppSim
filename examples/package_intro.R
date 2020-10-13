# Introduction to MathBioSim package
# Run code line by line to reproduce graphs
# 

# This installs library for ease of package management 
# p_load loads existing package or downloads it from CRAN and then loads it.
if (!require(pacman)){
  install.packages('pacman')
}
library(pacman)

# Takes lots of time on first run. 
# tidyverse is framework with many usefull packages (plots, pipes, functional tools)
# spatstat handles spatial moment calculation
# devtools allows installing MathBioSim from github
p_load(tidyverse,spatstat,devtools)


# Reqires Rtools installed on Windows
# See https://cran.r-project.org/bin/windows/Rtools/
if (!require(MathBioSim)){
  devtools::install_github('https://github.com/YegorGalkin/RcppSim/')
}
library(MathBioSim)

# Lets start with 2d simulation (pretty plots)
# First of all we need to initialise list of parameters that are passed to C++ code

sim_params <-
  list("area_length_x"=15, # The exact order of parameters is not important
       "area_length_y"=15, # First we specify area size, 15x15 here
       
       "cell_count_x"=100, # This parameter tells simulator how finely to split the area 
       "cell_count_y"=100, # for neighbour search speed up. 100 is arbitrary default value.
       
       "periodic"=TRUE,    # If false, new individuals will die if spawned outside of area,
                           # If true, they will appear on the other side
       
       "b"=1.5,            # birth rate for single individual
       "d"=0.5,            # death rate for single individual
       "dd"=0.01,          # competitive death rate 
       
       "seed"=1235,        # random number generator seed
       
       # We can specify exact location of every individual at start
       # We can also pass random distributions from R with runif function
       "initial_population_x"=c(0.05), 
       "initial_population_y"=c(0.05), 
       
       "population_limit" = 100000L, # This is a hard cap for population - simulation stops once reached
       
       # Here we specify death and birth kernels
       # They are evaluated on a grid with 1001 points on [0,5] interval
       # They are standart normal distributions - 0 mean and 1 standart deviation
       # dnorm calculates density for given grid with given normal distribution parameters
       
       "death_kernel_r"=5, 
       "death_kernel_y"=dnorm(x_grid_death <- seq(0,5,length.out = 1001), sd = 1),
       
       "birth_kernel_r"=5,
       "birth_kernel_y"=dnorm(x_grid_birth <- seq(0,5,length.out = 1001), sd = 1), 
       
       # Technical parameter, specifies internal simulation precision for kernels
       "spline_precision" = 1e-9
  )

# This creates new simulation object with binding to C++ code
sim<-new(poisson_2d,sim_params)
# Type this line by line to see what's inside simulator
print(sim$get_all_x_coordinates())
print(sim$get_all_y_coordinates())
print(sim$get_all_death_rates())
# Since we only initialized simulation, there is just one individual at 0.05,0.05 with 0 death rate
# Lets try making a single simulation step
sim$make_event()
sim$total_population
# The population changed to 2, since death rate is zero 
# We can now see how coordinates and death rates updated
print(sim$get_all_x_coordinates())
print(sim$get_all_y_coordinates())
print(sim$get_all_death_rates())
# You should get the following:
# X coordinates: 0.0500 14.4512 - new individual passed the area bound and turned up on the other side
# Y coordinates: 0.0500000 0.7166756
# Death rates: 0.5026702 0.5026702 - death rates are the same since there are only 2 individuals
# As we can see, extra death rate from competition is only ~0.0028, about 0.5% from intrinsic death rate
# We are far away from equilibrium

# Now we can do something more complicated, lets try making a movie of simulation steps
# Loading packages for video rendering and plots
p_load(gganimate,av,glue)

# Lets try 100 frames
n_times = 100

# Preparing empty lists for results
population_list_x=list()
population_list_y=list()

# Run simulator 100 times for 0.33 time units, save population coordinates to lists
# Takes a few minutes
for(j in 1:n_times){
  sim$run_for(0.33)
  population_list_x[[j]]=sim$get_all_x_coordinates()
  population_list_y[[j]]=sim$get_all_y_coordinates()
  print(glue('{j} frame done'))
}

# We have population lists indexed by integers corresponding to frames
# print(population_list_x)
# To get usefull table we first enframe it - we get a table with a frame number column and 
# another column with list of coordinates.
# Then we unnest it - the list column is converted to a simple numeric column
# The frame number is kept for every point
# Then we merge it with similar table for all y coordinates
# In the end we have table of a type frame - x_coord - y_coord for every individual
# print(data)
data <- enframe(population_list_x)%>%
  unnest(value)%>%
  rename(x_coordinate=value)%>%
  bind_cols(enframe(population_list_y)%>%
              unnest(value)%>%
              select(value))%>%
  rename(y_coordinate=value,time_point=name)

# This uses ggplot and gganimate packages to create a video out of several plots
# Each plot is a scatterplot - geom_point - and frame is used as a transition point
animation <- data%>%
  ggplot(aes(x=x_coordinate,y=y_coordinate)) +
  geom_point(size=0.5,alpha=0.5)+
  labs(x='X axis',y='Y axis')+
  transition_time(time_point)

# Then we render it and save to file
# Saves to working directory of R, change population_periodic_2d.mp4 to where you want to save it
animate(animation,renderer=av_renderer("population_periodic_2d.mp4"),fps=5)

# Now we can explore population growth and spatial moments in greater detail
# First lets draw a plot of population growth for our frames

data%>%
  group_by(time_point)%>%
  summarise(population=n())%>%
  ggplot(aes(x=time_point,y=population))+
  geom_line()+
  labs(x='Time (in frames)',y='Total population')

# As we can see, the equilibrium happens somewhere after 35-50 frames
# Lets check last frame distribution just in case
data%>%
  filter(time_point==100)%>%
  ggplot(aes(x=x_coordinate,y=y_coordinate)) +
  geom_point(size=0.5,alpha=0.5)+
  labs(x='X axis',y='Y axis')
# Seems ok, hard to say whether we have equilibrium or not from a single population plot

# Next we can calculate pair correlation function for the last frame and plot it as well
# First we initialise spatstat point process object with our simulation parameters
# We need to pass x coordinates, y coordinates, x area limits and y area limits
points<-ppp(data%>%filter(time_point==100)%>%pull(x_coordinate),
            data%>%filter(time_point==100)%>%pull(y_coordinate),
            c(0,15),
            c(0,15))

# We pass the points object, the grid on which we want to calculate the pcf and the PCF method
# In this case we estimate K function using border correction and take derivative of it
# See ?Kest

# Estimating derivatives of samples function, especially divided by a radius is a pain
# See ?pcf.fasp
K_estimate <- Kest(points,correction = 'Ripley')
PCF_estimate <- pcf(K_estimate, spar=1, method="b")

# Lets plot the PCF and see if there is anything interesting
# For completely spatially random point processes it should be equal to 1 everywhere

ggplot(data=data.frame(x=PCF_estimate$r,y=PCF_estimate$pcf),
       aes(x=x,y=y))+
  geom_point()

# Looks like we have an anomaly in zero - lets take a closer look

ggplot(data=data.frame(x=PCF_estimate$r,y=PCF_estimate$pcf),
       aes(x=x,y=y))+
  geom_point()+
  xlim(0,0.25)
# Similar results, with anomaly only at zero, are generally unreliable - since estimation
# at 0 of pair correlation function in 2d case is unstable
# We could try averaging over several frames to check if results persist.

pcfs <- data%>%
  filter(time_point>=40)%>%
  group_by(time_point)%>%
  summarise(points = list(ppp(x_coordinate,y_coordinate,c(0,15),c(0,15))))%>%
  mutate(K_estimate = map(points,~Kest(.x,r=seq(0,3,length.out = 501),correction = 'Ripley')))%>%
  mutate(PCF_estimate = map(K_estimate,~pcf(.x, spar=1, method="b")))

# We now have 60 PCF estimates
mean_pcfs <- pcfs %>% 
  mutate(PCF_estimate = map(PCF_estimate,as.data.frame))%>%
  pull(PCF_estimate)%>%
  bind_rows()%>%
  group_by(r)%>%
  summarise(mean_pcf=mean(pcf))

ggplot(data=mean_pcfs,
       aes(x=r,y=mean_pcf))+
  geom_point()
# Probably artefact of PCF calculation
# Anyway, it seems like in 2d case distribution is very similar to spatially random
# Lets try running full spatial randomness test

quadrat.test(points)
# P-value is 0.8471 - cant exclude CSP
dclf.test(points)
# P-value is 0.62 - cant exclude CSP
mad.test(points)
# P-value is 0.38 - cant exclude CSP
# Looks like its completely spatially random! It is unexpected.
# See https://github.com/YegorGalkin/matbio/blob/master/PCFs.ipynb 
# for proper PFC results (1d case)
