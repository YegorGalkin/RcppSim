# require(MathBioSim)
# library(MathBioSim)
# 
# get_seed <- function() {
#   return(floor(runif(1, min=0, max=2^32-1)))
# }
# 
# run_sim <- function(x, params, area_length) {
#   require(MathBioSim)
#   params["seed"] = seed = floor(runif(1, min=0, max=2^32-1))
#   cat("seed: ", seed, "\t");
#   sim = new(poisson_2spicies_1d, params)
#   COUNT = 400000
#   pop = matrix(nrow = 2, ncol = COUNT)
#   for (k in 1:COUNT){
#     sim$make_event()
#     pop[,k] = sim$total_population / area_length
#   }
#   return(
#     c(
#       mean(pop[1,COUNT - 10:COUNT]),
#       mean(pop[2,COUNT - 10:COUNT])
#     )
#   )
# }
# 
# 
# area_length = 90
# cell_count_x = 50
# x_grid_points = 1000
# 
# b = c(0.4, 0.4)
# d = c(0.2, 0.2)
# dd = matrix(
#   c(0.001, 0.001, 0.001, 0.001),
#   nrow = 2,
#   ncol = 2
# )
# 
# time_start<-Sys.time()
# 
# SM = c(0.04, 0.04)
# SW = matrix(
#   c(0.04, 0.04, 0.04, 0.04),
#   ncol = 2,
#   nrow = 2
# )
# x_length = max(SW) * 5
# x_grid = (0:x_grid_points) / x_grid_points * x_length
# 
# m1=dnorm(x_grid, sd = SM[1])
# m2=dnorm(x_grid, sd = SM[2])
# 
# w11=dnorm(x_grid, sd = SW[1, 1])
# w12=dnorm(x_grid, sd = SW[1, 2])
# w21=dnorm(x_grid, sd = SW[2, 1])
# w22=dnorm(x_grid, sd = SW[2, 2])
# 
# linespec = function(from, to, count) {
#   return((0:(count - 1)) / count * (to - from) + from)
# }
# 
# COUNTS = 1000;
# 
# counts = c();
# N1 = c();
# N2 = c();
# Y1 = c();
# Y2 = c();
# 
# for (i in 1:COUNTS) {
#   cat("\r ", i, "\t");
#   
#   params<-list(
#     "area_length_x"=area_length,
#     "cell_count_x"=cell_count_x,
#     "seed"=42,
#     
#     "b1"=b[1],
#     "b2"=b[2],
#     "d1"=d[1],
#     "d2"=d[2],
#     "dd11"=dd[1, 1],
#     "dd12"=dd[1, 2],
#     "dd21"=dd[2, 1],
#     "dd22"=dd[2, 2],
#     
#     "init_density1"=20,
#     "init_density2"=20,
#     
#     "death_kernel_y11"=w11,
#     "death_kernel_y12"=w12,
#     "death_kernel_y21"=w21,
#     "death_kernel_y22"=w22,
#     
#     "birth_kernel_y1"=m1,
#     "birth_kernel_y2"=m2,
#     
#     "death_kernel_r"=x_length,
#     "birth_kernel_r"=x_length,
#     
#     "spline_precision"=1e-6
#   )
#   a = run_sim(0, params = params, area_length = area_length);
#   N1[i] = a[1];
#   N2[i] = a[2];
#   counts[i] = i;
#   Y1[i] = mean(N1[1:i]);
#   Y2[i] = mean(N2[1:i]);
#   plot(counts[1:i], Y1[1:i], type = 'l', col = 'red', ylim = c(0, 200));
#   lines(counts[1:i], Y2[1:i], type = 'l', col = 'blue');
#   abline(h=198.158, col = 'red')
#   abline(h=0.000197447, col = 'blue')
# }
