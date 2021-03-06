require(MathBioSim)
library(foreach)
library(doParallel)



linespec = function(from, to, count) {
    return((0:(count)) / count * (to - from) + from)
}

file_name = "NEW_TEST"

area_length = 10
cell_count_x = 100
x_grid_points = 5000
epsilon = 1e-12

SPECIES_COUNT = 2;




b = c(0.4, 0.4)
d = c(0.2, 0.2)

dd = matrix(
    c(
        0.001, 0.001,
        0.001, 0.001
    ),
    nrow = SPECIES_COUNT,
    ncol = SPECIES_COUNT
)

sigma_m = c(0.04, 0.1)
sigma_w = matrix(
    c(
        0.04, 0.04,
        0.04, 0.04
    ),
    nrow = SPECIES_COUNT,
    ncol = SPECIES_COUNT
)

params<-list(
    "species_count"=SPECIES_COUNT,
    "area_length_x"=area_length,
    "cell_count_x"=cell_count_x,
    "seed"=42,
    "spline_precision"=1e-6,
    "periodic"=FALSE
)

popul_x = c()
popul_s = c()

for (i in 1:SPECIES_COUNT) {
    params[paste("b", i, sep = "_")] = b[i]
    params[paste("d", i, sep = "_")] = d[i]
    
    for (j in 1:(area_length * 100)) {
        popul_x = c(popul_x, c(runif(1, 0, area_length)))
        popul_s = c(popul_s, c(i - 1))
    }
    grid = linespec(0 + epsilon, 1 - epsilon, x_grid_points)
    sigma = sigma_m[i]
    
    params[[paste("birth_kernel", "y", i, sep = "_")]] = qnorm(grid, sd = sigma)
    
    for (j in 1:SPECIES_COUNT) {
        params[[paste("dd", i, j, sep = "_")]] = dd[i, j]
        
        sigma = sigma_w[i, j]
        length = sigma * 5
        grid = linespec(0, length, x_grid_points)
        
        params[[paste("death_kernel", "r", i, j, sep = "_")]] = length
        params[[paste("death_kernel", "y", i, j, sep = "_")]] = dnorm(grid, sd = sigma)
    }
}

params[[paste("initial_population", "x", sep = "_")]] = popul_x
params[["initial_population_species"]] = popul_s
sim = new(poisson_1d_n_species, params);

iter_count = 1000
events_coef = 1000
res = matrix(
    ncol = SPECIES_COUNT + 1,
    nrow = iter_count + 1
)

for (i in 1:iter_count) {
    res[i, 1] = iter_count * (i - 1)
    res[i, 2:(SPECIES_COUNT + 1)] = sim$total_population() / area_length
    sim$run_events(events_coef)
}
  
res[(iter_count+1), 1] = iter_count * (iter_count)
res[(iter_count+1), 2:(SPECIES_COUNT + 1)] = sim$total_population() / area_length 
