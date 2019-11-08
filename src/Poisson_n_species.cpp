// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>

#include "grid.h"

#ifndef POISSON_1D_N_SPECIES
#define POISSON_1D_N_SPECIES

RCPP_MODULE(poisson_1d_n_species_module)
{
  using namespace Rcpp;
  
  class_<class Grid>("poisson_1d_n_species")
    .constructor<List>()
    .field_readonly("area_length_x", &Grid::area_length_x)
    .field_readonly("cell_count_x", &Grid::cell_count_x)
  
  .field_readonly("b", &Grid::b)
  .field_readonly("d", &Grid::d)
  .field_readonly("dd", &Grid::dd)
  
  .field_readonly("seed", &Grid::seed)
  .field_readonly("initial_density", &Grid::initial_density)

  .field_readonly("death_cutoff_r", &Grid::death_cutoff_r)

  .field_readonly("cell_death_rates", &Grid::cell_death_rates)
  .field_readonly("cell_population", &Grid::cell_population)
  
  .method("get_all_coordinates", &Grid::get_all_coords)
  .method("GetAllCoordsForSpecies", &Grid::GetAllCoordsForSpecies)
  .method("get_all_death_rates", &Grid::get_all_death_rates)
  
  .method("get_x_coordinates_in_cell", &Grid::get_coords_at_cell)
  
  .method("make_event", &Grid::make_event)
  .method("run_events", &Grid::run_events)
  .method("run_for", &Grid::run_for)
  
  .field_readonly("total_population", &Grid::total_population)
  .field_readonly("total_death_rate", &Grid::total_death_rate)
  .field_readonly("events", &Grid::event_count)
  .field_readonly("time", &Grid::time);
}


#endif
