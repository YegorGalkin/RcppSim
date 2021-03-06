// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include "grid.h"

using namespace Rcpp;

RCPP_MODULE(poisson_n_species_module) {
    using namespace Rcpp;
    
    class_<Grid<1>>("poisson_1d_n_species")
        .constructor<List>()
        
        .method("make_event", &Grid<1>::MakeEvent)
        .method("run_events", &Grid<1>::RunEvents)
        .method("total_population", &Grid<1>::GetTotalPopulation);
        
    class_<Grid<2>>("poisson_2d_n_species")
        .constructor<List>()
        
        .method("make_event", &Grid<2>::MakeEvent)
        .method("run_events", &Grid<2>::RunEvents)
        .method("total_population", &Grid<2>::GetTotalPopulation);
        
    class_<Grid<3>>("poisson_3d_n_species")
        .constructor<List>()
        
        .method("make_event", &Grid<3>::MakeEvent)
        .method("run_events", &Grid<3>::RunEvents)
        .method("total_population", &Grid<3>::GetTotalPopulation);
}
