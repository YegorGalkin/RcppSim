#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <boost/random.hpp>

#include <Rcpp.h>

#ifndef GRID
#define GRID

#include "defines.h"
#include "cell.h"
#include "iterators.h"
#include "unit.h"

using boost::math::cubic_b_spline;

struct Grid {
  VEC<Cell> cells;
  MAT<double> cell_death_rates;
  MAT<int> cell_population;
  
  double area_length_x;
  
  int cell_count_x;
  int species_count;
  
  VEC<double> b, d;
  
  MAT<double> dd;
  int seed;
  boost::random::lagged_fibonacci2281 rng;
  
  VEC<double> initial_density;
  
  VEC<int> total_population;
  VEC<double> total_death_rate;
  
  double time;
  size_t event_count;
  
  MAT<double> death_cutoff_r;
  
  MAT<cubic_b_spline<double>> death_kernel_spline;
  
  VEC<cubic_b_spline<double>> birth_reverse_cdf_spline;
  
  int cull_x;
  
  std::tuple<double, int> last_event;
  
public:
  Cell &cell_at(int i);
  
  double &cell_death_rate_at(int s, int i);
  
  void chek();
  
  int get_all_population();
  
  VEC<DCoord> get_coords_at_cell(int i);
  
  Range RangeSpecies();
  
  UnitIterating RangeAllUnits();
  
  UnitIterating RangeLocalUnits(int i);
  
  double CalcInteraction(const Unit& cell1, const Unit& cell2);
  
  VEC<DCoord> get_all_coords();
  
  VEC<DCoord> GetAllCoordsForSpecies(int species);
  
  VEC<double> get_all_death_rates();
  
  void AddDeathRate(Unit& unit);
  void SubDeathRate(Unit& unit);
  
  void Initialize_death_rates();
  
  void AddInteraction(Unit& cell, double interaction);
  
  Unit GetUnit(int i, int j);
  
  double GetRandomDeathIndex(int s);
  
  void kill_random(int s);
  
  int GetRandomSpawnCell(int species);
  
  Unit GetRandomSpawnUnit(int cellIndex, int species);
  
  double GetNewCoord(const Unit& unit);
  
  double IsInArea(DCoord coord);
  
  void SetLastEvent(double d, int i);
  
  int GetNewCellIndex(DCoord x);
  
  Unit GetLastUnit(int cell);
  
  void IncrementPopulation(Unit& unit);
  void DecrementPopulation(Unit& unit);
  
  void spawn_random(int s);
  
  double GetRandomTime();
  
  void make_event();
  void run_events(int events);
  
  void run_for(double time);
  
  Grid(Rcpp::List params);
};

#endif
