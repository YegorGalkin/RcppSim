#include <boost/math/tools/roots.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>

#include "grid.h"

Cell &Grid::cell_at(int i)
{
  return cells[i];
}

double &Grid::cell_death_rate_at(int s, int i)
{
  return cell_death_rates[s][i];
}

void Grid::chek() {
  for (auto x : total_death_rate) {
    if (isnan(x)) {
      abort();
    }
  }
}

int Grid::get_all_population() {
  int res = 0;
  for (auto i : RangeSpecies()) {
    res += total_population[i];
  }
  return res;
}

VEC<DCoord> Grid::get_coords_at_cell(int i)
{
  return cells[i].coords_x;
}

Range Grid::RangeSpecies() {
  return Range(species_count);
}

UnitIterating Grid::RangeAllUnits() {
  return UnitIterating(
    cells,
    cell_death_rates,
    cell_population,
    0,
    cell_count_x
  );
}

UnitIterating Grid::RangeLocalUnits(int i) {
  return UnitIterating(
    cells,
    cell_death_rates,
    cell_population,
    std::max(0, i - cull_x),
    std::min(cell_count_x, i + cull_x + 1)
  );
}

double Grid::CalcInteraction(const Unit& cell1, const Unit& cell2) {
  double distance = Cell::Ro(cell1, cell2);
  
  if (distance > death_cutoff_r[cell1.Species()][cell2.Species()])
    return -1; //Too far to interact
  
  return dd[cell1.Species()][cell2.Species()] *
    death_kernel_spline[cell1.Species()][cell2.Species()](distance);
}

VEC<DCoord> Grid::get_all_coords()
{
  VEC<DCoord> result;
  for (auto cell : cells)
  {
    if (cell.coords_x.size() != 0)
      result.insert(result.end(), cell.coords_x.begin(), cell.coords_x.end());
  }
  return result;
}

VEC<DCoord> Grid::GetAllCoordsForSpecies(int species) {
  VEC<DCoord> res;
  res.reserve(total_population[species]);
  for (auto unit : RangeAllUnits()) {
    if (unit.Species() == species) {
      res.emplace_back(unit.Coord());
    }
  }
  return res;
}

VEC<double> Grid::get_all_death_rates()
{
  VEC<double> result;
  for (auto cell : cells)
  {
    if (cell.coords_x.size() != 0)
      result.insert(result.end(), cell.death_rates.begin(), cell.death_rates.end());
  }
  return result;
}

void Grid::AddDeathRate(Unit& unit) {
  unit.CellDeathRate() += d[unit.Species()];
  total_death_rate[unit.Species()] += d[unit.Species()];
}

void Grid::SubDeathRate(Unit& unit) {
  unit.CellDeathRate() -= d[unit.Species()];
  total_death_rate[unit.Species()] -= d[unit.Species()];
  if (unit.CellDeathRate() < 1e-10) {
    unit.CellDeathRate() = 0;
  }
}

void Grid::Initialize_death_rates() {
  for (auto &x : total_death_rate) {
    x = 0;
  }
  cell_death_rates = VEC<VEC<double>>(species_count, VEC<double>(cells.size(), 0));
  cell_population = VEC<VEC<int>>(species_count, VEC<int>(cell_count_x, 0));
  
  for (auto s : RangeSpecies()) {
    double x_coord;
    for (auto population : Range(ceil(area_length_x * initial_density[s]))) {
      x_coord = boost::uniform_real<>(0, area_length_x)(rng);
      
      auto i = GetNewCellIndex(x_coord);
      cell_at(i).Add(d[s], x_coord, s);
      auto unit = GetLastUnit(i);
      
      AddDeathRate(unit);
      IncrementPopulation(unit);
      chek();
    }
  }
  
  for (auto cell1 : RangeAllUnits()) {
    for (auto cell2 : RangeLocalUnits(cell1.CellNum())) {
      if (cell1 == cell2)
        continue; // same speciment
      
      double interaction = CalcInteraction(cell1, cell2);
      if (interaction < 0) {
        continue;
      }
      if (isnan(interaction)) {
        abort();
      }
      
      AddInteraction(cell1, interaction);
    }
  }
  chek();
}

void Grid::AddInteraction(Unit& cell, double interaction) {
  cell.DeathRate() += interaction;
  cell.CellDeathRate() += interaction;
  total_death_rate[cell.Species()] += interaction;
}

Unit Grid::GetUnit(int i, int j) {
  return Unit(*this, i, j);
}

double Grid::GetRandomDeathIndex(int s) {
  return boost::random::discrete_distribution<>(cell_death_rates[s])(rng);
}

void Grid::kill_random(int s)
{
  if (total_population[s] == 0) {
    return;
  }
  int cellDeathIndex = GetRandomDeathIndex(s);
  if (cell_population[s][cellDeathIndex] == 0) {
    return;
  }
  
  auto &death_cell = cells[cellDeathIndex];
  
  VEC<double> tmp_vec;
  VEC<size_t> translate;
  tmp_vec.reserve(death_cell.death_rates.size());
  translate.reserve(death_cell.death_rates.size());
  for (size_t i = 0; i < death_cell.death_rates.size(); ++i) {
    if (death_cell.species[i] == s) {
      tmp_vec.push_back(death_cell.death_rates[i]);
      translate.push_back(i);
    }
  }
  if (tmp_vec.size() == 0) {
    abort();
  }
  int in_cellDeathIndex = boost::random::discrete_distribution<>(tmp_vec)(rng);
  in_cellDeathIndex = translate[in_cellDeathIndex];
  if(death_cell.species[in_cellDeathIndex] != s) {
    std::cout << "Not that species" << std::endl;
    abort();
  }
  
  SetLastEvent(death_cell.coords_x[in_cellDeathIndex], -1);
  
  auto cellKilled = GetUnit(cellDeathIndex, in_cellDeathIndex);
  
  for (auto cell : RangeLocalUnits(cellDeathIndex)) {
    if (cell == cellKilled)
      continue;
    
    double interaction = CalcInteraction(cellKilled, cell);
    if (interaction < 0) {
      continue;
    }
    
    AddInteraction(cell, -interaction);
    AddInteraction(cellKilled, -interaction);
  }
  //remove dead speciment
  SubDeathRate(cellKilled);
  DecrementPopulation(cellKilled);
  death_cell.SwapWithLast(in_cellDeathIndex);
  death_cell.Pop();
}

int Grid::GetRandomSpawnCell(int species) {
  return boost::random::discrete_distribution<>(cell_population[species])(rng);
}

Unit Grid::GetRandomSpawnUnit(int cellIndex, int species) {
  int eventIndex = boost::random::uniform_smallint<>(
    0,
    cell_population[species][cellIndex] - 1
  )(rng);
  ++eventIndex;
  for (int i = 0; ; ++i) {
    if (cell_at(cellIndex).species[i] == species) {
      --eventIndex;
    }
    if (eventIndex == 0) {
      eventIndex = i;
      break;
    }
  }
  
  return GetUnit(cellIndex, eventIndex);
}

double Grid::GetNewCoord(const Unit& unit) {
  return unit.Coord() +
    birth_reverse_cdf_spline[unit.Species()](
        boost::random::uniform_01<>()(rng)
    ) * (
        boost::random::bernoulli_distribution<>(0.5)(rng) * 2 - 1
    );
}

double Grid::IsInArea(double coord) {
  return coord >= 0 && coord <= area_length_x;
}

void Grid::SetLastEvent(double d, int i) {
  last_event = std::make_tuple(d, i);
}

int Grid::GetNewCellIndex(DCoord x) {
  int new_i = floor(x * cell_count_x / area_length_x);
  if (new_i >= cell_count_x)
    new_i = cell_count_x - 1;
  if (new_i < 0) {
    new_i = 0;
  }
  return new_i;
}

Unit Grid::GetLastUnit(int cell) {
  return GetUnit(cell, cell_at(cell).coords_x.size() - 1);
}

void Grid::IncrementPopulation(Unit& unit) {
  ++unit.CellPopulation();
  ++total_population[unit.Species()];
}

void Grid::DecrementPopulation(Unit& unit) {
  --unit.CellPopulation();
  assert(unit.CellPopulation() >= 0);
  --total_population[unit.Species()];
}

void Grid::spawn_random(int s) {
  int cellIndex = GetRandomSpawnCell(s);
  auto parentCell = GetRandomSpawnUnit(cellIndex, s);
  double coordNew = GetNewCoord(parentCell);
  
  if (!IsInArea(coordNew)) {
    SetLastEvent(coordNew, 0);
  } else {
    SetLastEvent(coordNew, 1);
    
    auto newCellIndex = GetNewCellIndex(coordNew);
    cell_at(newCellIndex).Add(d[s], coordNew, s);
    auto newCell = GetLastUnit(newCellIndex);
    
    AddDeathRate(newCell);
    IncrementPopulation(newCell);
    
    for (auto cell : RangeLocalUnits(newCellIndex)) {
      if (cell == newCell)
        continue;
      
      double interaction = CalcInteraction(newCell, cell);
      if (interaction < 0) {
        continue;
      }
      
      AddInteraction(cell, interaction);
      AddInteraction(newCell, interaction);
    }
  }
}

double Grid::GetRandomTime() {
  double rates = 0;
  for (auto s : RangeSpecies()) {
    rates += total_population[s] * b[s] + total_death_rate[s];
  }
  return boost::random::exponential_distribution<>(rates)(rng);
}

void Grid::make_event() {
  if (get_all_population() == 0)
    return;
  
  ++event_count;
  time += GetRandomTime();
  //Rolling event according to global birth \ death rate
  std::vector<double> dis(species_count * 2, 0);
  for (auto s : RangeSpecies()) {
    if (total_population[s] > 0) {
      dis[2 * s + 0] = total_death_rate[s];
      dis[2 * s + 1] = total_population[s] * b[s];
    }
  }
  
  auto t = boost::random::discrete_distribution<>(dis)(rng);
  int event = t % 2;
  int species = t / 2;
  if (event == 0) {
    kill_random(species);
  } else {
    spawn_random(species);
  }
  chek();
}
void Grid::run_events(int events) {
  if (events > 0) {
    for (int i = 0; i < events; i++) {
      make_event();
    }
  }
}

void Grid::run_for(double time)
{
  if (time > 0.0)
  {
    double time0 = this->time;
    while (this->time < time0 + time)
    {
      make_event();
      if (get_all_population() == 0)
        return;
    }
  }
}

Grid::Grid(Rcpp::List params) {
  using std::vector;
  using std::to_string;
  //Parse parameters
  
  Rcpp::Environment base("package:base");
  
  // Make function callable from C++
  Rcpp::Function print = base["print"];
  
  area_length_x = Rcpp::as<double>(params["area_length_x"]);
  cell_count_x = Rcpp::as<int>(params["cell_count_x"]);
  species_count = Rcpp::as<int>(params["species_count"]);
  seed = Rcpp::as<int>(params["seed"]);
  rng = boost::random::lagged_fibonacci2281(uint32_t(seed));
  
  cull_x = 3;
  event_count = 0;
  
  cells = VEC<Cell>(cell_count_x);
  
  b = VEC<double>(species_count);
  d = VEC<double>(species_count);
  initial_density = VEC<double>(species_count);
  dd = MakeMat<double>(species_count);
  death_cutoff_r = MakeMat<double>(species_count);
  death_kernel_spline = MakeMat<cubic_b_spline<double>>(species_count);
  birth_reverse_cdf_spline = VEC<cubic_b_spline<double>>(species_count);
  
  total_death_rate = VEC<double>(species_count);
  total_population = VEC<int>(species_count);
  for (auto i : RangeSpecies()) {
    auto n1 = to_string(i + 1);
    b[i] = Rcpp::as<double>(params["b_" + n1]);
    d[i] = Rcpp::as<double>(params["d_" + n1]);
    initial_density[i] = Rcpp::as<double>(params["init_density_" + n1]);
    auto birth_kernel_y = Rcpp::as<vector<double>>(params["birth_kernel_y_" + n1]);
    auto birth_cutoff_r = Rcpp::as<double>(params["birth_kernel_r_" + n1]);
    
    for (auto j : RangeSpecies()) {
      auto n2 = n1 + "_" + to_string(j + 1);
      dd[i][j] = Rcpp::as<double>(params["dd_" + n2]);
      auto death_kernel_y = Rcpp::as<vector<double>>(params["death_kernel_y_" + n2]);
      death_cutoff_r[i][j] = Rcpp::as<double>(params["death_kernel_r_" + n2]);
      
      int death_spline_nodes = death_kernel_y.size();
      double death_step = death_cutoff_r[i][j] / (death_spline_nodes - 1);
      death_kernel_spline[i][j] = cubic_b_spline<double>(
        death_kernel_y.begin(),
        death_kernel_y.end(),
        0,
        death_step,
        0,
        0
      );
      cull_x = std::max<int>(cull_x, ceil(death_cutoff_r[i][j] / (area_length_x / cell_count_x)));
    }
    
    auto birth_spline_nodes = birth_kernel_y.size();
    auto birth_step = birth_cutoff_r / (birth_spline_nodes - 1);
    
    //Build birth spline
    auto birth_kernel_spline = cubic_b_spline<double>(
      birth_kernel_y.begin(),
      birth_kernel_y.end(),
      0,
      birth_step,
      0,
      0
    );
    
    vector<double> x_quantile_1d_array(birth_spline_nodes);
    vector<double> y_quantile_1d_array(birth_spline_nodes);
    
    using boost::math::quadrature::trapezoidal;
    double approx_const = trapezoidal(
      [&](double y) {return birth_kernel_spline(y);},
      0.0,
      birth_cutoff_r
    );
    
    using boost::math::tools::newton_raphson_iterate;
    
    for (int j = 0; j < birth_spline_nodes; j++) {
      x_quantile_1d_array[j] = (double)j / (birth_spline_nodes - 1);
      y_quantile_1d_array[j] =
        newton_raphson_iterate(
          [&](double y) {
            return std::make_tuple(
              trapezoidal(
                [&](double z) {return birth_kernel_spline(z);},
                0.0,
                y
              ) / approx_const - x_quantile_1d_array[j],
                                                    birth_kernel_spline(y) / approx_const
            );
          },
          1e-10,
          0.0,
          birth_cutoff_r,
          std::numeric_limits<double>::digits
        );
    }
    int k = 0;
    while (y_quantile_1d_array[k] < 1e-300)
    {
      k++;
    }
    
    vector<double> y_quantile_1d_array_temp(birth_spline_nodes - k);
    
    for (auto j : Range(birth_spline_nodes - k)) {
      y_quantile_1d_array_temp[j] = y_quantile_1d_array[k + j];
    }
    
    auto birth_reverse_cdf_step = 1.0/(y_quantile_1d_array_temp.size()-1);
    // Extrapolate last quantile element
    double right_derivative = (y_quantile_1d_array_temp[y_quantile_1d_array_temp.size()-2] - y_quantile_1d_array_temp[y_quantile_1d_array_temp.size()-3])/birth_reverse_cdf_step;
    
    y_quantile_1d_array_temp[y_quantile_1d_array_temp.size() - 1] =
      y_quantile_1d_array_temp[y_quantile_1d_array_temp.size() - 2] +
      right_derivative * birth_reverse_cdf_step;
    
    //Ensure correct derivative at 0, equal to 1/2*birth_kernel(0)
    birth_reverse_cdf_spline[i] = cubic_b_spline<double>(
      y_quantile_1d_array_temp.begin(),
      y_quantile_1d_array_temp.end(),
      0,
      birth_reverse_cdf_step,
                            0.5 / birth_kernel_spline(0),
                            2 * right_derivative);
  }
  
  Initialize_death_rates();
}
