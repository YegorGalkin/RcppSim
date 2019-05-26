// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>

#include <array>
#include <stdio.h>
#include <functional>
#include <list>
#include <chrono>
#include <ctime>
#include <vector>
#include <string>
#include <fstream>
#include <math.h>
#include <tuple>
#include <iomanip>

#include <boost/random.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/math/tools/roots.hpp>
#include <boost/fusion/sequence.hpp>
#include <boost/fusion/include/sequence.hpp>
#include <boost/fusion/adapted/std_tuple.hpp>
#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>


using namespace std;

#ifndef POISSON_2species_1D_H
#define POISSON_2species_1D_H

constexpr int SPICIES_COUNT = 2;

using boost::math::cubic_b_spline;

struct Cell_multi_1d
{
  
  vector<double> coords_x;
  vector<double> death_rates;
  vector<int> spicies;
  
  Cell_multi_1d() {}
};

template <class T>
using VEC = std::array<T, SPICIES_COUNT>;
//using VEC = std::vector<T>;

template <class T>
using MAT = std::array<std::array<T, SPICIES_COUNT>, SPICIES_COUNT>;
//using MAT = std::vector<std::vector<T>>;

struct Grid_multy_1d
{
  std::vector<Cell_multi_1d> cells;
  VEC<std::vector<double>> cell_death_rates;
  VEC<std::vector<int>> cell_population;
  
  double area_length_x;
  
  int cell_count_x;
  
  VEC<double> b, d;
  
  MAT<double> dd;
  int seed;
  boost::random::lagged_fibonacci2281 rng;
  
  VEC<double> initial_density;
  
  VEC<int> total_population;
  VEC<double> total_death_rate;
  
  double time;
  int event_count;
  
  MAT<std::vector<double>> death_kernel_y;
  VEC<std::vector<double>> birth_kernel_y;
  
  double spline_precision;
  
  double death_cutoff_r;
  double birth_cutoff_r;
  
  MAT<double> death_step;
  VEC<double> birth_step;
  VEC<double> birth_reverse_cdf_step;
  
  MAT<int> death_spline_nodes;
  VEC<int> birth_spline_nodes;
  VEC<int> birth_reverse_cdf_nodes;
  
  MAT<boost::math::cubic_b_spline<double>> death_kernel_spline;
  VEC<boost::math::cubic_b_spline<double>> birth_kernel_spline;
  
  VEC<boost::math::cubic_b_spline<double>> birth_reverse_cdf_spline;
  
  int cull_x;
  
  std::tuple<double, int> last_event;
  
  Cell_multi_1d &cell_at(int i)
  {
    return cells[i];
  }
  
  double &cell_death_rate_at(int s, int i)
  {
    return cell_death_rates[s][i];
  }
  
  int &cell_population_at(int s, int i)
  {
    return cell_population[s][i];
  }
  
  int cell_population_all_at(int i) {
    int res = 0;
    for (int s = 0; s < SPICIES_COUNT; ++s) {
      res += cell_population_at(s, i);
    }
    return res;
  }
  
  int get_all_population() {
    int res = 0;
    for (int i = 0; i < SPICIES_COUNT; ++i) {
      res += total_population[i];
    }
    return res;
  }
  
  vector<double> get_coords_at_cell(int i)
  {
    return cells[i].coords_x;
  }
  
  vector<double> get_death_rates_at_cell(int i)
  {
    return cells[i].death_rates;
  }
  
  vector<double> get_all_coords()
  {
    vector<double> result;
    for (auto cell : cells)
    {
      if (cell.coords_x.size() != 0)
        result.insert(result.end(), cell.coords_x.begin(), cell.coords_x.end());
    }
    return result;
  }
  
  vector<double> get_all_death_rates()
  {
    vector<double> result;
    for (auto cell : cells)
    {
      if (cell.coords_x.size() != 0)
        result.insert(result.end(), cell.death_rates.begin(), cell.death_rates.end());
    }
    return result;
  }
  
  void Initialize_death_rates()
  {
    
    for (int i = 0; i < cell_count_x; i++)
    {
      cells.push_back(Cell_multi_1d());
      for (int s = 0; s < SPICIES_COUNT; ++s) {
        cell_death_rates[s].push_back(0);
        cell_population[s].push_back(0);
      }
    }
    
    //Spawn all speciments
    for (size_t i = 0; i < SPICIES_COUNT; ++i) {
      total_population[i] = static_cast<int>(ceil(area_length_x * initial_density[i])); //initial population at t=0
    }
    for (size_t s = 0; s < SPICIES_COUNT; ++s) {
      double x_coord;
      int i;
      for (int k = 0; k < total_population[s]; k++)
      {
        x_coord = boost::uniform_real<>(0, area_length_x)(rng);
        
        i = static_cast<int>(floor(x_coord * cell_count_x / area_length_x));
        
        if (i == cell_count_x)
          i--;
        
        cell_at(i).coords_x.push_back(x_coord);
        
        cell_at(i).death_rates.push_back(d[s]);
        cell_at(i).spicies.push_back(s);
        cell_death_rate_at(s, i) += d[s];
        total_death_rate[s] += d[s];
        
        cell_population_at(s, i)++;
      }
    }
    
    for (int i = 0; i < cell_count_x; i++)
    {
      for (int k = 0; k < cell_population_all_at(i); k++)
      {
        for (int n = max(0, i - cull_x); n < min(cell_count_x, i + cull_x + 1); n++)
        {
          for (int p = 0; p < cell_population_all_at(n); p++)
          {
            if (i == n && k == p)
              continue; // same speciment
            
            
            int type1 = cell_at(i).spicies[k];
            int type2 = cell_at(n).spicies[p];
            
            //Distance between k-th speciment in (i) cell and p-th speciment in (n) cell
            double distance = abs(cell_at(i).coords_x[k] - cell_at(n).coords_x[p]);
            
            if (distance > death_cutoff_r)
              continue; //Too far to interact
            
            double interaction = dd[type1][type2] * death_kernel_spline[type1][type2](distance);
            
            cell_at(i).death_rates[k] += interaction;
            cell_death_rate_at(type1, i) += interaction;
            total_death_rate[type1] += interaction;
          }
        }
      }
    }
  }
  
  void kill_random(int s)
  {
    if (total_population[s] == 0) {
      return;
    }
    int cell_death_index = boost::random::discrete_distribution<>(cell_death_rates[s])(rng);
    auto tmp_vec = cells[cell_death_index].death_rates;
    for (int i = 0; i < tmp_vec.size(); ++i) {
      if (cell_at(cell_death_index).spicies[i] != s) {
        tmp_vec[i] = 0;
      }
    }
    int in_cell_death_index = boost::random::discrete_distribution<>(tmp_vec)(rng);
    
    Cell_multi_1d &death_cell = cells[cell_death_index];
    
    last_event = make_tuple<>(death_cell.coords_x[in_cell_death_index], -1);
    
    int cell_death_x = cell_death_index;
    
    for (int i = max(0, cell_death_x - cull_x); i < min(cell_count_x, cell_death_x + cull_x + 1); i++)
    {
      for (int k = 0; k < cell_population_all_at(i); k++)
      {
        if (i == cell_death_x && k == in_cell_death_index)
          continue;
        
        int type2 = cell_at(i).spicies[k];
        double distance = abs(cell_at(i).coords_x[k] - death_cell.coords_x[in_cell_death_index]);
        
        if (distance > death_cutoff_r)
          continue; //Too far to interact
        
        double interaction = dd[s][type2] * death_kernel_spline[s][type2](distance);
        
        cell_at(i).death_rates[k] -= interaction;
        //ignore dying speciment death rates since it is to be deleted
        
        cell_death_rate_at(type2, i) -= interaction;
        cell_death_rate_at(s, cell_death_x) -= interaction;
        
        total_death_rate[s] -= interaction;
        total_death_rate[type2] -= interaction;
      }
    }
    //remove dead speciment
    cell_death_rates[s][cell_death_index] -= d[s];
    total_death_rate[s] -= d[s];
    
    if (abs(cell_death_rates[s][cell_death_index]) < 1e-10)
    {
      cell_death_rates[s][cell_death_index] = 0;
    }
    
    cell_population[s][cell_death_index]--;
    total_population[s]--;
    
    //swap dead and last
    death_cell.death_rates[in_cell_death_index] = death_cell.death_rates[death_cell.death_rates.size() - 1];
    death_cell.coords_x[in_cell_death_index] = death_cell.coords_x[death_cell.coords_x.size() - 1];
    
    death_cell.death_rates.erase(death_cell.death_rates.end() - 1);
    death_cell.coords_x.erase(death_cell.coords_x.end() - 1);
  }
  
  void spawn_random(int s)
  {
    int cell_index = boost::random::discrete_distribution<>(cell_population[s])(rng);
    
    int event_index = boost::random::uniform_smallint<>(0, cell_population[s][cell_index] - 1)(rng);
    ++event_index;
    for (int i = 0; ; ++i) {
      if (cell_at(cell_index).spicies[i] == s) {
        --event_index;
      }
      if (event_index == 0) {
        event_index = i;
        break;
      }
    }
    
    Cell_multi_1d &parent_cell = cells[cell_index];
    
    double x_coord_new = parent_cell.coords_x[event_index] +
      birth_reverse_cdf_spline[s](boost::random::uniform_01<>()(rng)) * (boost::random::bernoulli_distribution<>(0.5)(rng) * 2 - 1);
    
    if (x_coord_new < 0 || x_coord_new > area_length_x)
    {
      last_event = make_tuple<>(x_coord_new, 0);
      //Speciment failed to spawn and died outside area boundaries
    }
    else
    {
      last_event = make_tuple<>(x_coord_new, 1);
      
      int new_i = static_cast<int>(floor(x_coord_new * cell_count_x / area_length_x));
      
      if (new_i >= cell_count_x)
        new_i = cell_count_x - 1;
      if (new_i < 0) {
        new_i = 0;
      }
      //cout << "new_i: " << new_i << endl;
      //New speciment is added to the end of vector
      cout << "Size: " << cell_at(new_i).coords_x.size() << endl;
      cell_at(new_i).coords_x.push_back(x_coord_new);
      cell_at(new_i).death_rates.push_back(d[s]);
      cell_at(new_i).spicies.push_back(s);
      
      cell_death_rate_at(s, new_i) += d[s];
      total_death_rate[s] += d[s];
      
      cell_population_at(s, new_i)++;
      total_population[s]++;
      
      for (int i = max(0, new_i - cull_x); i < min(cell_count_x, new_i + cull_x + 1); i++)
      {
        for (int k = 0; k < cell_population_all_at(i); k++)
        {
          if (i == new_i && k == cell_population_all_at(new_i) - 1)
            continue;
          int type2 = cell_at(i).spicies[k];
          double distance = abs(cell_at(i).coords_x[k] - x_coord_new);
          
          if (distance > death_cutoff_r)
            continue; //Too far to interact
          
          double interaction = dd[s][type2] * death_kernel_spline[s][type2](distance);
          
          cell_at(i).death_rates[k] += interaction;
          cell_at(new_i).death_rates[cell_population_all_at(new_i) - 1] += interaction;
          
          cell_death_rate_at(s, i) += interaction;
          cell_death_rate_at(s, new_i) += interaction;
          
          total_death_rate[s] += interaction;
          total_death_rate[type2] += interaction;
        }
      }
    }
  }
  
  void make_event()
  {
    if (get_all_population() == 0)
      return;
    event_count++;
    time += boost::random::exponential_distribution<>(total_population[0] * b[0] + total_death_rate[0] + total_population[1] * b[1] + total_death_rate[1])(rng);
    //Rolling event according to global birth \ death rate
    std::vector<double> dis;
    dis.resize(4);
    dis[0] = 1 - total_population[0] * b[0] / (total_population[0] * b[0] + total_death_rate[0]);
    dis[1] = 1 - dis[0];
    dis[2] = 1 - total_population[1] * b[1] / (total_population[1] * b[1] + total_death_rate[1]);
    dis[3] = 1 - dis[2];
    if (total_population[0] == 0) {
      dis[0] = 0;
      dis[1] = 0;
    }
    if (total_population[1] == 0) {
      dis[2] = 0;
      dis[3] = 0;
    }
    //cout << "DR1: " << total_death_rate[0] << "\nDR2: " << total_death_rate[1] << endl;
    //cout << dis[0] << ' ' << dis[1] << ' ' << dis[2] << " " << dis[3] << endl;
    auto t = boost::random::discrete_distribution<>(dis)(rng);
    //cout << "t: " << t << endl;
    switch (t) {
      case 0:
        kill_random(0);
        break;
      case 1:
        spawn_random(0);
        break;
      case 2:
        kill_random(1);
        break;
      case 3:
        spawn_random(1);
        break;
      default:
        throw "Unknown result";
    }
  }
  void run_events(int events)
  {
    if (events > 0)
    {
      for (int i = 0; i < events; i++)
      {
        make_event();
      }
    }
  }
  
  void run_for(double time)
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
  
  Grid_multy_1d(Rcpp::List params) :
    cells(),
    cell_death_rates(),
    cell_population(),
    time(),
    event_count(), 
    death_kernel_spline(),
    birth_kernel_spline(),
    birth_reverse_cdf_spline()
  {
    
    //Parse parameters
    
    Rcpp::Environment base("package:base");
    
    // Make function callable from C++
    Rcpp::Function print = base["print"];
    
    area_length_x = Rcpp::as<double>(params["area_length_x"]);
    cell_count_x = Rcpp::as<int>(params["cell_count_x"]);
    
    
    b[0] = Rcpp::as<double>(params["b1"]);
    b[1] = Rcpp::as<double>(params["b2"]);
    d[0] = Rcpp::as<double>(params["d1"]);
    d[1] = Rcpp::as<double>(params["d2"]);
    dd[0][0] = Rcpp::as<double>(params["dd11"]);
    dd[0][1] = Rcpp::as<double>(params["dd12"]);
    dd[1][0] = Rcpp::as<double>(params["dd21"]);
    dd[1][1] = Rcpp::as<double>(params["dd22"]);
    
    seed = Rcpp::as<int>(params["seed"]);
    rng = boost::random::lagged_fibonacci2281(uint32_t(seed));
    
    initial_density[0] = Rcpp::as<double>(params["init_density1"]);
    initial_density[1] = Rcpp::as<double>(params["init_density1"]);
    
    death_kernel_y[0][0] = Rcpp::as<vector<double>>(params["death_kernel_y11"]);
    death_kernel_y[0][1] = Rcpp::as<vector<double>>(params["death_kernel_y12"]);
    death_kernel_y[1][0] = Rcpp::as<vector<double>>(params["death_kernel_y21"]);
    death_kernel_y[1][1] = Rcpp::as<vector<double>>(params["death_kernel_y22"]);
    birth_kernel_y[0] = Rcpp::as<vector<double>>(params["birth_kernel_y1"]);
    birth_kernel_y[1] = Rcpp::as<vector<double>>(params["birth_kernel_y2"]);
    death_cutoff_r = Rcpp::as<double>(params["death_kernel_r"]);
    birth_cutoff_r = Rcpp::as<double>(params["birth_kernel_r"]);
    for (int i = 0; i < SPICIES_COUNT; ++i) {
      for (size_t j = 0; j < SPICIES_COUNT; ++j) {
        death_spline_nodes[i][j] = death_kernel_y[i][j].size();
        death_step[i][j] = death_cutoff_r / (death_spline_nodes[i][j] - 1);
        
        //Build death spline, ensure 0 derivative at 0 (symmetric) and endpoint (expected no death interaction further)
        //trim_spline(death_kernel_x, death_kernel_y, death_kernel_spline, death_spline_nodes, death_cutoff_r, spline_precision);
        death_kernel_spline[i][j] = cubic_b_spline<double>(
          death_kernel_y[i][j].begin(),
          death_kernel_y[i][j].end(),
          0,
          death_step[i][j],
          0,
          0
        );
      }
      birth_spline_nodes[i] = birth_kernel_y[i].size();
      birth_step[i] = birth_cutoff_r / (birth_spline_nodes[i] - 1);
      
      //Build birth spline
      birth_kernel_spline[i] = cubic_b_spline<double>(
        birth_kernel_y[i].begin(),
        birth_kernel_y[i].end(),
        0,
        birth_step[i],
        0,
        0
      );
      
      vector<double> x_quantile_1d_array(birth_spline_nodes[i]);
      vector<double> y_quantile_1d_array(birth_spline_nodes[i]);
      
      using boost::math::quadrature::trapezoidal;
      double approx_const = trapezoidal([=](double y) {return birth_kernel_spline[i](y);}, 0.0, birth_cutoff_r);
      
      using boost::math::tools::newton_raphson_iterate;
      
      for (int j = 0; j < birth_spline_nodes[i]; j++)
      {
        x_quantile_1d_array[j] = (double)j / (birth_spline_nodes[i] - 1);
        y_quantile_1d_array[j] =
          newton_raphson_iterate(
            [=](double y) {
              return make_tuple(
                trapezoidal(
                  [=](double z) {return birth_kernel_spline[i](z);},
                  0.0,
                  y
                ) / approx_const - x_quantile_1d_array[j],
                birth_kernel_spline[i](y) / approx_const);
              },
              1e-10,
              0.0,
              birth_cutoff_r,
              numeric_limits<double>::digits
          );
      }
      int k = 0;
      while (y_quantile_1d_array[k] < 1e-300)
      {
        k++;
      }
      
      vector<double> y_quantile_1d_array_temp(birth_spline_nodes[i] - k);
      
      for (int j = 0; j < birth_spline_nodes[i] - k; j++)
      {
        y_quantile_1d_array_temp[j] = y_quantile_1d_array[k + j];
      }
      
      birth_reverse_cdf_step[i] = 1.0/(y_quantile_1d_array_temp.size()-1);
      birth_reverse_cdf_nodes[i] = y_quantile_1d_array_temp.size();
      // Extrapolate last quantile element
      double right_derivative = (y_quantile_1d_array_temp[y_quantile_1d_array_temp.size()-2] - y_quantile_1d_array_temp[y_quantile_1d_array_temp.size()-3])/birth_reverse_cdf_step[i];
      
      y_quantile_1d_array_temp[y_quantile_1d_array_temp.size() - 1] =
        y_quantile_1d_array_temp[y_quantile_1d_array_temp.size() - 2]
        + right_derivative * birth_reverse_cdf_step[i];
      
      //Ensure correct derivative at 0, equal to 1/2*birth_kernel(0)
      birth_reverse_cdf_spline[i] = cubic_b_spline<double>(
        y_quantile_1d_array_temp.begin(),
        y_quantile_1d_array_temp.end(),
        0,
        birth_reverse_cdf_step[i],
        0.5 / birth_kernel_spline[i](0),
        2 * right_derivative);
    }
    
    spline_precision = Rcpp::as<double>(params["spline_precision"]);
    

    
    //Calculate amount of cells to check around for death interaction
    cull_x = max(static_cast<int>(ceil(death_cutoff_r / (area_length_x / cell_count_x))), 3);
    

    //trim_spline(birth_kernel_x, birth_kernel_y, birth_kernel_spline, birth_spline_nodes, birth_cutoff_r, spline_precision);
    //Calculate reverse CDF for birth spline
    
    
    //Spawn speciments and calculate death rates
    Initialize_death_rates();
  }
};

RCPP_MODULE(poisson_2spicies_1d_module)
{
  using namespace Rcpp;
  
  class_<Grid_multy_1d>("poisson_2spicies_1d")
    .constructor<List>()
    .field_readonly("area_length_x", &Grid_multy_1d::area_length_x)
    .field_readonly("cell_count_x", &Grid_multy_1d::cell_count_x)
  
  .field_readonly("b", &Grid_multy_1d::b)
  .field_readonly("d", &Grid_multy_1d::d)
  .field_readonly("dd", &Grid_multy_1d::dd)
  
  .field_readonly("seed", &Grid_multy_1d::seed)
  .field_readonly("initial_density", &Grid_multy_1d::initial_density)
  
  .field_readonly("death_kernel_y", &Grid_multy_1d::death_kernel_y)
  .field_readonly("death_cutoff_r", &Grid_multy_1d::death_cutoff_r)
  .field_readonly("death_spline_nodes", &Grid_multy_1d::death_spline_nodes)
  .field_readonly("death_step", &Grid_multy_1d::death_step)
  
  .field_readonly("birth_kernel_y", &Grid_multy_1d::birth_kernel_y)
  .field_readonly("birth_cutoff_r", &Grid_multy_1d::birth_cutoff_r)
  .field_readonly("birth_spline_nodes", &Grid_multy_1d::birth_spline_nodes)
  .field_readonly("birth_step", &Grid_multy_1d::birth_step)
  
  .field_readonly("birth_reverse_cdf_nodes", &Grid_multy_1d::birth_reverse_cdf_nodes)
  .field_readonly("birth_reverse_cdf_step", &Grid_multy_1d::birth_reverse_cdf_step)
  
  .field_readonly("spline_precision", &Grid_multy_1d::spline_precision)
  
  .field_readonly("cell_death_rates", &Grid_multy_1d::cell_death_rates)
  .field_readonly("cell_population", &Grid_multy_1d::cell_population)
  
  .method("get_all_coordinates", &Grid_multy_1d::get_all_coords)
  .method("get_all_death_rates", &Grid_multy_1d::get_all_death_rates)
  
  .method("get_x_coordinates_in_cell", &Grid_multy_1d::get_coords_at_cell)
  .method("get_death_rates_in_cell", &Grid_multy_1d::get_death_rates_at_cell)
  
  .method("make_event", &Grid_multy_1d::make_event)
  .method("run_events", &Grid_multy_1d::run_events)
  .method("run_for", &Grid_multy_1d::run_for)
  
  .field_readonly("total_population", &Grid_multy_1d::total_population)
  .field_readonly("total_death_rate", &Grid_multy_1d::total_death_rate)
  .field_readonly("events", &Grid_multy_1d::event_count)
  .field_readonly("time", &Grid_multy_1d::time);
}

#endif
