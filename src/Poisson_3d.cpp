// [[Rcpp::depends(BH)]]
// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <stdio.h>
#include <chrono>
#include <vector>
#include <math.h>
#include <boost/random.hpp>
#include <boost/random/lagged_fibonacci.hpp>
#include <boost/random/exponential_distribution.hpp>
#include <boost/math/interpolators/cardinal_cubic_b_spline.hpp>


using namespace std;

#ifndef POISSON_3D_H
#define POISSON_3D_H
#define _USE_MATH_DEFINES

struct Cell_3d {
  
  vector < double > coords_x;
  vector < double > coords_y;
  vector < double > coords_z;
  vector < double > death_rates;
  
  Cell_3d() {}
};

struct Grid_3d {
  std::vector < Cell_3d > cells;
  std::vector < double > cell_death_rates;
  std::vector < int > cell_population;
  
  double area_length_x;
  double area_length_y;
  double area_length_z;

  bool periodic;
  
  int cell_count_x;
  int cell_count_y;
  int cell_count_z;
    
  int cull_x;
  int cull_y;
  int cull_z;

  double b, d, dd;
  int seed;
  boost::random::lagged_fibonacci2281 rng;
  
  std::vector < double > initial_population_x;
  std::vector < double > initial_population_y;
  std::vector < double > initial_population_z;

  int total_population;

  chrono::system_clock::time_point init_time;
  double realtime_limit;
  bool realtime_limit_reached;

  double total_death_rate;
  
  double time;
  int event_count;
  
  std::vector<double> death_y;
  double death_cutoff_r;
  double death_step;
  int death_spline_nodes;
  boost::math::interpolators::cardinal_cubic_b_spline<double> death_spline;

  std::vector<double> birth_inverse_rcdf_y;
  double birth_inverse_rcdf_step;
  int birth_inverse_rcdf_nodes;
  boost::math::interpolators::cardinal_cubic_b_spline<double> birth_inverse_rcdf_spline;
  
  Cell_3d & cell_at(int i,int j,int k) {
    if (periodic) {
      if (i < 0) i += cell_count_x;
      if (i >= cell_count_x) i -= cell_count_x;

      if (j < 0) j += cell_count_y;
      if (j >= cell_count_y) j -= cell_count_y;

      if (k < 0) k += cell_count_z;
      if (k >= cell_count_z) k -= cell_count_z;
    }
    return cells[i + cell_count_y * (j + cell_count_z * k)];
  }
  
  double & cell_death_rate_at(int i,int j, int k) {
    if (periodic) {
      if (i < 0) i += cell_count_x;
      if (i >= cell_count_x) i -= cell_count_x;

      if (j < 0) j += cell_count_y;
      if (j >= cell_count_y) j -= cell_count_y;

      if (k < 0) k += cell_count_z;
      if (k >= cell_count_z) k -= cell_count_z;
    }
    return cell_death_rates[i + cell_count_y * (j + cell_count_z * k)];
  }
  
  int & cell_population_at(int i,int j, int k) {
    if (periodic) {
      if (i < 0) i += cell_count_x;
      if (i >= cell_count_x) i -= cell_count_x;

      if (j < 0) j += cell_count_y;
      if (j >= cell_count_y) j -= cell_count_y;

      if (k < 0) k += cell_count_z;
      if (k >= cell_count_z) k -= cell_count_z;
    }
    return cell_population[i + cell_count_y * (j + cell_count_z * k)];
  }
  
  
  vector < double > get_x_coords_at_cell(int i,int j, int k)  {
    return cells[i + cell_count_y * (j + cell_count_z * k)].coords_x;
  }
  vector < double > get_y_coords_at_cell(int i,int j, int k)  {
    return cells[i + cell_count_y * (j + cell_count_z * k)].coords_y;
  }
  vector < double > get_z_coords_at_cell(int i,int j, int k)  {
    return cells[i + cell_count_y * (j + cell_count_z * k)].coords_z;
  }
  vector < double > get_death_rates_at_cell(int i,int j, int k)  {
    return cells[i + cell_count_y * (j + cell_count_z * k)].death_rates;
  }
  
  vector < double > get_all_x_coords() {
    vector < double > result;
    for (auto cell: cells) {
      if (cell.coords_x.size() != 0)
        result.insert(result.end(), cell.coords_x.begin(), cell.coords_x.end());
    }
    return result;
  }
  
  vector < double > get_all_y_coords() {
    vector < double > result;
    for (auto cell: cells) {
      if (cell.coords_y.size() != 0)
        result.insert(result.end(), cell.coords_y.begin(), cell.coords_y.end());
    }
    return result;
  }

  vector < double > get_all_z_coords() {
    vector < double > result;
    for (auto cell: cells) {
      if (cell.coords_z.size() != 0)
        result.insert(result.end(), cell.coords_z.begin(), cell.coords_z.end());
    }
    return result;
  }

  vector < double > get_all_death_rates() {
    vector < double > result;
    for (auto cell: cells) {
      if (cell.coords_x.size() != 0)
        result.insert(result.end(), cell.death_rates.begin(), cell.death_rates.end());
    }
    return result;
  }
  
  void Initialize_death_rates() {
    
    for (int i = 0; i < cell_count_x*cell_count_y*cell_count_z; i++) {
      cells.push_back(Cell_3d());
      cell_death_rates.push_back(0);
      cell_population.push_back(0);
    }
    
    //Spawn all speciments
    {
      for (int sp_index=0;sp_index< initial_population_x.size();sp_index++) {
        double x_coord=initial_population_x[sp_index];
        double y_coord=initial_population_y[sp_index];
        double z_coord=initial_population_z[sp_index];

        if (x_coord < 0 || x_coord > area_length_x) continue;
        if (y_coord < 0 || y_coord > area_length_y) continue;
        if (z_coord < 0 || z_coord > area_length_z) continue;

        int i = static_cast < int > (floor(x_coord * cell_count_x / area_length_x));
        int j = static_cast < int > (floor(y_coord * cell_count_y / area_length_y));
        int k = static_cast < int > (floor(z_coord * cell_count_z / area_length_z));

        if (i == cell_count_x) i--;
        if (j == cell_count_y) j--;
        if (k == cell_count_z) k--;

        cell_at(i,j,k).coords_x.push_back(x_coord);
        cell_at(i,j,k).coords_y.push_back(y_coord);
        cell_at(i,j,k).coords_z.push_back(z_coord);

        cell_at(i,j,k).death_rates.push_back(d);
        cell_death_rate_at(i,j,k) += d;
        total_death_rate += d;
        
        cell_population_at(i,j,k) ++;
        total_population++;
      }
    }
    
    for (int i = 0; i < cell_count_x; i++) {
      for (int j = 0; j < cell_count_y; j++) {
        for (int k = 0; k < cell_count_z; k++){
          for (int w = 0; w < cell_population_at(i,j,k); w++) {

            for (int n = i - cull_x; n < i + cull_x + 1; n++) {
              if (!periodic && (n < 0 || n >= cell_count_x)) continue;

              for (int m = j - cull_y; m< j + cull_y + 1; m++){
                if (!periodic && (m < 0 || m >= cell_count_y)) continue; 

                for (int p = k - cull_z; p< k + cull_z + 1; p++){
                  if (!periodic && (p < 0 || p >= cell_count_z)) continue; 
                  
                  for (int q = 0; q < cell_population_at(n,m,p); q++) {
                    if (i == n && j == m && k==p && w == q) continue; // same speciment

                    double delta_x, delta_y, delta_z, distance;
            
                    if (periodic) {
                      if (n < 0) {
                        delta_x = cell_at(i,j,k).coords_x[w] - cell_at(n,m,p).coords_x[q] + area_length_x;
                      } else if (n >= cell_count_x) {
                        delta_x = cell_at(i,j,k).coords_x[w] - cell_at(n,m,p).coords_x[q] - area_length_x;
                      } else {
                        delta_x = cell_at(i,j,k).coords_x[w] - cell_at(n,m,p).coords_x[q];
                      }

                      if (m < 0) {
                        delta_y = cell_at(i,j,k).coords_y[w] - cell_at(n,m,p).coords_y[q] + area_length_y;
                      } else if (m >= cell_count_y) {
                        delta_y = cell_at(i,j,k).coords_y[w] - cell_at(n,m,p).coords_y[q] - area_length_y;
                      } else {
                        delta_y = cell_at(i,j,k).coords_y[w] - cell_at(n,m,p).coords_y[q];
                      }

                      if (p < 0) {
                        delta_z = cell_at(i,j,k).coords_z[w] - cell_at(n,m,p).coords_z[q] + area_length_z;
                      } else if (p >= cell_count_z) {
                        delta_z = cell_at(i,j,k).coords_z[w] - cell_at(n,m,p).coords_z[q] - area_length_z;
                      } else {
                        delta_z = cell_at(i,j,k).coords_z[w] - cell_at(n,m,p).coords_z[q];
                      }
                    }else {
                      delta_x = cell_at(i,j,k).coords_x[w] - cell_at(n,m,p).coords_x[q];
                      delta_y = cell_at(i,j,k).coords_y[w] - cell_at(n,m,p).coords_y[q];
                      delta_z = cell_at(i,j,k).coords_z[w] - cell_at(n,m,p).coords_z[q];
                    }
                    distance =sqrt(delta_x*delta_x+delta_y*delta_y+delta_z*delta_z);
                    if (distance > death_cutoff_r)
                      continue; //Too far to interact
            
                    double interaction = dd * death_spline(distance);
            
                    cell_at(i,j,k).death_rates[w] += interaction;
                    cell_death_rate_at(i,j,k) += interaction;
                    total_death_rate += interaction;
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  void kill_random() {
    
    if (total_population == 1) {
      total_population--;
      return;
    }
    
    int cell_death_index = boost::random::discrete_distribution < > (cell_death_rates)(rng);
    int in_cell_death_index = boost::random::discrete_distribution < > (cells[cell_death_index].death_rates)(rng);
    
    Cell_3d & death_cell = cells[cell_death_index];

    int cell_death_z = cell_death_index / (cell_count_x*cell_count_y);
    int cell_death_y = (cell_death_index - cell_death_z * cell_count_x * cell_count_y) / cell_count_x;
    int cell_death_x = (cell_death_index - cell_death_z * cell_count_x * cell_count_y) % cell_count_x;
    
    
    for (int i = cell_death_x - cull_x; i < cell_death_x + cull_x + 1; i++) {
      if (!periodic && (i < 0 || i >= cell_count_x)) continue;

      for (int j = cell_death_y - cull_y; j < cell_death_y + cull_y + 1; j++){
        if (!periodic && (j < 0 || j >= cell_count_y)) continue;  

        for (int k = cell_death_z - cull_z; k < cell_death_z + cull_z + 1; k++){
          if (!periodic && (k < 0 || k >= cell_count_z)) continue; 

          for (int w = 0; w < cell_population_at(i,j,k); w++) {
            if (i == cell_death_x && j == cell_death_y && k == cell_death_z && w == in_cell_death_index) continue;
        
            double delta_x, delta_y, delta_z, distance;
        
            if (periodic) {
              if (i < 0) {
                delta_x = cell_at(cell_death_x,cell_death_y,cell_death_z).coords_x[in_cell_death_index] - cell_at(i,j,k).coords_x[w] + area_length_x;
              } else if (i >= cell_count_x) {
                delta_x = cell_at(cell_death_x,cell_death_y,cell_death_z).coords_x[in_cell_death_index] - cell_at(i,j,k).coords_x[w] - area_length_x;
              } else {
                delta_x = cell_at(cell_death_x,cell_death_y,cell_death_z).coords_x[in_cell_death_index] - cell_at(i,j,k).coords_x[w];
              }

              if (j < 0) {
                delta_y = cell_at(cell_death_x,cell_death_y,cell_death_z).coords_y[in_cell_death_index] - cell_at(i,j,k).coords_y[w] + area_length_y;
              } else if (j >= cell_count_y) {
                delta_y = cell_at(cell_death_x,cell_death_y,cell_death_z).coords_y[in_cell_death_index] - cell_at(i,j,k).coords_y[w] - area_length_y;
              } else {
                delta_y = cell_at(cell_death_x,cell_death_y,cell_death_z).coords_y[in_cell_death_index] - cell_at(i,j,k).coords_y[w];
              }

              if (k < 0) {
                delta_z = cell_at(cell_death_x,cell_death_y,cell_death_z).coords_z[in_cell_death_index] - cell_at(i,j,k).coords_z[w] + area_length_z;
              } else if (k >= cell_count_z) {
                delta_z = cell_at(cell_death_x,cell_death_y,cell_death_z).coords_z[in_cell_death_index] - cell_at(i,j,k).coords_z[w] - area_length_z;
              } else {
                delta_z = cell_at(cell_death_x,cell_death_y,cell_death_z).coords_z[in_cell_death_index] - cell_at(i,j,k).coords_z[w];
              }
            } else {
              delta_x = cell_at(cell_death_x,cell_death_y,cell_death_z).coords_x[in_cell_death_index] - cell_at(i,j,k).coords_x[w];
              delta_y = cell_at(cell_death_x,cell_death_y,cell_death_z).coords_y[in_cell_death_index] - cell_at(i,j,k).coords_y[w];
              delta_z = cell_at(cell_death_x,cell_death_y,cell_death_z).coords_z[in_cell_death_index] - cell_at(i,j,k).coords_z[w];
            }
            distance =sqrt(delta_x*delta_x+delta_y*delta_y+delta_z*delta_z);
        
            if (distance > death_cutoff_r) continue; //Too far to interact
        
            double interaction = dd * death_spline(distance);
        
            cell_at(i,j,k).death_rates[w] -= interaction;
            //ignore dying speciment death rates since it is to be deleted
        
            cell_death_rate_at(i,j,k) -= interaction;
            cell_death_rate_at(cell_death_x,cell_death_y,cell_death_z) -= interaction;
        
            total_death_rate -= 2 * interaction;
          }
        }
      }
    }
    //remove dead speciment
    cell_death_rates[cell_death_index] -= d;
    total_death_rate -= d;
    
    if (abs(cell_death_rates[cell_death_index]) < 1e-10) {
      cell_death_rates[cell_death_index] = 0;
    }
    
    cell_population[cell_death_index]--;
    total_population--;
    
    //swap dead and last
    death_cell.death_rates[in_cell_death_index] = death_cell.death_rates[death_cell.death_rates.size() - 1];
    death_cell.coords_x[in_cell_death_index] = death_cell.coords_x[death_cell.coords_x.size() - 1];
    death_cell.coords_y[in_cell_death_index] = death_cell.coords_y[death_cell.coords_y.size() - 1];
    death_cell.coords_z[in_cell_death_index] = death_cell.coords_z[death_cell.coords_z.size() - 1];

    death_cell.death_rates.erase(death_cell.death_rates.end() - 1);
    death_cell.coords_x.erase(death_cell.coords_x.end() - 1);
    death_cell.coords_y.erase(death_cell.coords_y.end() - 1);
    death_cell.coords_z.erase(death_cell.coords_z.end() - 1);
  }
  
  void spawn_random() {
    int cell_index = boost::random::discrete_distribution < > (cell_population)(rng);
    int event_index = boost::random::uniform_smallint < > (0, cell_population[cell_index] - 1)(rng);
    Cell_3d & parent_cell = cells[cell_index];
    
    // Generates x, y and z coordinates uniformly distributed on a sphere
    // See https://stackoverflow.com/questions/38038776/computationally-picking-a-random-point-on-a-n-sphere
    double x_normal = boost::random::normal_distribution<>()(rng);
    double y_normal = boost::random::normal_distribution<>()(rng);
    double z_normal = boost::random::normal_distribution<>()(rng);

    // Uses radius distribution to create spherical simmetrical function
    double r_dispacement = birth_inverse_rcdf_spline(boost::random::uniform_01 < > ()(rng));

    double x_dispacement = x_normal / sqrt(x_normal*x_normal + y_normal*y_normal + z_normal*z_normal) * r_dispacement;
    double y_dispacement = y_normal / sqrt(x_normal*x_normal + y_normal*y_normal + z_normal*z_normal) * r_dispacement;
    double z_dispacement = z_normal / sqrt(x_normal*x_normal + y_normal*y_normal + z_normal*z_normal) * r_dispacement;

    double x_coord_new = parent_cell.coords_x[event_index] + x_dispacement;
    double y_coord_new = parent_cell.coords_y[event_index] + y_dispacement;
    double z_coord_new = parent_cell.coords_z[event_index] + z_dispacement;

    if (x_coord_new < 0 || x_coord_new > area_length_x || y_coord_new <0 || y_coord_new > area_length_y || z_coord_new <0 || z_coord_new > area_length_z) {
      if (!periodic) {
        return;
      } else {
        if (x_coord_new < 0)             x_coord_new += area_length_x;
        if (x_coord_new > area_length_x) x_coord_new -= area_length_x;

        if (y_coord_new < 0)             y_coord_new += area_length_y;
        if (y_coord_new > area_length_y) y_coord_new -= area_length_y;

        if (z_coord_new < 0)             z_coord_new += area_length_z;
        if (z_coord_new > area_length_z) z_coord_new -= area_length_z;
      }
    }
    
    
    int new_i = static_cast<int>(floor(x_coord_new*cell_count_x / area_length_x));
    int new_j = static_cast<int>(floor(y_coord_new*cell_count_y / area_length_y));
    int new_k = static_cast<int>(floor(z_coord_new*cell_count_z / area_length_z));

    if (new_i == cell_count_x) new_i--;
    if (new_j == cell_count_y) new_j--;
    if (new_k == cell_count_z) new_k--;

    //New speciment is added to the end of vector
    
    cell_at(new_i,new_j,new_k).coords_x.push_back(x_coord_new);
    cell_at(new_i,new_j,new_k).coords_y.push_back(y_coord_new);
    cell_at(new_i,new_j,new_k).coords_z.push_back(z_coord_new);
    cell_at(new_i,new_j,new_k).death_rates.push_back(d);
    
    cell_death_rate_at(new_i,new_j,new_k) += d;
    total_death_rate += d;
    
    cell_population_at(new_i,new_j,new_k) ++;
    total_population++;
    
    for (int i = new_i - cull_x; i < new_i + cull_x + 1; i++) {
      if (!periodic && (i < 0 || i >= cell_count_x)) continue;

      for (int j = new_j - cull_y; j < new_j + cull_y + 1; j++) {
        if (!periodic && (j < 0 ||j >= cell_count_y)) continue;

        for (int k = new_k - cull_z; k < new_k + cull_z + 1; k++) {
          if (!periodic && (k < 0 ||k >= cell_count_z)) continue;

          for (int w = 0; w < cell_population_at(i,j,k); w++) {
            if (i == new_i && j == new_j && k == new_k && w == cell_population_at(new_i, new_j, new_k) - 1) continue;
        
            double delta_x,delta_y,delta_z,distance;
        
            if (periodic) {
              if (i < 0) {
                delta_x = cell_at(new_i,new_j,new_k).coords_x[cell_population_at(new_i,new_j,new_k) - 1] - cell_at(i,j,k).coords_x[w] + area_length_x;
              } else if (i >= cell_count_x) {
                delta_x = cell_at(new_i,new_j,new_k).coords_x[cell_population_at(new_i,new_j,new_k) - 1] - cell_at(i,j,k).coords_x[w] - area_length_x;
              } else {
                delta_x = cell_at(new_i,new_j,new_k).coords_x[cell_population_at(new_i,new_j,new_k) - 1] - cell_at(i,j,k).coords_x[w];
              }

              if (j < 0) {
                delta_y = cell_at(new_i,new_j,new_k).coords_y[cell_population_at(new_i,new_j,new_k) - 1] - cell_at(i,j,k).coords_y[w] + area_length_y;
              } else if (j >= cell_count_y) {
                delta_y = cell_at(new_i,new_j,new_k).coords_y[cell_population_at(new_i,new_j,new_k) - 1] - cell_at(i,j,k).coords_y[w] - area_length_y;
              } else {
                delta_y = cell_at(new_i,new_j,new_k).coords_y[cell_population_at(new_i,new_j,new_k) - 1] - cell_at(i,j,k).coords_y[w];
              }

              if (k < 0) {
                delta_z = cell_at(new_i,new_j,new_k).coords_z[cell_population_at(new_i,new_j,new_k) - 1] - cell_at(i,j,k).coords_z[w] + area_length_z;
              } else if (k >= cell_count_z) {
                delta_z = cell_at(new_i,new_j,new_k).coords_z[cell_population_at(new_i,new_j,new_k) - 1] - cell_at(i,j,k).coords_z[w] - area_length_z;
              } else {
                delta_z = cell_at(new_i,new_j,new_k).coords_z[cell_population_at(new_i,new_j,new_k) - 1] - cell_at(i,j,k).coords_z[w];
              }
            } else {
              delta_x = cell_at(new_i,new_j,new_k).coords_x[cell_population_at(new_i,new_j,new_k) - 1] - cell_at(i,j,k).coords_x[w];
              delta_y = cell_at(new_i,new_j,new_k).coords_y[cell_population_at(new_i,new_j,new_k) - 1] - cell_at(i,j,k).coords_y[w];
              delta_z = cell_at(new_i,new_j,new_k).coords_z[cell_population_at(new_i,new_j,new_k) - 1] - cell_at(i,j,k).coords_z[w];

            }

            distance =sqrt(delta_x*delta_x+delta_y*delta_y+delta_z*delta_z);
            if (distance > death_cutoff_r) continue; //Too far to interact
        
            double interaction = dd * death_spline(distance);
        
            cell_at(i,j,k).death_rates[w] += interaction;
            cell_at(new_i,new_j,new_k).death_rates[cell_population_at(new_i,new_j,new_k) - 1] += interaction;
        
            cell_death_rate_at(i,j,k) += interaction;
            cell_death_rate_at(new_i,new_j,new_k) += interaction;
        
            total_death_rate += 2 * interaction;
          }
        }
      }
    }
    
  }
  
  void make_event() {
    if (total_population == 0)
      return;
    event_count++;
    time += boost::random::exponential_distribution < > (total_population * b + total_death_rate)(rng);
    //Rolling event according to global birth \ death rate
    if (boost::random::bernoulli_distribution < > (total_population * b / (total_population * b + total_death_rate))(rng) == 0) {
      kill_random();
    } else {
      spawn_random();
    }
  }


void run_events(int events)
  {
    if (events > 0)
    {
      for (int i = 0; i < events; i++)
      {
        if (chrono::system_clock::now() > init_time + chrono::duration<double>(realtime_limit))
        {
          realtime_limit_reached = true;
          return;
        }
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
        if (chrono::system_clock::now() > init_time + chrono::duration<double>(realtime_limit))
        {
          realtime_limit_reached = true;
          return;
        }
        make_event();
      }
    }
  }
  
  double get_death_spline_value(double at) {
    return death_spline(at);
  }
  
  double get_birth_inverse_rcdf_spline_value(double at) {
    return birth_inverse_rcdf_spline(at);
  }
  
  Grid_3d(Rcpp::List params): time(), event_count(), 
  cells(), cell_death_rates(), cell_population(),
  death_spline(), birth_inverse_rcdf_spline(), 
  total_population(), realtime_limit_reached(false) {
    
    //Parse parameters
    
    area_length_x = Rcpp::as < double > (params["area_length_x"]);
    area_length_y = Rcpp::as < double > (params["area_length_y"]);
    area_length_z = Rcpp::as < double > (params["area_length_z"]);

    cell_count_x = Rcpp::as < int > (params["cell_count_x"]);
    cell_count_y = Rcpp::as < int > (params["cell_count_y"]);
    cell_count_z = Rcpp::as < int > (params["cell_count_z"]);

    b = Rcpp::as < double > (params["b"]);
    d = Rcpp::as < double > (params["d"]);
    dd = Rcpp::as < double > (params["dd"]);
    
    seed = Rcpp::as < int > (params["seed"]);
    rng = boost::random::lagged_fibonacci2281(uint32_t(seed));
    
    initial_population_x = Rcpp::as < vector < double >> (params["initial_population_x"]);
    initial_population_y = Rcpp::as < vector < double >> (params["initial_population_y"]);
    initial_population_z = Rcpp::as < vector < double >> (params["initial_population_z"]);

    death_y = Rcpp::as < vector < double >> (params["death_y"]);
    death_cutoff_r = Rcpp::as < double > (params["death_r"]);
    death_spline_nodes = death_y.size();
    death_step = death_cutoff_r / (death_spline_nodes - 1);
    
    birth_inverse_rcdf_y = Rcpp::as<vector<double>>(params["birth_ircdf_y"]);
    birth_inverse_rcdf_nodes = birth_inverse_rcdf_y.size();
    birth_inverse_rcdf_step = 1.0 / (birth_inverse_rcdf_nodes - 1);
    
    periodic = Rcpp::as < bool > (params["periodic"]);

    init_time = chrono::system_clock::now();
    realtime_limit = Rcpp::as<double>(params["realtime_limit"]);
    
    using boost::math::interpolators::cardinal_cubic_b_spline;
    //Build death spline, ensure 0 derivative at 0 (symmetric) and endpoint (expected no death interaction further)
    death_spline = cardinal_cubic_b_spline < double > (death_y.begin(), death_y.end(), 0, death_step, 0, 0);
     
    //Calculate amount of cells to check around for death interaction
    
    cull_x = max(static_cast < int > (ceil(death_cutoff_r / (area_length_x / cell_count_x))), 3);
    cull_y = max(static_cast < int > (ceil(death_cutoff_r / (area_length_y / cell_count_y))), 3);
    cull_z = max(static_cast < int > (ceil(death_cutoff_r / (area_length_z / cell_count_z))), 3);

    //Build birth inverse rcdf spline, endpoint derivatives not specified
    birth_inverse_rcdf_spline = cardinal_cubic_b_spline<double>(birth_inverse_rcdf_y.begin(), birth_inverse_rcdf_y.end(), 0, birth_inverse_rcdf_step);

    //Spawn speciments and calculate death rates
    Initialize_death_rates();
    total_death_rate = accumulate(cell_death_rates.begin(), cell_death_rates.end(), 0.0);
  }
};

RCPP_EXPOSED_CLASS(poisson_3d)
RCPP_MODULE(poisson_3d_module) {
  using namespace Rcpp;
  
  class_ < Grid_3d > ("poisson_3d")
  .constructor < List > ("Creates an instance of 3d simulator")
  .field_readonly("area_length_x", & Grid_3d::area_length_x)
  .field_readonly("area_length_y", & Grid_3d::area_length_y)
  .field_readonly("area_length_z", & Grid_3d::area_length_z)
  .field_readonly("cell_count_x", & Grid_3d::cell_count_x)
  .field_readonly("cell_count_y", & Grid_3d::cell_count_y)
  .field_readonly("cell_count_z", & Grid_3d::cell_count_z)

  .field_readonly("cull_x", & Grid_3d::cull_x)
  .field_readonly("cull_y", & Grid_3d::cull_y)
  .field_readonly("cull_z", & Grid_3d::cull_z)
  .field_readonly("periodic", & Grid_3d::periodic)
  
  .field_readonly("b", & Grid_3d::b)
  .field_readonly("d", & Grid_3d::d)
  .field_readonly("dd", & Grid_3d::dd)
  
  .field_readonly("seed", & Grid_3d::seed)
  .field_readonly("initial_population_x", & Grid_3d::initial_population_x)
  .field_readonly("initial_population_y", & Grid_3d::initial_population_y)
  .field_readonly("initial_population_z", & Grid_3d::initial_population_y)

  .field_readonly("death_y", & Grid_3d::death_y)
  .field_readonly("death_cutoff_r", & Grid_3d::death_cutoff_r)
  .field_readonly("death_spline_nodes", & Grid_3d::death_spline_nodes)
  .field_readonly("death_step", & Grid_3d::death_step)
  
  .field_readonly("birth_inverse_rcdf_y", &Grid_3d::birth_inverse_rcdf_y)
  .field_readonly("birth_inverse_rcdf_nodes", &Grid_3d::birth_inverse_rcdf_nodes)
  .field_readonly("birth_inverse_rcdf_step", &Grid_3d::birth_inverse_rcdf_step)
  
  .field_readonly("cell_death_rates", & Grid_3d::cell_death_rates)
  .field_readonly("cell_population", & Grid_3d::cell_population)
  
  .method("get_all_x_coordinates", & Grid_3d::get_all_x_coords)
  .method("get_all_y_coordinates", & Grid_3d::get_all_y_coords)
  .method("get_all_z_coordinates", & Grid_3d::get_all_z_coords)
  .method("get_all_death_rates", & Grid_3d::get_all_death_rates)
  
  .method("get_x_coordinates_in_cell", & Grid_3d::get_x_coords_at_cell)
  .method("get_y_coordinates_in_cell", & Grid_3d::get_y_coords_at_cell)
  .method("get_z_coordinates_in_cell", & Grid_3d::get_z_coords_at_cell)

  .method("get_death_rates_in_cell", & Grid_3d::get_death_rates_at_cell)
  
  .method("death_spline_at", & Grid_3d::get_death_spline_value)
  .method("birth_inverse_rcdf_spline_at", & Grid_3d::get_birth_inverse_rcdf_spline_value)
  
  .method("make_event", & Grid_3d::make_event)
  .method("run_events", & Grid_3d::run_events)
  .method("run_for", & Grid_3d::run_for)
  
  .field_readonly("total_population", & Grid_3d::total_population)
  .field_readonly("total_death_rate", & Grid_3d::total_death_rate)
  .field_readonly("events", & Grid_3d::event_count)
  .field_readonly("time", & Grid_3d::time)
  
  .field_readonly("realtime_limit", &Grid_3d::realtime_limit)
  .field_readonly("realtime_limit_reached", &Grid_3d::realtime_limit_reached);

}

# endif
