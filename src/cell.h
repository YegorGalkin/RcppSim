#ifndef CELL
#define CELL

#include "defines.h"
#include "unit.h"

struct Cell {
  std::vector<DCoord> coords_x;
  std::vector<double> death_rates;
  std::vector<int> species;
  
  void SwapWithLast(int i);
  void Pop();
  void Add(double deathRate, DCoord x, int species);
  static DCoord Ro(const Unit& a, const Unit& b);
};

#endif
