#ifndef UNIT
#define UNIT

#include "defines.h"
#include "cell.h"
#include "grid.h"

class Unit {
  Cell& cell;
  double& _CellDeathRate;
  int& _CellPopulation;
  const int _CellNum;
  const int I;
  
public:
  Unit(
    struct Cell& cell,
    MAT<double>& cellsDeathRates,
    MAT<int>& cellsPopulation,
    int cellNum,
    int i
  );
  
  Unit(class Grid& grid, int i, int j);
  
  DCoord Coord() const;
  
  double& DeathRate();
  
  int Species() const;
  
  double& CellDeathRate();
  
  int& CellPopulation();
  
  bool operator==(const Unit& a) const;
  
  int CellNum() const;
};

#endif
