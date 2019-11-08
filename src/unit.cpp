#include "unit.h"


Unit::Unit(
  struct Cell& cell,
  MAT<double>& cellsDeathRates,
  MAT<int>& cellsPopulation,
  int cellNum,
  int i
)
  : cell(cell)
, _CellDeathRate(cellsDeathRates[cell.species[i]][cellNum])
, _CellPopulation(cellsPopulation[cell.species[i]][cellNum])
, _CellNum(cellNum)
, I(i)
{ }

Unit::Unit(Grid& grid, int i, int j)
  : cell(grid.cells[i])
  , _CellDeathRate(grid.cell_death_rates[cell.species[j]][i])
  , _CellPopulation(grid.cell_population[cell.species[j]][i])
  , _CellNum(i)
  , I(j)
{
  assert(_CellPopulation >= 0);
  }

DCoord Unit::Coord() const {
  return cell.coords_x[I];
}

double& Unit::DeathRate() {
  return cell.death_rates[I];
}

int Unit::Species() const {
  return cell.species[I];
}

double& Unit::CellDeathRate() {
  return _CellDeathRate;
}

int& Unit::CellPopulation() {
  return _CellPopulation;
}

bool Unit::operator==(const Unit& a) const {
  return &cell == &a.cell && I == a.I;
}

int Unit::CellNum() const {
  return _CellNum;
}
