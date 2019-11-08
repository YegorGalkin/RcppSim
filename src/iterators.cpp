#include "iterators.h"

Iterator::Iterator() : I(0) {}
Iterator::Iterator(int i) : I(i) {}

int Iterator::operator*() const {
  return I;
}

Iterator& Iterator::operator++() {
  ++I;
  return *this;
}

bool Iterator::operator!=(const Iterator& a) const {
  return I != a.I;
}
  
Range::Range(int end)
  : Begin(0)
  , End(end)
{ }

Range::Range(int begin, int end)
  : Begin(begin)
  , End(end)
{}

Iterator Range::begin() const {
  return Begin;
}

Iterator Range::end() const {
  return End;
}

UnitIterator::UnitIterator(
  VEC<Cell>& cells,
  MAT<double>& cellsDeathRate,
  MAT<int>& cellsPopulation,
  int i,
  bool isEnd
)
  : Cells(cells)
  , CellsDeathRate(cellsDeathRate)
  , CellsPopulation(cellsPopulation)
  , I(i)
  , J(0)
  , IsEnd(isEnd)
{
  if (isEnd) {
    return;
  }
  while (Cells[I].coords_x.empty() && I != Cells.size()) {
    ++I;
  }
}

Unit UnitIterator::operator*() {
  return Unit(Cells[I], CellsDeathRate, CellsPopulation, I, J);
}

void UnitIterator::operator++() {
  ++J;
  if (J == Cells[I].coords_x.size()) {
    ++I;
    while (Cells[I].coords_x.empty() && I != Cells.size()) {
      ++I;
    }
    J = 0;
  }
}

bool UnitIterator::operator==(const UnitIterator& b) const {
  if (IsEnd == b.IsEnd) {
    if (IsEnd) {
      return true;
    } else {
      return I == b.I && J == b.J;
    }
  } else {
    if (!IsEnd) {
      return I >= b.I;
    } else {
      return I <= b.I;
    }
  }
}

bool UnitIterator::operator!=(const UnitIterator& b) const {
  return !(*this == b);
}

UnitIterating::UnitIterating(
  VEC<Cell>& cells,
  MAT<double>& cellsDeathRate,
  MAT<int>& cellsPopulation,
  int minCell,
  int maxCell
)
  : Cells(cells)
, CellsDeathRate(cellsDeathRate)
, CellsPopulation(cellsPopulation)
, MinCell(minCell)
, MaxCell(maxCell)
{}

UnitIterator UnitIterating::begin() {
  return UnitIterator(
    Cells,
    CellsDeathRate,
    CellsPopulation,
    MinCell,
    false
  );
}

UnitIterator UnitIterating::end() {
  return UnitIterator(
    Cells,
    CellsDeathRate,
    CellsPopulation,
    MaxCell,
    true
  );
}
