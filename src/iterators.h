#ifndef ITERATORS
#define ITERATORS

#include "defines.h"
#include "cell.h"
#include "unit.h"

class Iterator {
  int I;
  
public:
  Iterator();
  Iterator(int i);
  
  int operator*() const;
  Iterator& operator++();
  bool operator!=(const Iterator& a) const;
};

struct Range {
  const int Begin;
  const int End;
  
  Range(int end);
  Range(int begin, int end);
  
  Iterator begin() const;
  Iterator end() const;
};

class UnitIterator {
  VEC<Cell>& Cells;
  MAT<double>& CellsDeathRate;
  MAT<int>& CellsPopulation;
  int I, J;
  bool IsEnd;
  
public:
  UnitIterator(
    VEC<Cell>& cells,
    MAT<double>& cellsDeathRate,
    MAT<int>& cellsPopulation,
    int i,
    bool isEnd
  );
  
  Unit operator*();
  void operator++();
  bool operator==(const UnitIterator& b) const;
  bool operator!=(const UnitIterator& b) const;
};

class UnitIterating {
  VEC<Cell>& Cells;
  MAT<double>& CellsDeathRate;
  MAT<int>& CellsPopulation;
  const int MinCell, MaxCell;
  
public:
  UnitIterating(
    VEC<Cell>& cells,
    MAT<double>& cellsDeathRate,
    MAT<int>& cellsPopulation,
    int minCell,
    int maxCell
  );
  
  UnitIterator begin();
  UnitIterator end();
};

#endif
