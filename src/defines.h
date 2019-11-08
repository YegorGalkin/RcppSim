#include <vector>
#include <cassert>

#ifndef DEFINES
#define DEFINES

class Unit;
struct Cell;
class UnitIterating;
class Range;
class UnitIteration;

template <class T>
using VEC = std::vector<T>;

template <class T>
using MAT = VEC<VEC<T>>;

template<class T>
MAT<T> MakeMat(int i) {
  return MAT<T>(i, VEC<T>(i));
}

using DCoord = double;

#endif
