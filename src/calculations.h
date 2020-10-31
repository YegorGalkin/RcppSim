#pragma once

#include "defines.h"

#include "unit.h"
#include "grid.h"

template <size_t dim>
double Ro(const Grid<dim> &grid, const Unit<dim> &a, const Unit<dim> &b);

template<size_t dim>
double Ro(const Coord<dim>& a, const Coord<dim>& b);
