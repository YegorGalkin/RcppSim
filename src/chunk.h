#pragma once

#include <vector>

#include "defines.h"

template <size_t dim>
struct Chunk {
    std::vector<Coord<dim>> coords;
    std::vector<double> deathRate;
    std::vector<size_t> species;
};
