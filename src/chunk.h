#pragma once

#include <vector>

#include "defines.h"

template <size_t dim>
class Chunk {
    std::vector<Coord<dim>> Coords;
    std::vector<double> DeathRate;
    std::vector<size_t> Species;

public:
    Unit<dim>& AddUnit(Grid<dim>& grid, Position<dim> chunkPosition, Coord<dim>& coord, size_t species);
    size_t GetPopulation() const;
    double& GetDeathRate(size_t i);
    size_t GetSpecies(size_t i) const;
    const Coord<dim>& GetCoord(size_t i) const;
    void RemoveUnit(size_t i);
};
