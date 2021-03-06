#include "chunk.h"

#include "unit.h"

template <size_t dim>
Unit<dim> Chunk<dim>::AddUnit(Grid<dim>& grid, Position<dim> chunkPosition, Coord<dim>& coord, size_t species) {
    const auto i = Coords.size();
    
    Coords.push_back(coord);
    DeathRate.push_back(0.0);
    Species.push_back(species);   
    return Unit<dim>(
        grid,
        chunkPosition,
        i
    );
}

template <size_t dim>
size_t Chunk<dim>::GetPopulation() const {
    return Coords.size();
}

template <size_t dim>
size_t Chunk<dim>::GetSpecies(size_t i) const {
    return Species[i];
}

template <size_t dim>
double& Chunk<dim>::GetDeathRate(size_t i) {
    return DeathRate[i];
}

template <size_t dim>
const Coord<dim>& Chunk<dim>::GetCoord(size_t i) const {
    return Coords[i];
}

template <size_t dim>
void Chunk<dim>::RemoveUnit(size_t i) {
    if (i >= GetPopulation()) {
        throw std::runtime_error("Invalid Deletion unit i: " + std::to_string(i) + " : Population: " + std::to_string(GetPopulation()));
    }
    auto lastIndex = Coords.size() - 1;
    if (i != lastIndex) {
        std::swap(Coords[i], Coords[lastIndex]);
        std::swap(DeathRate[i], DeathRate[lastIndex]);
        std::swap(Species[i], Species[lastIndex]);
    }
    Coords.pop_back();
    DeathRate.pop_back();
    Species.pop_back();
}

template <size_t dim>
const std::vector<double>& Chunk<dim>::GetDeathRates() const {
    return DeathRate;
}

template class Chunk<1>;
template class Chunk<2>;
template class Chunk<3>;
