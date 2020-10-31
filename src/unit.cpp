#include "unit.h"

template <size_t dim>
const Coord<dim>& Unit<dim>::Coord() const {
    return chunk.coords[i];
}

template <size_t dim>
double& Unit<dim>::DeathRate() {
    return chunk.deathRate[i];
}

template <size_t dim>
size_t Unit<dim>::Species() const {
    return chunk.species[i];
}

template <size_t dim>
double& Unit<dim>::ChunkDeathRate() {
    return chunkDeathRate;
}

template <size_t dim>
size_t& Unit<dim>::ChunkPopulation() {
    return chunkPopulation;
}

template <size_t dim>
const Position<dim>& Unit<dim>::ChunkPosition() const {
    return chunkPos;
}

template class Unit<1>;
template class Unit<2>;
template class Unit<3>;
