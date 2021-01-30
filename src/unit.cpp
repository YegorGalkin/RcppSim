#include "unit.h"

using boost::irange;

template <size_t dim>
const Coord<dim>& Unit<dim>::Coord() const {
    return chunk.GetCoord(i);
}

template <size_t dim>
double& Unit<dim>::DeathRate() {
    return chunk.GetDeathRate(i);
}

template <size_t dim>
size_t Unit<dim>::Species() const {
    return chunk.GetSpecies(i);
}

template <size_t dim>
double& Unit<dim>::ChunkDeathRate() {
    return chunkDeathRate;
}

template <size_t dim>
size_t Unit<dim>::ChunkPopulation() const {
    return chunk.GetPopulation();
}

template <size_t dim>
bool Unit<dim>::operator==(const Unit<dim>& other) const {
    if (i != other.i) {
        return false;
    }
    for (auto i : irange(dim)) {
        if (chunkPos[i] != other.chunkPos[i]) {
            return false;
        }
    }
    return true;
}

template <size_t dim>
bool Unit<dim>::operator!=(const Unit<dim>& other) const {
    return !((*this) == other);
}

template <size_t dim>
const Position<dim>& Unit<dim>::ChunkPosition() const {
    return chunkPos;
}

template <size_t dim>
void Unit<dim>::RemoveUnit() {
    chunk.RemoveUnit(i);
}

template class Unit<1>;
template class Unit<2>;
template class Unit<3>;
