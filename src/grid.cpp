#include "grid.h"

template<>
size_t Grid<1>::GetOffset(const Position<1>& pos) const {
    return pos[0];
}

template<>
size_t Grid<2>::GetOffset(const Position<2>& pos) const {
    return pos[0] * cellCounts[1] + pos[1];
}

template<>
size_t Grid<3>::GetOffset(const Position<3>& pos) const {
    return (pos[0] * cellCounts[1] + pos[1]) * cellCounts[2] + pos[2];
}

template <size_t dim>
size_t Grid<dim>::GetOffset(const Position<dim>& pos, size_t species) const {
    return GetOffset(pos) * speciesCount + species;
}

template <size_t dim>
Chunk<dim>& Grid<dim>::GetChunk(const Position<dim>& chunkPos) {
    return chunks[GetOffset(chunkPos)];
}

template <size_t dim>
double& Grid<dim>::GetChunkDeathRate(const Position<dim>& chunkPos, size_t species) {
    return chunkDeathRate[GetOffset(chunkPos, species)];
}

template <size_t dim>
size_t& Grid<dim>::GetChunkPopulation(const Position<dim>& chunkPos, size_t species) {
    return chunkPopulation[GetOffset(chunkPos, species)];
}

template <size_t dim>
size_t Grid<dim>::GetChunkPopulation(const Position<dim>& chunkPos) const {
    return chunks[GetOffset(chunkPos)].coords.size();
}

template <size_t dim>
size_t Grid<dim>::GetCellCount(size_t i) const {
    return cellCounts[i];
}

template <size_t dim>
void Grid<dim>::AddInteraction(Unit<dim>& unit, double interaction) {
    unit.DeathRate() += interaction;
    unit.ChunkDeathRate() += interaction;
    TotalDeathRate[modelParameters.SpeciesCount] += interaction;
}

template class Grid<1>;
template class Grid<2>;
template class Grid<3>;
