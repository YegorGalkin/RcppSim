#include "grid.h"

#include "calculations.h"

template<>
size_t Grid<1>::GetOffset(const Position<1>& pos) const {
    return pos[0];
}

template<>
size_t Grid<2>::GetOffset(const Position<2>& pos) const {
    return pos[0] * CellCounts[1] + pos[1];
}

template<>
size_t Grid<3>::GetOffset(const Position<3>& pos) const {
    return (pos[0] * CellCounts[1] + pos[1]) * CellCounts[2] + pos[2];
}

template <size_t dim>
size_t Grid<dim>::GetOffset(const Position<dim>& pos, size_t species) const {
    return GetOffset(pos) * ModelParameters.SpeciesCount + species;
}

template <size_t dim>
Chunk<dim>& Grid<dim>::GetChunk(const Position<dim>& chunkPos) {
    return Chunks[GetOffset(chunkPos)];
}

template <size_t dim>
double& Grid<dim>::GetChunkDeathRate(const Position<dim>& chunkPos, size_t species) {
    return ChunkDeathRate[GetOffset(chunkPos, species)];
}

template <size_t dim>
size_t& Grid<dim>::GetChunkPopulation(const Position<dim>& chunkPos, size_t species) {
    return ChunkPopulation[GetOffset(chunkPos, species)];
}

template <size_t dim>
size_t Grid<dim>::GetChunkPopulation(const Position<dim>& chunkPos) const {
    return Chunks[GetOffset(chunkPos)].coords.size();
}

template <size_t dim>
size_t Grid<dim>::GetCellCount(size_t i) const {
    return CellCounts[i];
}

template <size_t dim>
void Grid<dim>::AddInteraction(Unit<dim>& unit, double interaction) {
    unit.DeathRate() += interaction;
    unit.ChunkDeathRate() += interaction;
    TotalDeathRate[ModelParameters.SpeciesCount] += interaction;
}

template <size_t dim>
void Grid<dim>::AddInteraction(Unit<dim>& a, Unit<dim>& b) {
    auto interaction = ModelParameters.GetInteraction(
        a.Species(),
        b.Species(),
        Ro(*this, a, b)
    );
    if (interaction < 0) {
        return;
    }
    AddInteraction(a, interaction);
    AddInteraction(b, interaction);
}

template <size_t dim>
void Grid<dim>::AddDeathRate(Unit<dim>& a) {
    AddInteraction(a, ModelParameters.GetD(a.Species()));
}

template <size_t dim>
void Grid<dim>::SubDeathRate(Unit<dim>& a) {
    AddInteraction(a, -ModelParameters.GetD(a.Species()));
}

template class Grid<1>;
template class Grid<2>;
template class Grid<3>;
