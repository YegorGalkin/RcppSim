#pragma once

#include <vector>

#include "defines.h"
#include "chunk.h"
#include "model_parameters.h"

template <size_t dim>
class Grid {
    std::vector<Chunk<dim>> chunks;
    std::vector<double> chunkDeathRate;
    std::vector<size_t> chunkPopulation;
    
    Position<dim> cellCounts;

public:
    const size_t speciesCount;
    const bool isPeriodic;
    const Coord<dim> areaLength;
    const ModelParameters modelParameters;

private:
    size_t GetOffset(const Position<dim>& pos) const;
    size_t GetOffset(const Position<dim>& pos, size_t species) const;
    
public:
    Chunk<dim>& GetChunk(const Position<dim>& chunkPos);
    double& GetChunkDeathRate(const Position<dim>& chunkPos, size_t species);
    size_t& GetChunkPopulation(const Position<dim>& chunkPos, size_t species);
    size_t GetChunkPopulation(const Position<dim>& chunkPos) const;
    size_t GetCellCount(size_t i) const;
};
