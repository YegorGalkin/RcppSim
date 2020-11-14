#pragma once

#include <vector>

#include "defines.h"
#include "area.h"
#include "chunk.h"
#include "model_parameters.h"
#include "unit.h"

template <size_t dim>
class Grid {
    std::vector<Chunk<dim>> Chunks;
    std::vector<double> ChunkDeathRate;
    std::vector<size_t> ChunkPopulation;
    
    Position<dim> CellCounts;
    
    std::vector<double> TotalDeathRate;
    std::vector<size_t> TotalPopulation;

public:
    const Area<dim> Area;
    const ModelParameters ModelParameters;

private:
    size_t GetOffset(const Position<dim>& pos) const;
    size_t GetOffset(const Position<dim>& pos, size_t species) const;
    
    void AddInteraction(Unit<dim>& a, double interaction);
    void AddInteraction(Unit<dim>& a, Unit<dim>& b);
    
    void AddDeathRate(Unit<dim>& a);
    void SubDeathRate(Unit<dim>& a);
    
    void IncrementPopulation(Unit<dim>& a);
    void DecrementPopulation(Unit<dim>& a);
    
public:
    Chunk<dim>& GetChunk(const Position<dim>& chunkPos);
    double& GetChunkDeathRate(const Position<dim>& chunkPos, size_t species);
    size_t& GetChunkPopulation(const Position<dim>& chunkPos, size_t species);
    size_t GetChunkPopulation(const Position<dim>& chunkPos) const;
    size_t GetCellCount(size_t i) const;
};
