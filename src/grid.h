#pragma once

#include <vector>

#include <boost/random/lagged_fibonacci.hpp>
#include <boost/math/interpolators/cubic_b_spline.hpp>

#include "defines.h"
#include "area.h"
#include "chunk.h"
#include "model_parameters.h"
#include "unit.h"
#include "iterator.h"

template <size_t dim>
class Grid {
    std::vector<Chunk<dim>> Chunks;
    std::vector<std::vector<double>> ChunkDeathRate;
    std::vector<std::vector<size_t>> ChunkPopulation;
    
    Position<dim> CellCounts;
    Position<dim> LocalRadius;
    
    std::vector<double> TotalDeathRate;
    std::vector<size_t> TotalPopulation;
    
    boost::random::lagged_fibonacci2281 Rnd;
    
    std::vector<boost::math::cubic_b_spline<double>> BirthReverseCdfSpline;
    
    size_t EventCount;

public:
    const Area<dim> Area;
    const ModelParameters ModelParameters;

private:
    size_t GetOffset(const Position<dim>& pos) const;
    size_t GetOffset(const Position<dim>& pos, size_t species) const;
    
    template<class T>
    Position<dim> GetRandomDistribution(const std::vector<T>& dist);
    
    Position<dim> GetPositionByOffset(size_t offset) const;
    
    void AddInteraction(Unit<dim>& a, double interaction);
    void AddInteraction(Unit<dim>& a, Unit<dim>& b, bool isSub);
    
    void AddDeathRate(Unit<dim>& a);
    void SubDeathRate(Unit<dim>& a);
    
    Range<UnitIterator<dim>> GetLocalUnits(const Unit<dim>& unit);
    
public:
    bool AddUnit(Coord<dim> coord, size_t species); // return if unit added
    void RemoveUnit(Unit<dim>& unit);
    
    void KillRandom(size_t species);
    void SpawnRandom(size_t species);
    
    void MakeEvent();

    Chunk<dim>& GetChunk(const Position<dim>& chunkPos);
    double& GetChunkDeathRate(const Position<dim>& chunkPos, size_t species);
    size_t& GetChunkPopulation(const Position<dim>& chunkPos, size_t species);
    size_t GetChunkPopulation(const Position<dim>& chunkPos) const;
    size_t GetCellCount(size_t i) const;
    size_t GetAllPopulation() const;
};
