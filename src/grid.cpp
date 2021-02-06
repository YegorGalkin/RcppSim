#include "grid.h"

#include <algorithm>
#include <cmath>

#include <boost/random.hpp>

using boost::irange;

template<>
size_t Grid<1>::GetOffset(const Position<1>& pos) const {
    return pos[0];
}

template<>
Position<1> Grid<1>::GetPositionByOffset(size_t offset) const {
    return {
        offset
    };
}

template<>
size_t Grid<2>::GetOffset(const Position<2>& pos) const {
    return pos[0] * CellCounts[1] + pos[1];
}

template<>
Position<2> Grid<2>::GetPositionByOffset(size_t offset) const {
    return {
        offset / CellCounts[1],
        offset % CellCounts[1]
    };
}

template<>
size_t Grid<3>::GetOffset(const Position<3>& pos) const {
    return (pos[0] * CellCounts[1] + pos[1]) * CellCounts[2] + pos[2];
}

template<>
Position<3> Grid<3>::GetPositionByOffset(size_t offset) const {
    return {
        (offset / CellCounts[2]) / CellCounts[1],
        (offset / CellCounts[2]) % CellCounts[1],
        offset % CellCounts[2]
    };
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
    return ChunkDeathRate[species][GetOffset(chunkPos)];
}

template <size_t dim>
size_t& Grid<dim>::GetChunkPopulation(const Position<dim>& chunkPos, size_t species) {
    return ChunkPopulation[species][GetOffset(chunkPos)];
}

template <size_t dim>
size_t Grid<dim>::GetChunkPopulation(const Position<dim>& chunkPos) const {
    return Chunks[GetOffset(chunkPos)].GetPopulation();
}

template <size_t dim>
size_t Grid<dim>::GetCellCount(size_t i) const {
    return CellCounts[i];
}

void FixDoubleZero(double& x, double epsilon = 1e-10) {
    if (std::abs(x) < epsilon) {
        x = 0;
    }
}

template <size_t dim>
void Grid<dim>::AddInteraction(Unit<dim>& unit, double interaction) {
    unit.DeathRate() += interaction;
    FixDoubleZero(unit.DeathRate());
    
    unit.ChunkDeathRate() += interaction;
    FixDoubleZero(unit.ChunkDeathRate());
    
    TotalDeathRate[unit.Species()] += interaction;
    FixDoubleZero(TotalDeathRate[unit.Species()]);
}

template <size_t dim>
void Grid<dim>::AddInteraction(Unit<dim>& a, Unit<dim>& b, bool isSub) {
    auto ro = Area.Ro(a, b);
    auto aS = a.Species();
    auto bS = b.Species();
    if (auto interaction = ModelParameters.GetInteraction(aS, bS, ro) >= 0) {
        interaction = isSub ? -interaction : interaction;
        AddInteraction(b, interaction);
    }
    if (auto interaction = ModelParameters.GetInteraction(bS, aS, ro) >= 0) {
        interaction = isSub ? -interaction : interaction;
        AddInteraction(a, interaction);
    }
}

template <size_t dim>
void Grid<dim>::AddDeathRate(Unit<dim>& a) {
    AddInteraction(a, ModelParameters.GetD(a.Species()));
}

template <size_t dim>
void Grid<dim>::SubDeathRate(Unit<dim>& a) {
    AddInteraction(a, -ModelParameters.GetD(a.Species()));
}

template <size_t dim>
bool Grid<dim>::AddUnit(Coord<dim> coord, size_t species) {
    if (!Area.IsInArea(coord)) {
        return false;
    }
    Area.FixCoord(coord);
    
    auto pos = Area.GetCellIndex(CellCounts, coord);
    auto& chunk = GetChunk(pos);
    auto& newUnit = chunk.AddUnit(*this, pos, coord, species);
    
    ++TotalPopulation[species];
    AddDeathRate(newUnit);
    
    for (auto unit : GetLocalUnits(newUnit)) {
        if (newUnit != unit) {
            AddInteraction(newUnit, unit, false);
        }
    }
    return true;
}

template <size_t dim>
void Grid<dim>::RemoveUnit(Unit<dim>& unit) {
    --TotalPopulation[unit.Species()];
    SubDeathRate(unit);
    
    for (auto other : GetLocalUnits(unit)) {
        if (unit != other) {
            AddInteraction(unit, other, true);
        }
    }
    
    unit.RemoveUnit();
}

template <size_t dim>
Range<UnitIterator<dim>> Grid<dim>::GetLocalUnits(const Unit<dim>& unit) {
    Position<dim> startPos, endPos;
    const auto& unitCoord = unit.Coord();
    for (auto i : irange(dim)) {
        startPos[i] = unitCoord[i] - LocalRadius[i];
        endPos[i] = unitCoord[i] + LocalRadius[i];
    }
    return {
        UnitIterator<dim>(
            *this,
            startPos,
            endPos,
            Area.GetIsPeriodic()
        ),
        UnitIterator<dim>(
            *this,
            /* isEnd = */ true
        )
    };
}

template <size_t dim>
template <class T>
Position<dim> Grid<dim>::GetRandomDistribution(const std::vector<T>& dist) {
    size_t index = boost::random::discrete_distribution<>(dist)(Rnd);
    return GetPositionByOffset(index);
}


template <size_t dim>
void Grid<dim>::KillRandom(size_t species) {
    auto chunkPosition = GetRandomDistribution(ChunkDeathRate[species]);
    auto& chunk = GetChunk(chunkPosition);
    std::vector<double> chunkDeathRateTmp = chunk.GetDeathRates();
    for (size_t i = 0; i < chunk.GetPopulation(); ++i) {
        if (chunk.GetSpecies(i) != species) {
            chunkDeathRateTmp[i] = 0;
        }
    }
    size_t unitIndex = boost::random::discrete_distribution<>(chunkDeathRateTmp)(Rnd);
    auto unit = Unit<dim>(*this, chunkPosition, unitIndex);
    
    RemoveUnit(unit);
}

template <size_t dim>
void Grid<dim>::SpawnRandom(size_t species) {
    auto chunkPosition = GetRandomDistribution(ChunkPopulation[species]);
    auto chunk = GetChunk(chunkPosition);
    
    size_t eventIndex = boost::random::uniform_smallint<>(
        1,
        chunk.GetPopulation()
    )(Rnd);
    for (size_t i = 0; ; ++i) {
        if (chunk.GetSpecies(i) == species) {
            --eventIndex;
        }
        if (eventIndex == 0) {
            eventIndex = i;
            break;
        }
    }
    
    auto parentUnit = Unit<dim>(*this, chunkPosition, eventIndex);
    const auto oldCoord = parentUnit.Coord();
    Coord<dim> newCoord;
    for (size_t i = 0; i < dim; ++i) {
        newCoord[i] = oldCoord[i] + BirthReverseCdfSpline[species](
            boost::random::uniform_01<>()(Rnd)
        ) * (
            boost::random::bernoulli_distribution<>(0.5)(Rnd) * 2 - 1
        );
    }
    AddUnit(newCoord, species);
}

template <size_t dim>
void Grid<dim>::MakeEvent() {
    if (GetAllPopulation()) {
        return;
    }
    
    ++EventCount;
    
    //Rolling event according to global birth \ death rate
    std::vector<double> dis(ModelParameters.SpeciesCount * 2, 0);
    for (auto s : ModelParameters.IterSpecies()) {
        if (TotalPopulation[s] > 0) {
            dis[2 * s + 0] = TotalDeathRate[s];
            dis[2 * s + 1] = TotalPopulation[s] * ModelParameters.GetD(s);
        }
    }
    auto t = boost::random::discrete_distribution<>(dis)(Rnd);
    size_t event = t % 2;
    size_t species = t / 2;
    if (event == 0) {
        KillRandom(species);
    } else {
        SpawnRandom(species);
    }
}

template <size_t dim>
size_t Grid<dim>::GetAllPopulation() const {
    size_t res = 0;
    for (auto x : TotalPopulation) {
        res += x;
    }
    return res;
}

template class Grid<1>;
template class Grid<2>;
template class Grid<3>;
