#pragma once

#include "defines.h"

#include "unit.h"
#include "grid.h"

template <size_t dim>
class UnitIterator {
    Grid<dim>& grid;
    Position<dim> pos;
    Position<dim> end;
    size_t i;
    bool isPeriodic;
    bool isEnd;

private:
    void FixPositions();
    Position<dim> GetPosition() const;
    void Next();

public:
    UnitIterator<dim>(
        Grid<dim>& grid,
        bool isEnd = false
    )
        : grid(grid)
        , i(0)
        , isPeriodic(false)
        , isEnd(isEnd)
    {
        for (size_t j = 0; j < dim; ++j) {
            pos[j] = 0;
            end[j] = grid.GetCellCount(j);
        }
        if (!IsEnd() && grid.GetChunkPopulation(pos) == 0) {
            this->operator++();
        }
    }
    UnitIterator<dim>(
        Grid<dim>& grid,
        Position<dim> pos,
        Position<dim> end,
        bool isPeriodic,
        bool isEnd = false
    )
        : grid(grid)
        , pos(pos)
        , end(end)
        , i(0)
        , isPeriodic(isPeriodic)
        , isEnd(isEnd)
    {
        FixPositions();
        if (!IsEnd() && grid.GetChunkPopulation(GetPosition()) == 0) {
            this->operator++();
        }
    }
    
    bool IsEnd() const;
    
    Unit<dim> operator*() const;
    
    void operator++();
    bool operator!=(const UnitIterator<dim>& other) const;
    bool operator==(const UnitIterator<dim>& other) const;
};

// TODO figure out why UnitIterator not work with boost::iterator_range
template <class It>
class Range {
    It Begin;
    It End;
    
public:
    Range(It begin, It end)
        : Begin(begin)
        , End(end)
    {}
    
    It begin() const {
        return Begin;
    }
    
    It end() const {
        return End;
    } 
};
