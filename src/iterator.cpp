#include "iterator.h"

#include <algorithm>
#include <cstdlib>

template <size_t dim>
bool UnitIterator<dim>::IsEnd() const {
    if (isEnd) {
        return true;
    }
    return pos.back() >= end.back();
}

template <size_t dim>
void UnitIterator<dim>::FixPositions() {
    for (size_t j = 0; j < dim; ++j) {
        const pos_t count = grid.GetCellCount(j);
        if (isPeriodic) {
            if (pos[j] < 0) {
                pos[j] = std::div(pos[j], count).rem;
            }
            if (end[j] >= count) {
                end[j] = count + std::div(end[j], count).rem;
            }
            if (end[j] - pos[j] > count) {
                pos[j] = 0;
                end[j] = count;
                isPeriodic = false;
            }
        } else {
            pos[j] = std::max<pos_t>(pos[j], 0);
            end[j] = std::min<pos_t>(end[j], count);
        }
    }
}

template <size_t dim>
Position<dim> UnitIterator<dim>::GetPosition() const {
    if (!isPeriodic) {
        return pos;
    }
    auto res = pos;
    for (size_t j = 0; j < dim; ++j) {
        const pos_t count = grid.GetCellCount(j);
        if (res[j] < 0) {
            res[j] += count;
        } else if (res[j] >= count) {
            res[j] -= count;
        }
    }
    return res;
}

template <size_t dim>
void UnitIterator<dim>::Next() {
    if (isEnd) {
        return;
    }
    ++i;
    if (i >= grid.GetChunkPopulation(GetPosition())) {
        i = 0;
        for (size_t j = 0; j != dim; ++j) {
            ++pos[j];
            if (pos[j] < end[j]) {
                break;
            }
            if (j != dim - 1) {
                pos[j] = 0;
            }
        }
    }
    if (IsEnd()) {
        isEnd = true;
    }
}

template <size_t dim>
void UnitIterator<dim>::operator++() {
    if (isEnd) {
        return;
    }
    do {
        Next();
    } while (!isEnd && grid.GetChunkPopulation(GetPosition()) == 0);
}

template <size_t dim>
bool UnitIterator<dim>::operator==(const UnitIterator<dim>& other) const {
    if (isEnd == other.isEnd) {
        if (isEnd) {
            return true;
        }
        if (i != other.i) {
            return false;
        }
        auto a = GetPosition();
        auto b = other.GetPosition();
        for (size_t j = 0; j < dim; ++j) {
            if (a[j] != b[j]) {
                return false;
            }
        }
        return true;
    }
    return false;
}

template <size_t dim>
bool UnitIterator<dim>::operator!=(const UnitIterator<dim>& other) const {
    return !((*this) == other);
}

template <size_t dim>
Unit<dim> UnitIterator<dim>::operator*() const {
    if (isEnd) {
        throw std::runtime_error("Invalid unpacking unitIterator");
    }
    return Unit<dim>(
        grid,
        pos,
        i
    );
}

template class UnitIterator<1>;
template class UnitIterator<2>;
template class UnitIterator<3>;
