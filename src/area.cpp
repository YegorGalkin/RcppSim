#include <algorithm>
#include <cmath>

#include <boost/range/irange.hpp>

#include "area.h"

using boost::irange;

const std::string AREA_RENAME[] = {"x", "y", "z"};

std::string GetAreaName(const std::string& name, size_t i) {
    return name + "_" + AREA_RENAME[i];
}

template <size_t dim>
double Area<dim>::Ro(const Unit<dim> &a, const Unit<dim> &b) const {
    if (IsPeriodic) {
        double calc = 0;
        for (auto i : irange(dim)) {
            auto dist = std::abs(a.Coord()[i] - b.Coord()[i]);
            dist = std::min<double>(
                dist,
                AreaLength[i] - dist
            );
            calc += dist * dist;
        }
        return std::sqrt(calc);
        
    }
    return Ro(a.Coord(), b.Coord());
}

template <size_t dim>
bool Area<dim>::IsInArea(const Coord<dim>& coord) const {
    if (IsPeriodic) {
        return true;
    }
    for (auto i : irange(dim)) {
        if (coord[i] < 0 || coord[i] > AreaLength[i]) {
            return false;
        }
    }
    return true;
}

template <size_t dim>
void Area<dim>::FixCoord(Coord<dim>& coord) const {
    for (auto i : irange(dim)) {
        const auto len = AreaLength[i];
        auto& x = coord[i];
        
        while (x < 0) {
            x += len;
        }
        while (x > len) {
            x -= len;
        }
    }
}

template <size_t dim>
bool Area<dim>::GetIsPeriodic() const {
    return IsPeriodic;
}

template <size_t dim>
Position<dim> Area<dim>::GetCellIndex(const Position<dim>& cellCounts, const Coord<dim>& coord) const {
    Position<dim> res;
    for (auto i : irange(dim)) {
        const auto x = coord[i];
        const auto len = AreaLength[i];
        const auto counts = cellCounts[i];
        auto& pos = res[i];
        
        pos = std::floor(x * counts / len);
        if (pos >= counts) {
            pos = counts - 1;
        }
        if (pos <= 0) {
            pos = 0;
        }
    }
    return res;
}

template<>
double Area<1>::Ro(const Coord<1>& a, const Coord<1>& b) {
    return std::abs(a[0] - b[0]);
}

template<>
double Area<2>::Ro(const Coord<2>& a, const Coord<2>& b) {
    return std::hypot(a[0] - b[0], a[1] - b[1]);
}

template<>
double Area<3>::Ro(const Coord<3>& a, const Coord<3>& b) {
    return std::sqrt(
        std::pow(a[0] - b[0], 2) +
        std::pow(a[1] - b[1], 2) +
        std::pow(a[2] - b[2], 2)
    );
}

template<size_t dim>
const Coord<dim>& Area<dim>::GetAreaLength() const {
    return AreaLength;
}

template class Area<1>;
template class Area<2>;
template class Area<3>;
