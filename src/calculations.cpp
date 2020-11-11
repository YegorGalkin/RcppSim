#include <algorithm>
#include <cmath>

#include "calculations.h"

template <size_t dim>
double Ro(const Grid<dim> &grid, const Unit<dim> &a, const Unit<dim> &b) {
    if (grid.IsPeriodic) {
        double calc = 0;
        for (size_t i = 0; i < dim; ++i) {
            auto dist = std::abs(a.Coord()[i] - b.Coord()[i]);
            dist = std::min(
                dist,
                grid.AreaLength[i] - dist
            );
            calc += dist * dist;
        }
        return std::sqrt(calc);
        
    }
    return Ro(a.Coord(), b.Coord());
}

template double Ro<1>(const Grid<1> &grid, const Unit<1> &a, const Unit<1> &b);
template double Ro<2>(const Grid<2> &grid, const Unit<2> &a, const Unit<2> &b);
template double Ro<3>(const Grid<3> &grid, const Unit<3> &a, const Unit<3> &b);

template<>
double Ro(const Coord<1>& a, const Coord<1>& b) {
    return std::abs(a[0] - b[0]);
}

template<>
double Ro(const Coord<2>& a, const Coord<2>& b) {
    return std::hypot(a[0] - b[0], a[1] - b[1]);
}

template<>
double Ro(const Coord<3>& a, const Coord<3>& b) {
    return std::sqrt(
        std::pow(a[0] - b[0], 2) +
        std::pow(a[1] - b[1], 2) +
        std::pow(a[2] - b[2], 2)
    );
}
