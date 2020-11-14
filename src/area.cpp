#include <algorithm>
#include <cmath>

#include "area.h"

template<> std::string Area<1>::GetName(const std::string& name, size_t i) {
    return name;
}

const std::string AREA_RENAME[] = {"x", "y", "z"};

template <size_t dim>
std::string Area<dim>::GetName(const std::string& name, size_t i) {
    return name + "_" + AREA_RENAME[i];
}

template <size_t dim>
double Area<dim>::Ro(const Unit<dim> &a, const Unit<dim> &b) const {
    if (IsPeriodic) {
        double calc = 0;
        for (size_t i = 0; i < dim; ++i) {
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

template class Area<1>;
template class Area<2>;
template class Area<3>;
