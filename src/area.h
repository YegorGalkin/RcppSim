#pragma once

#include <Rcpp.h>

#include "defines.h"

#include "unit.h"
#include "grid.h"

template <size_t dim>
class Area {
    Coord<dim> AreaLength;
    bool IsPeriodic;
    
public:
    Area(const Rcpp::List& params)
        : IsPeriodic(Rcpp::as<bool>(params["periodic"]))
    {
        for (size_t i = 0; i < dim; ++i) {
            AreaLength[i] = Rcpp::as<double>(params[GetName("area_length", i)]);
        }
    }

public:
    static double Ro(const Coord<dim>& a, const Coord<dim>& b);
    double Ro(const Unit<dim> &a, const Unit<dim> &b) const;
    
private:
    std::string GetName(const std::string& name, size_t i);
};
