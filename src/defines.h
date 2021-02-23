#pragma once

#include <Rcpp.h>

#include <array>
#include <iostream>

using pos_t = int32_t;

template<size_t dim>
struct Area;

template<size_t dim>
struct Chunk;

template <size_t dim>
struct Grid;

template <class It>
class Range;

template<size_t dim>
class Unit;

template <size_t dim>
class UnitIterator;

template <size_t dim>
using Coord = std::array<double, dim>;

template <size_t dim>
using Position = std::array<pos_t, dim>;

template<size_t dim>
std::string ToString(const Position<dim>& pos) {
    std::stringstream a;
    a << "(";
    for (size_t i = 0; i < dim; ++i) {
        a << pos[i];
        if (i != 0) {
            a << ", ";
        }
    }
    a << ")";
    return a.str();
}

template <typename T>
T GetParameter(const Rcpp::List& params, const std::string& name, T def = T()) {
    if (params.containsElementNamed(name.c_str())) {
        return Rcpp::as<T>(params[name]);
    }
    return def;
}

static Rcpp::Environment base("package:base");
static Rcpp::Function print = base["print"];
