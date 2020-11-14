#pragma once

#include <array>

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
