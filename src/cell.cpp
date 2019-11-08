#include "cell.h"

void Cell::SwapWithLast(int i) {
  death_rates[i] = death_rates.back();
  coords_x[i] = coords_x.back();
  species[i] = species.back();
}

void Cell::Pop() {
  death_rates.pop_back();
  coords_x.pop_back();
  species.pop_back();
}

void Cell::Add(double deathRate, DCoord x, int species) {
  death_rates.emplace_back(deathRate);
  coords_x.emplace_back(std::move(x));
  this->species.emplace_back(species);
}

DCoord Cell::Ro(const Unit& a, const Unit& b) {
  return std::abs(a.Coord() - b.Coord());
}
