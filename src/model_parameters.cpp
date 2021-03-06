#include "model_parameters.h"

#include <boost/math/tools/roots.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>
#include "defines.h"

using boost::irange;

using Spline = boost::math::cubic_b_spline<double>;

constexpr auto SPECIES_FIELD = "species_count";

Spline GetSpline(const std::vector<double>& func, double cutoff) {
    auto step = cutoff / (func.size() - 1);
    
    return Spline{
        func.begin(),
        func.end(),
        0,
        step,
        0,
        0
    };
}

double GetTrapezoidal(const std::function<double(double)>& func, double x) {
    return boost::math::quadrature::trapezoidal(
        func,
        0.0,
        x
    );
}

ModelParameters::ModelParameters(const Rcpp::List& params) : SpeciesCount(GetParameter<size_t>(params, SPECIES_FIELD, 1)) {
    B.resize(SpeciesCount);
    D.resize(SpeciesCount);
    DD.resize(SpeciesCount * SpeciesCount);

    BirthReverseKernel.resize(SpeciesCount);
    DeathCutoff.resize(SpeciesCount * SpeciesCount);
    DeathKernel.resize(SpeciesCount * SpeciesCount);

    for (auto i : IterSpecies()) {
        B[i] = Rcpp::as<double>(params[GetName("b", i)]);
        D[i] = Rcpp::as<double>(params[GetName("d", i)]);
        
        for (auto j : IterSpecies()) {
            auto offset = GetOffset(i, j);
            DD[offset] = Rcpp::as<double>(params[GetName("dd", i, j)]);
            
            auto deathKernelY = Rcpp::as<std::vector<double>>(params[GetName("death_kernel_y", i, j)]);
            auto deathKernelCutoff = Rcpp::as<double>(params[GetName("death_kernel_r", i, j)]);
            DeathKernel[offset] = GetSpline(deathKernelY, deathKernelCutoff);
            DeathCutoff[offset] = deathKernelCutoff;
        }
        
        BirthReverseKernel[i] = GetSpline(Rcpp::as<std::vector<double>>(params[GetName("birth_kernel_y", i)]), 1);
    }
}

boost::integer_range<size_t> ModelParameters::IterSpecies() const {
    return boost::irange(SpeciesCount);
}

inline size_t ModelParameters::GetOffset(size_t i, size_t j) const {
    if (i >= SpeciesCount || j >= SpeciesCount) {
        throw std::runtime_error(std::string(__PRETTY_FUNCTION__) + ":" + std::to_string(__LINE__) + "Invalid index: (" + std::to_string(i) + " || " + std::to_string(j) + ") >= " + std::to_string(SpeciesCount));
    }
    return SpeciesCount * i + j;
}

inline double ModelParameters::GetCutoff(size_t a, size_t b) const {
    return DeathCutoff[GetOffset(a, b)];
}

inline double ModelParameters::GetDeathKernel(size_t a, size_t b, double distance) const {
    return DeathKernel[GetOffset(a, b)](distance);
}

double ModelParameters::GetDD(size_t a, size_t b) const {
    return DD[GetOffset(a, b)];
}

double ModelParameters::GetInteraction(size_t a, size_t b, double distance) const {
    if (distance > GetCutoff(a, b)) {
        return -1;
    }
    return GetDD(a, b) * GetDeathKernel(a, b, distance);
}

double ModelParameters::GetD(size_t species) const {
    return D[species];
}

const boost::math::cubic_b_spline<double>& ModelParameters::GetBirthSpline(size_t species) const {
    return BirthReverseKernel[species];
}

std::string ModelParameters::GetName(const std::string& name, size_t i) const {
    if (SpeciesCount == 1) {
        return name;
    }
    return name + "_" + std::to_string(i + 1);
}

std::string ModelParameters::GetName(const std::string& name, size_t i, size_t j) const {
    if (SpeciesCount == 1) {
        return name;
    }
    return name + "_" + std::to_string(i + 1) + "_" + std::to_string(j + 1);
}

double ModelParameters::GetMaximumCutoff() const {
    double cutoff = 0;
    for (auto i : IterSpecies()) {
        for (auto j : IterSpecies()) {
            cutoff = std::max(cutoff, GetCutoff(i, j));
        }
    }
    return cutoff;
}
