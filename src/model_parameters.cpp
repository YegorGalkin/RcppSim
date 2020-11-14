#include "model_parameters.h"

#include <boost/math/tools/roots.hpp>
#include <boost/math/quadrature/trapezoidal.hpp>

using boost::irange;

using Spline = boost::math::cubic_b_spline<double>;

constexpr auto SPECIES_FIELD = "species_count";

Spline GetSpline(std::vector<double> func, double cutoff) {
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

double GetTrapezoidal(std::function<double(double)> func, double x) {
    return boost::math::quadrature::trapezoidal(
        func,
        0.0,
        x
    );
}

Spline GetReverseSpline(std::vector<double> func, double cutoff) {
    auto nodes = func.size();
    auto step = cutoff / (nodes - 1);
    auto spline = Spline{
        func.begin(),
        func.end(),
        0,
        step,
        0,
        0
    };
    
    double approxConst = GetTrapezoidal(spline, cutoff);
    
    std::vector<double> quantile(nodes);
    
    for (auto i : irange(nodes)) {
        quantile[i] = boost::math::tools::newton_raphson_iterate(
            [&](double y) {
                return std::make_tuple(
                    GetTrapezoidal(spline, y) / approxConst - (double)i / (nodes - 1),
                    spline(y) / approxConst
                );
            },
            1e-10,
            0.0,
            cutoff,
            std::numeric_limits<double>::digits
        );
    }
    
    auto tempSize = std::find_if_not(
        quantile.begin(),
        quantile.end(),
        [](double x){return x < 1e-300;}
    ) - quantile.begin();
    
    std::vector<double> quantileTemp(quantile.begin() + tempSize, quantile.end());
    
    auto reverseStep = 1.0 / (quantileTemp.size() - 1);
    
    auto rightDeretive = (*(quantileTemp.rbegin()+ 2) - *(quantileTemp.rbegin() + 3)) / reverseStep;
    
    quantileTemp.back() = *(quantileTemp.rbegin() + 2) + rightDeretive * reverseStep;
    
    return Spline(
        quantileTemp.begin(),
        quantileTemp.end(),
        0,
        reverseStep,
        0.5 / spline(0),
        2.0 * rightDeretive
    );
}

ModelParameters::ModelParameters(const Rcpp::List& params) : SpeciesCount(1) {
    if (params.hasSlot(SPECIES_FIELD)) {
        const_cast<size_t&>(SpeciesCount) =  Rcpp::as<size_t>(params[SPECIES_FIELD]);
    }
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
            
            auto deathKernelY = Rcpp::as<std::vector<double>>(params[GetName("death_kernel_y", i)]);
            auto deathKernelCutoff = Rcpp::as<double>(params[GetName("death_kernel_r", i)]);
            DeathKernel[offset] = GetSpline(deathKernelY, deathKernelCutoff);
            DeathCutoff[offset] = deathKernelCutoff;
        }
        
        auto birthKernelY = Rcpp::as<std::vector<double>>(params[GetName("birth_kernel_y", i)]);
        auto birthKernelCutoff = Rcpp::as<double>(params[GetName("birth_kernel_r", i)]);
        BirthReverseKernel[i] = GetReverseSpline(birthKernelY, birthKernelCutoff);
    }
}

boost::integer_range<size_t> ModelParameters::IterSpecies() const {
    return boost::irange(SpeciesCount);
}

inline size_t ModelParameters::GetOffset(size_t i, size_t j) const {
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
