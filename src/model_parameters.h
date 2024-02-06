#pragma once

#include <Rcpp.h>

#include <boost/math/interpolators/cubic_b_spline.hpp>
#include <boost/range/irange.hpp>
#include <vector>

class ModelParameters {
public:
    const size_t SpeciesCount;
    
private:
    std::vector<double> B;
    std::vector<double> D;
    std::vector<double> DD;
    
    std::vector<double> DeathCutoff;
    
    std::vector<boost::math::cubic_b_spline<double>> DeathKernel;
    std::vector<boost::math::cubic_b_spline<double>> BirthReverseKernel;

public:
    ModelParameters(const Rcpp::List& params);
    double GetDD(size_t a, size_t b) const;
    double GetDeathKernel(size_t a, size_t b, double distance) const;
    boost::integer_range<size_t> IterSpecies() const;
    
    
private:
    std::string GetName(const std::string& name, size_t i) const;
    std::string GetName(const std::string& name, size_t i, size_t j) const;
    inline size_t GetOffset(size_t i, size_t j) const;
    double GetCutoff(size_t a, size_t b) const;
    
public:
    double GetInteraction(size_t speciesA, size_t speciesB, double distance) const;
    double GetD(size_t species) const;
    double GetB(size_t species) const;
    const boost::math::cubic_b_spline<double>& GetBirthSpline(size_t species) const;
    double GetMaximumCutoff() const;
}; 
