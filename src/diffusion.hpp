// Diffusion coefficients and step sizes (calculate.py). Units: radius in nm, viscosity in Pa*s,
// temperature in K, time in ns.
#pragma once

#include <cmath>
#include <stdexcept>
#include <vector>

namespace diffusion {

// Written as in calculate.py (1.380649 * 1e-23, not 1.380649e-23) so results match to the last digit.
inline const double kBoltzmann = 1.380649 * 1e-23; // J/K
constexpr double kPi = 3.14159265358979323846;

// calc.calculate_diff_trans: Stokes-Einstein, D = kT / (6 pi eta r), returned in nm^2/ns.
inline double trans(double radius, double viscosity, double temperature) {
    if (viscosity == 0) throw std::invalid_argument("viscosity cannot be zero!");
    if (radius == 0) throw std::invalid_argument("radius cannot be zero!");
    const double d = kBoltzmann * temperature / (6 * kPi * viscosity * 1e-9 * radius);
    return d * 1e9; // m^2/s -> nm^2/ns
}

// calc.calculate_diff_rotation: Stokes-Einstein-Debye, D = kT / (8 pi eta r^3), in rad^2/ns.
// CORRECTED relative to the Python code: kT / (8 pi eta r^3) is in rad^2/s, and 1/s = 1e-9 /ns,
// so the conversion is * 1e-9. calculate.py multiplies by 1e9, giving values 1e18 times too large.
inline double rot(double radius, double viscosity, double temperature) {
    if (viscosity == 0) throw std::invalid_argument("viscosity cannot be zero!");
    if (radius == 0) throw std::invalid_argument("radius cannot be zero!");
    const double d = kBoltzmann * temperature / (8 * kPi * viscosity * std::pow(1e-9 * radius, 3));
    return d * 1e-9; // rad^2/s -> rad^2/ns
}

// calc.diff_to_sigma: standard deviation of the step in one time step, sqrt(2 D dt).
inline double to_sigma(double diff, double time_step) { return std::sqrt(2 * diff * time_step); }

// calc.diffusion_cluster: D of a cluster of n particles = smallest D among its particles,
// divided by sqrt(n) for translation (D ~ N^(-1/2)) and, FIX 4 (user's choice), by sqrt(n^3) for
// rotation (D_rot ~ N^(-3/2)); the user's latest Python divides by sqrt(n) for both.
// The starting value 1e14 also acts as an upper limit on D, as in the Python code.
// `diffs` holds the per-particle D values of the cluster's particles.
inline double cluster(const std::vector<double>& diffs, bool rotation) {
    const double n = static_cast<double>(diffs.size());
    double d = 1e14;
    for (double v : diffs)
        if (v < d) d = v;
    return rotation ? d / std::sqrt(n * n * n) : d / std::sqrt(n);
}

} // namespace diffusion
