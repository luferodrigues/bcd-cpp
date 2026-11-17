// Parsing of the three .par input files (the "import" half of utils.py).
#pragma once

#include <array>
#include <cstdint>
#include <iosfwd>
#include <optional>
#include <string>
#include <vector>

#include "geometry.hpp" // Vec3

// Row-major rectangular matrix (rows x cols).
template <typename T>
struct Matrix {
    int rows = 0;
    int cols = 0;
    std::vector<T> data;

    bool empty() const { return data.empty(); }
    T& operator()(int i, int j) { return data[static_cast<std::size_t>(i) * cols + j]; }
    const T& operator()(int i, int j) const { return data[static_cast<std::size_t>(i) * cols + j]; }
};

// simulation_params.par
struct SimParams {
    double box_limits = 0.0;     // coordinates live in [-box_limits, box_limits]
    double temperature = 0.0;
    double viscosity = 0.0;
    double time_step = 0.0;      // Python casts this to int; see params.cpp
    std::int64_t n_steps = 0;
    std::int64_t n_interval = 0; // sampling stride
    int n_dimensions = 3;
};

// One block of particle_params.par
struct ParticleType {
    int type = 0;
    std::string name;
    int number = 0;
    double radius = 0.0;
    bool interacting = false;
    std::vector<double> patch_radius;
    std::vector<Vec3> patch_positions;   // body frame, as written in the file (not normalized)
    std::vector<double> patch_alpha_deg; // Kern-Frenkel half-angle, degrees
    std::vector<int> patch_bonds;        // max bonds per patch
    // Optional overrides; if absent, computed from radius, viscosity and temperature
    // (Stokes-Einstein, as in simulation.py).
    std::optional<double> diff_trans;
    std::optional<double> diff_rot;

    int n_patches() const { return static_cast<int>(patch_positions.size()); }
};

// One "type a,b" block of interaction_params.par. Matrices are n_patches(a) x n_patches(b);
// any that are absent from the file are left empty.
struct InteractionPair {
    int type_a = 0;
    int type_b = 0;
    Matrix<int> interact;
    Matrix<double> distances;
    Matrix<double> energies;
    Matrix<double> p_ass;
    Matrix<double> p_diss;
};

struct Params {
    SimParams sim;
    std::vector<ParticleType> types;          // indexed by type id
    std::vector<InteractionPair> interactions;

    int total_particles() const;
    // Interaction block for a type pair in either order, or nullptr if they don't interact.
    const InteractionPair* find_interaction(int a, int b) const;
};

SimParams load_sim_params(const std::string& path);
std::vector<ParticleType> load_particle_params(const std::string& path);
std::vector<InteractionPair> load_interaction_params(const std::string& path);

// Loads all three files and cross-checks them (patch counts vs matrix sizes, known types, ...).
// Throws std::runtime_error with a file:line message on any problem.
Params load_params(const std::string& sim_path, const std::string& parts_path,
                   const std::string& inters_path);

// utils.add_probs_from_energies (called for every pair in run.py): where p_ass / p_diss are not given,
// compute them from the energies (in kT): p_ass = exp(-E), p_diss = exp(E), with values below 1e-8 set
// to 0, above 1 - 1e-8 set to 1, and p_ass = 0 where E == 0. Throws if neither is available.
void add_probs_from_energies(Params& p);

// Canonical text dump, matched line-for-line by tools/dump_params.py for validation.
void dump_params(const Params& p, std::ostream& os);
