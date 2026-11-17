// One simulation step (steps.py, plus simulation_step from simulation.py).
//
// Port of the user's latest Python code (verified step by step against it, see tools/test_steps.py
// in archive bcd-cpp-phase5d), followed by the fixes the user approved. Every deviation from the
// Python code is marked "FIX <number>" (numbers as in CLAUDE.md, "Decisions"). The pair order and
// the decision rules otherwise follow the Python code.
#pragma once

#include <cstdint>
#include <vector>

#include "params.hpp"
#include "rng.hpp"
#include "system.hpp"

// Counters for one or more steps (not in the Python code; for information).
struct StepStats {
    std::int64_t single_moves = 0, single_accepted = 0;   // particles moving alone
    std::int64_t cluster_moves = 0, cluster_accepted = 0; // clusters moving as rigid bodies
    std::int64_t association_tests = 0, associations = 0; // unbonded pairs in bonding position
    std::int64_t dissociation_tests = 0, dissociations = 0;
    std::int64_t bond_steps = 0; // sum over steps of the number of bonds (bond_steps / dissociations
                                 // = average bond lifetime in steps, in a steady state)

    StepStats& operator+=(const StepStats& o);
};

// steps.generate_step_coords: random displacement for diffusion coefficient `diff`, as a length
// and a direction (FIX 2: length = that of a 3D Gaussian step with per-axis sigma sqrt(2 diff dt);
// direction uniform over the sphere). Each of x, y, z then has variance 2 diff dt.
Vec3 step_displacement(double diff, double dt, Rng& rng);

// steps.generate_next_step: every particle or cluster, in random order, tries one move: clusters
// rotate + translate as one rigid move; single interacting particles try a translation and,
// separately, a rotation (FIX 9); crowders only translate. A move that would overlap another particle
// is rejected. Accepted positions are wrapped into the box at once (FIX 1A).
void move_particles(System& s, const Params& p, Rng& rng, StepStats& stats);

// utils.find_bonds_all: number of bonds of every patch, counted from the bonds (flat, like
// System::n_bonds, which update_bonds keeps equal to this).
std::vector<int> count_bonds_per_patch(const System& s);

// steps.update_interactions, with FIXES 5a, 5b, 5c, 6, 8: first every existing bond may break
// (p_diss of its own patch pair); then every nearby unbonded pair of interacting particles bonds with
// probability p_ass through the first patch pair that is allowed (interact = 1), has room on both
// patches, faces each other and is within reach. Keeps s.n_bonds current. Rebuilds s.cells.
void update_bonds(System& s, const Params& p, Rng& rng, StepStats& stats);

// Diagnostic (not in Python): number of pairs of particles that overlap, and the deepest overlap
// (sum of radii minus distance). Uses the cell list, so s.cells must be up to date.
struct OverlapReport {
    int pairs = 0;
    double deepest = 0.0;
};
OverlapReport find_overlaps(const System& s);

// The bonds, clusters and bond counts of the starting configuration (end of the set-up in
// run_simulation).
void initial_bonds(System& s, const Params& p, Rng& rng, StepStats& stats);

// simulation.simulation_step: moves, bonds, clusters.
void simulation_step(System& s, const Params& p, Rng& rng, StepStats& stats);
