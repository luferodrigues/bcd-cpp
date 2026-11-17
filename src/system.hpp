// Simulation state (the `simulation` dict of simulation.py) and its initialization
// (initialize.py plus the set-up part of run_simulation).
#pragma once

#include <algorithm>
#include <cstdint>
#include <vector>

#include "cells.hpp"
#include "geometry.hpp"
#include "params.hpp"
#include "rng.hpp"

// Diffusion of one particle type (the diff_trans / sigma_trans / diff_rot / sigma_rot entries
// that run_simulation adds to part_params).
struct TypeMobility {
    double diff_trans = 0.0;  // nm^2/ns: from the file, or Stokes-Einstein
    double sigma_trans = 0.0; // sqrt(2 D dt)
    bool rotates = false;     // only interacting types rotate
    double diff_rot = 0.0;    // only set if rotates
    double sigma_rot = 0.0;
};

// One bond, seen from one of its two particles (simulation['interactions'][i][j] = [patch,
// partner_patch] in Python). Each bond is stored twice: in the list of each particle.
struct Bond {
    int partner;       // the other particle
    int patch;         // this particle's patch
    int partner_patch; // the other particle's patch
};

struct System {
    double box_limits = 0.0; // the box is [-box_limits, box_limits]^3
    int n_particles = 0;

    // One entry per particle, indexed 0 .. n_particles-1.
    std::vector<int> type;      // simulation['particles']
    std::vector<double> radius; // simulation['radii']
    std::vector<Vec3> pos;      // simulation['coordinates']

    // simulation['patches']: unit patch vectors in the world (box) frame.
    // Stored flat instead of as a list of lists: the patches of particle i are
    //   patch_dir[patch_begin[i]], ..., patch_dir[patch_begin[i+1] - 1]
    // Use n_patches(i) and patch(i, k) rather than indexing patch_dir directly.
    // Particles of non-interacting types have no patches (as in initialize_patches).
    std::vector<int> patch_begin; // n_particles + 1 entries
    std::vector<Vec3> patch_dir;

    // One entry per particle type.
    std::vector<TypeMobility> mobility;

    // Neighbour-search settings (threshold, n_cells_1d, cell_size in simulation.py) and the
    // cell list itself, kept up to date with `pos`.
    double threshold = 0.0;
    int n_cells_1d = 0;
    double cell_size = 0.0;
    CellList cells;

    // simulation['interactions'], stored per particle: bonds[i] lists the bonds of particle i.
    // (Python keeps an N x N x 2 table; this list form needs memory only for existing bonds.)
    std::vector<std::vector<Bond>> bonds;

    // simulation['clusters']: groups of >= 2 bonded particles, in the order the Python
    // update_clusters produces them (single particles are not listed, as in Python).
    std::vector<std::vector<int>> clusters;
    std::vector<int> cluster_of; // index into `clusters` for each particle, -1 if in none

    // simulation['n_bonds']: number of bonds per patch, flat like patch_dir (patch k of particle i
    // is n_bonds[patch_begin[i] + k]). As in Python this table is ONE STEP OLD when the bonds are
    // updated: simulation_step stores the counts of the state at the START of the step.
    std::vector<int> n_bonds;

    // Diagnostics: fraction of the box volume filled by particles, and how many times a
    // particle was re-placed because it overlapped another.
    double volume_fraction = 0.0;
    std::int64_t init_redraws = 0;

    int n_patches(int i) const { return patch_begin[i + 1] - patch_begin[i]; }
    const Vec3& patch(int i, int k) const { return patch_dir[patch_begin[i] + k]; }
    Vec3& patch(int i, int k) { return patch_dir[patch_begin[i] + k]; }

    // The bond between i and j as seen from i, or nullptr if they are not bonded.
    const Bond* find_bond(int i, int j) const {
        for (const Bond& b : bonds[i])
            if (b.partner == j) return &b;
        return nullptr;
    }
    void add_bond(int i, int j, int patch_i, int patch_j) {
        bonds[i].push_back({j, patch_i, patch_j});
        bonds[j].push_back({i, patch_j, patch_i});
    }
    void remove_bond(int i, int j) {
        auto drop = [](std::vector<Bond>& list, int partner) {
            list.erase(std::remove_if(list.begin(), list.end(), [&](const Bond& b) { return b.partner == partner; }),
                       list.end());
        };
        drop(bonds[i], j);
        drop(bonds[j], i);
    }
};

// Everything that does not need random numbers: types, radii, patch layout, diffusion
// coefficients, neighbour-cell settings, empty bonds. Positions and patch directions are zero.
System empty_system(const Params& p);

// Builds the initial state, drawing random numbers in the same order as the Python code:
//   1. all coordinates, uniform in the box             (initialize_coordinates)
//   2. one random orientation per interacting particle (initialize_patches)
//   3. particle by particle, redraw the position until it overlaps no other particle
//      (the loop over update_distances in run_simulation)
System build_system(const Params& p, Rng& rng);
