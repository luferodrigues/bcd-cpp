#include "system.hpp"

#include "diffusion.hpp"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <string>

namespace {
// Give up placing a particle after this many tries: the box is then too crowded for random
// insertion (the Python loop would run forever instead).
constexpr std::int64_t kMaxPlacementAttempts = 10'000'000;

// Random insertion of spheres gets stuck ("jams") at a volume fraction of about 0.38.
constexpr double kWarnVolumeFraction = 0.30;
constexpr double kPi = 3.14159265358979323846;
} // namespace

namespace {
// A uniformly random rotation (Shoemake's method: a uniformly random unit quaternion, i.e. a random
// point on the 4D unit sphere, built from three uniform numbers).
Mat3 random_rotation(Rng& rng) {
    const double u1 = rng.uniform(0.0, 1.0), u2 = rng.uniform(0.0, 1.0), u3 = rng.uniform(0.0, 1.0);
    const double a = std::sqrt(1.0 - u1), b = std::sqrt(u1);
    const double x = a * std::sin(2 * kPi * u2), y = a * std::cos(2 * kPi * u2);
    const double z = b * std::sin(2 * kPi * u3), w = b * std::cos(2 * kPi * u3);
    return rotation_from_quaternion(w, x, y, z);
}
} // namespace

System empty_system(const Params& p) {
    System s;
    const double L = p.sim.box_limits;
    s.box_limits = L;
    s.n_particles = p.total_particles();
    if (s.n_particles == 0) throw std::runtime_error("no particles: all 'number' values are 0");

    // Particle types in file order: all particles of type 0, then type 1, ... (initialize_particles)
    for (const ParticleType& t : p.types)
        for (int k = 0; k < t.number; ++k) {
            s.type.push_back(t.type);
            s.radius.push_back(t.radius);
        }
    s.pos.assign(s.n_particles, Vec3{0.0, 0.0, 0.0});

    // Patch layout: interacting particles have their type's patches, others none.
    s.patch_begin.push_back(0);
    for (int i = 0; i < s.n_particles; ++i) {
        const ParticleType& t = p.types[s.type[i]];
        const int n = t.interacting ? t.n_patches() : 0;
        s.patch_begin.push_back(s.patch_begin.back() + n);
    }
    s.patch_dir.assign(s.patch_begin.back(), Vec3{0.0, 0.0, 0.0});

    // Diffusion per type (run_simulation): values from the file if given, else Stokes-Einstein.
    for (const ParticleType& t : p.types) {
        TypeMobility m;
        m.diff_trans = t.diff_trans ? *t.diff_trans : diffusion::trans(t.radius, p.sim.viscosity, p.sim.temperature);
        m.sigma_trans = diffusion::to_sigma(m.diff_trans, p.sim.time_step);
        m.rotates = t.interacting;
        if (m.rotates) {
            m.diff_rot = t.diff_rot ? *t.diff_rot : diffusion::rot(t.radius, p.sim.viscosity, p.sim.temperature);
            m.sigma_rot = diffusion::to_sigma(m.diff_rot, p.sim.time_step);
        }
        s.mobility.push_back(m);
    }

    // Neighbour-cell settings (utils.max_distance with manual_value = 0, then run_simulation).
    double max_hs = 0.0, max_patch = 0.0;
    for (const ParticleType& t : p.types) {
        max_hs = std::max(max_hs, t.radius);
        for (double r : t.patch_radius) max_patch = std::max(max_patch, r);
    }
    s.threshold = 2.0 * max_hs + 2.0 * max_patch;
    s.n_cells_1d = static_cast<int>(std::floor(2.0 * L / s.threshold));
    if (s.n_cells_1d < 1)
        throw std::runtime_error("box is smaller than the interaction threshold (" + std::to_string(s.threshold) +
                                 "); increase box_limits");
    s.cell_size = 2.0 * L / s.n_cells_1d;
    s.cells = CellList(L, s.cell_size, s.n_cells_1d);

    double particle_volume = 0.0;
    for (double r : s.radius) particle_volume += 4.0 / 3.0 * kPi * r * r * r;
    s.volume_fraction = particle_volume / std::pow(2.0 * L, 3);

    s.bonds.assign(s.n_particles, {});
    s.cluster_of.assign(s.n_particles, -1);
    s.n_bonds.assign(s.patch_dir.size(), 0); // initialize.initialize_n_bonds
    return s;
}

System build_system(const Params& p, Rng& rng) {
    System s = empty_system(p);
    const double L = s.box_limits;
    if (s.n_cells_1d < 3)
        std::cerr << "warning: only " << s.n_cells_1d << " neighbour cell(s) per side; the cell list will "
                  << "visit some cells more than once\n";

    // Fail fast if the particles cannot fit, instead of searching for a free spot forever.
    if (s.volume_fraction >= 1.0)
        throw std::runtime_error("the particles fill " + std::to_string(s.volume_fraction) +
                                 " times the box volume; they cannot fit");
    if (s.volume_fraction > kWarnVolumeFraction)
        std::cerr << "warning: volume fraction " << s.volume_fraction << " is high; random placement "
                  << "may be very slow or fail (it jams near 0.38)\n";

    // 1. Uniform random coordinates in [-L, L) (np.random.uniform(-L, L, size=(N, 3)))
    for (Vec3& r : s.pos) r = {rng.uniform(-L, L), rng.uniform(-L, L), rng.uniform(-L, L)};

    // 2. Patches: rotate the body-frame vectors by a random rotation, then normalize
    //    (initialize_patches). FIX 10: a uniformly random rotation (every orientation equally
    //    likely). Python used Rotation.from_rotvec(uniform(0, 1, 3)), which only reaches rotations
    //    about axes in one octant, by at most ~1.7 rad.
    for (int i = 0; i < s.n_particles; ++i) {
        if (s.n_patches(i) == 0) continue;
        const ParticleType& t = p.types[s.type[i]];
        const Mat3 R = random_rotation(rng);
        for (int k = 0; k < s.n_patches(i); ++k) {
            const Vec3 v = R * t.patch_positions[k];
            s.patch(i, k) = v / norm(v);
        }
    }

    // 3. Remove overlaps: each particle in turn is moved to a new random position until it
    //    overlaps nobody. Later particles are checked against all others, so the final
    //    configuration has no overlaps. The cell list only looks at nearby particles; this works
    //    because a cell is at least `threshold` = 2*max radius + 2*max patch radius wide.
    s.cells.build(s.pos);
    for (int i = 0; i < s.n_particles; ++i) {
        const std::vector<int> skip{i}; // don't compare i with itself
        std::int64_t attempts = 0;
        while (s.cells.clash(i, s.pos[i], s.pos, s.radius, skip)) {
            if (++attempts > kMaxPlacementAttempts)
                throw std::runtime_error("could not place particle " + std::to_string(i) + " without overlap after " +
                                         std::to_string(kMaxPlacementAttempts) + " attempts; the box is too crowded");
            s.pos[i] = {rng.uniform(-L, L), rng.uniform(-L, L), rng.uniform(-L, L)};
            s.cells.move(i, s.pos[i]);
            ++s.init_redraws;
        }
    }
    return s;
}
