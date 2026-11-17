#include "steps.hpp"

#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>
#include <string>

#include "clusters.hpp"
#include "diffusion.hpp"
#include "geometry.hpp"

namespace {
constexpr double kPi = 3.14159265358979323846;

// Rotate the patch vectors of particle i by R (rotation.apply in rotate_patches).
void rotate_patches(System& s, int i, const Mat3& R) {
    for (int k = 0; k < s.n_patches(i); ++k) s.patch(i, k) = R * s.patch(i, k);
}

// Three normal numbers as a rotation vector (np.random.normal(0, sigma, size=3)).
Vec3 random_rot_vec(double sigma, Rng& rng) {
    const double x = rng.normal(0.0, sigma);
    const double y = rng.normal(0.0, sigma);
    const double z = rng.normal(0.0, sigma);
    return {x, y, z};
}
} // namespace

StepStats& StepStats::operator+=(const StepStats& o) {
    single_moves += o.single_moves;
    single_accepted += o.single_accepted;
    cluster_moves += o.cluster_moves;
    cluster_accepted += o.cluster_accepted;
    association_tests += o.association_tests;
    associations += o.associations;
    dissociation_tests += o.dissociation_tests;
    dissociations += o.dissociations;
    bond_steps += o.bond_steps;
    return *this;
}

Vec3 step_displacement(double diff, double dt, Rng& rng) {
    const double sigma = std::sqrt(2 * diff * dt); // per-axis standard deviation of a Brownian step
    // FIX 2 (keeps the Python "length + direction" form, with the right distributions):
    // Length: the length of a 3D Gaussian step whose x, y, z each have standard deviation sigma
    // (Maxwell distribution). Python used |N(0, sigma)|, a 1D length, which is too short on average.
    const double g1 = rng.normal(0.0, 1.0), g2 = rng.normal(0.0, 1.0), g3 = rng.normal(0.0, 1.0);
    const double r = sigma * std::sqrt(g1 * g1 + g2 * g2 + g3 * g3);
    // Direction: uniform over the sphere. That needs cos(theta) uniform in [-1, 1]; Python drew
    // theta itself uniformly in [0, pi], which crowds directions towards the poles (+z and -z).
    const double cos_theta = rng.uniform(-1.0, 1.0);
    const double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
    const double phi = rng.uniform(0.0, 2 * kPi);
    return {r * sin_theta * std::cos(phi), r * sin_theta * std::sin(phi), r * cos_theta};
}

void move_particles(System& s, const Params& p, Rng& rng, StepStats& stats) {
    const int n = s.n_particles;
    const double dt = p.sim.time_step;

    // Random order, so that no particle has priority (np.random.shuffle(indexes)).
    std::vector<int> order(n);
    std::iota(order.begin(), order.end(), 0); // 0, 1, ..., n-1  (np.arange)
    rng.shuffle(order);

    std::vector<char> processed(n, 0);
    std::vector<double> d_trans, d_rot;
    for (int i : order) {
        if (processed[i]) continue; // already moved as part of its cluster

        const int type = s.type[i];
        const bool interacting = p.types[type].interacting;
        double diff = s.mobility[type].diff_trans;
        double diff_rot = 0.0;

        // An interacting particle in a cluster moves together with the whole cluster.
        std::vector<int> to_process{i};
        if (interacting && s.cluster_of[i] >= 0) to_process = s.clusters[s.cluster_of[i]];
        if (to_process.size() > 1) { // calc.diffusion_cluster
            d_trans.clear();
            d_rot.clear();
            for (int m : to_process) {
                d_trans.push_back(s.mobility[s.type[m]].diff_trans);
                d_rot.push_back(s.mobility[s.type[m]].diff_rot);
            }
            diff = diffusion::cluster(d_trans, false);
            diff_rot = diffusion::cluster(d_rot, true);
        }

        const Vec3 step = step_displacement(diff, dt, rng);

        if (interacting && to_process.size() > 1) {
            // Cluster: rotate about its centre, then translate; reject if any member would overlap
            // a particle outside the cluster.
            ++stats.cluster_moves;
            const double sigma_rot = diffusion::to_sigma(diff_rot, dt);
            const Mat3 R = rotation_from_rotvec(random_rot_vec(sigma_rot, rng));
            std::vector<Vec3> new_pos = rotate_cluster(s, to_process, R); // FIX 1C inside
            for (Vec3& r : new_pos) r = wrap_into_box(r + step, s.box_limits); // FIX 1A

            bool clash = false;
            for (std::size_t k = 0; k < to_process.size(); ++k)
                if (s.cells.clash(to_process[k], new_pos[k], s.pos, s.radius, to_process)) {
                    clash = true;
                    break;
                }
            if (!clash) {
                ++stats.cluster_accepted;
                for (std::size_t k = 0; k < to_process.size(); ++k) {
                    const int m = to_process[k];
                    rotate_patches(s, m, R);
                    s.pos[m] = new_pos[k];
                    s.cells.move(m, new_pos[k]);
                }
            }
        } else if (interacting) {
            // Single interacting particle.
            ++stats.single_moves;
            // FIX 3: rotation step with sigma_rot = sqrt(2 D_rot dt). Python used D_rot itself as the
            // standard deviation (steps.py line 104).
            const Mat3 R = rotation_from_rotvec(random_rot_vec(s.mobility[type].sigma_rot, rng));
            // FIX 9: translation and rotation are separate attempts. The translation is rejected if
            // it would overlap another particle; the rotation of a sphere about its own centre cannot
            // overlap anything, so it is always accepted (Python rejected it with the translation).
            const Vec3 new_position = wrap_into_box(s.pos[i] + step, s.box_limits); // FIX 1A
            if (!s.cells.clash(i, new_position, s.pos, s.radius, to_process)) {
                ++stats.single_accepted;
                s.pos[i] = new_position;
                s.cells.move(i, new_position);
            }
            rotate_patches(s, i, R);
        } else {
            // Non-interacting particle (crowder): translate only.
            ++stats.single_moves;
            const Vec3 new_position = wrap_into_box(s.pos[i] + step, s.box_limits); // FIX 1A
            if (!s.cells.clash(i, new_position, s.pos, s.radius, to_process)) {
                ++stats.single_accepted;
                s.pos[i] = new_position;
                s.cells.move(i, new_position);
            }
        }
        for (int m : to_process) processed[m] = 1;
    }
}

std::vector<int> count_bonds_per_patch(const System& s) {
    std::vector<int> counts(s.patch_dir.size(), 0);
    for (int i = 0; i < s.n_particles; ++i)
        for (const Bond& b : s.bonds[i]) ++counts[s.patch_begin[i] + b.patch];
    return counts;
}

void update_bonds(System& s, const Params& p, Rng& rng, StepStats& stats) {
    const double L = s.box_limits;
    s.cells.build(s.pos);

    // --- 1. Dissociation. Every existing bond may break with the p_diss of its patch pair.
    // FIX 5a: a bond is tested whether or not its patches are full (Python skipped the test when a
    //         patch was full, so such bonds could never break).
    // FIX 8:  the bond is tested with, and keeps, the patches it formed with (Python re-assigned a
    //         surviving bond to the first facing patch pair found).
    // Bonded particles move together rigidly, so a bond stays in its bonding position.
    std::vector<std::pair<int, int>> broken; // pairs (i < j) whose bond broke in this step
    for (int i = 0; i < s.n_particles; ++i) {
        for (const Bond& b : s.bonds[i]) {
            if (b.partner < i) continue; // each bond once, from its lower-index particle
            // Particles are ordered by type, so type(i) <= type(partner): rows of p_diss are i's patches.
            const InteractionPair* ip = p.find_interaction(s.type[i], s.type[b.partner]);
            ++stats.dissociation_tests;
            if (rng.uniform(0.0, 1.0) < ip->p_diss(b.patch, b.partner_patch)) broken.emplace_back(i, b.partner);
        }
    }
    for (const auto& [i, j] : broken) { // "structured binding": i, j = pair  (Python: for i, j in broken)
        const Bond* b = s.find_bond(i, j);
        --s.n_bonds[s.patch_begin[i] + b->patch];
        --s.n_bonds[s.patch_begin[j] + b->partner_patch];
        s.remove_bond(i, j);
    }
    stats.dissociations += static_cast<std::int64_t>(broken.size());

    // --- 2. Association. Every nearby unbonded pair of interacting particles (in the Python pair
    // order) looks for the first patch pair (a, b) that
    //   - is allowed to bond: interact(a, b) == 1                         (FIX 6: Python ignored it)
    //   - has room on both patches: n_bonds < patches_bonds               (FIX 5a/5c: Python used ==)
    //   - faces each other (Kern-Frenkel) and is within reach;
    // then bonds with probability p_ass(a, b). Bond counts are updated at once, so a patch cannot
    // exceed its limit within one step (FIX 5b: Python used counts from the previous step).
    // A pair whose bond broke in this step is not tested again until the next step (as in Python,
    // each pair is tested once per step).
    std::vector<int> near;
    for (int i = 0; i < s.n_particles; ++i) {
        const ParticleType& ti = p.types[s.type[i]];
        if (!ti.interacting) continue;
        s.cells.neighbours_above(i, near);
        for (int j : near) {
            if (s.find_bond(i, j) != nullptr) continue;
            if (std::find(broken.begin(), broken.end(), std::make_pair(i, j)) != broken.end()) continue;
            const double distance12 = std::sqrt(distance_sq_pbc(L, s.pos[i], s.pos[j]));
            if (!(distance12 <= s.threshold)) continue;
            const ParticleType& tj = p.types[s.type[j]];
            if (!tj.interacting) continue;
            const InteractionPair* ip = p.find_interaction(ti.type, tj.type);
            if (ip == nullptr) continue;

            const Vec3 d_vector = distance_vector_pbc(L, s.pos[j], s.pos[i]); // r_i - r_j
            const Vec3 d_norm = d_vector / norm(d_vector);
            int a_found = -1, b_found = -1;
            for (int a = 0; a < s.n_patches(i) && a_found < 0; ++a) {
                if (s.n_bonds[s.patch_begin[i] + a] >= ti.patch_bonds[a]) continue; // patch a full
                for (int b = 0; b < s.n_patches(j); ++b) {
                    if (ip->interact(a, b) != 1) continue;
                    if (s.n_bonds[s.patch_begin[j] + b] >= tj.patch_bonds[b]) continue; // patch b full
                    if (!patches_aligned(s.patch(i, a), s.patch(j, b), d_norm, ti.patch_alpha_deg[a],
                                         tj.patch_alpha_deg[b]))
                        continue;
                    // Reach: the pair's `distances` entry if the file has one, else the two patch radii.
                    const double r_int =
                        ip->distances.empty() ? ti.patch_radius[a] + tj.patch_radius[b] : ip->distances(a, b);
                    if (!(distance12 < ti.radius + tj.radius + r_int)) continue;
                    a_found = a;
                    b_found = b;
                    break;
                }
            }
            if (a_found < 0) continue;

            ++stats.association_tests;
            if (rng.uniform(0.0, 1.0) < ip->p_ass(a_found, b_found)) {
                ++stats.associations;
                s.add_bond(i, j, a_found, b_found);
                ++s.n_bonds[s.patch_begin[i] + a_found];
                ++s.n_bonds[s.patch_begin[j] + b_found];
            }
        }
    }
    // Each bond is listed under both of its particles, hence the division by 2.
    std::int64_t listed = 0;
    for (int i = 0; i < s.n_particles; ++i) listed += static_cast<std::int64_t>(s.bonds[i].size());
    stats.bond_steps += listed / 2;
}

OverlapReport find_overlaps(const System& s) {
    OverlapReport rep;
    std::vector<int> near;
    for (int i = 0; i < s.n_particles; ++i) {
        s.cells.neighbours_above(i, near);
        for (int j : near) {
            const double limit = s.radius[i] + s.radius[j];
            const double d2 = distance_sq_pbc(s.box_limits, s.pos[i], s.pos[j]);
            if (d2 < limit * limit) {
                ++rep.pairs;
                rep.deepest = std::max(rep.deepest, limit - std::sqrt(d2));
            }
        }
    }
    return rep;
}

void initial_bonds(System& s, const Params& p, Rng& rng, StepStats& stats) {
    s.n_bonds.assign(s.patch_dir.size(), 0);
    update_bonds(s, p, rng, stats);
    update_clusters(s);
}

void simulation_step(System& s, const Params& p, Rng& rng, StepStats& stats) {
    move_particles(s, p, rng, stats); // positions stay inside the box (FIX 1A)
    update_bonds(s, p, rng, stats);   // keeps s.n_bonds current (FIX 5b)
    update_clusters(s);
}
