#include "transitions.hpp"

#include <map>

namespace {
// cluster index of every particle, -1 for particles in no cluster (utils.particle_cluster_indexing)
std::vector<int> cluster_index(const std::vector<std::vector<int>>& clusters, int n_particles) {
    std::vector<int> idx(n_particles, -1);
    for (int c = 0; c < static_cast<int>(clusters.size()); ++c)
        for (int m : clusters[c]) idx[m] = c;
    return idx;
}

// The pieces of one cluster: how many of its particles are in each cluster of the other step
// (std::map keeps them in cluster order), then one piece of size 1 per single particle.
std::vector<int> pieces_of(const std::vector<int>& members, const std::vector<int>& other_index) {
    std::map<int, int> per_cluster; // other cluster index -> count (like a dict sorted by key)
    int singles = 0;
    for (int m : members) {
        if (other_index[m] >= 0) ++per_cluster[other_index[m]];
        else ++singles;
    }
    std::vector<int> pieces;
    for (const auto& [cluster, count] : per_cluster) pieces.push_back(count);
    pieces.insert(pieces.end(), singles, 1); // `singles` pieces of size 1
    return pieces;
}
} // namespace

std::vector<Transition> find_transitions(const std::vector<std::vector<int>>& previous,
                                         const std::vector<std::vector<int>>& current, int n_particles) {
    const std::vector<int> prev_index = cluster_index(previous, n_particles);
    const std::vector<int> curr_index = cluster_index(current, n_particles);
    std::vector<Transition> out;

    // Dissociations (assign_transitions_row): a previous cluster that fell apart.
    for (const auto& cluster : previous) {
        const std::vector<int> pieces = pieces_of(cluster, curr_index);
        int size = static_cast<int>(cluster.size());
        for (std::size_t k = 1; k < pieces.size(); ++k) { // pieces[0] is what remains
            out.push_back({{size}, {size - pieces[k], pieces[k]}});
            size -= pieces[k];
        }
    }
    // Associations (assign_transitions_column): a current cluster that was put together.
    for (const auto& cluster : current) {
        const std::vector<int> pieces = pieces_of(cluster, prev_index);
        int size = pieces.empty() ? 0 : pieces[0];
        for (std::size_t k = 1; k < pieces.size(); ++k) {
            out.push_back({{size, pieces[k]}, {size + pieces[k]}});
            size += pieces[k];
        }
    }
    return out;
}
