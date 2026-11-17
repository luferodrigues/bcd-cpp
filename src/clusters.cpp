#include "clusters.hpp"

#include <algorithm>
#include <stdexcept>
#include <unordered_map>

void update_clusters(System& s) {
    const int n = s.n_particles;
    std::vector<char> visited(n, 0); // std::vector<char> as a list of yes/no flags
    std::vector<int> stack, cluster, neighbours;
    s.clusters.clear();
    s.cluster_of.assign(n, -1);

    for (int i = 0; i < n; ++i) {
        if (visited[i]) continue;
        stack.assign(1, i); // stack = [i]
        cluster.clear();
        while (!stack.empty()) {
            const int p = stack.back(); // p = stack.pop()
            stack.pop_back();
            if (visited[p]) continue;
            visited[p] = 1;
            cluster.push_back(p);
            // np.where(bonded[p])[0]: the bonded partners of p in ascending order
            neighbours.clear();
            for (const Bond& b : s.bonds[p]) neighbours.push_back(b.partner);
            std::sort(neighbours.begin(), neighbours.end());
            for (int nb : neighbours)
                if (!visited[nb]) stack.push_back(nb);
        }
        if (cluster.size() > 1) {
            const int index = static_cast<int>(s.clusters.size());
            for (int m : cluster) s.cluster_of[m] = index;
            s.clusters.push_back(cluster);
        }
    }
}

Vec3 pbc_neighbor(const Vec3& v, const Vec3& ref, double box_limits) {
    Vec3 out;
    for (int k = 0; k < 3; ++k) {
        const double cands[3] = {v[k], v[k] + 2 * box_limits, v[k] - 2 * box_limits};
        int best = 0; // np.argmin: the first of equal minima
        double best_d2 = (ref[k] - cands[0]) * (ref[k] - cands[0]);
        for (int c = 1; c < 3; ++c) {
            const double d2 = (ref[k] - cands[c]) * (ref[k] - cands[c]);
            if (d2 < best_d2) { best = c; best_d2 = d2; }
        }
        out[k] = cands[best];
    }
    return out;
}

std::vector<Vec3> rotate_cluster(const System& s, const std::vector<int>& members, const Mat3& R) {
    const std::size_t n = members.size();

    // slot[particle] = its position in `members` (a dict in Python terms).
    std::unordered_map<int, std::size_t> slot;
    for (std::size_t k = 0; k < n; ++k) slot[members[k]] = k;

    // FIX 1C: make the cluster whole by walking along its bonds (breadth-first search), starting
    // from the first member. Each newly reached member is put at its periodic image closest to the
    // already placed neighbour it is bonded to. Bonds are short, so this is always the right image.
    std::vector<Vec3> whole(n);
    std::vector<char> placed(n, 0);
    std::vector<std::size_t> queue{0};
    whole[0] = s.pos[members[0]];
    placed[0] = 1;
    for (std::size_t q = 0; q < queue.size(); ++q) {
        const std::size_t k = queue[q];
        for (const Bond& b : s.bonds[members[k]]) {
            const auto found = slot.find(b.partner); // like slot.get(b.partner) in Python
            if (found == slot.end()) continue;       // partner not in this cluster (cannot happen)
            const std::size_t kb = found->second;
            if (placed[kb]) continue;
            whole[kb] = pbc_neighbor(s.pos[b.partner], whole[k], s.box_limits);
            placed[kb] = 1;
            queue.push_back(kb);
        }
    }
    if (queue.size() != n) throw std::logic_error("rotate_cluster: cluster members are not all connected by bonds");

    Vec3 com{0.0, 0.0, 0.0};
    for (const Vec3& r : whole) com = com + r;
    com = com / static_cast<double>(n);

    std::vector<Vec3> out(n);
    for (std::size_t p = 0; p < n; ++p) out[p] = R * (whole[p] - com) + com;
    return out;
}
