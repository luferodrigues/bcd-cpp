// Clusters of bonded particles and rigid-body rotation (update_clusters, get_pbc_neighbor and
// rotate_cluster in steps.py).
#pragma once

#include <vector>

#include "geometry.hpp"
#include "system.hpp"

// steps.update_clusters: find groups of bonded particles with a depth-first search, in exactly
// the order of the Python code (start from the lowest unvisited index; neighbours pushed in
// ascending order and taken from the end of the stack). Fills s.clusters and s.cluster_of.
void update_clusters(System& s);

// steps.get_pbc_neighbor: the periodic image of v (v, v + 2L or v - 2L, per coordinate) closest
// to ref.
Vec3 pbc_neighbor(const Vec3& v, const Vec3& ref, double box_limits);

// steps.rotate_cluster: new positions of `members` after rotating the cluster by R about its
// centre. The cluster is first made whole across the box edge, so that its centre is right even
// when it crosses the edge. FIX 1C: this follows the bonds (each member is placed at the image
// closest to the bonded neighbour it is reached from), not the distance to the first member as
// in Python, which fails for clusters reaching about half a box from their first member.
// The returned positions are not wrapped back into the box (the caller does that).
std::vector<Vec3> rotate_cluster(const System& s, const std::vector<int>& members, const Mat3& R);
