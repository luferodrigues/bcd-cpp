// Cluster transitions between two consecutive steps (transitions.csv): the utils.py functions
// clusters_origins_reference ... assign_transitions, with FIX 7.
#pragma once

#include <vector>

// One line of transitions.csv: sizes before "-" and after it.
//   association:  before = {a, b},  after = {a + b}      -> "frame,a,b,-,a+b"
//   dissociation: before = {a + b}, after = {a, b}       -> "frame,a+b,-,a,b"
struct Transition {
    std::vector<int> before;
    std::vector<int> after;
};

// Compares the clusters of the previous step with the current ones (both as in System::clusters:
// only groups of >= 2 particles; `n_particles` in total). Same method as the Python code:
//   - For every previous cluster, its particles now form "pieces": parts of current clusters, in
//     current-cluster order, then single particles. With more than one piece, the first piece is
//     what remains and the others split off one by one (dissociations).
//   - For every current cluster, its particles came from pieces: parts of previous clusters, in
//     previous-cluster order, then former single particles. With more than one piece, they join the
//     first piece one by one (associations).
// All dissociations come first (previous-cluster order), then all associations.
// FIX 7 (vs Python): (a) each event is listed once (Python repeated an event's lines once per line);
// (b) single particles always join or leave one at a time (Python lumped them together whenever
// other clusters existed); (c) works when the current step has no clusters (Python skipped it).
std::vector<Transition> find_transitions(const std::vector<std::vector<int>>& previous,
                                         const std::vector<std::vector<int>>& current, int n_particles);
