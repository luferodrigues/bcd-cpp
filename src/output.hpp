// Output files (the export half of utils.py and the writing part of run_simulation).
// Formats are identical to the Python writers (checked by tools/test_outputs.py).
#pragma once

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <ostream>
#include <string>
#include <vector>

#include "params.hpp"
#include "system.hpp"
#include "transitions.hpp"

// A double formatted exactly as Python prints it: repr(x), str(x) and f"{x}" for floats and
// NumPy float64 all give the same text, e.g. 70.0, 0.30000000000000004, 1e-05, -3.25.
std::string py_repr(double v);

// Folder the output files go into: os.path.join(folder, prefix) in simulation.py.
std::filesystem::path output_dir(const std::string& folder, const std::string& prefix);

// Atom names for visualization.xyz and vmd_init.tcl (utils.build_label_dictionaries):
// centre of a type-t particle = labels_com[t]; patches are numbered across all types in file order
// (type 0's patches first, then type 1's, ...), including non-interacting types.
struct AtomLabels {
    std::vector<std::string> centre;             // per type
    std::vector<std::vector<std::string>> patch; // per type, per patch
};
AtomLabels build_labels(const Params& p);

// ---- Single writers (one call = what one call of the Python writer writes) ----
void write_particle_types(std::ostream& os, const System& s);                    // particles.csv
void write_tcl(std::ostream& os, const Params& p, const AtomLabels& labels);       // vmd_init.tcl
void write_trajectory_frame(std::ostream& os, const System& s, std::int64_t frame, std::int64_t n_steps);
void write_xyz_frame(std::ostream& os, const System& s, const AtomLabels& labels, std::int64_t frame,
                     std::int64_t n_steps);                                          // visualization.xyz
void write_clusters_frame(std::ostream& os, const std::vector<std::vector<int>>& clusters,
                          std::int64_t frame);                                       // clusters.csv
void write_cluster_sizes(std::ostream& os, const std::vector<std::vector<int>>& clusters); // cluster_agg.csv
void write_int_line(std::ostream& os, std::int64_t value); // cluster_n.csv, cluster_max.csv
void write_n_bonds_frame(std::ostream& os, const System& s, std::int64_t frame); // n_bonds.csv
void write_transitions(std::ostream& os, std::int64_t frame, const std::vector<Transition>& events);

// Largest cluster size, 0 if there are no clusters (cluster_max.csv).
std::int64_t largest_cluster(const std::vector<std::vector<int>>& clusters);

// ---- All files of one run ----
struct OutputOptions {
    bool transitions = true;         // -nt / --no-transitions turns this off
    bool clusters_every_step = true; // -cs / --no-clusters-steps: clusters.csv only at saved frames
};

class OutputFiles {
public:
    // Creates (overwrites) all files in `dir` and writes particles.csv and vmd_init.tcl.
    OutputFiles(const std::filesystem::path& dir, const Params& p, const System& s, OutputOptions options);

    // Call after every step (1, 2, ..., n_steps). `previous_clusters` = the clusters before the step.
    // Saved frames: step 1 and every n_interval steps, as in run_simulation.
    void after_step(std::int64_t step, const System& s, const std::vector<std::vector<int>>& previous_clusters);

private:
    std::ofstream open(const std::string& name);

    std::filesystem::path dir_;
    AtomLabels labels_;
    OutputOptions options_;
    std::int64_t n_steps_, n_interval_;
    std::ofstream trajectory_, xyz_, clusters_, cluster_agg_, cluster_n_, cluster_max_, n_bonds_, transitions_;
};
