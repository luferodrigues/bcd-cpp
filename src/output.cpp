#include "output.hpp"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <string_view>

namespace fs = std::filesystem;

std::string py_repr(double v) {
    if (std::isnan(v)) return "nan";
    if (std::isinf(v)) return v > 0 ? "inf" : "-inf";
    if (v == 0.0) return std::signbit(v) ? "-0.0" : "0.0";

    // std::to_chars gives the shortest digits that read back to exactly v (the same digits
    // Python chooses), here in scientific form such as "-1.2345e+02" or "5e-05".
    char buf[64];
    const auto res = std::to_chars(buf, buf + sizeof buf, v, std::chars_format::scientific);
    std::string_view sci(buf, static_cast<std::size_t>(res.ptr - buf));

    std::string out;
    if (sci.front() == '-') { out = "-"; sci.remove_prefix(1); }
    const std::size_t e = sci.find('e');
    std::string digits; // significant digits without the decimal point, e.g. "12345"
    for (char c : sci.substr(0, e))
        if (c != '.') digits += c;
    const int exp = std::stoi(std::string(sci.substr(e + 1)));
    const int nd = static_cast<int>(digits.size());

    // Python uses plain decimal notation when -4 <= exponent < 16, scientific otherwise.
    if (exp >= 16 || exp < -4) {
        out += digits[0];
        if (nd > 1) out += "." + digits.substr(1);
        out += exp < 0 ? "e-" : "e+";
        const std::string ed = std::to_string(std::abs(exp));
        if (ed.size() < 2) out += '0'; // at least two exponent digits: 1e-05, 1e+16
        out += ed;
    } else if (exp >= 0) {
        if (nd <= exp + 1) out += digits + std::string(exp + 1 - nd, '0') + ".0"; // 70.0, 1000.0
        else out += digits.substr(0, exp + 1) + "." + digits.substr(exp + 1);   // 3.25
    } else {
        out += "0." + std::string(-exp - 1, '0') + digits; // 0.0072
    }
    return out;
}

fs::path output_dir(const std::string& folder, const std::string& prefix) {
    return prefix.empty() ? fs::path(folder) : fs::path(folder) / prefix;
}

// ---------------------------------------------------------------------------
// Atom names (labels_com and labels_patch in simulation.run_simulation)
// ---------------------------------------------------------------------------
namespace {
const std::vector<std::string> kLabelsCom = {
    "C",  "N",  "O",  "F",  "B",  "Ne", "Al", "Si", "P",  "S",  "Cl", "Ga",  "Ge", "As",  "Se", "Br", "Kr", "In",
    "Sn", "Sb", "Te", "I",  "Xe", "Tl", "Pb", "Bi", "Po", "At", "Rn", "UUt", "Fl", "Uup", "Lv", "Uus", "Uuo"};
const std::vector<std::string> kLabelsPatch = {
    "H",  "He", "Li", "Be", "Na", "Mg", "K",  "Ca", "Sc", "Ti", "V",  "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn",
    "Rb", "Sr", "Y",  "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "Cs", "Ba", "Hf", "Ta", "W",  "Re",
    "Os", "Ir", "Pt", "Au", "Hg", "Fr", "Ra", "Rf", "Db", "Sg", "Bh", "Hs", "Mt", "Ds", "Rg", "Cn"};
} // namespace

AtomLabels build_labels(const Params& p) {
    AtomLabels labels;
    std::size_t patch_counter = 0;
    for (const ParticleType& t : p.types) {
        if (t.type >= static_cast<int>(kLabelsCom.size()))
            throw std::runtime_error("too many particle types for the atom names (" + std::to_string(kLabelsCom.size()) +
                                     " at most)");
        labels.centre.push_back(kLabelsCom[t.type]);
        labels.patch.emplace_back();
        for (int k = 0; k < static_cast<int>(t.patch_radius.size()); ++k) {
            if (patch_counter >= kLabelsPatch.size())
                throw std::runtime_error("too many patches in total for the atom names (" +
                                         std::to_string(kLabelsPatch.size()) + " at most)");
            labels.patch.back().push_back(kLabelsPatch[patch_counter++]);
        }
    }
    return labels;
}

// ---------------------------------------------------------------------------
// Writers. Comments name the Python function each one reproduces.
// ---------------------------------------------------------------------------
void write_particle_types(std::ostream& os, const System& s) { // utils.write_particle_types
    for (int t : s.type) os << t << '\n';
}

void write_tcl(std::ostream& os, const Params& p, const AtomLabels& labels) { // utils.write_tcl
    const double box_length = 2 * p.sim.box_limits;
    os << "package require pbctools\n\n"
       << "mol modstyle 0 top VDW\n"
       << "foreach {elem rad} {\n";
    for (const ParticleType& t : p.types) {
        os << "    " << labels.centre[t.type] << ' ' << py_repr(t.radius) << '\n';
        for (std::size_t k = 0; k < t.patch_radius.size(); ++k)
            os << "    " << labels.patch[t.type][k] << ' ' << py_repr(t.patch_radius[k]) << '\n';
    }
    os << "} {\n"
       << "    set sel [atomselect top \"element $elem\"]\n"
       << "    $sel set radius $rad\n"
       << "}\n"
       << "mol reanalyze top\n"
       << "pbc set {" << py_repr(box_length) << ' ' << py_repr(box_length) << ' ' << py_repr(box_length) << "} -all\n"
       << "pbc box -center origin\n"
       << "display projection orthographic"; // no newline at the end, as in Python
}

void write_trajectory_frame(std::ostream& os, const System& s, std::int64_t frame, std::int64_t n_steps) {
    // utils.write_coords_csv
    os << frame << ',' << n_steps << '\n';
    for (const Vec3& c : s.pos) os << py_repr(c[0]) << ',' << py_repr(c[1]) << ',' << py_repr(c[2]) << '\n';
}

void write_xyz_frame(std::ostream& os, const System& s, const AtomLabels& labels, std::int64_t frame,
                     std::int64_t n_steps) { // utils.write_coords_xyz
    os << s.n_particles + static_cast<int>(s.patch_dir.size()) << '\n';
    os << "Trajectory: Frame " << frame << '/' << n_steps << '\n';
    for (int i = 0; i < s.n_particles; ++i) {
        const Vec3& c = s.pos[i];
        os << labels.centre[s.type[i]] << "    " << py_repr(c[0]) << "    " << py_repr(c[1]) << "    "
           << py_repr(c[2]) << '\n';
        for (int k = 0; k < s.n_patches(i); ++k) {
            // Patch site on the particle surface: (x * radius) + coords, in the same order of
            // operations as the Python code so the numbers are identical.
            const Vec3 site = s.patch(i, k) * s.radius[i] + c;
            os << labels.patch[s.type[i]][k] << "    " << py_repr(site[0]) << "    " << py_repr(site[1]) << "    "
               << py_repr(site[2]) << '\n';
        }
    }
}

void write_clusters_frame(std::ostream& os, const std::vector<std::vector<int>>& clusters, std::int64_t frame) {
    // utils.write_clusters: "Frame f", one line of particle indices per cluster, an empty line
    os << "Frame " << frame << '\n';
    for (const auto& c : clusters) {
        for (std::size_t k = 0; k < c.size(); ++k) os << (k ? "," : "") << c[k];
        os << '\n';
    }
    os << '\n';
}

void write_cluster_sizes(std::ostream& os, const std::vector<std::vector<int>>& clusters) {
    // utils.write_clusters_feats with the sizes from obtain_cluster_features: every value followed by a
    // comma ("5,3,2,"); "0," when there are no clusters
    if (clusters.empty()) os << "0,";
    for (const auto& c : clusters) os << c.size() << ',';
    os << '\n';
}

void write_int_line(std::ostream& os, std::int64_t value) { os << value << '\n'; } // utils.write_clusters_int

std::int64_t largest_cluster(const std::vector<std::vector<int>>& clusters) {
    std::size_t largest = 0;
    for (const auto& c : clusters) largest = std::max(largest, c.size());
    return static_cast<std::int64_t>(largest);
}

void write_n_bonds_frame(std::ostream& os, const System& s, std::int64_t frame) {
    // utils.write_n_bonds for every particle of every cluster (run_simulation):
    // "frame,particle,b1,b2,...," with the bonds per patch. FIX 11: current counts (Python wrote the
    // counts of the previous step).
    for (const auto& c : s.clusters)
        for (int i : c) {
            os << frame << ',' << i << ',';
            for (int k = 0; k < s.n_patches(i); ++k) os << s.n_bonds[s.patch_begin[i] + k] << ',';
            os << '\n';
        }
}

void write_transitions(std::ostream& os, std::int64_t frame, const std::vector<Transition>& events) {
    // utils.write_transitions: "frame,a,b,-,c" or "frame,a,-,b,c"
    for (const Transition& t : events) {
        os << frame;
        for (int v : t.before) os << ',' << v;
        os << ",-";
        for (int v : t.after) os << ',' << v;
        os << '\n';
    }
}

// ---------------------------------------------------------------------------
// OutputFiles
// ---------------------------------------------------------------------------
std::ofstream OutputFiles::open(const std::string& name) {
    std::ofstream f(dir_ / name); // overwrites an existing file
    if (!f) throw std::runtime_error("cannot write " + (dir_ / name).string());
    return f;
}

OutputFiles::OutputFiles(const fs::path& dir, const Params& p, const System& s, OutputOptions options)
    : dir_(dir), labels_(build_labels(p)), options_(options), n_steps_(p.sim.n_steps), n_interval_(p.sim.n_interval) {
    {
        std::ofstream f = open("particles.csv");
        write_particle_types(f, s);
    }
    {
        std::ofstream f = open("vmd_init.tcl");
        write_tcl(f, p, labels_);
    }
    trajectory_ = open("trajectory.csv");
    xyz_ = open("visualization.xyz");
    clusters_ = open("clusters.csv");
    cluster_agg_ = open("cluster_agg.csv");
    cluster_n_ = open("cluster_n.csv");
    cluster_max_ = open("cluster_max.csv");
    n_bonds_ = open("n_bonds.csv");
    if (options_.transitions) transitions_ = open("transitions.csv");
}

void OutputFiles::after_step(std::int64_t step, const System& s, const std::vector<std::vector<int>>& previous) {
    // Every step (FIX 7: also step 1, compared with the starting clusters, and steps without clusters).
    if (options_.transitions) write_transitions(transitions_, step, find_transitions(previous, s.clusters, s.n_particles));
    if (options_.clusters_every_step) write_clusters_frame(clusters_, s.clusters, step);

    // Saved frames: step 1 and every n_interval steps.
    if (step == 1 || step % n_interval_ == 0) {
        write_trajectory_frame(trajectory_, s, step, n_steps_);
        write_xyz_frame(xyz_, s, labels_, step, n_steps_);
        if (!options_.clusters_every_step) write_clusters_frame(clusters_, s.clusters, step);
        write_cluster_sizes(cluster_agg_, s.clusters);
        write_int_line(cluster_n_, static_cast<std::int64_t>(s.clusters.size()));
        write_int_line(cluster_max_, largest_cluster(s.clusters));
        write_n_bonds_frame(n_bonds_, s, step);
        // Push everything to disk at each saved frame, so the files are usable during a long run.
        for (std::ofstream* f : {&trajectory_, &xyz_, &clusters_, &cluster_agg_, &cluster_n_, &cluster_max_,
                                 &n_bonds_, &transitions_})
            if (f->is_open()) f->flush();
    }
}
