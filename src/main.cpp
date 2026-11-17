// Entry point (port of run.py and the main loop of simulation.run_simulation).
//
//   patchy_bcd SIM.par PARTS.par INTERS.par [-f FOLDER] [-p PREFIX] [-nt] [-cs] [-s SEED] ...

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "output.hpp"
#include "params.hpp"
#include "rng.hpp"
#include "steps.hpp"
#include "system.hpp"

namespace fs = std::filesystem;

namespace {

struct Args {
    std::string sim_file, parts_file, inters_file;
    std::string folder = "outputs";
    std::string prefix;
    bool transitions = true;    // -nt turns off
    bool clusters_steps = true; // -cs turns off
    bool has_seed = false;
    std::uint64_t seed = 0;
    bool dump_params = false;
    bool init_only = false;
    bool check = false;
};

void usage(std::ostream& os) {
    os << "usage: patchy_bcd params parts inters [-f FOLDER] [-p PREFIX] [-nt] [-cs] [-s SEED]\n"
          "                  [--check] [--dump-params] [--init-only]\n\n"
          "MC simulation of patchy particles\n\n"
          "positional arguments:\n"
          "  params                    Simulation parameters (.par) file\n"
          "  parts                     Particle parameters (.par) file\n"
          "  inters                    Interaction parameters (.par) file\n\n"
          "options:\n"
          "  -h, --help                show this help message and exit\n"
          "  -f, --folder FOLDER       Output folder (Default: outputs)\n"
          "  -p, --prefix PREFIX       Prefix for output files (Default: None)\n"
          "  -nt, --no-transitions     Flag for disabling transitions calculations\n"
          "  -cs, --no-clusters-steps  Flag for disabling storing clusters for every step\n"
          "  -s, --seed SEED           Random seed (Default: random, printed at start)\n"
          "  --check                   At every saved frame, check that no particles overlap and that\n"
          "                            the bond counts are right and within the patch limits\n"
          "  --dump-params             Print the parsed parameters and exit (for validation)\n"
          "  --init-only               Only build the starting configuration: write particles.csv and\n"
          "                            visualization.xyz (frame 0), then exit (for validation)\n";
}

Args parse_args(int argc, char** argv) {
    Args a;
    std::vector<std::string> positional;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        std::string value;
        bool has_inline = false;
        if (arg.rfind("--", 0) == 0) {
            auto eq = arg.find('=');
            if (eq != std::string::npos) { value = arg.substr(eq + 1); arg = arg.substr(0, eq); has_inline = true; }
        }
        auto next = [&](const std::string& name) -> std::string {
            if (has_inline) return value;
            if (i + 1 >= argc) throw std::invalid_argument("argument " + name + ": expected one argument");
            return argv[++i];
        };

        if (arg == "-h" || arg == "--help") { usage(std::cout); std::exit(0); }
        else if (arg == "-f" || arg == "--folder") a.folder = next(arg);
        else if (arg == "-p" || arg == "--prefix") a.prefix = next(arg);
        else if (arg == "-nt" || arg == "--no-transitions") a.transitions = false;
        else if (arg == "-cs" || arg == "--no-clusters-steps") a.clusters_steps = false;
        else if (arg == "-s" || arg == "--seed") {
            std::string v = next(arg);
            std::size_t pos = 0;
            try { a.seed = std::stoull(v, &pos); } catch (...) { pos = 0; }
            if (pos == 0 || pos != v.size() || v[0] == '-')
                throw std::invalid_argument("argument " + arg + ": invalid seed '" + v + "'");
            a.has_seed = true;
        }
        else if (arg == "--check") a.check = true;
        else if (arg == "--dump-params") a.dump_params = true;
        else if (arg == "--init-only") a.init_only = true;
        else if (arg.size() > 1 && arg[0] == '-') throw std::invalid_argument("unrecognized argument: " + arg);
        else positional.push_back(arg);
    }
    if (positional.size() != 3)
        throw std::invalid_argument("expected 3 parameter files (params, parts, inters), got " +
                                    std::to_string(positional.size()));
    a.sim_file = positional[0];
    a.parts_file = positional[1];
    a.inters_file = positional[2];
    return a;
}

// Same as shutil.copy(file, folder): keeps the file name, overwrites.
void copy_into(const std::string& file, const fs::path& folder) {
    fs::path dest = folder / fs::path(file).filename();
    if (fs::exists(dest) && fs::equivalent(file, dest)) return;
    fs::copy_file(file, dest, fs::copy_options::overwrite_existing);
}

// --check: overlaps, and bond counts equal to the bonds and within each patch's limit.
// Returns a description of the problems found ("" if none).
std::string self_check(const System& s, const Params& p) {
    std::string problems;
    const OverlapReport rep = find_overlaps(s);
    if (rep.pairs > 0)
        problems += std::to_string(rep.pairs) + " overlapping pair(s), deepest " + std::to_string(rep.deepest) + "; ";
    if (count_bonds_per_patch(s) != s.n_bonds) problems += "bond counts differ from the bonds; ";
    for (int i = 0; i < s.n_particles; ++i)
        for (int k = 0; k < s.n_patches(i); ++k)
            if (s.n_bonds[s.patch_begin[i] + k] > p.types[s.type[i]].patch_bonds[k])
                problems += "particle " + std::to_string(i) + " patch " + std::to_string(k) + " over its limit; ";
    return problems;
}

// The main loop of run_simulation: initial bonds, then n_steps steps, writing the output files.
void run(System& sys, const Params& params, Rng& rng, const fs::path& out, const Args& args) {
    StepStats total;
    initial_bonds(sys, params, rng, total);
    OutputFiles files(out, params, sys, {args.transitions, args.clusters_steps}); // particles.csv, vmd_init.tcl
    std::cout << "Initial bonds: " << total.associations << ", clusters: " << sys.clusters.size() << '\n'
              << "----------------\n";

    const std::int64_t n_steps = params.sim.n_steps;
    const std::int64_t every_percent = std::max<std::int64_t>(1, n_steps / 100);
    const auto t0 = std::chrono::steady_clock::now();
    std::int64_t checks = 0, failed_checks = 0;
    std::vector<std::vector<int>> previous_clusters;
    for (std::int64_t step = 1; step <= n_steps; ++step) {
        previous_clusters = sys.clusters;
        simulation_step(sys, params, rng, total);
        files.after_step(step, sys, previous_clusters);
        if (args.check && (step == 1 || step % params.sim.n_interval == 0)) {
            ++checks;
            const std::string problems = self_check(sys, params);
            if (!problems.empty()) {
                ++failed_checks;
                std::cout << "\nWARNING: step " << step << ": " << problems << '\n';
            }
        }
        if (step % every_percent == 0)
            std::cout << "\rProgress: " << step * 100 / n_steps << "%" << std::flush;
    }
    const std::chrono::duration<double> dt = std::chrono::steady_clock::now() - t0;
    std::cout << "\n----------------\nDone!\n";

    auto pct = [](std::int64_t a, std::int64_t b) { return b ? 100.0 * a / b : 0.0; };
    std::cout << "Final state: " << sys.clusters.size() << " clusters, largest " << largest_cluster(sys.clusters)
              << " particles\n"
              << "Moves accepted: single " << pct(total.single_accepted, total.single_moves) << "%, cluster "
              << pct(total.cluster_accepted, total.cluster_moves) << "% (" << total.cluster_moves << " cluster moves)\n"
              << "Bonds formed " << total.associations << " (of " << total.association_tests << " tests), broken "
              << total.dissociations << " (of " << total.dissociation_tests << " tests)";
    if (total.dissociations > 0)
        std::cout << "; average bond lifetime " << static_cast<double>(total.bond_steps) / total.dissociations
                  << " steps";
    std::cout << "\nTime per step: " << (n_steps ? 1e6 * dt.count() / n_steps : 0.0) << " microseconds\n";
    if (args.check)
        std::cout << "Check: " << failed_checks << " of " << checks << " saved frames had problems\n";
}

} // namespace

int main(int argc, char** argv) {
    const auto start = std::chrono::steady_clock::now();

    Args args;
    try {
        args = parse_args(argc, argv);
    } catch (const std::exception& e) {
        usage(std::cerr);
        std::cerr << "patchy_bcd: error: " << e.what() << '\n';
        return 2;
    }

    try {
        // Parse before touching the output folder, so a bad input file leaves nothing behind.
        Params params = load_params(args.sim_file, args.parts_file, args.inters_file);
        if (args.dump_params) {
            dump_params(params, std::cout);
            return 0;
        }
        add_probs_from_energies(params); // as run.py; fails before anything is written

        fs::create_directories(args.folder);
        for (const auto& f : {args.sim_file, args.parts_file, args.inters_file}) copy_into(f, args.folder);
        // Output files go into folder/prefix (os.path.join(folder, prefix, name) in simulation.py).
        // Unlike the Python code, the prefix subfolder is created here.
        const fs::path out = output_dir(args.folder, args.prefix);
        fs::create_directories(out);

        if (!args.has_seed) {
            std::random_device rd;
            args.seed = (static_cast<std::uint64_t>(rd()) << 32) ^ rd();
        }
        std::cout << "Seed: " << args.seed << '\n';
        MtRng rng(args.seed);

        std::cout << "Loaded " << params.types.size() << " particle type(s), " << params.total_particles()
                  << " particles, " << params.interactions.size() << " interacting type pair(s)\n";

        System sys = build_system(params, rng);
        std::cout << "Initialized " << sys.n_particles << " particles with " << sys.patch_dir.size()
                  << " patches (volume fraction " << sys.volume_fraction << "); " << sys.init_redraws
                  << " re-placement(s) to remove overlaps\n"
                  << "Neighbour cells: threshold " << sys.threshold << ", " << sys.n_cells_1d
                  << " cells per side of " << sys.cell_size << '\n';

        if (args.init_only) {
            const AtomLabels labels = build_labels(params);
            std::ofstream types(out / "particles.csv");
            write_particle_types(types, sys);
            std::ofstream xyz(out / "visualization.xyz");
            write_xyz_frame(xyz, sys, labels, 0, params.sim.n_steps);
        } else {
            run(sys, params, rng, out, args);
        }
    } catch (const std::exception& e) {
        std::cerr << "patchy_bcd: error: " << e.what() << '\n';
        return 1;
    }

    // As run.py: print(f"Finished in {round(minutes, 1)} minutes")
    const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start;
    std::printf("Finished in %.1f minutes\n", elapsed.count() / 60.0);
    return 0;
}
