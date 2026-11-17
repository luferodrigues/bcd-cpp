#include "params.hpp"

#include <algorithm>
#include <cctype>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <type_traits>

// ---------------------------------------------------------------------------
// Generic .par reader
//
// All three files share one layout:
//   key                 <- line without a comma
//   v1,v2,v3,           <- one or more value rows (trailing comma optional)
//   <blank line>        <- separates particle / interaction blocks
// After a key, the next line is always its first value row; further rows are
// consumed while lines contain a comma (this is how the Python parsers detect
// the end of multi-line parameters such as patches_positions or interact).
// ---------------------------------------------------------------------------
namespace {

struct Row {
    int line = 0;
    std::vector<std::string> tokens;
};

struct Entry {
    std::string key;
    int line = 0;
    int block = 0; // index of the blank-line-separated block this entry is in
    std::vector<Row> rows;
};

[[noreturn]] void fail(const std::string& file, int line, const std::string& msg) {
    std::ostringstream os;
    os << file;
    if (line > 0) os << ':' << line;
    os << ": " << msg;
    throw std::runtime_error(os.str());
}

void warn(const std::string& file, int line, const std::string& msg) {
    std::cerr << "warning: " << file;
    if (line > 0) std::cerr << ':' << line;
    std::cerr << ": " << msg << '\n';
}

std::string trim(const std::string& s) {
    std::size_t b = 0, e = s.size();
    while (b < e && std::isspace(static_cast<unsigned char>(s[b]))) ++b;
    while (e > b && std::isspace(static_cast<unsigned char>(s[e - 1]))) --e;
    return s.substr(b, e - b);
}

std::vector<std::string> split_row(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    std::istringstream is(line);
    while (std::getline(is, cur, ',')) out.push_back(trim(cur));
    // getline drops a final empty field after a trailing comma; also drop an
    // explicit empty last field so "1,2," and "1,2" both give {1,2}.
    while (!out.empty() && out.back().empty()) out.pop_back();
    return out;
}

std::vector<Entry> read_entries(const std::string& path) {
    std::ifstream in(path);
    if (!in) fail(path, 0, "cannot open file");

    std::vector<std::string> lines;
    for (std::string l; std::getline(in, l);) lines.push_back(trim(l));

    std::vector<Entry> entries;
    int block = 0;
    bool block_has_entries = false;
    std::size_t i = 0;
    while (i < lines.size()) {
        const std::string& l = lines[i];
        if (l.empty()) {
            if (block_has_entries) { ++block; block_has_entries = false; }
            ++i;
            continue;
        }
        if (l.find(',') != std::string::npos)
            fail(path, static_cast<int>(i) + 1, "expected a parameter name, found values: '" + l + "'");

        Entry e;
        e.key = l;
        e.line = static_cast<int>(i) + 1;
        e.block = block;
        ++i;
        if (i >= lines.size() || lines[i].empty())
            fail(path, e.line, "parameter '" + e.key + "' has no value");
        e.rows.push_back({static_cast<int>(i) + 1, split_row(lines[i])});
        ++i;
        while (i < lines.size() && !lines[i].empty() && lines[i].find(',') != std::string::npos) {
            e.rows.push_back({static_cast<int>(i) + 1, split_row(lines[i])});
            ++i;
        }
        entries.push_back(std::move(e));
        block_has_entries = true;
    }
    return entries;
}

double to_double(const std::string& file, int line, const std::string& tok) {
    if (tok.empty()) fail(file, line, "empty value");
    errno = 0;
    char* end = nullptr;
    double v = std::strtod(tok.c_str(), &end);
    if (end != tok.c_str() + tok.size() || errno == ERANGE)
        fail(file, line, "not a number: '" + tok + "'");
    return v;
}

// Python uses int(float(x)) or int(x); either way the value must be integral here.
std::int64_t to_int(const std::string& file, int line, const std::string& tok) {
    double v = to_double(file, line, tok);
    if (v != std::floor(v) || std::fabs(v) > 9.0e15)
        fail(file, line, "expected an integer, found '" + tok + "'");
    return static_cast<std::int64_t>(v);
}

const std::string& single_token(const std::string& file, const Entry& e) {
    if (e.rows.size() != 1 || e.rows[0].tokens.empty())
        fail(file, e.line, "parameter '" + e.key + "' expects a single value");
    if (e.rows[0].tokens.size() > 1)
        warn(file, e.rows[0].line, "parameter '" + e.key + "': only the first value is used");
    return e.rows[0].tokens[0];
}

std::vector<double> double_list(const std::string& file, const Entry& e) {
    if (e.rows.size() != 1) fail(file, e.line, "parameter '" + e.key + "' expects a single row");
    std::vector<double> v;
    for (const auto& t : e.rows[0].tokens) v.push_back(to_double(file, e.rows[0].line, t));
    return v;
}

template <typename T>
Matrix<T> matrix(const std::string& file, const Entry& e) {
    Matrix<T> m;
    m.rows = static_cast<int>(e.rows.size());
    m.cols = static_cast<int>(e.rows[0].tokens.size());
    for (const Row& r : e.rows) {
        if (static_cast<int>(r.tokens.size()) != m.cols)
            fail(file, r.line, "parameter '" + e.key + "': rows have different lengths");
        for (const auto& t : r.tokens) {
            if constexpr (std::is_integral_v<T>)
                m.data.push_back(static_cast<T>(to_int(file, r.line, t)));
            else
                m.data.push_back(to_double(file, r.line, t));
        }
    }
    return m;
}

} // namespace

// ---------------------------------------------------------------------------
// simulation_params.par
// ---------------------------------------------------------------------------
SimParams load_sim_params(const std::string& path) {
    SimParams s;
    std::set<std::string> seen;
    for (const Entry& e : read_entries(path)) {
        if (!seen.insert(e.key).second) fail(path, e.line, "duplicate parameter '" + e.key + "'");
        const std::string& tok = single_token(path, e);
        const int ln = e.rows[0].line;
        if (e.key == "box_limits") s.box_limits = to_double(path, ln, tok);
        else if (e.key == "temperature") s.temperature = to_double(path, ln, tok);
        else if (e.key == "viscosity") s.viscosity = to_double(path, ln, tok);
        else if (e.key == "time_step") {
            s.time_step = to_double(path, ln, tok);
            if (s.time_step != std::floor(s.time_step))
                warn(path, ln, "time_step is not an integer; the Python code truncates it with int(), "
                               "this port keeps the full value");
        }
        else if (e.key == "n_steps") s.n_steps = to_int(path, ln, tok);
        else if (e.key == "n_interval") s.n_interval = to_int(path, ln, tok);
        else if (e.key == "n_dimensions") s.n_dimensions = static_cast<int>(to_int(path, ln, tok));
        else warn(path, e.line, "unknown parameter '" + e.key + "' ignored");
    }
    for (std::string_view k : {"box_limits", "temperature", "viscosity", "time_step", "n_steps", "n_interval"})
        if (!seen.count(std::string(k))) fail(path, 0, "missing parameter '" + std::string(k) + "'");

    if (s.box_limits <= 0) fail(path, 0, "box_limits must be positive");
    if (s.n_steps < 0) fail(path, 0, "n_steps must be non-negative");
    if (s.n_interval <= 0) fail(path, 0, "n_interval must be positive");
    if (s.n_dimensions != 3) fail(path, 0, "only n_dimensions = 3 is supported");
    return s;
}

// ---------------------------------------------------------------------------
// particle_params.par
// ---------------------------------------------------------------------------
std::vector<ParticleType> load_particle_params(const std::string& path) {
    std::vector<Entry> entries = read_entries(path);
    std::vector<ParticleType> types;
    std::vector<std::set<std::string>> seen;
    std::vector<int> block_line;

    for (const Entry& e : entries) {
        if (e.block >= static_cast<int>(types.size())) {
            types.emplace_back();
            seen.emplace_back();
            block_line.push_back(e.line);
        }
        ParticleType& t = types[e.block];
        if (!seen[e.block].insert(e.key).second) fail(path, e.line, "duplicate parameter '" + e.key + "'");
        const int ln = e.rows[0].line;

        if (e.key == "type") t.type = static_cast<int>(to_int(path, ln, single_token(path, e)));
        else if (e.key == "name") t.name = single_token(path, e);
        else if (e.key == "number") t.number = static_cast<int>(to_int(path, ln, single_token(path, e)));
        else if (e.key == "radius") t.radius = to_double(path, ln, single_token(path, e));
        else if (e.key == "diff_trans" || e.key == "diff_rot") {
            double d = to_double(path, ln, single_token(path, e));
            if (d < 0) fail(path, ln, e.key + " must be non-negative");
            (e.key == "diff_trans" ? t.diff_trans : t.diff_rot) = d;
        }
        else if (e.key == "interacting") {
            std::string v = single_token(path, e);
            std::transform(v.begin(), v.end(), v.begin(), [](unsigned char c) { return std::tolower(c); });
            if (v == "yes" || v == "true") t.interacting = true;
            else if (v == "no" || v == "false") t.interacting = false;
            else fail(path, ln, "interacting must be yes/no/true/false, found '" + v + "'");
        }
        else if (e.key == "patches_radius") t.patch_radius = double_list(path, e);
        else if (e.key == "patches_alphas") t.patch_alpha_deg = double_list(path, e);
        else if (e.key == "patches_bonds") {
            for (double b : double_list(path, e)) {
                if (b != std::floor(b) || b < 0) fail(path, ln, "patches_bonds must be non-negative integers");
                t.patch_bonds.push_back(static_cast<int>(b));
            }
        }
        else if (e.key == "patches_positions") {
            for (const Row& r : e.rows) {
                if (r.tokens.size() != 3) fail(path, r.line, "patch positions need 3 coordinates");
                t.patch_positions.push_back({to_double(path, r.line, r.tokens[0]),
                                             to_double(path, r.line, r.tokens[1]),
                                             to_double(path, r.line, r.tokens[2])});
            }
        }
        else warn(path, e.line, "unknown parameter '" + e.key + "' ignored");
    }

    if (types.empty()) fail(path, 0, "no particle types defined");
    for (std::size_t k = 0; k < types.size(); ++k) {
        ParticleType& t = types[k];
        const int ln = block_line[k];
        for (std::string_view key : {"type", "number", "radius", "interacting", "patches_radius",
                                "patches_positions", "patches_alphas", "patches_bonds"})
            if (!seen[k].count(std::string(key)))
                fail(path, ln, "particle block is missing '" + std::string(key) + "'");

        // The Python code indexes types by block order, so the declared type must match it.
        if (t.type != static_cast<int>(k))
            fail(path, ln, "type " + std::to_string(t.type) + " declared in block " + std::to_string(k) +
                               "; types must be numbered 0,1,2,... in file order");
        if (t.number < 0) fail(path, ln, "number must be non-negative");
        if (t.radius <= 0) fail(path, ln, "radius must be positive");

        const std::size_t np = t.patch_positions.size();
        if (np == 0) fail(path, ln, "a particle type needs at least one patch");
        if (t.patch_radius.size() != np) fail(path, ln, "patches_radius count differs from patches_positions");
        if (t.patch_alpha_deg.size() != np) fail(path, ln, "patches_alphas count differs from patches_positions");
        if (t.patch_bonds.size() == 1 && np > 1) t.patch_bonds.assign(np, t.patch_bonds[0]);
        if (t.patch_bonds.size() != np) fail(path, ln, "patches_bonds count differs from patches_positions");
        for (const Vec3& n : t.patch_positions)
            if (std::sqrt(n[0] * n[0] + n[1] * n[1] + n[2] * n[2]) == 0.0)
                fail(path, ln, "patch position vector of length zero");
    }
    return types;
}

// ---------------------------------------------------------------------------
// interaction_params.par
// ---------------------------------------------------------------------------
std::vector<InteractionPair> load_interaction_params(const std::string& path) {
    std::vector<InteractionPair> pairs;
    std::vector<std::set<std::string>> seen;
    std::vector<int> pair_line;
    for (const Entry& e : read_entries(path)) {
        if (e.key == "type") {
            const Row& r = e.rows[0];
            if (e.rows.size() != 1 || r.tokens.size() != 2) fail(path, r.line, "type needs two values, e.g. '0,1,'");
            InteractionPair p;
            p.type_a = static_cast<int>(to_int(path, r.line, r.tokens[0]));
            p.type_b = static_cast<int>(to_int(path, r.line, r.tokens[1]));
            if (p.type_a > p.type_b)
                fail(path, r.line, "write the lower type first (the Python code looks pairs up as (min, max))");
            for (const auto& q : pairs)
                if (q.type_a == p.type_a && q.type_b == p.type_b) fail(path, r.line, "duplicate type pair");
            pairs.push_back(std::move(p));
            seen.emplace_back();
            pair_line.push_back(e.line);
            continue;
        }
        if (pairs.empty()) fail(path, e.line, "'" + e.key + "' appears before any 'type' line");
        InteractionPair& p = pairs.back();
        if (!seen.back().insert(e.key).second) fail(path, e.line, "duplicate parameter '" + e.key + "'");
        if (e.key == "interact") p.interact = matrix<int>(path, e);
        else if (e.key == "distances") p.distances = matrix<double>(path, e);
        else if (e.key == "energies") p.energies = matrix<double>(path, e);
        else if (e.key == "p_ass") p.p_ass = matrix<double>(path, e);
        else if (e.key == "p_diss") p.p_diss = matrix<double>(path, e);
        else warn(path, e.line, "unknown parameter '" + e.key + "' ignored");
    }
    for (std::size_t k = 0; k < pairs.size(); ++k) {
        if (pairs[k].interact.empty()) fail(path, pair_line[k], "type pair has no 'interact' matrix");
        for (int v : pairs[k].interact.data)
            if (v != 0 && v != 1) fail(path, pair_line[k], "interact entries must be 0 or 1");
    }
    return pairs;
}

// ---------------------------------------------------------------------------
// Cross-file checks
// ---------------------------------------------------------------------------
int Params::total_particles() const {
    int n = 0;
    for (const auto& t : types) n += t.number;
    return n;
}

const InteractionPair* Params::find_interaction(int a, int b) const {
    if (a > b) std::swap(a, b);
    for (const auto& p : interactions)
        if (p.type_a == a && p.type_b == b) return &p;
    return nullptr;
}

Params load_params(const std::string& sim_path, const std::string& parts_path,
                   const std::string& inters_path) {
    Params p;
    p.sim = load_sim_params(sim_path);
    p.types = load_particle_params(parts_path);
    p.interactions = load_interaction_params(inters_path);

    const int ntypes = static_cast<int>(p.types.size());
    for (const auto& ip : p.interactions) {
        const std::string pair = "(" + std::to_string(ip.type_a) + "," + std::to_string(ip.type_b) + ")";
        if (ip.type_a < 0 || ip.type_b >= ntypes)
            fail(inters_path, 0, "pair " + pair + " refers to an undefined particle type");
        const int ra = p.types[ip.type_a].n_patches();
        const int cb = p.types[ip.type_b].n_patches();
        auto check = [&](std::string_view name, int rows, int cols) {
            if (rows == 0) return; // absent
            if (rows != ra || cols != cb)
                fail(inters_path, 0, std::string(name) + " for pair " + pair + " is " + std::to_string(rows) + "x" +
                                         std::to_string(cols) + ", expected " + std::to_string(ra) + "x" +
                                         std::to_string(cb) + " (patches of type a x patches of type b)");
        };
        check("interact", ip.interact.rows, ip.interact.cols);
        check("distances", ip.distances.rows, ip.distances.cols);
        check("energies", ip.energies.rows, ip.energies.cols);
        check("p_ass", ip.p_ass.rows, ip.p_ass.cols);
        check("p_diss", ip.p_diss.rows, ip.p_diss.cols);

        if (ip.type_a == ip.type_b)
            for (int i = 0; i < ra; ++i)
                for (int j = 0; j < i; ++j)
                    if (ip.interact(i, j) != ip.interact(j, i))
                        warn(inters_path, 0, "interact matrix for same-type pair " + pair + " is not symmetric");
        if (!p.types[ip.type_a].interacting || !p.types[ip.type_b].interacting)
            warn(inters_path, 0, "pair " + pair + " involves a type marked 'interacting: no'");
    }
    return p;
}

void add_probs_from_energies(Params& p) {
    constexpr double tol = 1e-8;
    auto prob = [](double delta_energy) { return std::exp(-delta_energy); }; // calc.calculate_probability_from_energy
    for (InteractionPair& ip : p.interactions) {
        if (!ip.p_ass.empty() && !ip.p_diss.empty()) continue;
        if (ip.energies.empty())
            throw std::runtime_error("interaction pair (" + std::to_string(ip.type_a) + "," + std::to_string(ip.type_b) +
                                     ") has no p_ass/p_diss and no energies to compute them from");
        const Matrix<double>& E = ip.energies;
        if (ip.p_ass.empty()) {
            ip.p_ass = E;
            for (std::size_t k = 0; k < E.data.size(); ++k) {
                double v = prob(E.data[k]);
                if (v < tol || E.data[k] == 0) v = 0;
                else if (v > 1.0 - tol) v = 1;
                ip.p_ass.data[k] = v;
            }
        }
        if (ip.p_diss.empty()) {
            ip.p_diss = E;
            for (std::size_t k = 0; k < E.data.size(); ++k) {
                double v = prob(-E.data[k]);
                if (v < tol) v = 0;
                else if (v > 1.0 - tol) v = 1;
                ip.p_diss.data[k] = v;
            }
        }
    }
}

// ---------------------------------------------------------------------------
// Canonical dump (must stay in sync with tools/dump_params.py)
// ---------------------------------------------------------------------------
namespace {
std::string num(double v) {
    char buf[64];
    std::snprintf(buf, sizeof buf, "%.17g", v);
    return buf;
}

template <typename T>
void dump_matrix(std::ostream& os, std::string_view name, const Matrix<T>& m) {
    if (m.empty()) { os << name << " absent\n"; return; }
    os << name << '\n';
    for (int i = 0; i < m.rows; ++i) {
        for (int j = 0; j < m.cols; ++j) os << (j ? " " : "") << num(static_cast<double>(m(i, j)));
        os << '\n';
    }
}
} // namespace

void dump_params(const Params& p, std::ostream& os) {
    const SimParams& s = p.sim;
    os << "[simulation]\n"
       << "box_limits " << num(s.box_limits) << '\n'
       << "temperature " << num(s.temperature) << '\n'
       << "viscosity " << num(s.viscosity) << '\n'
       << "time_step " << num(s.time_step) << '\n'
       << "n_steps " << num(static_cast<double>(s.n_steps)) << '\n'
       << "n_interval " << num(static_cast<double>(s.n_interval)) << '\n'
       << "n_dimensions " << s.n_dimensions << '\n';
    for (const auto& t : p.types) {
        os << "[type " << t.type << "]\n"
           << "name " << t.name << '\n'
           << "number " << t.number << '\n'
           << "radius " << num(t.radius) << '\n'
           << "interacting " << (t.interacting ? "true" : "false") << '\n'
           << "diff_trans " << (t.diff_trans ? num(*t.diff_trans) : "none") << '\n'
           << "diff_rot " << (t.diff_rot ? num(*t.diff_rot) : "none") << '\n';
        for (int k = 0; k < t.n_patches(); ++k) {
            const Vec3& n = t.patch_positions[k];
            os << "patch " << k << " radius " << num(t.patch_radius[k]) << " alpha " << num(t.patch_alpha_deg[k])
               << " bonds " << t.patch_bonds[k] << " position " << num(n[0]) << ' ' << num(n[1]) << ' '
               << num(n[2]) << '\n';
        }
    }
    for (const auto& ip : p.interactions) {
        os << "[interaction " << ip.type_a << ' ' << ip.type_b << "]\n";
        dump_matrix(os, "interact", ip.interact);
        dump_matrix(os, "distances", ip.distances);
        dump_matrix(os, "energies", ip.energies);
        dump_matrix(os, "p_ass", ip.p_ass);
        dump_matrix(os, "p_diss", ip.p_diss);
    }
}
