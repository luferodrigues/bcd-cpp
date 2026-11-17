// Small vector/matrix helpers (part of calculate.py). Header-only: these are tiny functions the
// compiler should inline into the hot loops.
#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <vector>

// A 3D vector: v[0], v[1], v[2] are x, y, z (like a NumPy array of shape (3,)).
using Vec3 = std::array<double, 3>;

// Arithmetic on Vec3, so `a + b`, `a - b`, `2.0 * a`, `a / s` work as they do on NumPy arrays.
inline Vec3 operator+(const Vec3& a, const Vec3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
inline Vec3 operator-(const Vec3& a, const Vec3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
inline Vec3 operator*(double s, const Vec3& a) { return {s * a[0], s * a[1], s * a[2]}; }
inline Vec3 operator*(const Vec3& a, double s) { return {a[0] * s, a[1] * s, a[2] * s}; }
inline Vec3 operator/(const Vec3& a, double s) { return {a[0] / s, a[1] / s, a[2] / s}; }

inline double dot(const Vec3& a, const Vec3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
inline double norm(const Vec3& a) { return std::sqrt(dot(a, a)); }

// 3x3 matrix, row-major: m[3*i + j] is row i, column j.
struct Mat3 {
    std::array<double, 9> m{};

    Vec3 operator*(const Vec3& v) const {
        return {m[0] * v[0] + m[1] * v[1] + m[2] * v[2],
                m[3] * v[0] + m[4] * v[1] + m[5] * v[2],
                m[6] * v[0] + m[7] * v[1] + m[8] * v[2]};
    }
};

// Rotation matrix for a rotation vector (axis times angle in radians); the same rotation as
// scipy.spatial.transform.Rotation.from_rotvec(v). Rodrigues' formula: R = I + A*K + B*K^2,
// where K is the cross-product matrix of v, A = sin(t)/t, B = (1 - cos(t))/t^2 and t = |v|.
inline Mat3 rotation_from_rotvec(const Vec3& v) {
    const double t2 = dot(v, v);
    const double t = std::sqrt(t2);
    double A, B;
    if (t < 1e-6) { // Taylor series, avoids 0/0 for tiny angles
        A = 1.0 - t2 / 6.0;
        B = 0.5 - t2 / 24.0;
    } else {
        A = std::sin(t) / t;
        B = (1.0 - std::cos(t)) / t2;
    }
    const double x = v[0], y = v[1], z = v[2];
    Mat3 R;
    R.m = {1.0 - B * (y * y + z * z), -A * z + B * x * y,        A * y + B * x * z,
           A * z + B * x * y,         1.0 - B * (x * x + z * z), -A * x + B * y * z,
           -A * y + B * x * z,        A * x + B * y * z,         1.0 - B * (x * x + y * y)};
    return R;
}

// Rotation matrix of a unit quaternion (w, x, y, z).
inline Mat3 rotation_from_quaternion(double w, double x, double y, double z) {
    Mat3 R;
    R.m = {1 - 2 * (y * y + z * z), 2 * (x * y - w * z),     2 * (x * z + w * y),
           2 * (x * y + w * z),     1 - 2 * (x * x + z * z), 2 * (y * z - w * x),
           2 * (x * z - w * y),     2 * (y * z + w * x),     1 - 2 * (x * x + y * y)};
    return R;
}

// ---------------------------------------------------------------------------
// Periodic box [-box_limits, box_limits]^3 (calculate.py). Written exactly like the Python
// functions, so results agree to the last digit.
// ---------------------------------------------------------------------------

// FIX 1A: shift a position back into the box (what steps.apply_periodic_boundaries does for every
// particle at the end of the move phase), used right after each accepted move.
inline Vec3 wrap_into_box(Vec3 r, double box_limits) {
    for (double& c : r)
        if (std::abs(c) > box_limits) c = c > 0 ? c - 2 * box_limits : c + 2 * box_limits;
    return r;
}

// calc.calculate_distance_sq: squared distance; per coordinate the smallest of the direct
// difference and the differences to the two periodic images of a.
inline double distance_sq_pbc(double box_limits, const Vec3& a, const Vec3& b) {
    double dist_squared = 0.0;
    for (int k = 0; k < 3; ++k) {
        const double d1 = b[k] - a[k];
        const double d2 = b[k] - (a[k] + 2 * box_limits);
        const double d3 = b[k] - (a[k] - 2 * box_limits);
        dist_squared += std::min({d1 * d1, d2 * d2, d3 * d3});
    }
    return dist_squared;
}

// calc.calculate_distance_vector_pbc: vector from a to b (b - a), shifted by one box length
// where a component is longer than half the box.
inline Vec3 distance_vector_pbc(double box_limits, const Vec3& a, const Vec3& b) {
    const double l = 2.0 * box_limits;
    Vec3 d = b - a;
    for (double& c : d) {
        if (c > box_limits) c -= l;
        else if (c < -box_limits) c += l;
    }
    return d;
}

// calc.center_of_mass: plain average of the coordinates. NOT periodic-aware: for a cluster that
// crosses the box edge this is not its real centre (see CLAUDE.md, phase 4).
inline Vec3 center_of_mass(const std::vector<Vec3>& coords) {
    Vec3 sum{0.0, 0.0, 0.0};
    for (const Vec3& r : coords)
        for (int k = 0; k < 3; ++k) sum[k] += r[k];
    return sum / static_cast<double>(coords.size());
}

// ---------------------------------------------------------------------------
// Kern-Frenkel orientation check (utils.check_alignment / check_alignment_all)
// ---------------------------------------------------------------------------

// utils.check_alignment: patch 1 (unit vector n1, half-angle alpha1_deg) and patch 2 (n2,
// alpha2_deg) face each other along the unit vector d_norm if
//     n1 . (-d_norm) >= cos(alpha1)   and   n2 . (+d_norm) >= cos(alpha2).
// For "patch 1 points at particle 2 and patch 2 points at particle 1", d must therefore point
// from particle 2 to particle 1 (d = r1 - r2). Which way steps.py computes d: check in phase 5.
inline bool patches_aligned(const Vec3& n1, const Vec3& n2, const Vec3& d_norm, double alpha1_deg,
                            double alpha2_deg) {
    constexpr double kDegToRad = 3.14159265358979323846 / 180.0; // np.pi / 180.0
    const Vec3 minus_d = -1.0 * d_norm;
    return dot(n1, minus_d) >= std::cos(alpha1_deg * kDegToRad) &&
           dot(n2, d_norm) >= std::cos(alpha2_deg * kDegToRad);
}

struct PatchPair {
    int patch1 = -1; // index among the patches of particle 1, or -1 if none is aligned
    int patch2 = -1;
};

// utils.check_alignment_all: the FIRST aligned pair, looping over the patches of particle 1
// (outer loop) and particle 2 (inner loop). d_vector does not need to be normalized.
// Each list is passed as a "slice": a pointer to its first element plus a count (like the
// char* text earlier, but with the length given explicitly). patches1/alphas1 are the n1 patch
// vectors/half-angles of particle 1; patches2/alphas2 the n2 of particle 2.
inline PatchPair first_aligned_pair(const Vec3* patches1, const double* alphas1, int n1, const Vec3* patches2,
                                    const double* alphas2, int n2, const Vec3& d_vector) {
    const Vec3 d_norm = d_vector / norm(d_vector);
    for (int i = 0; i < n1; ++i)
        for (int j = 0; j < n2; ++j)
            if (patches_aligned(patches1[i], patches2[j], d_norm, alphas1[i], alphas2[j])) return {i, j};
    return {};
}
