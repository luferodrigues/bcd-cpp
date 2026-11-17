#include "cells.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace {
// Python's % always gives a result in [0, n), C++'s % can be negative: (-1) % 12 is -1 in C++
// but 11 in Python. This helper behaves like Python.
int py_mod(int a, int n) {
    const int r = a % n;
    return r < 0 ? r + n : r;
}
} // namespace

CellList::CellList(double box_limits, double cell_size, int n_cells_1d)
    : box_limits_(box_limits), cell_size_(cell_size), n_(n_cells_1d) {
    if (n_ < 1) throw std::invalid_argument("CellList needs at least one cell per side");
    const int n_cells = n_ * n_ * n_;
    head_.assign(n_cells, -1);

    // The neighbour cells of every cell, in the order the Python loops visit them (dx, dy, dz from
    // -1 to 1). With fewer than 3 cells per side the wrap reaches the same cell more than once; each
    // cell is kept only once.
    neighbour_begin_.push_back(0);
    for (int c = 0; c < n_cells; ++c) {
        const int ix = c % n_, iy = (c / n_) % n_, iz = c / (n_ * n_);
        const std::size_t first = neighbour_cells_.size();
        for (int dx = -1; dx <= 1; ++dx)
            for (int dy = -1; dy <= 1; ++dy)
                for (int dz = -1; dz <= 1; ++dz) {
                    const int nc = py_mod(ix + dx, n_) + n_ * (py_mod(iy + dy, n_) + n_ * py_mod(iz + dz, n_));
                    if (std::find(neighbour_cells_.begin() + first, neighbour_cells_.end(), nc) == neighbour_cells_.end())
                        neighbour_cells_.push_back(nc);
                }
        neighbour_begin_.push_back(static_cast<int>(neighbour_cells_.size()));
    }
}

int CellList::cell_index(const Vec3& r) const {
    int idx[3];
    for (int k = 0; k < 3; ++k) {
        // int() in Python and static_cast<int> both round towards zero.
        const int c = static_cast<int>((r[k] + box_limits_) / cell_size_);
        idx[k] = std::max(0, std::min(c, n_ - 1)); // clip, as in update_cell_list
    }
    return idx[0] + n_ * (idx[1] + n_ * idx[2]);
}

void CellList::build(const std::vector<Vec3>& pos) {
    std::fill(head_.begin(), head_.end(), -1);
    next_.assign(pos.size(), -1);
    cell_.assign(pos.size(), -1);
    for (int i = 0; i < static_cast<int>(pos.size()); ++i) {
        const int c = cell_index(pos[i]);
        next_[i] = head_[c]; // put i at the front of the chain of cell c
        head_[c] = i;
        cell_[i] = c;
    }
}

void CellList::move(int i, const Vec3& new_pos) {
    const int c_new = cell_index(new_pos);
    const int c_old = cell_[i];
    if (c_new == c_old) return;
    // Unlink i from the chain of its old cell...
    if (head_[c_old] == i) {
        head_[c_old] = next_[i];
    } else {
        int j = head_[c_old];
        while (next_[j] != i) j = next_[j];
        next_[j] = next_[i];
    }
    // ...and put it at the front of the new one.
    next_[i] = head_[c_new];
    head_[c_new] = i;
    cell_[i] = c_new;
}

void CellList::neighbours_above(int i, std::vector<int>& out) const {
    out.clear();
    const int c = cell_[i];
    for (int k = neighbour_begin_[c]; k < neighbour_begin_[c + 1]; ++k)
        for (int j = head_[neighbour_cells_[k]]; j != -1; j = next_[j])
            if (j > i) out.push_back(j);
    std::sort(out.begin(), out.end()); // each cell is listed once, so there are no duplicates
}

bool CellList::clash(int i, const Vec3& new_pos, const std::vector<Vec3>& pos, const std::vector<double>& radius,
                     const std::vector<int>& skip) const {
    // Cell of the new position, not clipped, so a position just outside the box still finds the
    // right neighbours through the periodic wrap below.
    // FIX 1B: std::floor, not int() as in Python. int() rounds towards zero, so for a position below
    // -box_limits it gave a cell index one too high and one layer of neighbour cells was not searched.
    const int ix = static_cast<int>(std::floor((new_pos[0] + box_limits_) / cell_size_));
    const int iy = static_cast<int>(std::floor((new_pos[1] + box_limits_) / cell_size_));
    const int iz = static_cast<int>(std::floor((new_pos[2] + box_limits_) / cell_size_));
    const int c = py_mod(ix, n_) + n_ * (py_mod(iy, n_) + n_ * py_mod(iz, n_));
    for (int k = neighbour_begin_[c]; k < neighbour_begin_[c + 1]; ++k)
        for (int j = head_[neighbour_cells_[k]]; j != -1; j = next_[j]) {
            if (std::find(skip.begin(), skip.end(), j) != skip.end()) continue; // j in skip
            const double limit = radius[i] + radius[j];
            if (distance_sq_pbc(box_limits_, new_pos, pos[j]) < limit * limit) return true;
        }
    return false;
}
