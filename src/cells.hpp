// Cell list for neighbour search (utils.update_cell_list / check_steric_clash_cell).
//
// The box is divided into n x n x n cells of side cell_size. Each cell keeps a chain of the
// particles inside it, exactly like the Python `head` / `linked_list` arrays:
//   head[c]  = one particle in cell c (-1 if empty)
//   next[i]  = the next particle in the same cell as particle i (-1 at the end of the chain)
// To visit cell c: j = head[c]; while (j != -1) { ...; j = next[j]; }
#pragma once

#include <vector>

#include "geometry.hpp"

class CellList {
public:
    CellList() = default;
    CellList(double box_limits, double cell_size, int n_cells_1d);

    // utils.update_cell_list: put every particle in its cell (same order as Python, so `head`
    // and `next` come out identical).
    void build(const std::vector<Vec3>& pos);

    // Move particle i to the cell of new_pos (call after changing its position).
    void move(int i, const Vec3& new_pos);

    // utils.check_steric_clash_cell: would particle i at new_pos overlap a particle in its own or
    // the 26 neighbouring cells? Particles listed in `skip` are ignored (normally the cluster
    // being moved, which must include i itself).
    bool clash(int i, const Vec3& new_pos, const std::vector<Vec3>& pos, const std::vector<double>& radius,
               const std::vector<int>& skip) const;

    // All particles j > i in particle i's cell and the 26 around it, sorted ascending.
    // Written into `out` (cleared first) so the vector's memory can be reused between calls.
    void neighbours_above(int i, std::vector<int>& out) const;

    int cell_index(const Vec3& r) const; // with the clipping of update_cell_list

    const std::vector<int>& head() const { return head_; }
    const std::vector<int>& next() const { return next_; }

private:
    double box_limits_ = 0.0;
    double cell_size_ = 0.0;
    int n_ = 0;               // cells per side
    std::vector<int> head_;   // n^3 entries
    std::vector<int> next_;   // one per particle
    std::vector<int> cell_;   // current cell of each particle (needed by move)

    // Neighbour cells, worked out once (speed-up, phase 7): the cells around cell c (itself and the
    // 26 around it, each listed once, with the periodic wrap) are
    //   neighbour_cells_[neighbour_begin_[c]] ... neighbour_cells_[neighbour_begin_[c + 1] - 1]
    std::vector<int> neighbour_begin_;
    std::vector<int> neighbour_cells_;
};
