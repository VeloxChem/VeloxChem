#ifndef PrecisionCut_hpp
#define PrecisionCut_hpp

#include <vector>     // std::vector
#include <cstdint>    // uint32_t

// Compute the maximum absolute value within each tile (per-tile |max|)
std::vector<double>
tile_max_abs(const std::vector<double>& arr, int tile);

// For a descending-sorted array:
// return the index of the first element <= thr
// (equivalent to counting how many elements satisfy desc[i] > thr)
uint32_t
cut_descending(const std::vector<double>& desc, double thr);

// scheme A: per ij-tile cut
std::vector<uint32_t>
build_cut_ij_tile(
    const std::vector<double>& Q_ij_local,
    const std::vector<double>& Q_kl,
    const std::vector<double>& D_kl,
    uint32_t ij_count_local,
    uint32_t kl_count,
    int tile,
    double tau);

std::vector<uint32_t>
build_cut_ij_tile_dd(
    const std::vector<double>& Q_ij_local,
    const std::vector<double>& Q_kl,
    const std::vector<double>& D_kl,
    uint32_t ij_count_local,
    uint32_t kl_count,
    int ij_tile_dim,
    int kl_tile_dim,
    double tau);

#endif
