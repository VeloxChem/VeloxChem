#ifndef PrecisionCut_hpp
#define PrecisionCut_hpp

#include <vector>     // std::vector
#include <cstdint>    // uint32_t

#include "GpuRuntime.hpp"
#include "GpuWrapper.hpp"

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

struct ExchangeCuts {
    std::vector<uint32_t> prec_cut_flat;    // optional host-side flat cut array
    std::vector<uint32_t> screen_cut_flat;  // optional host-side flat cut array
    std::vector<uint32_t> displ_cuts;       // [n_ik] offset into flat arrays
    std::vector<uint32_t> cut_weights;      // optional host-only m-tile weight per cut entry
    uint32_t total_cut_entries = 0;         // number of flat cut entries
};

ExchangeCuts
build_exchange_cuts(
    const std::vector<uint32_t>& pair_inds_i,
    const std::vector<uint32_t>& pair_inds_k,
    const std::vector<double>&   Q_K_AB,
    const std::vector<double>&   Q_K_CD,
    const std::vector<uint32_t>& pair_displs_AB,
    const std::vector<uint32_t>& pair_displs_CD,
    const std::vector<uint32_t>& pair_counts_AB,
    const std::vector<uint32_t>& pair_counts_CD,
    int      tile_dim_y,
    int      tile_dim_x,
    double   max_D,
    double   tau,
    double   eri_threshold);


ExchangeCuts
build_exchange_cuts_strided(
    const std::vector<uint32_t>& pair_inds_i,
    const std::vector<uint32_t>& pair_inds_k,
    const std::vector<double>&   Q_K_AB,
    const std::vector<double>&   Q_K_CD,
    const std::vector<uint32_t>& pair_displs_AB,
    const std::vector<uint32_t>& pair_displs_CD,
    const std::vector<uint32_t>& pair_counts_AB,
    const std::vector<uint32_t>& pair_counts_CD,
    int      tile_dim_y,
    int      tile_dim_x,
    double   max_D,
    double   tau,
    double   eri_threshold);


ExchangeCuts
build_exchange_cuts_cached_k(
    const std::vector<uint32_t>& pair_inds_i,
    const std::vector<uint32_t>& pair_inds_k,
    const std::vector<double>&   Q_K_AB,
    const std::vector<double>&   Q_K_CD,
    const std::vector<uint32_t>& pair_displs_AB,
    const std::vector<uint32_t>& pair_displs_CD,
    const std::vector<uint32_t>& pair_counts_AB,
    const std::vector<uint32_t>& pair_counts_CD,
    int      tile_dim_y,
    int      tile_dim_x,
    double   max_D,
    double   tau,
    double   eri_threshold);

ExchangeCuts
build_exchange_cuts_grouped_m_cached_k(
    const std::vector<uint32_t>& pair_inds_i,
    const std::vector<uint32_t>& pair_inds_k,
    const std::vector<double>&   Q_K_AB,
    const std::vector<double>&   Q_K_CD,
    const std::vector<uint32_t>& pair_displs_AB,
    const std::vector<uint32_t>& pair_displs_CD,
    const std::vector<uint32_t>& pair_counts_AB,
    const std::vector<uint32_t>& pair_counts_CD,
    int      tile_dim_y,
    int      tile_dim_x,
    uint32_t m_group_size,
    double   max_D,
    double   tau,
    double   eri_threshold);

ExchangeCuts
build_exchange_cut_layout_grouped_m(
    const std::vector<uint32_t>& pair_inds_i,
    const std::vector<uint32_t>& pair_counts_AB,
    int      tile_dim_y,
    uint32_t m_group_size);

ExchangeCuts
build_exchange_cut_layout(
    const std::vector<uint32_t>& pair_inds_i,
    const std::vector<uint32_t>& pair_counts_AB,
    int      tile_dim_y);

namespace gpu {  // gpu namespace

void
build_exchange_cuts_grouped_m_device(
    uint32_t*       d_prec_cut_flat,
    uint32_t*       d_screen_cut_flat,
    const uint32_t* d_displ_cuts,
    const uint32_t* d_pair_inds_i,
    const uint32_t* d_pair_inds_k,
    const double*   d_Q_K_AB,
    const double*   d_Q_K_CD,
    const uint32_t* d_pair_displs_AB,
    const uint32_t* d_pair_displs_CD,
    const uint32_t* d_pair_counts_AB,
    const uint32_t* d_pair_counts_CD,
    uint32_t        n_ik,
    uint32_t        tile_dim_y,
    uint32_t        tile_dim_x,
    uint32_t        m_group_size,
    double          max_D,
    double          tau,
    double          eri_threshold,
    gpuStream_t     stream);

void
build_exchange_cuts_device(
    uint32_t*       d_prec_cut_flat,
    uint32_t*       d_screen_cut_flat,
    const uint32_t* d_displ_cuts,
    const uint32_t* d_pair_inds_i,
    const uint32_t* d_pair_inds_k,
    const double*   d_Q_K_AB,
    const double*   d_Q_K_CD,
    const uint32_t* d_pair_displs_AB,
    const uint32_t* d_pair_displs_CD,
    const uint32_t* d_pair_counts_AB,
    const uint32_t* d_pair_counts_CD,
    uint32_t        n_ik,
    uint32_t        tile_dim_y,
    uint32_t        tile_dim_x,
    double          max_D,
    double          tau,
    double          eri_threshold,
    gpuStream_t     stream);

}  // namespace gpu

#endif
