#include "PrecisionCut.hpp"

#include <algorithm>
#include <cmath>
#include <functional>

#if defined(USE_CUDA) || defined(USE_HIP)
#include "GpuRuntime.hpp"
#include "GpuWrapper.hpp"
#endif

// Return the max value in each tile
std::vector<double> tile_max_abs(const std::vector<double>& arr, int tile_dim) {
    const size_t n = arr.size();
    const size_t n_tiles = (n + tile_dim - 1) / tile_dim;
    std::vector<double> res(n_tiles, 0.0);
    for(size_t i = 0; i < n_tiles; ++i) {
        size_t beg = i * tile_dim;
        size_t end = std::min(beg + tile_dim, n);
        double max_value = 0.0;
        for(size_t j = beg; j < end; ++j) {
            // Only the magnitude matters; the sign is irrelevant
            max_value = std::max(max_value, std::abs(arr[j]));
        }
        res[i] = max_value;
    }
    return res;
}

// Reutrn the index of the first <= thr (i.e. count(>thr))
uint32_t cut_descending(const std::vector<double>& desc, double thr) {
    uint32_t left = 0, right = (uint32_t)desc.size();
    while(left < right) {
        uint32_t mid = left + (right - left) / 2;
        if(desc[mid] > thr) left = mid + 1;
        else right = mid;
    }
    return left;
}

// Build cut_ij_tile at tile granularity using scheme A / per ij-tile cuts.
std::vector<uint32_t> build_cut_ij_tile(
    const std::vector<double>& Q_ij_local,   // len >= ij_count_local
    const std::vector<double>& Q_kl,         // len >= kl_count
    const std::vector<double>& D_kl,         // len >= kl_count
    uint32_t ij_count_local,
    uint32_t kl_count,
    int tile_dim,
    double tau) 
{
    // Q_tile_max along ij
    std::vector<double> Q_local(Q_ij_local.begin(), Q_ij_local.begin() + ij_count_local);
    // auto Q_tile_max = tile_max_abs(Q_local, tile_dim);
    // const int32_t nij_tiles = (int32_t)Q_tile_max.size();
    const uint32_t nij_tiles =
    (ij_count_local + tile_dim - 1) / tile_dim;

    // QD_tile_max along kl
    //    QD_kl = |Q_kl| * |D_kl|
    std::vector<double> QD_kl(kl_count);
    for (uint32_t kl = 0; kl < kl_count; ++kl) {
        QD_kl[kl] = std::abs(Q_kl[kl]) * std::abs(D_kl[kl]);
    }
    // auto QD_tile_max = tile_max_abs(QD_kl, tile_dim);
    const uint32_t nkl_tiles = (kl_count + tile_dim - 1) / tile_dim;
    std::vector<double> QD_tile_max(nkl_tiles);

    for(uint32_t t = 0; t < nkl_tiles; ++t)
    {
        QD_tile_max[t] = QD_kl[t * tile_dim];
    }

    // Build cut_ij_tile at tile granularity using scheme A / per ij-tile cuts.
    std::vector<uint32_t> cut(nij_tiles, 0);
    for(uint32_t t = 0; t < nij_tiles; ++t) {
        // const double qmax = Q_tile_max[t];
        const double qmax = Q_ij_local[t*tile_dim];
        if(qmax <= 0.0) {
            // This ij-tile has no contribution.
            cut[t] = 0; 
            continue;
        }

        const double thr = tau/qmax;
        // cut[t] = cut_descending(QD_tile_max, thr);
        auto it = std::lower_bound(
            QD_tile_max.begin(),
            QD_tile_max.end(),
            thr,
            std::greater<double>()
        );
        cut[t] = static_cast<uint32_t>(it - QD_tile_max.begin());
    }
    return cut;
}

std::vector<uint32_t> build_cut_ij_tile_dd(
    const std::vector<double>& Q_ij_local,
    const std::vector<double>& Q_kl,
    const std::vector<double>& D_kl,
    uint32_t ij_count_local,
    uint32_t kl_count,
    int ij_tile_dim,
    int kl_tile_dim,
    double tau)
{
    const uint32_t nij_tiles = (ij_count_local + ij_tile_dim - 1) / ij_tile_dim;

    std::vector<double> QD_kl(kl_count);
    for (uint32_t kl = 0; kl < kl_count; ++kl) {
        QD_kl[kl] = std::abs(Q_kl[kl]) * std::abs(D_kl[kl]);
    }

    const uint32_t nkl_tiles = (kl_count + kl_tile_dim - 1) / kl_tile_dim;
    std::vector<double> QD_tile_max(nkl_tiles);

    for (uint32_t t = 0; t < nkl_tiles; ++t)
    {
        QD_tile_max[t] = QD_kl[t * kl_tile_dim];
    }

    std::vector<uint32_t> cut(nij_tiles, 0);
    for (uint32_t t = 0; t < nij_tiles; ++t) {
        const double qmax = Q_ij_local[t * ij_tile_dim];
        if (qmax <= 0.0) {
            cut[t] = 0;
            continue;
        }

        const double thr = tau / qmax;
        auto it = std::lower_bound(
            QD_tile_max.begin(),
            QD_tile_max.end(),
            thr,
            std::greater<double>()
        );
        cut[t] = static_cast<uint32_t>(it - QD_tile_max.begin());
    }
    return cut;
}

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
    double   eri_threshold)
{
    const uint32_t n_ik = static_cast<uint32_t>(pair_inds_i.size());

    // Pass 1: compute displ_cuts (prefix sum of n_m per ik) and total flat size.
    // n_m[ik] = number of j-tiles for this ik block, determined by the i-primitive's pair count.
    std::vector<uint32_t> displ_cuts(n_ik);
    uint32_t total = 0;
    for (uint32_t ik = 0; ik < n_ik; ++ik) {
        displ_cuts[ik] = total;
        const uint32_t n_m = (pair_counts_AB[pair_inds_i[ik]] + tile_dim_y - 1) / tile_dim_y;
        total += n_m;
    }

    std::vector<uint32_t> prec_cut_flat(total, 0);
    std::vector<uint32_t> screen_cut_flat(total, 0);

    // Pass 2: for each (ik, m-tile), binary search in n-tiles to find:
    //   prec_cut  [ik][m] = first n s.t. Q_ij_m * Q_kl_n * max_D <= tau,
    //                         where FP32 is sufficient
    //   screen_cut[ik][m] = first n s.t. Q_ij_m * Q_kl_n * max_D <= eri_threshold,
    //                         where the tile is screened out
    // Stored flat: index = displ_cuts[ik] + m
    for (uint32_t ik = 0; ik < n_ik; ik++) {
        const auto i = pair_inds_i[ik];
        const auto k = pair_inds_k[ik];

        const uint32_t displ_i = pair_displs_AB[i];
        const uint32_t count_i = pair_counts_AB[i];
        const uint32_t displ_k = pair_displs_CD[k];
        const uint32_t count_k = pair_counts_CD[k];

        const uint32_t n_m = (count_i + tile_dim_y - 1) / tile_dim_y;
        const uint32_t n_n = (count_k + tile_dim_x - 1) / tile_dim_x;

        // Representative Q for each n-tile: first element (largest, since Q_K_CD is descending).
        // Shared across all m-tiles for this ik block.
        std::vector<double> Q_tile(n_n);
        for (uint32_t n = 0; n < n_n; ++n) {
            Q_tile[n] = Q_K_CD[displ_k + n * tile_dim_x];
        }

        for (uint32_t m = 0; m < n_m; m++) {
            // Representative Q_ij for this j-tile: first element (largest).
            const double Q_ij_m = Q_K_AB[displ_i + m * tile_dim_y];

            if (Q_ij_m <= 0.0) {
                // No contribution from this j-tile; both cuts stay 0.
                continue;
            }

            // Rearrange threshold: Q_kl_n > tau / (Q_ij_m * max_D), so FP64 is needed.
            const double thr_prec   = tau           / (Q_ij_m * max_D);
            const double thr_screen = eri_threshold / (Q_ij_m * max_D);

            // lower_bound on descending Q_tile finds first n where Q_tile[n] <= thr
            auto it_prec   = std::lower_bound(Q_tile.begin(), Q_tile.end(), thr_prec,   std::greater<double>());
            auto it_screen = std::lower_bound(Q_tile.begin(), Q_tile.end(), thr_screen, std::greater<double>());

            prec_cut_flat  [displ_cuts[ik] + m] = static_cast<uint32_t>(it_prec   - Q_tile.begin());
            screen_cut_flat[displ_cuts[ik] + m] = static_cast<uint32_t>(it_screen - Q_tile.begin());
        }
    }

    return {prec_cut_flat, screen_cut_flat, displ_cuts};
}


static uint32_t
lower_bound_desc_strided(const std::vector<double>& q,
                         uint32_t displ,
                         uint32_t n_tiles,
                         uint32_t tile_dim,
                         double thr)
{
    uint32_t left = 0;
    uint32_t right = n_tiles;

    while (left < right) {
        const uint32_t mid = left + (right - left) / 2;
        const double val = q[displ + mid * tile_dim];

        if (val > thr) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }

    return left;
}

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
    double   eri_threshold)
{
    const uint32_t n_ik = static_cast<uint32_t>(pair_inds_i.size());

    std::vector<uint32_t> displ_cuts(n_ik);
    uint32_t total = 0;
    for (uint32_t ik = 0; ik < n_ik; ++ik) {
        displ_cuts[ik] = total;
        const uint32_t n_m = (pair_counts_AB[pair_inds_i[ik]] + tile_dim_y - 1) / tile_dim_y;
        total += n_m;
    }

    std::vector<uint32_t> prec_cut_flat(total, 0);
    std::vector<uint32_t> screen_cut_flat(total, 0);

    for (uint32_t ik = 0; ik < n_ik; ++ik) {
        const auto i = pair_inds_i[ik];
        const auto k = pair_inds_k[ik];

        const uint32_t displ_i = pair_displs_AB[i];
        const uint32_t count_i = pair_counts_AB[i];
        const uint32_t displ_k = pair_displs_CD[k];
        const uint32_t count_k = pair_counts_CD[k];

        const uint32_t n_m = (count_i + tile_dim_y - 1) / tile_dim_y;
        const uint32_t n_n = (count_k + tile_dim_x - 1) / tile_dim_x;

        for (uint32_t m = 0; m < n_m; ++m) {
            const double Q_ij_m = Q_K_AB[displ_i + m * tile_dim_y];

            if (Q_ij_m <= 0.0) {
                continue;
            }

            const double thr_prec   = tau           / (Q_ij_m * max_D);
            const double thr_screen = eri_threshold / (Q_ij_m * max_D);

            prec_cut_flat[displ_cuts[ik] + m] =
                lower_bound_desc_strided(Q_K_CD, displ_k, n_n, tile_dim_x, thr_prec);
            screen_cut_flat[displ_cuts[ik] + m] =
                lower_bound_desc_strided(Q_K_CD, displ_k, n_n, tile_dim_x, thr_screen);
        }
    }

    return {prec_cut_flat, screen_cut_flat, displ_cuts};
}


static uint32_t
lower_bound_desc(const std::vector<double>& q_tile,
                        double thr)
{
    const auto it = std::lower_bound(q_tile.begin(), q_tile.end(), thr, std::greater<double>());
    return static_cast<uint32_t>(it - q_tile.begin());
}

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
    double   eri_threshold)
{
    const uint32_t n_ik = static_cast<uint32_t>(pair_inds_i.size());

    std::vector<uint32_t> displ_cuts(n_ik);
    uint32_t total = 0;
    for (uint32_t ik = 0; ik < n_ik; ++ik) {
        displ_cuts[ik] = total;
        const uint32_t n_m = (pair_counts_AB[pair_inds_i[ik]] + tile_dim_y - 1) / tile_dim_y;
        total += n_m;
    }

    std::vector<uint32_t> prec_cut_flat(total, 0);
    std::vector<uint32_t> screen_cut_flat(total, 0);

    std::vector<std::vector<double>> q_tile_cache(pair_counts_CD.size());
    std::vector<unsigned char> q_tile_cached(pair_counts_CD.size(), 0);

    for (uint32_t ik = 0; ik < n_ik; ++ik) {
        const auto i = pair_inds_i[ik];
        const auto k = pair_inds_k[ik];

        const uint32_t displ_i = pair_displs_AB[i];
        const uint32_t count_i = pair_counts_AB[i];
        const uint32_t displ_k = pair_displs_CD[k];
        const uint32_t count_k = pair_counts_CD[k];

        const uint32_t n_m = (count_i + tile_dim_y - 1) / tile_dim_y;
        const uint32_t n_n = (count_k + tile_dim_x - 1) / tile_dim_x;

        if (!q_tile_cached[k]) {
            auto& q_tile = q_tile_cache[k];
            q_tile.resize(n_n);
            for (uint32_t n = 0; n < n_n; ++n) {
                q_tile[n] = Q_K_CD[displ_k + n * tile_dim_x];
            }
            q_tile_cached[k] = 1;
        }

        const auto& q_tile = q_tile_cache[k];

        for (uint32_t m = 0; m < n_m; ++m) {
            const double Q_ij_m = Q_K_AB[displ_i + m * tile_dim_y];

            if (Q_ij_m <= 0.0) {
                continue;
            }

            const double inv_qd = 1.0 / (Q_ij_m * max_D);
            const double thr_prec = tau * inv_qd;
            const double thr_screen = eri_threshold * inv_qd;

            prec_cut_flat[displ_cuts[ik] + m] = lower_bound_desc(q_tile, thr_prec);
            screen_cut_flat[displ_cuts[ik] + m] = lower_bound_desc(q_tile, thr_screen);
        }
    }

    return {prec_cut_flat, screen_cut_flat, displ_cuts};
}


ExchangeCuts
build_exchange_cut_layout_grouped_m(
    const std::vector<uint32_t>& pair_inds_i,
    const std::vector<uint32_t>& pair_counts_AB,
    int      tile_dim_y,
    uint32_t m_group_size)
{
    if (m_group_size == 0) {
        m_group_size = 1;
    }

    const uint32_t n_ik = static_cast<uint32_t>(pair_inds_i.size());

    std::vector<uint32_t> displ_cuts(n_ik);
    uint32_t total = 0;
    for (uint32_t ik = 0; ik < n_ik; ++ik) {
        displ_cuts[ik] = total;
        const uint32_t n_m = (pair_counts_AB[pair_inds_i[ik]] + tile_dim_y - 1) / tile_dim_y;
        total += (n_m + m_group_size - 1) / m_group_size;
    }

    return {{}, {}, displ_cuts, {}, total};
}

ExchangeCuts
build_exchange_cut_layout(
    const std::vector<uint32_t>& pair_inds_i,
    const std::vector<uint32_t>& pair_counts_AB,
    int      tile_dim_y)
{
    const uint32_t n_ik = static_cast<uint32_t>(pair_inds_i.size());

    std::vector<uint32_t> displ_cuts(n_ik);
    uint32_t total = 0;
    for (uint32_t ik = 0; ik < n_ik; ++ik) {
        displ_cuts[ik] = total;
        const uint32_t n_m = (pair_counts_AB[pair_inds_i[ik]] + tile_dim_y - 1) / tile_dim_y;
        total += n_m;
    }

    return {{}, {}, displ_cuts, {}, total};
}

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
    double   eri_threshold)
{
    if (m_group_size == 0) {
        m_group_size = 1;
    }

    const uint32_t n_ik = static_cast<uint32_t>(pair_inds_i.size());

    std::vector<uint32_t> displ_cuts(n_ik);
    uint32_t total = 0;
    for (uint32_t ik = 0; ik < n_ik; ++ik) {
        displ_cuts[ik] = total;
        const uint32_t n_m = (pair_counts_AB[pair_inds_i[ik]] + tile_dim_y - 1) / tile_dim_y;
        total += (n_m + m_group_size - 1) / m_group_size;
    }

    std::vector<uint32_t> prec_cut_flat(total, 0);
    std::vector<uint32_t> screen_cut_flat(total, 0);
    std::vector<uint32_t> cut_weights(total, 1);

    std::vector<std::vector<double>> q_tile_cache(pair_counts_CD.size());
    std::vector<unsigned char> q_tile_cached(pair_counts_CD.size(), 0);

    for (uint32_t ik = 0; ik < n_ik; ++ik) {
        const auto i = pair_inds_i[ik];
        const auto k = pair_inds_k[ik];

        const uint32_t displ_i = pair_displs_AB[i];
        const uint32_t count_i = pair_counts_AB[i];
        const uint32_t displ_k = pair_displs_CD[k];
        const uint32_t count_k = pair_counts_CD[k];

        const uint32_t n_m = (count_i + tile_dim_y - 1) / tile_dim_y;
        const uint32_t n_n = (count_k + tile_dim_x - 1) / tile_dim_x;
        const uint32_t n_m_groups = (n_m + m_group_size - 1) / m_group_size;

        if (!q_tile_cached[k]) {
            auto& q_tile = q_tile_cache[k];
            q_tile.resize(n_n);
            for (uint32_t n = 0; n < n_n; ++n) {
                q_tile[n] = Q_K_CD[displ_k + n * tile_dim_x];
            }
            q_tile_cached[k] = 1;
        }

        const auto& q_tile = q_tile_cache[k];

        for (uint32_t g = 0; g < n_m_groups; ++g) {
            const uint32_t m_begin = g * m_group_size;
            const uint32_t m_end = std::min(n_m, m_begin + m_group_size);

            cut_weights[displ_cuts[ik] + g] = m_end - m_begin;

            // double Q_ij_group = 0.0;
            // for (uint32_t m = m_begin; m < m_end; ++m) {
            //     Q_ij_group = std::max(Q_ij_group, Q_K_AB[displ_i + m * tile_dim_y]);
            // }
            const double Q_ij_group = Q_K_AB[displ_i + m_begin * tile_dim_y];

            if (Q_ij_group <= 0.0) {
                continue;
            }

            const double inv_qd = 1.0 / (Q_ij_group * max_D);
            const double thr_prec = tau * inv_qd;
            const double thr_screen = eri_threshold * inv_qd;

            prec_cut_flat[displ_cuts[ik] + g] = lower_bound_desc(q_tile, thr_prec);
            screen_cut_flat[displ_cuts[ik] + g] = lower_bound_desc(q_tile, thr_screen);
        }
    }

    return {prec_cut_flat, screen_cut_flat, displ_cuts, cut_weights};
}

#if defined(USE_CUDA) || defined(USE_HIP)

namespace {

__device__ uint32_t
lower_bound_desc_device(const double* q,
                        uint32_t      displ,
                        uint32_t      n_tiles,
                        uint32_t      tile_dim,
                        double        thr)
{
    uint32_t left = 0;
    uint32_t right = n_tiles;

    while (left < right) {
        const uint32_t mid = left + (right - left) / 2;
        const double val = q[displ + mid * tile_dim];

        if (val > thr) {
            left = mid + 1;
        } else {
            right = mid;
        }
    }

    return left;
}

__global__ void
build_exchange_cuts_kernel(
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
    double          eri_threshold)
{
    const uint32_t ik = blockIdx.x;
    if (ik >= n_ik) {
        return;
    }

    const uint32_t i = d_pair_inds_i[ik];
    const uint32_t k = d_pair_inds_k[ik];

    const uint32_t displ_i = d_pair_displs_AB[i];
    const uint32_t displ_k = d_pair_displs_CD[k];
    const uint32_t count_i = d_pair_counts_AB[i];
    const uint32_t count_k = d_pair_counts_CD[k];

    const uint32_t n_m = (count_i + tile_dim_y - 1) / tile_dim_y;
    const uint32_t n_n = (count_k + tile_dim_x - 1) / tile_dim_x;

    for (uint32_t m = threadIdx.x; m < n_m; m += blockDim.x) {
        const uint32_t entry = d_displ_cuts[ik] + m;
        const double Q_ij_m = d_Q_K_AB[displ_i + m * tile_dim_y];

        if (Q_ij_m <= 0.0) {
            d_prec_cut_flat[entry] = 0;
            d_screen_cut_flat[entry] = 0;
            continue;
        }

        const double inv_qd = 1.0 / (Q_ij_m * max_D);
        const double thr_prec = tau * inv_qd;
        const double thr_screen = eri_threshold * inv_qd;

        d_prec_cut_flat[entry] = lower_bound_desc_device(d_Q_K_CD, displ_k, n_n, tile_dim_x, thr_prec);
        d_screen_cut_flat[entry] = lower_bound_desc_device(d_Q_K_CD, displ_k, n_n, tile_dim_x, thr_screen);
    }
}

// Grouped-m cut builder kept for benchmarking/reference. The active PPPP MP path
// uses build_exchange_cuts_kernel because testing showed m_group_size=1 is fastest.
__global__ void
build_exchange_cuts_grouped_m_kernel(
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
    double          eri_threshold)
{
    const uint32_t ik = blockIdx.x;
    if (ik >= n_ik) {
        return;
    }

    if (m_group_size == 0) {
        m_group_size = 1;
    }

    const uint32_t i = d_pair_inds_i[ik];
    const uint32_t k = d_pair_inds_k[ik];

    const uint32_t displ_i = d_pair_displs_AB[i];
    const uint32_t displ_k = d_pair_displs_CD[k];
    const uint32_t count_i = d_pair_counts_AB[i];
    const uint32_t count_k = d_pair_counts_CD[k];

    const uint32_t n_m = (count_i + tile_dim_y - 1) / tile_dim_y;
    const uint32_t n_n = (count_k + tile_dim_x - 1) / tile_dim_x;
    const uint32_t n_m_groups = (n_m + m_group_size - 1) / m_group_size;

    for (uint32_t g = threadIdx.x; g < n_m_groups; g += blockDim.x) {
        const uint32_t m_begin = g * m_group_size;
        const uint32_t entry = d_displ_cuts[ik] + g;
        const double Q_ij_group = d_Q_K_AB[displ_i + m_begin * tile_dim_y];

        if (Q_ij_group <= 0.0) {
            d_prec_cut_flat[entry] = 0;
            d_screen_cut_flat[entry] = 0;
            continue;
        }

        const double inv_qd = 1.0 / (Q_ij_group * max_D);
        const double thr_prec = tau * inv_qd;
        const double thr_screen = eri_threshold * inv_qd;

        d_prec_cut_flat[entry] = lower_bound_desc_device(d_Q_K_CD, displ_k, n_n, tile_dim_x, thr_prec);
        d_screen_cut_flat[entry] = lower_bound_desc_device(d_Q_K_CD, displ_k, n_n, tile_dim_x, thr_screen);
    }
}

} // namespace

namespace gpu {  // gpu namespace

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
    gpuStream_t     stream)
{
    if (n_ik == 0) {
        return;
    }

    constexpr uint32_t threads_per_block = 128;
    build_exchange_cuts_kernel<<<n_ik, threads_per_block, 0, stream>>>(
        d_prec_cut_flat,
        d_screen_cut_flat,
        d_displ_cuts,
        d_pair_inds_i,
        d_pair_inds_k,
        d_Q_K_AB,
        d_Q_K_CD,
        d_pair_displs_AB,
        d_pair_displs_CD,
        d_pair_counts_AB,
        d_pair_counts_CD,
        n_ik,
        tile_dim_y,
        tile_dim_x,
        max_D,
        tau,
        eri_threshold);
}

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
    gpuStream_t     stream)
{
    if (n_ik == 0) {
        return;
    }

    constexpr uint32_t threads_per_block = 128;
    build_exchange_cuts_grouped_m_kernel<<<n_ik, threads_per_block, 0, stream>>>(
        d_prec_cut_flat,
        d_screen_cut_flat,
        d_displ_cuts,
        d_pair_inds_i,
        d_pair_inds_k,
        d_Q_K_AB,
        d_Q_K_CD,
        d_pair_displs_AB,
        d_pair_displs_CD,
        d_pair_counts_AB,
        d_pair_counts_CD,
        n_ik,
        tile_dim_y,
        tile_dim_x,
        m_group_size,
        max_D,
        tau,
        eri_threshold);
}

}  // namespace gpu

#endif  // USE_CUDA || USE_HIP
