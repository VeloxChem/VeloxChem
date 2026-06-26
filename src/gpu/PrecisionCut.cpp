#include "PrecisionCut.hpp"

#include <algorithm>
#include <cmath>
#include <functional>

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
