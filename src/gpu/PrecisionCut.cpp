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

// build cut_ij_tile (粒度为 tile，scheme A / per ij-tile cut)
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

    // build cut_ij_tile (粒度为 tile，scheme A / per ij-tile cut)
    std::vector<uint32_t> cut(nij_tiles, 0);
    for(uint32_t t = 0; t < nij_tiles; ++t) {
        // const double qmax = Q_tile_max[t];
        const double qmax = Q_ij_local[t*tile_dim];
        if(qmax <= 0.0) {
            // 这个ij-tile没贡献
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