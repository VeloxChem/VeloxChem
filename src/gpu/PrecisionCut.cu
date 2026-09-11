#include "PrecisionCut.hpp"

#include <algorithm>
#include <cmath>
#include <functional>

#if defined(USE_CUDA) || defined(USE_HIP)
#include "GpuRuntime.hpp"
#include "GpuWrapper.hpp"
#endif

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
    const uint32_t nij_tiles =
    (ij_count_local + tile_dim - 1) / tile_dim;

    // QD_tile_max along kl
    //    QD_kl = |Q_kl| * |D_kl|
    std::vector<double> QD_kl(kl_count);
    for (uint32_t kl = 0; kl < kl_count; ++kl) {
        QD_kl[kl] = std::abs(Q_kl[kl]) * std::abs(D_kl[kl]);
    }
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

    return {{}, {}, displ_cuts, total};
}

#if defined(USE_CUDA) || defined(USE_HIP)

namespace {

// Search descending values at stride tile_dim directly, without building a q_tile array.
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

__global__ void
accumulate_exchange_cut_work_kernel(
    const uint32_t* d_prec_cut_flat,
    const uint32_t* d_screen_cut_flat,
    const uint32_t* d_displ_cuts,
    const uint32_t* d_pair_inds_i,
    const uint32_t* d_pair_inds_k,
    const uint32_t* d_pair_counts_AB,
    const uint32_t* d_pair_counts_CD,
    uint32_t        n_ik,
    uint32_t        tile_dim_y,
    uint32_t        tile_dim_x,
    unsigned long long* d_work_counts)
{
    const uint32_t ik = blockIdx.x;
    if (ik >= n_ik) {
        return;
    }

    const uint32_t i = d_pair_inds_i[ik];
    const uint32_t k = d_pair_inds_k[ik];
    const uint32_t n_m = (d_pair_counts_AB[i] + tile_dim_y - 1) / tile_dim_y;
    const uint32_t n_n = (d_pair_counts_CD[k] + tile_dim_x - 1) / tile_dim_x;

    unsigned long long fp64 = 0;
    unsigned long long fp32 = 0;
    unsigned long long screened = 0;

    for (uint32_t m = threadIdx.x; m < n_m; m += blockDim.x) {
        const uint32_t entry = d_displ_cuts[ik] + m;
        const uint32_t prec = d_prec_cut_flat[entry];
        const uint32_t screen = d_screen_cut_flat[entry];

        fp64 += prec;
        fp32 += screen - prec;
        screened += n_n - screen;
    }

    __shared__ unsigned long long block_counts[3][128];
    block_counts[0][threadIdx.x] = fp64;
    block_counts[1][threadIdx.x] = fp32;
    block_counts[2][threadIdx.x] = screened;
    __syncthreads();

    for (uint32_t stride = blockDim.x / 2; stride > 0; stride /= 2) {
        if (threadIdx.x < stride) {
            block_counts[0][threadIdx.x] += block_counts[0][threadIdx.x + stride];
            block_counts[1][threadIdx.x] += block_counts[1][threadIdx.x + stride];
            block_counts[2][threadIdx.x] += block_counts[2][threadIdx.x + stride];
        }
        __syncthreads();
    }

    if (threadIdx.x == 0) {
        atomicAdd(&d_work_counts[0], block_counts[0][0]);
        atomicAdd(&d_work_counts[1], block_counts[1][0]);
        atomicAdd(&d_work_counts[2], block_counts[2][0]);
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
    gpuStream_t     stream,
    unsigned long long* d_work_counts)
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

    // Check the host-held pointer, not the contents of device memory.
    // FockDriverGPU allocates this buffer only when VLX_EXCHANGE_FRACTION_STATS=1;
    // otherwise it passes nullptr to skip work statistics.
    if (d_work_counts != nullptr) {
        accumulate_exchange_cut_work_kernel<<<n_ik, threads_per_block, 0, stream>>>(
            d_prec_cut_flat,
            d_screen_cut_flat,
            d_displ_cuts,
            d_pair_inds_i,
            d_pair_inds_k,
            d_pair_counts_AB,
            d_pair_counts_CD,
            n_ik,
            tile_dim_y,
            tile_dim_x,
            d_work_counts);
    }
}

}  // namespace gpu

#endif  // USE_CUDA || USE_HIP
