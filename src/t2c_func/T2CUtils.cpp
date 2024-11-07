#include "T2CUtils.hpp"

#include <algorithm>
#include <ranges>

#include "CustomViews.hpp"
#include "TensorComponents.hpp"

namespace t2cfunc {  // t2cfunc namespace

auto
comp_distances_ab(CSimdArray<double>& buffer, const size_t index_ab, const size_t index_b, const TPoint<double>& r_a) -> void
{
    // set up R(AB) distances

    auto ab_x = buffer.data(index_ab);

    auto ab_y = buffer.data(index_ab + 1);

    auto ab_z = buffer.data(index_ab + 2);

    // set up Cartesian B coordinates

    auto b_x = buffer.data(index_b);

    auto b_y = buffer.data(index_b + 1);

    auto b_z = buffer.data(index_b + 2);

    // set up Cartesian A coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // compute R(AB) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(ab_x, ab_y, ab_z, b_x, b_y, b_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        ab_x[i] = a_x - b_x[i];

        ab_y[i] = a_y - b_y[i];

        ab_z[i] = a_z - b_z[i];
    }
}

auto
comp_coordinates_p(CSimdArray<double>& buffer, const size_t index_p, const size_t index_b, const TPoint<double>& r_a, const double a_exp) -> void
{
    // Set up exponents

    auto b_exps = buffer.data(0);

    // set up Cartesian P coordinates

    auto p_x = buffer.data(index_p);

    auto p_y = buffer.data(index_p + 1);

    auto p_z = buffer.data(index_p + 2);

    // set up Cartesian B coordinates

    auto b_x = buffer.data(index_b);

    auto b_y = buffer.data(index_b + 1);

    auto b_z = buffer.data(index_b + 2);

    // set up Cartesian A coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // compute R(AB) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(p_x, p_y, p_z, b_x, b_y, b_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        double fact = 1.0 / (a_exp + b_exps[i]);

        p_x[i] = fact * (a_x * a_exp + b_x[i] * b_exps[i]);

        p_y[i] = fact * (a_y * a_exp + b_y[i] * b_exps[i]);

        p_z[i] = fact * (a_z * a_exp + b_z[i] * b_exps[i]);
    }
}

auto
comp_distances_pb(CSimdArray<double>& buffer, const size_t index_pb, const size_t index_ab, const double a_exp) -> void
{
    // Set up exponents

    auto b_exps = buffer.data(0);

    // set up R(PB) distances

    auto pb_x = buffer.data(index_pb);

    auto pb_y = buffer.data(index_pb + 1);

    auto pb_z = buffer.data(index_pb + 2);

    // set up R(AB) distances

    auto ab_x = buffer.data(index_ab);

    auto ab_y = buffer.data(index_ab + 1);

    auto ab_z = buffer.data(index_ab + 2);

    // compute R(PB) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(pb_x, pb_y, pb_z, ab_x, ab_y, ab_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        double fact = a_exp / (a_exp + b_exps[i]);

        pb_x[i] = fact * ab_x[i];

        pb_y[i] = fact * ab_y[i];

        pb_z[i] = fact * ab_z[i];
    }
}

auto
comp_distances_pa(CSimdArray<double>& buffer, const size_t index_pa, const size_t index_ab, const double a_exp) -> void
{
    // Set up exponents

    auto b_exps = buffer.data(0);

    // set up R(PA) distances

    auto pa_x = buffer.data(index_pa);

    auto pa_y = buffer.data(index_pa + 1);

    auto pa_z = buffer.data(index_pa + 2);

    // set up R(AB) distances

    auto ab_x = buffer.data(index_ab);

    auto ab_y = buffer.data(index_ab + 1);

    auto ab_z = buffer.data(index_ab + 2);

    // compute R(PA) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(pa_x, pa_y, pa_z, ab_x, ab_y, ab_z, b_exps : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        double fact = -b_exps[i] / (a_exp + b_exps[i]);

        pa_x[i] = fact * ab_x[i];

        pa_y[i] = fact * ab_y[i];

        pa_z[i] = fact * ab_z[i];
    }
}

auto
comp_distances_pb_from_p(CSimdArray<double>& buffer, const size_t index_pb, const size_t index_p, const size_t index_b) -> void
{
    // set up R(PB) distances

    auto pb_x = buffer.data(index_pb);

    auto pb_y = buffer.data(index_pb + 1);

    auto pb_z = buffer.data(index_pb + 2);

    // set up Cartesian P coordinates

    auto p_x = buffer.data(index_p);

    auto p_y = buffer.data(index_p + 1);

    auto p_z = buffer.data(index_p + 2);

    // set up Cartesian B coordinates

    auto b_x = buffer.data(index_b);

    auto b_y = buffer.data(index_b + 1);

    auto b_z = buffer.data(index_b + 2);

    // compute R(PB) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(pb_x, pb_y, pb_z, p_x, p_y, p_z, b_x, b_y, b_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        pb_x[i] = p_x[i] - b_x[i];

        pb_y[i] = p_y[i] - b_y[i];

        pb_z[i] = p_z[i] - b_z[i];
    }
}

auto
comp_distances_pa_from_p(CSimdArray<double>& buffer, const size_t index_pa, const size_t index_p, const TPoint<double>& r_a) -> void
{
    // set up R(PA) distances

    auto pa_x = buffer.data(index_pa);

    auto pa_y = buffer.data(index_pa + 1);

    auto pa_z = buffer.data(index_pa + 2);

    // set up Cartesian P coordinates

    auto p_x = buffer.data(index_p);

    auto p_y = buffer.data(index_p + 1);

    auto p_z = buffer.data(index_p + 2);

    // set up Cartesian A coordinates

    const auto xyz = r_a.coordinates();

    const auto a_x = xyz[0];

    const auto a_y = xyz[1];

    const auto a_z = xyz[2];

    // compute R(PA) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(pa_x, pa_y, pa_z, p_x, p_y, p_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        pa_x[i] = p_x[i] - a_x;

        pa_y[i] = p_y[i] - a_y;

        pa_z[i] = p_z[i] - a_z;
    }
}

auto
comp_distances_pc(CSimdArray<double>& buffer, const size_t index_pc, const size_t index_p, const TPoint<double>& r_c) -> void
{
    t2cfunc::comp_distances_pa_from_p(buffer, index_pc, index_p, r_c);
}

void
comp_boys_args(CSimdArray<double>& bf_data, const size_t index_args, CSimdArray<double>& buffer, const size_t index_pc, const double a_exp)
{
    // Set up exponents

    auto b_exps = buffer.data(0);

    // set up Boys function arguments

    auto bargs = bf_data.data(index_args);

    // set up R(PC) distances

    auto pc_x = buffer.data(index_pc);

    auto pc_y = buffer.data(index_pc + 1);

    auto pc_z = buffer.data(index_pc + 2);

    // compute Boys function arguments

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(bargs, pc_x, pc_y, pc_z, b_exps : 64)
    for (int i = 0; i < nelems; i++)
    {
        bargs[i] = (a_exp + b_exps[i]) * (pc_x[i] * pc_x[i] + pc_y[i] * pc_y[i] + pc_z[i] * pc_z[i]);
    }
}

void
comp_boys_args(CSimdArray<double>& bf_data,
               const size_t        index_args,
               CSimdArray<double>& buffer,
               const size_t        index_pc,
               const double        a_exp,
               const double        omega)
{
    // Set up exponents

    auto b_exps = buffer.data(0);

    // set up Boys function arguments

    auto bargs = bf_data.data(index_args);

    // set up R(PC) distances

    auto pc_x = buffer.data(index_pc);

    auto pc_y = buffer.data(index_pc + 1);

    auto pc_z = buffer.data(index_pc + 2);

    // compute Boys function arguments

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(bargs, pc_x, pc_y, pc_z, b_exps : 64)
    for (int i = 0; i < nelems; i++)
    {
        const double frho = a_exp + b_exps[i];

        bargs[i] = omega * omega / (omega * omega + frho);

        bargs[i] *= frho * (pc_x[i] * pc_x[i] + pc_y[i] * pc_y[i] + pc_z[i] * pc_z[i]);
    }
}

auto
reduce(CSimdArray<double>& cbuffer, CSimdArray<double>& pbuffer, const size_t position, const size_t ndims, const size_t nblocks) -> void
{
    if (const auto nrows = cbuffer.number_of_rows(); nrows > 0)
    {
        std::ranges::for_each(views::rectangular(nrows, nblocks), [&](const auto& index) {
            const auto [i, j] = index;
            auto pdata        = pbuffer.data(position + i);
            auto cdata        = cbuffer.data(i);
            auto pvals        = &pdata[j * ndims];
#pragma omp simd
            for (size_t k = 0; k < ndims; k++)
            {
                cdata[k] += pvals[k];
            }
        });
    }
}

auto
reduce(CSimdArray<double>& cbuffer, const size_t cposition, CSimdArray<double>& pbuffer, const size_t pposition, const size_t nrows, const size_t ndims, const size_t nblocks) -> void
{
    if (nrows > 0)
    {
        std::ranges::for_each(views::rectangular(nrows, nblocks), [&](const auto& index) {
            const auto [i, j] = index;
            auto pdata        = pbuffer.data(pposition + i);
            auto cdata        = cbuffer.data(cposition + i);
            auto pvals        = &pdata[j * ndims];
#pragma omp simd
            for (size_t k = 0; k < ndims; k++)
            {
                cdata[k] += pvals[k];
            }
        });
    }
}

auto
reduce(CSimdArray<double>& cbuffer, CSimdArray<double>& pbuffer, const size_t position, const double factor, const size_t ndims, const size_t nblocks)
    -> void
{
    if (const auto nrows = cbuffer.number_of_rows(); nrows > 0)
    {
        std::ranges::for_each(views::rectangular(nrows, nblocks), [&](const auto& index) {
            const auto [i, j] = index;
            auto pdata        = pbuffer.data(position + i);
            auto cdata        = cbuffer.data(i);
            auto pvals        = &pdata[j * ndims];
#pragma omp simd
            for (size_t k = 0; k < ndims; k++)
            {
                cdata[k] += factor * pvals[k];
            }
        });
    }
}

auto
reduce(CSimdArray<double>&        cbuffer,
       CSimdArray<double>&        pbuffer,
       const size_t               position,
       const std::vector<double>& facts,
       const size_t               nfacts,
       const size_t               index,
       const size_t               ndims,
       const size_t               nblocks) -> void
{
    if (const auto nrows = cbuffer.number_of_rows(); nrows > 0)
    {
        const auto loc_rows = nrows / nfacts;

        const auto ncenters = facts.size() / nfacts;

        std::ranges::for_each(std::views::iota(size_t{0}, nfacts), [&](const size_t idx) {
            const auto factor = facts[idx * ncenters + index];
            std::ranges::for_each(views::rectangular(loc_rows, nblocks), [&](const auto& rindex) {
                const auto [i, j] = rindex;
                auto pdata        = pbuffer.data(idx * loc_rows + position + i);
                auto cdata        = cbuffer.data(idx * loc_rows + i);
                auto pvals        = &pdata[j * ndims];
#pragma omp simd
                for (size_t k = 0; k < ndims; k++)
                {
                    cdata[k] += factor * pvals[k];
                }
            });
        });
    }
}

auto
distribute(CSubMatrix*                      matrix,
           const double*                    buffer,
           const std::vector<size_t>&       indices,
           const int                        bra_comp,
           const int                        ket_comp,
           const size_t                     bra_igto,
           const std::pair<size_t, size_t>& ket_range) -> void
{
    const auto bra_idx = indices[0] * bra_comp + indices[bra_igto + 1];

    std::ranges::for_each(std::views::iota(ket_range.first, ket_range.second),
                          [&](const auto i) { matrix->operator[]({bra_idx, indices[0] * ket_comp + indices[i + 1]}) = buffer[i - ket_range.first]; });
}

auto
distribute(CSubMatrix*                      matrix,
           const double*                    buffer,
           const std::vector<size_t>&       bra_indices,
           const std::vector<size_t>&       ket_indices,
           const int                        bra_comp,
           const int                        ket_comp,
           const size_t                     bra_igto,
           const std::pair<size_t, size_t>& ket_range,
           const mat_t                      mat_type) -> void
{
    const auto bra_idx = bra_indices[0] * bra_comp + bra_indices[bra_igto + 1];

    std::ranges::for_each(std::views::iota(ket_range.first, ket_range.second), [&](const auto i) {
        const auto ket_idx                     = ket_indices[0] * ket_comp + ket_indices[i + 1];
        const auto fval                        = buffer[i - ket_range.first];
        matrix->operator[]({bra_idx, ket_idx}) = fval;
        if (mat_type == mat_t::symmetric)
        {
            matrix->operator[]({ket_idx, bra_idx}) = fval;
        }
        if (mat_type == mat_t::antisymmetric)
        {
            matrix->operator[]({ket_idx, bra_idx}) = -fval;
        }
    });
}

auto
distribute(CSubMatrix*                      matrix,
           const double*                    buffer,
           const std::vector<size_t>&       bra_indices,
           const std::vector<size_t>&       ket_indices,
           const int                        bra_comp,
           const int                        ket_comp,
           const size_t                     bra_igto,
           const std::pair<size_t, size_t>& ket_range,
           const bool                       ang_order) -> void
{
    const auto bra_idx = bra_indices[0] * bra_comp + bra_indices[bra_igto + 1];

    std::ranges::for_each(std::views::iota(ket_range.first, ket_range.second), [&](const auto i) {
        const auto ket_idx = ket_indices[0] * ket_comp + ket_indices[i + 1];
        if (ang_order)
        {
            matrix->operator[]({bra_idx, ket_idx}) = buffer[i - ket_range.first];
        }
        else
        {
            matrix->operator[]({ket_idx, bra_idx}) = buffer[i - ket_range.first];
        }
    });
}

auto
distribute(CSubMatrix*                      matrix,
           const CSimdArray<double>&        buffer,
           const size_t                     offset,
           const std::vector<size_t>&       indices,
           const int                        bra_angmom,
           const size_t                     bra_igto,
           const std::pair<size_t, size_t>& ket_range) -> void
{
    const auto bra_comps = tensor::number_of_spherical_components(std::array<int, 1>{bra_angmom});

    std::ranges::for_each(views::rectangular(bra_comps, bra_comps), [&](const auto& index) {
        const auto [i, j] = index;
        t2cfunc::distribute(matrix, buffer.data(offset + i * bra_comps + j), indices, i, j, bra_igto, ket_range);
    });
}

auto
distribute(CSubMatrix*                      matrix,
           const CSimdArray<double>&        buffer,
           const size_t                     offset,
           const std::vector<size_t>&       bra_indices,
           const std::vector<size_t>&       ket_indices,
           const int                        bra_angmom,
           const int                        ket_angmom,
           const size_t                     bra_igto,
           const std::pair<size_t, size_t>& ket_range,
           const mat_t                      mat_type) -> void
{
    const auto bra_comps = tensor::number_of_spherical_components(std::array<int, 1>{bra_angmom});

    const auto ket_comps = tensor::number_of_spherical_components(std::array<int, 1>{ket_angmom});

    std::ranges::for_each(views::rectangular(bra_comps, ket_comps), [&](const auto& index) {
        const auto [i, j] = index;
        t2cfunc::distribute(matrix, buffer.data(offset + i * ket_comps + j), bra_indices, ket_indices, i, j, bra_igto, ket_range, mat_type);
    });
}

auto
distribute(CSubMatrix*                      matrix,
           const CSimdArray<double>&        buffer,
           const size_t                     offset,
           const std::vector<size_t>&       bra_indices,
           const std::vector<size_t>&       ket_indices,
           const int                        bra_angmom,
           const int                        ket_angmom,
           const size_t                     bra_igto,
           const std::pair<size_t, size_t>& ket_range,
           const bool                       ang_order) -> void
{
    const auto bra_comps = tensor::number_of_spherical_components(std::array<int, 1>{bra_angmom});

    const auto ket_comps = tensor::number_of_spherical_components(std::array<int, 1>{ket_angmom});

    std::ranges::for_each(views::rectangular(bra_comps, ket_comps), [&](const auto& index) {
        const auto [i, j] = index;
        t2cfunc::distribute(matrix, buffer.data(offset + i * ket_comps + j), bra_indices, ket_indices, i, j, bra_igto, ket_range, ang_order);
    });
}

}  // namespace t2cfunc
