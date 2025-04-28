#include "T2CUtils.hpp"

#include <algorithm>
#include <ranges>

#include "CustomViews.hpp"
#include "TensorComponents.hpp"

#include <iostream>

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
comp_distances_gb(CSimdArray<double>& buffer, const size_t index_gb, const size_t index_p, const size_t index_b, const TPoint<double>& r_c, const double a_exp, const double c_exp) -> void
{
    // Set up exponents

    auto b_exps = buffer.data(0);
    
    // set up R(GB) distances

    auto gb_x = buffer.data(index_gb);

    auto gb_y = buffer.data(index_gb + 1);

    auto gb_z = buffer.data(index_gb + 2);

    // set up Cartesian P coordinates

    auto p_x = buffer.data(index_p);

    auto p_y = buffer.data(index_p + 1);

    auto p_z = buffer.data(index_p + 2);

    // set up Cartesian B coordinates

    auto b_x = buffer.data(index_b);

    auto b_y = buffer.data(index_b + 1);

    auto b_z = buffer.data(index_b + 2);
    
    // set up Cartesian C coordinates

    const auto xyz = r_c.coordinates();

    const auto c_x = xyz[0];

    const auto c_y = xyz[1];

    const auto c_z = xyz[2];

    // compute R(PB) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(gb_x, gb_y, gb_z, p_x, p_y, p_z, b_x, b_y, b_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const auto fact = a_exp + b_exps[i];
        
        const auto finv = 1.0 / (fact + c_exp);
        
        gb_x[i] = (fact * p_x[i] + c_exp * c_x) * finv - b_x[i];

        gb_y[i] = (fact * p_y[i] + c_exp * c_y) * finv - b_y[i];

        gb_z[i] = (fact * p_z[i] + c_exp * c_z) * finv - b_z[i];
    }
}

auto
comp_distances_ga(CSimdArray<double>& buffer, const size_t index_ga, const size_t index_p, const TPoint<double>& r_a, const TPoint<double>& r_c, const double a_exp, const double c_exp) -> void
{
    // Set up exponents

    auto b_exps = buffer.data(0);
    
    // set up R(GA) distances

    auto ga_x = buffer.data(index_ga);

    auto ga_y = buffer.data(index_ga + 1);

    auto ga_z = buffer.data(index_ga + 2);

    // set up Cartesian P coordinates

    auto p_x = buffer.data(index_p);

    auto p_y = buffer.data(index_p + 1);

    auto p_z = buffer.data(index_p + 2);

    // set up Cartesian A coordinates

    const auto a_xyz = r_a.coordinates();

    const auto a_x = a_xyz[0];

    const auto a_y = a_xyz[1];

    const auto a_z = a_xyz[2];
    
    // set up Cartesian C coordinates

    const auto c_xyz = r_c.coordinates();

    const auto c_x = c_xyz[0];

    const auto c_y = c_xyz[1];

    const auto c_z = c_xyz[2];

    // compute R(PB) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(ga_x, ga_y, ga_z, p_x, p_y, p_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const auto fact = a_exp + b_exps[i];
        
        const auto finv = 1.0 / (fact + c_exp);
        
        ga_x[i] = (fact * p_x[i] + c_exp * c_x) * finv - a_x;

        ga_y[i] = (fact * p_y[i] + c_exp * c_y) * finv - a_y;

        ga_z[i] = (fact * p_z[i] + c_exp * c_z) * finv - a_z;
    }
}

auto
comp_distances_gc(CSimdArray<double>& buffer, const size_t index_gc, const size_t index_p, const TPoint<double>& r_c, const double a_exp, const double c_exp) -> void
{
    // Set up exponents

    auto b_exps = buffer.data(0);
    
    // set up R(GB) distances

    auto gc_x = buffer.data(index_gc);

    auto gc_y = buffer.data(index_gc + 1);

    auto gc_z = buffer.data(index_gc + 2);

    // set up Cartesian P coordinates

    auto p_x = buffer.data(index_p);

    auto p_y = buffer.data(index_p + 1);

    auto p_z = buffer.data(index_p + 2);
    
    // set up Cartesian C coordinates

    const auto xyz = r_c.coordinates();

    const auto c_x = xyz[0];

    const auto c_y = xyz[1];

    const auto c_z = xyz[2];

    // compute R(PB) distances

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(gc_x, gc_y, gc_z, p_x, p_y, p_z : 64)
    for (size_t i = 0; i < nelems; i++)
    {
        const auto fact = a_exp + b_exps[i];
        
        const auto finv = 1.0 / (fact + c_exp);
        
        gc_x[i] = (fact * p_x[i] + c_exp * c_x) * finv - c_x;

        gc_y[i] = (fact * p_y[i] + c_exp * c_y) * finv - c_y;

        gc_z[i] = (fact * p_z[i] + c_exp * c_z) * finv - c_z;
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

void
comp_boys_args_with_rho(CSimdArray<double>& bf_data, const size_t index_args, CSimdArray<double>& buffer, const size_t index_ab, const double a_exp)
{
    // Set up exponents

    auto b_exps = buffer.data(0);

    // set up Boys function arguments

    auto bargs = bf_data.data(index_args);

    // set up R(AB) distances

    auto ab_x = buffer.data(index_ab);

    auto ab_y = buffer.data(index_ab + 1);

    auto ab_z = buffer.data(index_ab + 2);

    // compute Boys function arguments

    const auto nelems = buffer.number_of_active_elements();

#pragma omp simd aligned(bargs, ab_x, ab_y, ab_z, b_exps : 64)
    for (int i = 0; i < nelems; i++)
    {
        bargs[i] = a_exp * b_exps[i] * (ab_x[i] * ab_x[i] + ab_y[i] * ab_y[i] + ab_z[i] * ab_z[i]) / (a_exp + b_exps[i]);
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

auto
comp_distances_pc(      CSubMatrix&          buffer,
                  const size_t               index_pc,
                  const std::vector<double>& gcoords_x,
                  const std::vector<double>& gcoords_y,
                  const std::vector<double>& gcoords_z,
                  const double               p_x,
                  const double               p_y,
                  const double               p_z) -> void
{
    // set up number of grid points
    
    const auto npoints = buffer.number_of_columns();
    
    // set up R(PC) = P - C distances
    
    auto pc_x = &(buffer.data()[npoints * index_pc]);
    
    auto pc_y = &(buffer.data()[npoints * (index_pc + 1)]);
    
    auto pc_z = &(buffer.data()[npoints * (index_pc + 2)]);
    
    // set up grid point coordinates
    
    auto g_x = gcoords_x.data();
    
    auto g_y = gcoords_y.data();
    
    auto g_z = gcoords_z.data();
    
    // compute R(PC) distances on grid points
    
    #pragma omp simd
    for (size_t i = 0; i < npoints; i++)
    {
        pc_x[i] = p_x - g_x[i];
        
        pc_y[i] = p_y - g_y[i];
        
        pc_z[i] = p_z - g_z[i];
    }
}

auto
comp_boys_args(     CSubMatrix& buffer,
               const size_t     index_args,
               const size_t     index_pc,
               const double     factor) -> void
{
    // set up number of grid points
    
    const auto npoints = buffer.number_of_columns();
    
    // set up R(PC) = P - C distances
    
    auto pc_x = &(buffer.data()[npoints * index_pc]);
    
    auto pc_y = &(buffer.data()[npoints * (index_pc + 1)]);
    
    auto pc_z = &(buffer.data()[npoints * (index_pc + 2)]);
    
    // set up Boys function arguments
    
    auto bargs = &(buffer.data()[npoints * index_args]);
    
    #pragma omp simd
    for (size_t i = 0; i < npoints; i++)
    {
        bargs[i] = factor * (pc_x[i] * pc_x[i] + pc_y[i] * pc_y[i] + pc_z[i] * pc_z[i]);
    }
}

auto
reduce(CSubMatrix& buffer, const size_t index_contr, const size_t index_prim, const size_t ndims) -> void
{
    // set up number of grid points
    
    const auto npoints = buffer.number_of_columns();
    
    for (size_t i = 0; i < ndims; i++)
    {
        // set up buffers
        
        auto cvals = &(buffer.data()[npoints * (index_contr + i)]);
        
        auto pvals = &(buffer.data()[npoints * (index_prim + i)]);
        
        #pragma omp simd
        for (size_t j = 0; j < npoints; j++)
        {
            cvals[j] += pvals[j];
        }
    }
}

auto
distribute(      CDenseMatrix&             gmatrix,
           const CSubMatrix&               buffer,
           const size_t                    offset,
           const CDenseMatrix&             fmatrix,
           const std::vector<double>&      weights, 
           const std::map<size_t, size_t>& ao_mask,
           const std::vector<size_t>&      bra_indices,
           const std::vector<size_t>&      ket_indices,
           const int                       bra_angmom,
           const int                       ket_angmom,
           const size_t                    bra_igto,
           const size_t                    ket_igto,
           const bool                      bra_eq_ket) -> void
{
    // reference indexes on bra and ket sides

    const size_t refp = bra_indices[bra_igto + 1];
    
    const size_t refq = ket_indices[ket_igto + 1];
    
    // dimensions of bra and ket orbital indexes

    const auto adim = bra_indices[0];

    const auto bdim = ket_indices[0];
    
    // set up angular components

    const auto acomps = tensor::number_of_spherical_components(std::array<int, 1>{bra_angmom});
    
    const auto bcomps = tensor::number_of_spherical_components(std::array<int, 1>{ket_angmom});
    
    // set up pointer to G matrix
    
    auto gmat = gmatrix.values();
    
    // set up pointer to F matrix
    
    auto fmat = fmatrix.values();
    
    // set up numnber of AOs
    
    auto naos = gmatrix.getNumberOfColumns();
    
    // set up number of grid points
    
    auto npoints = fmatrix.getNumberOfColumns();
    
    for (int i = 0; i < acomps; i++)
    {
        const auto p = ao_mask.at(i * adim + refp);
        
        for (int j = 0; j < bcomps; j++)
        {
            const size_t idx_ij = offset + static_cast<size_t>(i * bcomps + j);
            
            const auto q = ao_mask.at(j * bdim + refq);
                
            if (bra_eq_ket)
            {
                for (size_t m = 0; m < npoints; m++)
                {
                    gmat[m * naos + p] += weights[m] * buffer.at({idx_ij, m}) * fmat[q * npoints + m];
                }
            }
            else
            {
                for (size_t m = 0; m < npoints; m++)
                {
                    gmat[m * naos + p] += weights[m] * buffer.at({idx_ij, m}) * fmat[q * npoints + m];
                    
                    gmat[m * naos + q] += weights[m] * buffer.at({idx_ij, m}) * fmat[p * npoints + m];
                }
            }
        }
    }
}



}  // namespace t2cfunc
