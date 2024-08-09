#include "T4CUtils.hpp"

#include <cmath>
#include <set>

#include "MathConst.hpp"

#include <iostream>

namespace t4cfunc {  // t2cfunc namespace

auto
masked_indices(const std::vector<int>& indices) -> std::vector<int>
{
    std::vector<int> loc_indices;
    
    std::set<int> unique_indices(std::next(indices.cbegin()), indices.cend());
    
    if (!unique_indices.empty())
    {
        loc_indices.push_back(static_cast<int>(unique_indices.size()));
        
        for (const auto index : indices)
        {
            int position = 0;
            
            for (const auto unique_index : unique_indices)
            {
                if (unique_index == index)
                {
                    loc_indices.push_back(position);
                    
                    break;
                }
                
                position++;
            }
        }
    }
    
    return loc_indices;
}

auto
comp_coordinates_q(double*       q_x,
                   double*       q_y,
                   double*       q_z,
                   const double* c_x,
                   const double* c_y,
                   const double* c_z,
                   const double* d_x,
                   const double* d_y,
                   const double* d_z,
                   const double* c_exps,
                   const double* d_exps,
                   const int     ndims) -> void
{
#pragma omp simd aligned(q_x, q_y, q_z, c_x, c_y, c_z, d_x, d_y, d_z, c_exps, d_exps : 64)
    for (int i = 0; i < ndims; i++)
    {
        double fact = 1.0 / (c_exps[i] + d_exps[i]);
        
        q_x[i] = fact * (c_x[i] * c_exps[i] + d_x[i] * d_exps[i]);
        
        q_y[i] = fact * (c_y[i] * c_exps[i] + d_y[i] * d_exps[i]);
        
        q_z[i] = fact * (c_z[i] * c_exps[i] + d_z[i] * d_exps[i]);
    }
}

auto
comp_coordinates_w(double*       w_x,
                   double*       w_y,
                   double*       w_z,
                   const double  p_x,
                   const double  p_y,
                   const double  p_z,
                   const double* q_x,
                   const double* q_y,
                   const double* q_z,
                   const double  a_exp,
                   const double  b_exp,
                   const double* c_exps,
                   const double* d_exps,
                   const int     ndims) -> void
{
#pragma omp simd aligned(w_x, w_y, w_z, q_x, q_y, q_z, c_exps, d_exps : 64)
    for (int i = 0; i < ndims; i++)
    {
        double ab_exp = a_exp + b_exp;
        
        double cd_exp = c_exps[i] + d_exps[i];
        
        double fact = 1.0 / (ab_exp + cd_exp);
        
        w_x[i] = fact * (p_x * ab_exp + q_x[i] * cd_exp);
        
        w_y[i] = fact * (p_y * ab_exp + q_y[i] * cd_exp);
        
        w_z[i] = fact * (p_z * ab_exp + q_z[i] * cd_exp);
    }
}

auto
comp_distances_pq(double*       pq_x,
                  double*       pq_y,
                  double*       pq_z,
                  const double  p_x,
                  const double  p_y,
                  const double  p_z,
                  const double* q_x,
                  const double* q_y,
                  const double* q_z,
                  const int     ndims) -> void
{
#pragma omp simd aligned(pq_x, pq_y, pq_z, q_x, q_y, q_z : 64)
    for (int i = 0; i < ndims; i++)
    {
        pq_x[i] = p_x - q_x[i];
        
        pq_y[i] = p_y - q_y[i];
        
        pq_z[i] = p_z - q_z[i];
    }
}

auto comp_distances_wq(double*       wq_x,
                       double*       wq_y,
                       double*       wq_z,
                       const double* w_x,
                       const double* w_y,
                       const double* w_z,
                       const double* q_x,
                       const double* q_y,
                       const double* q_z,
                       const int     ndims) -> void
{
#pragma omp simd aligned(wq_x, wq_y, wq_z, w_x, w_y, w_z, q_x, q_y, q_z : 64)
    for (int i = 0; i < ndims; i++)
    {
        wq_x[i] = w_x[i] - q_x[i];
        
        wq_y[i] = w_y[i] - q_y[i];
        
        wq_z[i] = w_z[i] - q_z[i];
    }
}

auto comp_distances_wp(double*       wp_x,
                       double*       wp_y,
                       double*       wp_z,
                       const double* w_x,
                       const double* w_y,
                       const double* w_z,
                       const double  p_x,
                       const double  p_y,
                       const double  p_z,
                       const int     ndims) -> void
{
#pragma omp simd aligned(wp_x, wp_y, wp_z, w_x, w_y, w_z : 64)
    for (int i = 0; i < ndims; i++)
    {
        wp_x[i] = w_x[i] - p_x;
        
        wp_y[i] = w_y[i] - p_y;
        
        wp_z[i] = w_z[i] - p_z;
    }
}

auto comp_distances_qd(double*       qd_x,
                       double*       qd_y,
                       double*       qd_z,
                       const double* q_x,
                       const double* q_y,
                       const double* q_z,
                       const double* d_x,
                       const double* d_y,
                       const double* d_z,
                       const int     ndims) -> void
{
#pragma omp simd aligned(qd_x, qd_y, qd_z, q_x, q_y, q_z, d_x, d_y, d_z : 64)
    for (int i = 0; i < ndims; i++)
    {
        qd_x[i] = q_x[i] - d_x[i];
        
        qd_y[i] = q_y[i] - d_y[i];
        
        qd_z[i] = q_z[i] - d_z[i];
    }
}

auto comp_distances_qc(double*       qc_x,
                       double*       qc_y,
                       double*       qc_z,
                       const double* q_x,
                       const double* q_y,
                       const double* q_z,
                       const double* c_x,
                       const double* c_y,
                       const double* c_z,
                       const int     ndims) -> void
{
#pragma omp simd aligned(qc_x, qc_y, qc_z, q_x, q_y, q_z, c_x, c_y, c_z : 64)
    for (int i = 0; i < ndims; i++)
    {
        qc_x[i] = q_x[i] - c_x[i];
        
        qc_y[i] = q_y[i] - c_y[i];
        
        qc_z[i] = q_z[i] - c_z[i];
    }
}

auto comp_distances_cd(double*       cd_x,
                       double*       cd_y,
                       double*       cd_z,
                       const double* c_x,
                       const double* c_y,
                       const double* c_z,
                       const double* d_x,
                       const double* d_y,
                       const double* d_z,
                       const int     ndims) -> void
{
#pragma omp simd aligned(cd_x, cd_y, cd_z, c_x, c_y, c_z, d_x, d_y, d_z : 64)
    for (int i = 0; i < ndims; i++)
    {
        cd_x[i] = c_x[i] - d_x[i];
        
        cd_y[i] = c_y[i] - d_y[i];
        
        cd_z[i] = c_z[i] - d_z[i];
    }
}

auto
comp_boys_args(CSimdArray<double>& bf_args,
               const double*       pq_x,
               const double*       pq_y,
               const double*       pq_z,
               const double        a_exp,
               const double        b_exp,
               const double*       c_exps,
               const double*       d_exps) -> void
{
    const auto ndims = bf_args.number_of_columns();

    auto bargs = bf_args[0];

#pragma omp simd aligned(bargs, pq_x, pq_y, pq_z, c_exps, d_exps : 64)
    for (int i = 0; i < ndims; i++)
    {
        double ab_exp = a_exp + b_exp;
        
        double cd_exp = c_exps[i] + d_exps[i];
        
        bargs[i] = ab_exp * cd_exp * (pq_x[i] * pq_x[i] + pq_y[i] * pq_y[i] + pq_z[i] * pq_z[i]) / (ab_exp + cd_exp);
    }
}

auto
comp_ovl_factors(CSimdArray<double>& fss_abcd,
                      const double        bra_ovl,
                      const double*       ket_ovls,
                      const double        bra_norm,
                      const double*       ket_norms,
                      const double        a_exp,
                      const double        b_exp,
                      const double*       c_exps,
                      const double*       d_exps) -> void
{
    const auto ndims = fss_abcd.number_of_columns();

    auto fss = fss_abcd[0];
    
    const auto invfpi = 1.0 / mathconst::pi_value();
    
#pragma omp simd aligned(fss, ket_ovls, ket_norms, c_exps, d_exps : 64)
    for (int i = 0; i < ndims; i++)
    {
        double ab_exp = a_exp + b_exp;
        
        double cd_exp = c_exps[i] + d_exps[i];
        
        fss[i] = 2.0 * bra_norm * ket_norms[i] * bra_ovl * ket_ovls[i] * std::sqrt(invfpi * ab_exp * cd_exp / (ab_exp + cd_exp));
    }
}

auto
update_max_values(      std::vector<double>& max_values,
                  const CSimdArray<double>&  buffer,
                  const int                  index) -> void
{
    for (int i = 0; i < buffer.number_of_rows(); i++)
    {
        if (max_values[index] < buffer[i][0]) max_values[index] = buffer[i][0];
    }
}

auto store_values(     CSimdArray<double>&   buffer,
                  const CSimdArray<double>&  values,
                  const int                  offset) -> void
{
    const auto ndims = values.number_of_columns();
    
    for (int i = 0; i < values.number_of_rows(); i++)
    {
        auto sdata = values[i];
        
        auto ddata = buffer[offset + i];
        
        #pragma omp simd aligned(sdata, ddata : 64)
        for (int j = 0; j < ndims; j++)
        {
            ddata[j] = sdata[j];
        }
    }
}

auto accumulate(      CSubMatrix* glob_matrix,
                const CSubMatrix* loc_matrix,
                const std::vector<int>& bra_loc_indices,
                const std::vector<int>& ket_loc_indices,
                const std::vector<int>& bra_glob_indices,
                const std::vector<int>& ket_glob_indices,
                const int               bra_comps,
                const int               ket_comps,
                const bool              ang_order) -> void
{
    for (int i = 0; i < bra_comps; i++)
    {
        const auto bra_loff = i * bra_loc_indices[0];
        
        const auto bra_goff = i * bra_glob_indices[0];
        
        for (int j = 0; j < ket_comps; j ++)
        {
            const auto ket_loff = j * ket_loc_indices[0];
            
            const auto ket_goff = j * ket_glob_indices[0];
            
            for (int k = 0; k < bra_loc_indices[0]; k++)
            {
                const auto kg = bra_goff + bra_glob_indices[k + 1];
                
                for (int l = 0; l < ket_loc_indices[0]; l++)
                {
                    const auto lg = ket_goff + ket_glob_indices[l + 1];
                    
                   // std::cout << "accum : (kg,lg) = (" << kg << "," << lg << ") (l, k) = " << bra_loff + k << " , " << ket_loff + l << std::endl;
                    
                    if (ang_order)
                    {
                        glob_matrix->operator[]({kg, lg}) += loc_matrix->at({bra_loff + k, ket_loff + l});
                    }
                    else
                    {
                        glob_matrix->operator[]({lg, kg}) += loc_matrix->at({bra_loff + k, ket_loff + l});
                    }
                }
            }
        }
    }
}

}  // namespace t4cfunc
