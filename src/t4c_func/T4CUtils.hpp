#ifndef T4CUtils_hpp
#define T4CUtils_hpp

#include <vector>
#include <set>

#include "SubMatrix.hpp"
#include "SimdArray.hpp"
#include "SphericalMomentum.hpp"
#include "TensorComponents.hpp"
#include "CustomConstrains.hpp"

namespace t4cfunc {  // t2cfunc namespace

/// @brief Generates vector of masked local indices for given vector of indices.
/// @param indices The vector of indices to generate masked indices.
/// @return The vector of masked indices.
template <Integral T>
inline auto
masked_indices(const std::vector<T>& indices) -> std::vector<T>
{
    std::vector<T> loc_indices;
        
    auto unique_indices = std::set<T>(std::next(indices.cbegin()), indices.cend());
        
    if (!unique_indices.empty())
    {
        loc_indices.push_back(static_cast<T>(unique_indices.size()));
            
        for (const auto index : indices)
        {
            auto position = T{0};
                
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

/// @brief Computes Cartesian Q center coordinates.
/// @param buffer The SIMD array containing factors data.
/// @param index_q The primary row index of Cartesian Q points coordinates in SIMD array.
/// @param index_c The primary row index of  Cartesian C points coordinates in SIMD array.
/// @param index_d The primary row index of  Cartesian D points coordinates in SIMD array.
auto comp_coordinates_q(CSimdArray<double>& buffer,
                        const size_t        index_q,
                        const size_t        index_c,
                        const size_t        index_d) -> void;

/// @brief Computes Cartesian W center coordinates.
/// @param buffer The SIMD array containing factors data.
/// @param index_w The primary row index of Cartesian Q points coordinates in SIMD array.
/// @param index_q The primary row index of Cartesian Q points coordinates in SIMD array.
/// @param r_p The Cartesian P point coordinates.
/// @param a_exp The exponent on A center.
/// @param b_exp The exponent on B center.
auto comp_coordinates_w(CSimdArray<double>& buffer,
                        const size_t        index_w,
                        const size_t        index_q,
                        const TPoint<double>& r_p,
                        const double          a_exp,
                        const double          b_exp) -> void;

/// @brief Computes R(PQ) = P - Q distances.
/// @param buffer The SIMD array containing factors data.
/// @param index_pq The primary row index of R(PQ) distances in SIMD array.
/// @param index_q The primary row index of  Cartesian Q points coordinates in SIMD array.
/// @param r_p The Cartesian P point coordinates.
auto comp_distances_pq(CSimdArray<double>& buffer,
                       const size_t        index_pq,
                       const size_t        index_q,
                       const TPoint<double>& r_p) -> void;

/// @brief Computes R(WQ) = W - Q distances.
/// @param buffer The SIMD array containing factors data.
/// @param index_wq The primary row index of R(WQ) distances in SIMD array.
/// @param index_w The primary row index of  Cartesian W points coordinates in SIMD array.
/// @param index_q The primary row index of  Cartesian Q points coordinates in SIMD array.
auto comp_distances_wq(CSimdArray<double>& buffer,
                       const size_t        index_wq,
                       const size_t        index_w,
                       const size_t        index_q) -> void;

/// @brief Computes R(WP) = W - P distances.
/// @param buffer The SIMD array containing factors data.
/// @param index_wp The primary row index of R(WP) distances in SIMD array.
/// @param index_w The primary row index of  Cartesian W points coordinates in SIMD array.
/// @param r_p The Cartesian P point coordinates.
auto comp_distances_wp(CSimdArray<double>& buffer,
                       const size_t        index_wp,
                       const size_t        index_w,
                       const TPoint<double>& r_p) -> void;

/// @brief Computes R(CD) = C - D distances.
/// @param buffer The SIMD array containing factors data.
/// @param index_cd The primary row index of R(CD) distances in SIMD array.
/// @param index_c The primary row index of  Cartesian C points coordinates in SIMD array.
/// @param index_d The primary row index of  Cartesian D points coordinates in SIMD array.
auto comp_distances_cd(CSimdArray<double>& buffer,
                       const size_t        index_cd,
                       const size_t        index_c,
                       const size_t        index_d) -> void;

/// @brief Computes R(QD) = Q - D distances.
/// @param buffer The SIMD array containing factors data.
/// @param index_qd The primary row index of R(QD) distances in SIMD array.
/// @param index_q The primary row index of  Cartesian Q points coordinates in SIMD array.
/// @param index_d The primary row index of  Cartesian D points coordinates in SIMD array.
auto comp_distances_qd(CSimdArray<double>& buffer,
                       const size_t        index_qd,
                       const size_t        index_q,
                       const size_t        index_d) -> void;

/// @brief Computes R(QC) = Q - C distances.
/// @param buffer The SIMD array containing factors data.
/// @param index_qc The primary row index of R(QC) distances in SIMD array.
/// @param index_q The primary row index of  Cartesian Q points coordinates in SIMD array.
/// @param index_c The primary row index of  Cartesian C points coordinates in SIMD array.
auto comp_distances_qc(CSimdArray<double>& buffer,
                       const size_t        index_qc,
                       const size_t        index_q,
                       const size_t        index_c) -> void;

/// @brief Computes Boys function arguments.
/// @param bf_data The Boys function data.
/// @param index_args The primary row index of arguments in Boys function data.
/// @param buffer The SIMD array containing factors data.
/// @param index_pq The primary row index of R(PQ) distances in SIMD array.
/// @param a_exp The exponent on A center.
/// @param b_exp The exponent on B center.
auto comp_boys_args(CSimdArray<double>& bf_data,
                    const size_t index_args,
                    const CSimdArray<double>& buffer,
                    const size_t        index_pq,
                    const double        a_exp,
                    const double        b_exp) -> void;

/// @brief Computes combined overlap factors.
/// @param buffer The SIMD array containing factors data.
/// @param index_ovl The primary row index of combined overlap in SIMD array.
/// @param index_ket_ovl The primary row index of ket overlap in SIMD array.
/// @param index_ket_norm The primary row index of ket overlap in SIMD array.
/// @param bra_ovl The overlap on bra side.
/// @param bra_norm The normalization factor on bra side.
/// @param a_exp The exponent on A center.
/// @param b_exp The exponent on B center.
auto comp_ovl_factors(const CSimdArray<double>& buffer,
                      const size_t        index_ovl,
                      const size_t        index_ket_ovl,
                      const size_t        index_ket_norm,
                      const double        bra_ovl,
                      const double        bra_norm,
                      const double        a_exp,
                      const double        b_exp) -> void;

/// Updates maximum values .
/// - Parameter max_values : the vector of maximum values of integral shells.
/// - Parameter buffer: the buffer of integrals.
/// - Parameter index: the index in vector of maximum values.
auto update_max_values(      std::vector<double>& max_values,
                       const CSimdArray<double>&  buffer,
                       const int                  index) -> void;

/// Updates maximum values.
/// - Parameter buffer: the buffer to store integrals.
/// - Parameter values : the buffer of integrals to store.
/// - Parameter offset: the offset in storage buffer.
auto store_values(    CSimdArray<double>&    buffer,
                  const CSimdArray<double>&  values,
                  const int                  offset) -> void;

/// Accumulates local matrix into global matrix using dual indices space.
/// - Parameter glob_matrix: the global matrix to add contributions.
/// - Parameter loc_matrix: the local matrix to be added.
/// - Parameter bra_loc_indices : the local indices on bra side.
/// - Parameter ket_loc_indices : the local indices on ket side.
/// - Parameter bra_glob_indices : the global indices on bra side.
/// - Parameter ket_glob_indices : the global indices on ket side.
/// - Parameter bra_comps : the number of angular components on bra side.
/// - Parameter ket_comps : the number of angular components on bra side.
/// - Parameter ang_order : the angular order of local/global matrix.
auto accumulate(      CSubMatrix* glob_matrix,
                const CSubMatrix* loc_matrix,
                const std::vector<int>& bra_loc_indices,
                const std::vector<int>& ket_loc_indices,
                const std::vector<int>& bra_glob_indices,
                const std::vector<int>& ket_glob_indices,
                const int               bra_comps,
                const int               ket_comps,
                const bool              ang_order) -> void;

/// Transforms Cartesian integrals buffer to half-transformed integrals buffer.
/// - Parameter sbuffer: the spherical  integrals array.
/// - Parameter cbuffer: the Cartesian integrals array.
template <int N, int M>
inline auto
ket_transform(CSimdArray<double>& sbuffer, const CSimdArray<double>& cbuffer, const int a_angmom, const int b_angmom) -> void
{
    const auto ndims = sbuffer.number_of_active_elements();

    const auto c_spher_comps = tensor::number_of_spherical_components(std::array<int, 1>{N});

    const auto d_spher_comps = tensor::number_of_spherical_components(std::array<int, 1>{M});
    
    const auto c_cart_comps = tensor::number_of_cartesian_components(std::array<int, 1>{N});
    
    const auto d_cart_comps = tensor::number_of_cartesian_components(std::array<int, 1>{M});

    const auto a_cart_comps = tensor::number_of_cartesian_components(std::array<int, 1>{a_angmom});

    const auto b_cart_comps = tensor::number_of_cartesian_components(std::array<int, 1>{b_angmom});

    for (int i = 0; i < a_cart_comps; i++)
    {
        for (int j = 0; j < b_cart_comps; j++)
        {
            const auto cart_off  = (i * b_cart_comps + j) * c_cart_comps * d_cart_comps;
            
            const auto spher_off = (i * b_cart_comps + j) * c_spher_comps * d_spher_comps;
            
            for (int k = 0; k < c_spher_comps; k++)
            {
                for (int l = 0; l < d_spher_comps; l++)
                {
                    auto dst_ptr = sbuffer.data(spher_off + k * d_spher_comps + l);
                    
                    for (const auto& [c_idx, cfact] : spher_mom::transformation_factors<N>(k))
                    {
                        for (const auto& [d_idx, dfact] : spher_mom::transformation_factors<M>(l))
                        {
                            auto src_ptr = cbuffer.data(cart_off + c_idx * d_cart_comps + d_idx);
                            
                            const double fact = cfact * dfact;

                            #pragma omp simd aligned(dst_ptr, src_ptr : 64)
                            for (size_t n = 0; n < ndims; n++)
                            {
                                dst_ptr[n] += fact * src_ptr[n];
                            }
                        }
                    }
                }
            }
        }
    }
}

/// Transforms half-transformed integrals buffer to spherical integrals buffer.
/// - Parameter sbuffer: the spherical  integrals array.
/// - Parameter cbuffer: the Cartesian integrals array.
template <int N, int M>
inline auto
bra_transform(CSimdArray<double>& sbuffer, const CSimdArray<double>& cbuffer, const int c_angmom, const int d_angmom) -> void
{
    const auto ndims = sbuffer.number_of_active_elements();

    const auto c_spher_comps = tensor::number_of_spherical_components(std::array<int, 1>{c_angmom});

    const auto d_spher_comps = tensor::number_of_spherical_components(std::array<int, 1>{d_angmom});

    const auto a_spher_comps = tensor::number_of_spherical_components(std::array<int, 1>{N});

    const auto b_spher_comps = tensor::number_of_spherical_components(std::array<int, 1>{M});
    
    const auto b_cart_comps = tensor::number_of_cartesian_components(std::array<int, 1>{M});

    for (int i = 0; i < c_spher_comps; i++)
    {
        for (int j = 0; j < d_spher_comps; j++)
        {
            const auto ket_off  = i * d_spher_comps + j;
            
            for (int k = 0; k < a_spher_comps; k++)
            {
                for (int l = 0; l < b_spher_comps; l++)
                {
                    auto dst_ptr = sbuffer.data((k * b_spher_comps + l) * c_spher_comps * d_spher_comps + ket_off);
                    
                    for (const auto& [a_idx, afact] : spher_mom::transformation_factors<N>(k))
                    {
                        for (const auto& [b_idx, bfact] : spher_mom::transformation_factors<M>(l))
                        {
                            auto src_ptr = cbuffer.data((a_idx * b_cart_comps + b_idx) * c_spher_comps * d_spher_comps + ket_off);
                            
                            const double fact = afact * bfact;

                            #pragma omp simd aligned(dst_ptr, src_ptr : 64)
                            for (size_t n = 0; n < ndims; n++)
                            {
                                dst_ptr[n] += fact * src_ptr[n];
                            }
                        }
                    }
                }
            }
        }
    }
}

}  // namespace t4cfunc

#endif /* T4CUtils_hpp */
