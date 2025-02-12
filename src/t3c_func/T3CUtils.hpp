#ifndef T3CUtils_hpp
#define T3CUtils_hpp

#include <set>
#include <vector>

#include "CustomConstrains.hpp"
#include "Matrices.hpp"
#include "SimdArray.hpp"
#include "SphericalMomentum.hpp"
#include "SubMatrix.hpp"
#include "TensorComponents.hpp"
#include "GtoBlock.hpp"

namespace t3cfunc {  // t3cfunc namespace

/// @brief Generates unique linear orbital indices vector for the given basis functions blocks.
/// @param gto_blocks The vector of basis functions blocks.
/// @return The vector of unique orbital indices.
auto unique_indices(const std::vector<CGtoBlock>& gto_blocks) -> std::vector<size_t>;

/// @brief Computes Cartesian W center coordinates.
/// @param buffer The SIMD array containing factors data.
/// @param index_w The primary row index of Cartesian Q points coordinates in SIMD array.
/// @param index_q The primary row index of Cartesian Q points coordinates in SIMD array.
/// @param r_a The Cartesian A point coordinates.
/// @param a_exp The exponent on A center.
auto comp_coordinates_w(CSimdArray<double>&   buffer,
                        const size_t          index_w,
                        const size_t          index_q,
                        const TPoint<double>& r_a,
                        const double          a_exp) -> void;

/// @brief Computes R(AQ) = A - Q distances.
/// @param buffer The SIMD array containing factors data.
/// @param index_aq The primary row index of R(AQ) distances in SIMD array.
/// @param index_q The primary row index of  Cartesian Q points coordinates in SIMD array.
/// @param r_a The Cartesian A point coordinates.
auto comp_distances_aq(CSimdArray<double>&   buffer,
                       const size_t          index_aq,
                       const size_t          index_q,
                       const TPoint<double>& r_a) -> void;

/// @brief Computes R(WA) = W - A distances.
/// @param buffer The SIMD array containing factors data.
/// @param index_wa The primary row index of R(WA) distances in SIMD array.
/// @param index_w The primary row index of  Cartesian W points coordinates in SIMD array.
/// @param r_a The Cartesian A point coordinates.
auto comp_distances_wa(CSimdArray<double>&   buffer,
                       const size_t          index_wa,
                       const size_t          index_w,
                       const TPoint<double>& r_a) -> void;

/// @brief Computes Boys function arguments.
/// @param bf_data The Boys function data.
/// @param index_args The primary row index of arguments in Boys function data.
/// @param buffer The SIMD array containing factors data.
/// @param index_aq The primary row index of R(AQ) distances in SIMD array.
/// @param a_exp The exponent on A center.
auto comp_boys_args(CSimdArray<double>&       bf_data,
                    const size_t              index_args,
                    const CSimdArray<double>& buffer,
                    const size_t              index_aq,
                    const double              a_exp) -> void;

/// @brief Computes combined overlap factors.
/// @param buffer The SIMD array containing factors data.
/// @param index_ovl The primary row index of combined overlap in SIMD array.
/// @param index_ket_ovl The primary row index of ket overlap in SIMD array.
/// @param index_ket_norm The primary row index of ket overlap in SIMD array.
/// @param a_norm The normalization factor on A center.
/// @param a_exp The exponent on A center.
auto comp_ovl_factors(CSimdArray<double>& buffer,
                      const size_t        index_ovl,
                      const size_t        index_ket_ovl,
                      const size_t        index_ket_norm,
                      const double        a_norm,
                      const double        a_exp) -> void;

/// Transforms half-transformed integrals buffer to spherical integrals buffer.
/// - Parameter sbuffer: the spherical  integrals array.
/// - Parameter cbuffer: the Cartesian integrals array.
template <int N>
inline auto
bra_transform(CSimdArray<double>&       sbuffer,
              const size_t              sposition,
              const CSimdArray<double>& cbuffer,
              const size_t              cposition,
              const int                 c_angmom,
              const int                 d_angmom) -> void
{
    const auto ndims = sbuffer.number_of_active_elements();

    const auto c_cart_comps = tensor::number_of_cartesian_components(std::array<int, 1>{c_angmom});

    const auto d_cart_comps = tensor::number_of_cartesian_components(std::array<int, 1>{d_angmom});

    const auto a_spher_comps = tensor::number_of_spherical_components(std::array<int, 1>{N});

    for (int i = 0; i < c_cart_comps; i++)
    {
        for (int j = 0; j < d_cart_comps; j++)
        {
            const auto ket_off = i * d_cart_comps + j;

            for (int k = 0; k < a_spher_comps; k++)
            {
                auto dst_ptr = sbuffer.data(sposition + k * c_cart_comps * d_cart_comps + ket_off);

                for (const auto& [a_idx, afact] : spher_mom::transformation_factors<N>(k))
                {
                    auto src_ptr = cbuffer.data(cposition + a_idx * c_cart_comps * d_cart_comps + ket_off);

                    const double fact = afact;

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

/// Transforms Cartesian integrals buffer to half-transformed integrals buffer.
/// - Parameter sbuffer: the spherical  integrals array.
/// - Parameter cbuffer: the Cartesian integrals array.
template <int N, int M>
inline auto
ket_transform(CSimdArray<double>&       sbuffer,
              const size_t              sposition,
              const CSimdArray<double>& cbuffer,
              const size_t              cposition,
              const int                 a_angmom) -> void
{
    const auto ndims = sbuffer.number_of_active_elements();

    const auto c_spher_comps = tensor::number_of_spherical_components(std::array<int, 1>{N});

    const auto d_spher_comps = tensor::number_of_spherical_components(std::array<int, 1>{M});

    const auto c_cart_comps = tensor::number_of_cartesian_components(std::array<int, 1>{N});

    const auto d_cart_comps = tensor::number_of_cartesian_components(std::array<int, 1>{M});

    const auto a_spher_comps = tensor::number_of_spherical_components(std::array<int, 1>{a_angmom});

    for (int i = 0; i < a_spher_comps; i++)
    {
        const auto cart_off = i * c_cart_comps * d_cart_comps;

        const auto spher_off = i * c_spher_comps * d_spher_comps;

        for (int k = 0; k < c_spher_comps; k++)
        {
            for (int l = 0; l < d_spher_comps; l++)
            {
                auto dst_ptr = sbuffer.data(sposition + spher_off + k * d_spher_comps + l);

                for (const auto& [c_idx, cfact] : spher_mom::transformation_factors<N>(k))
                {
                    for (const auto& [d_idx, dfact] : spher_mom::transformation_factors<M>(l))
                    {
                        auto src_ptr = cbuffer.data(cposition + cart_off + c_idx * d_cart_comps + d_idx);

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

}  // namespace t3cfunc

#endif /* T3CUtils_hpp */
