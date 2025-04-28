#ifndef T2CTransform_hpp
#define T2CTransform_hpp

#include <algorithm>
#include <array>
#include <ranges>

#include "CustomViews.hpp"
#include "SphericalMomentum.hpp"
#include "TensorComponents.hpp"
#include "SubMatrix.hpp"

namespace t2cfunc {  // t2cfunc namespace

/// @brief Transforms Cartesian integrals buffer to spherical integrals buffer.
/// @tparam N The order of angular momentum tensor on bra side.
/// @tparam M The order of angular momentum tensor on ket side.
/// @param sbuffer The spherical integrals buffer.
/// @param cbuffer The Cartesian integrals array.
template <int N, int M>
inline auto
transform(CSimdArray<double>& sbuffer, const CSimdArray<double>& cbuffer) -> void
{
    const auto ndims = sbuffer.number_of_active_elements();

    const auto bra_spher_comps = tensor::number_of_spherical_components(std::array<int, 1>{N});

    const auto ket_spher_comps = tensor::number_of_spherical_components(std::array<int, 1>{M});

    const auto bra_cart_comps = tensor::number_of_cartesian_components(std::array<int, 1>{N});

    const auto ket_cart_comps = tensor::number_of_cartesian_components(std::array<int, 1>{M});

    const auto nblocks = sbuffer.number_of_rows() / (bra_spher_comps * ket_spher_comps);

    std::ranges::for_each(std::views::iota(size_t{0}, nblocks), [&](const auto n) {
        const auto cart_off  = n * bra_cart_comps * ket_cart_comps;
        const auto spher_off = n * bra_spher_comps * ket_spher_comps;
        std::ranges::for_each(views::rectangular(bra_spher_comps, ket_spher_comps), [&](const auto& index) {
            const auto [i, j] = index;
            auto dst_ptr      = sbuffer.data(spher_off + i * ket_spher_comps + j);
            for (const auto& [bra_idx, bra_fact] : spher_mom::transformation_factors<N>(i))
            {
                for (const auto& [ket_idx, ket_fact] : spher_mom::transformation_factors<M>(j))
                {
                    auto         src_ptr = cbuffer.data(cart_off + bra_idx * ket_cart_comps + ket_idx);
                    const double fact    = bra_fact * ket_fact;
#pragma omp simd aligned(dst_ptr, src_ptr : 64)
                    for (size_t k = 0; k < ndims; k++)
                    {
                        dst_ptr[k] += fact * src_ptr[k];
                    }
                }
            }
        });
    });
}

/// @brief Transforms Cartesian integrals buffer to spherical integrals buffer.
/// @tparam N The order of angular momentum tensor on bra side.
/// @tparam M The order of angular momentum tensor on ket side.
/// @param sbuffer The spherical integrals buffer.
/// @param cbuffer The Cartesian integrals array.
template <int N, int M>
inline auto
transform(CSubMatrix& sbuffer, const CSubMatrix& cbuffer, const size_t position) -> void
{
    const auto ndims = sbuffer.number_of_columns();

    const auto bra_spher_comps = tensor::number_of_spherical_components(std::array<int, 1>{N});

    const auto ket_spher_comps = tensor::number_of_spherical_components(std::array<int, 1>{M});

    const auto bra_cart_comps = tensor::number_of_cartesian_components(std::array<int, 1>{N});

    const auto ket_cart_comps = tensor::number_of_cartesian_components(std::array<int, 1>{M});

    const auto nblocks = sbuffer.number_of_rows() / (bra_spher_comps * ket_spher_comps);

    std::ranges::for_each(std::views::iota(size_t{0}, nblocks), [&](const auto n) {
        const auto cart_off  = position + n * bra_cart_comps * ket_cart_comps;
        const auto spher_off = n * bra_spher_comps * ket_spher_comps;
        std::ranges::for_each(views::rectangular(bra_spher_comps, ket_spher_comps), [&](const auto& index) {
            const auto [i, j] = index;
            auto dst_ptr      = &(sbuffer.data()[ndims * (spher_off + i * ket_spher_comps + j)]);
            for (const auto& [bra_idx, bra_fact] : spher_mom::transformation_factors<N>(i))
            {
                for (const auto& [ket_idx, ket_fact] : spher_mom::transformation_factors<M>(j))
                {
                    auto         src_ptr = &(cbuffer.data()[ndims * (cart_off + bra_idx * ket_cart_comps + ket_idx)]);
                    const double fact    = bra_fact * ket_fact;
#pragma omp simd 
                    for (size_t k = 0; k < ndims; k++)
                    {
                        dst_ptr[k] += fact * src_ptr[k];
                    }
                }
            }
        });
    });
}


}  // namespace t2cfunc

#endif /* T2CTransform_hpp */
