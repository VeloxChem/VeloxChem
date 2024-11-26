#ifndef T4CUtils_hpp
#define T4CUtils_hpp

#include <set>
#include <vector>

#include "CustomConstrains.hpp"
#include "Matrices.hpp"
#include "SimdArray.hpp"
#include "SphericalMomentum.hpp"
#include "SubMatrix.hpp"
#include "TensorComponents.hpp"

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

/// @brief Generates vector of global compresed indices for given vector of indices.
/// @param indices The vector of indices to generate global compresed indices.
/// @return The vector of global compresed indices.
template <Integral T>
inline auto
compresed_indices(const std::vector<T>& indices) -> std::vector<T>
{
    std::set<T> unique_indices(std::next(indices.cbegin()), indices.cend());

    std::vector<T> glob_indices;

    glob_indices.push_back(indices[0]);

    glob_indices.insert(glob_indices.end(), unique_indices.begin(), unique_indices.end());

    return glob_indices;
}

/// @brief Computes Cartesian Q center coordinates.
/// @param buffer The SIMD array containing factors data.
/// @param index_q The primary row index of Cartesian Q points coordinates in SIMD array.
/// @param index_c The primary row index of  Cartesian C points coordinates in SIMD array.
/// @param index_d The primary row index of  Cartesian D points coordinates in SIMD array.
auto comp_coordinates_q(CSimdArray<double>& buffer, const size_t index_q, const size_t index_c, const size_t index_d) -> void;

/// @brief Computes Cartesian W center coordinates.
/// @param buffer The SIMD array containing factors data.
/// @param index_w The primary row index of Cartesian Q points coordinates in SIMD array.
/// @param index_q The primary row index of Cartesian Q points coordinates in SIMD array.
/// @param r_p The Cartesian P point coordinates.
/// @param a_exp The exponent on A center.
/// @param b_exp The exponent on B center.
auto comp_coordinates_w(CSimdArray<double>&   buffer,
                        const size_t          index_w,
                        const size_t          index_q,
                        const TPoint<double>& r_p,
                        const double          a_exp,
                        const double          b_exp) -> void;

/// @brief Computes R(PQ) = P - Q distances.
/// @param buffer The SIMD array containing factors data.
/// @param index_pq The primary row index of R(PQ) distances in SIMD array.
/// @param index_q The primary row index of  Cartesian Q points coordinates in SIMD array.
/// @param r_p The Cartesian P point coordinates.
auto comp_distances_pq(CSimdArray<double>& buffer, const size_t index_pq, const size_t index_q, const TPoint<double>& r_p) -> void;

/// @brief Computes R(WQ) = W - Q distances.
/// @param buffer The SIMD array containing factors data.
/// @param index_wq The primary row index of R(WQ) distances in SIMD array.
/// @param index_w The primary row index of  Cartesian W points coordinates in SIMD array.
/// @param index_q The primary row index of  Cartesian Q points coordinates in SIMD array.
auto comp_distances_wq(CSimdArray<double>& buffer, const size_t index_wq, const size_t index_w, const size_t index_q) -> void;

/// @brief Computes R(WP) = W - P distances.
/// @param buffer The SIMD array containing factors data.
/// @param index_wp The primary row index of R(WP) distances in SIMD array.
/// @param index_w The primary row index of  Cartesian W points coordinates in SIMD array.
/// @param r_p The Cartesian P point coordinates.
auto comp_distances_wp(CSimdArray<double>& buffer, const size_t index_wp, const size_t index_w, const TPoint<double>& r_p) -> void;

/// @brief Computes R(CD) = C - D distances.
/// @param buffer The SIMD array containing factors data.
/// @param index_cd The primary row index of R(CD) distances in SIMD array.
/// @param index_c The primary row index of  Cartesian C points coordinates in SIMD array.
/// @param index_d The primary row index of  Cartesian D points coordinates in SIMD array.
auto comp_distances_cd(CSimdArray<double>& buffer, const size_t index_cd, const size_t index_c, const size_t index_d) -> void;

/// @brief Computes R(QD) = Q - D distances.
/// @param buffer The SIMD array containing factors data.
/// @param index_qd The primary row index of R(QD) distances in SIMD array.
/// @param index_q The primary row index of  Cartesian Q points coordinates in SIMD array.
/// @param index_d The primary row index of  Cartesian D points coordinates in SIMD array.
auto comp_distances_qd(CSimdArray<double>& buffer, const size_t index_qd, const size_t index_q, const size_t index_d) -> void;

/// @brief Computes R(QC) = Q - C distances.
/// @param buffer The SIMD array containing factors data.
/// @param index_qc The primary row index of R(QC) distances in SIMD array.
/// @param index_q The primary row index of  Cartesian Q points coordinates in SIMD array.
/// @param index_c The primary row index of  Cartesian C points coordinates in SIMD array.
auto comp_distances_qc(CSimdArray<double>& buffer, const size_t index_qc, const size_t index_q, const size_t index_c) -> void;

/// @brief Computes Boys function arguments.
/// @param bf_data The Boys function data.
/// @param index_args The primary row index of arguments in Boys function data.
/// @param buffer The SIMD array containing factors data.
/// @param index_pq The primary row index of R(PQ) distances in SIMD array.
/// @param a_exp The exponent on A center.
/// @param b_exp The exponent on B center.
auto comp_boys_args(CSimdArray<double>&       bf_data,
                    const size_t              index_args,
                    const CSimdArray<double>& buffer,
                    const size_t              index_pq,
                    const double              a_exp,
                    const double              b_exp) -> void;

/// @brief Computes Boys function arguments.
/// @param bf_data The Boys function data.
/// @param index_args The primary row index of arguments in Boys function data.
/// @param buffer The SIMD array containing factors data.
/// @param index_pq The primary row index of R(PQ) distances in SIMD array.
/// @param a_exp The exponent on A center.
/// @param b_exp The exponent on B center.
/// @param omega The range separation parameter.
auto comp_boys_args(CSimdArray<double>&       bf_data,
                    const size_t              index_args,
                    const CSimdArray<double>& buffer,
                    const size_t              index_pq,
                    const double              a_exp,
                    const double              b_exp,
                    const double              omega) -> void;

/// @brief Computes combined overlap factors.
/// @param buffer The SIMD array containing factors data.
/// @param index_ovl The primary row index of combined overlap in SIMD array.
/// @param index_ket_ovl The primary row index of ket overlap in SIMD array.
/// @param index_ket_norm The primary row index of ket overlap in SIMD array.
/// @param bra_ovl The overlap on bra side.
/// @param bra_norm The normalization factor on bra side.
/// @param a_exp The exponent on A center.
/// @param b_exp The exponent on B center.
auto comp_ovl_factors(CSimdArray<double>& buffer,
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
auto update_max_values(std::vector<double>& max_values, const CSimdArray<double>& buffer, const size_t index) -> void;

/// Updates maximum values.
/// - Parameter buffer: the buffer to store integrals.
/// - Parameter values : the buffer of integrals to store.
/// - Parameter offset: the offset in storage buffer.
auto store_values(CSimdArray<double>& buffer, const CSimdArray<double>& values, const int offset) -> void;

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
auto accumulate(CSubMatrix*                glob_matrix,
                const CSubMatrix*          loc_matrix,
                const std::vector<size_t>& bra_loc_indices,
                const std::vector<size_t>& ket_loc_indices,
                const std::vector<size_t>& bra_glob_indices,
                const std::vector<size_t>& ket_glob_indices,
                const int                  bra_comps,
                const int                  ket_comps,
                const bool                 ang_order) -> void;

/// @brief Adds batch of  local matrices to matrices container.
/// @param matrices The matrices container to be added to.
/// @param label The identifier of Fock matrix type.
/// @param mtype The matrix type.
/// @param suffix The suffix of local matrices identifiers.
/// @param adims The dimensions along center A.
/// @param bdims The dimensions along center B.
/// @param cdims The dimensions along center C.
/// @param ddims The dimensions along center D.
auto add_local_matrices(CMatrices&         matrices,
                        const std::string& label,
                        const mat_t        mtype,
                        const std::string& suffix,
                        const size_t       adims,
                        const size_t       bdims,
                        const size_t       cdims,
                        const size_t       ddims) -> void;

/// @brief Distributes buffer of integrals into local Fock matrix.
/// @param focks  The local Fock matrices.
/// @param suffix The suffix of Fock matrix identifier.
/// @param density  The pointer to AO density matrix.
/// @param label The label of Fock matrix.
/// @param exchange_factor The exchange scaling factor.
/// @param buffer The integrals buffer.
/// @param offset  The intgeral buffer offset.
/// @param a_indices The compressed basis function indexes on center A.
/// @param b_indices The compressed basis function indexes on center B.
/// @param c_indices The compressed basis function indexes on center C.
/// @param d_indices The compressed basis function indexes on center D.
/// @param a_angmom The angular momentum of integrals buffer on center A.
/// @param b_angmom The angular momentum of integrals buffer on center B.
/// @param c_angmom The angular momentum of integrals buffer on center C.
/// @param d_angmom Tthe angular momentum of integrals buffer on center D.
/// @param bra_igto The index of basis function on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of basis functions on ket side.
/// @param diagonal The flag to indicate diagonal subblock of matrix.
auto local_distribute(CMatrices&                       focks,
                      const std::string&               suffix,
                      const CMatrix*                   density,
                      const std::string&               label,
                      const double                     exchange_factor,
                      const CSimdArray<double>&        buffer,
                      const size_t                     offset,
                      const std::vector<size_t>&       a_indices,
                      const std::vector<size_t>&       b_indices,
                      const std::vector<size_t>&       c_indices,
                      const std::vector<size_t>&       d_indices,
                      const std::vector<size_t>&       a_loc_indices,
                      const std::vector<size_t>&       b_loc_indices,
                      const std::vector<size_t>&       c_loc_indices,
                      const std::vector<size_t>&       d_loc_indices,
                      const int                        a_angmom,
                      const int                        b_angmom,
                      const int                        c_angmom,
                      const int                        d_angmom,
                      const size_t                     bra_igto,
                      const std::pair<size_t, size_t>& ket_range,
                      const bool                       diagonal) -> void;

/// @brief Distributes buffer of integrals into local Fock matrix.
/// @param focks  The local Fock matrices.
/// @param suffix The suffix of Fock matrix identifier.
/// @param density  The pointer to AO density matrix.
/// @param label The label of Fock matrix.
/// @param exchange_factor The exchange scaling factor.
/// @param buffer The integrals buffer.
/// @param offset  The intgeral buffer offset.
/// @param a_indices The compressed basis function indexes on center A.
/// @param b_indices The compressed basis function indexes on center B.
/// @param c_indices The compressed basis function indexes on center C.
/// @param d_indices The compressed basis function indexes on center D.
/// @param a_angmom The angular momentum of integrals buffer on center A.
/// @param b_angmom The angular momentum of integrals buffer on center B.
/// @param c_angmom The angular momentum of integrals buffer on center C.
/// @param d_angmom Tthe angular momentum of integrals buffer on center D.
/// @param bra_igto The index of basis function on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of basis functions on ket side.
auto local_distribute_geom_ket_symm(CMatrices&                       focks,
                                    const std::string&               suffix,
                                    const CMatrix*                   density,
                                    const std::string&               label,
                                    const double                     exchange_factor,
                                    const CSimdArray<double>&        buffer,
                                    const size_t                     offset,
                                    const std::vector<size_t>&       a_indices,
                                    const std::vector<size_t>&       b_indices,
                                    const std::vector<size_t>&       c_indices,
                                    const std::vector<size_t>&       d_indices,
                                    const std::vector<size_t>&       a_loc_indices,
                                    const std::vector<size_t>&       b_loc_indices,
                                    const std::vector<size_t>&       c_loc_indices,
                                    const std::vector<size_t>&       d_loc_indices,
                                    const int                        a_angmom,
                                    const int                        b_angmom,
                                    const int                        c_angmom,
                                    const int                        d_angmom,
                                    const size_t                     bra_igto,
                                    const std::pair<size_t, size_t>& ket_range) -> void;


/// @brief Distributes buffer of integrals into local Fock matrix.
/// @param focks  The local Fock matrices.
/// @param suffix The suffix of Fock matrix identifier.
/// @param density  The pointer to AO density matrix.
/// @param label The label of Fock matrix.
/// @param exchange_factor The exchange scaling factor.
/// @param buffer The integrals buffer.
/// @param offset  The intgeral buffer offset.
/// @param a_indices The compressed basis function indexes on center A.
/// @param b_indices The compressed basis function indexes on center B.
/// @param c_indices The compressed basis function indexes on center C.
/// @param d_indices The compressed basis function indexes on center D.
/// @param a_angmom The angular momentum of integrals buffer on center A.
/// @param b_angmom The angular momentum of integrals buffer on center B.
/// @param c_angmom The angular momentum of integrals buffer on center C.
/// @param d_angmom Tthe angular momentum of integrals buffer on center D.
/// @param bra_igto The index of basis function on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of basis functions on ket side.
auto local_distribute_geom_ket_no_symm(CMatrices&                       focks,
                                       const std::string&               suffix,
                                       const CMatrix*                   density,
                                       const std::string&               label,
                                       const double                     exchange_factor,
                                       const CSimdArray<double>&        buffer,
                                       const size_t                     offset,
                                       const std::vector<size_t>&       a_indices,
                                       const std::vector<size_t>&       b_indices,
                                       const std::vector<size_t>&       c_indices,
                                       const std::vector<size_t>&       d_indices,
                                       const std::vector<size_t>&       a_loc_indices,
                                       const std::vector<size_t>&       b_loc_indices,
                                       const std::vector<size_t>&       c_loc_indices,
                                       const std::vector<size_t>&       d_loc_indices,
                                       const int                        a_angmom,
                                       const int                        b_angmom,
                                       const int                        c_angmom,
                                       const int                        d_angmom,
                                       const size_t                     bra_igto,
                                       const std::pair<size_t, size_t>& ket_range) -> void;

/// Transforms Cartesian integrals buffer to half-transformed integrals buffer.
/// - Parameter sbuffer: the spherical  integrals array.
/// - Parameter cbuffer: the Cartesian integrals array.
template <int N, int M>
inline auto
ket_transform(CSimdArray<double>&       sbuffer,
              const size_t              sposition,
              const CSimdArray<double>& cbuffer,
              const size_t              cposition,
              const int                 a_angmom,
              const int                 b_angmom) -> void
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
            const auto cart_off = (i * b_cart_comps + j) * c_cart_comps * d_cart_comps;

            const auto spher_off = (i * b_cart_comps + j) * c_spher_comps * d_spher_comps;

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
}

/// Transforms half-transformed integrals buffer to spherical integrals buffer.
/// - Parameter sbuffer: the spherical  integrals array.
/// - Parameter cbuffer: the Cartesian integrals array.
template <int N, int M>
inline auto
bra_transform(CSimdArray<double>&       sbuffer,
              const size_t              sposition,
              const CSimdArray<double>& cbuffer,
              const size_t              cposition,
              const int                 c_angmom,
              const int                 d_angmom) -> void
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
            const auto ket_off = i * d_spher_comps + j;

            for (int k = 0; k < a_spher_comps; k++)
            {
                for (int l = 0; l < b_spher_comps; l++)
                {
                    auto dst_ptr = sbuffer.data(sposition + (k * b_spher_comps + l) * c_spher_comps * d_spher_comps + ket_off);

                    for (const auto& [a_idx, afact] : spher_mom::transformation_factors<N>(k))
                    {
                        for (const auto& [b_idx, bfact] : spher_mom::transformation_factors<M>(l))
                        {
                            auto src_ptr = cbuffer.data(cposition + (a_idx * b_cart_comps + b_idx) * c_spher_comps * d_spher_comps + ket_off);

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
