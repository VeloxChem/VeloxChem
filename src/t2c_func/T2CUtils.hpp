#ifndef T2CUtils_hpp
#define T2CUtils_hpp

#include <cstddef>
#include <utility>
#include <vector>
#include <map>

#include "Matrix.hpp"
#include "Point.hpp"
#include "SimdArray.hpp"
#include "SubMatrix.hpp"
#include "DenseMatrix.hpp"

namespace t2cfunc {  // t2cfunc namespace

/// @brief Computes R(AB) = A -B distances.
/// @param buffer The SIMD array containing R(AB) distances and Cartesian A, B coordinates.
/// @param index_ab The primary row index of R(AB) distances in SIMD array.
/// @param index_b  The primary row index of  Cartesian B points coordinates in SIMD array.
/// @param r_a The Cartesian A point coordinates.
auto comp_distances_ab(CSimdArray<double>& buffer, const size_t index_ab, const size_t index_b, const TPoint<double>& r_a) -> void;

/// @brief Computes P center coordinates.
/// @param buffer The SIMD array containing R(AB) distances and Cartesian A, B coordinates.
/// @param index_p The primary row index of  Cartesian P points coordinates in SIMD array.
/// @param index_b  The primary row index of  Cartesian B points coordinates in SIMD array.
/// @param r_a The Cartesian A point coordinates.
/// @param a_exp The exponent on A center.
auto comp_coordinates_p(CSimdArray<double>& buffer, const size_t index_p, const size_t index_b, const TPoint<double>& r_a, const double a_exp)
    -> void;

/// @brief Computes R(PB) = P - B distances.
/// @param buffer The SIMD array containing R(PB) distances.
/// @param index_pb The primary row index of R(PB) distances in SIMD array.
/// @param index_ab The primary row index of R(AB) distances in SIMD array.
/// @param a_exp The exponent on A center.
auto comp_distances_pb(CSimdArray<double>& buffer, const size_t index_pb, const size_t index_ab, const double a_exp) -> void;

/// @brief Computes R(PA) = P - A distances.
/// @param buffer The SIMD array containing R(PA) distances.
/// @param index_pa The primary row index of R(PA) distances in SIMD array.
/// @param index_ab The primary row index of R(AB) distances in SIMD array.
/// @param a_exp The exponent on A center.
auto comp_distances_pa(CSimdArray<double>& buffer, const size_t index_pa, const size_t index_ab, const double a_exp) -> void;

/// @brief Computes R(PB) = P - B distances.
/// @param buffer The SIMD array containing R(PB) distances.
/// @param index_pb The primary row index of R(PB) distances in SIMD array.
/// @param index_p The primary row index of  Cartesian P points coordinates in SIMD array.
/// @param index_b  The primary row index of  Cartesian B points coordinates in SIMD array.
auto comp_distances_pb_from_p(CSimdArray<double>& buffer, const size_t index_pb, const size_t index_p, const size_t index_b) -> void;

/// @brief Computes R(PA) = P - A distances.
/// @param buffer The SIMD array containing R(PA) distances.
/// @param index_pa The primary row index of R(PA) distances in SIMD array.
/// @param index_p The primary row index of  Cartesian P points coordinates in SIMD array.
/// @param r_a The Cartesian A point coordinates.
auto comp_distances_pa_from_p(CSimdArray<double>& buffer, const size_t index_pa, const size_t index_p, const TPoint<double>& r_a) -> void;

/// @brief Computes R(GB) = G - B distances.
/// @param buffer The SIMD array containing R(GB) distances.
/// @param index_gb The primary row index of R(GB) distances in SIMD array.
/// @param index_p The primary row index of  Cartesian P points coordinates in SIMD array.
/// @param index_b  The primary row index of  Cartesian B points coordinates in SIMD array.
/// @param r_c The Cartesian C point coordinates.
/// @param a_exp The exponent on A center.
/// @param c_exp The exponent on C center.
auto comp_distances_gb(CSimdArray<double>& buffer, const size_t index_gb, const size_t index_p, const size_t index_b, const TPoint<double>& r_c, const double a_exp, const double c_exp) -> void;

/// @brief Computes R(GA) = G - A distances.
/// @param buffer The SIMD array containing R(GA) distances.
/// @param index_ga The primary row index of R(GA) distances in SIMD array.
/// @param index_p The primary row index of  Cartesian P points coordinates in SIMD array.
/// @param r_a The Cartesian A point coordinates.
/// @param r_c The Cartesian C point coordinates.
/// @param a_exp The exponent on A center.
/// @param c_exp The exponent on C center.
auto comp_distances_ga(CSimdArray<double>& buffer, const size_t index_ga, const size_t index_p, const TPoint<double>& r_a, const TPoint<double>& r_c, const double a_exp, const double c_exp) -> void;

/// @brief Computes R(GC) = G -CA distances.
/// @param buffer The SIMD array containing R(GC) distances.
/// @param index_gc The primary row index of R(GC) distances in SIMD array.
/// @param index_p The primary row index of  Cartesian P points coordinates in SIMD array.
/// @param r_c The Cartesian C point coordinates.
/// @param a_exp The exponent on A center.
/// @param c_exp The exponent on C center.
auto comp_distances_gc(CSimdArray<double>& buffer, const size_t index_gc, const size_t index_p, const TPoint<double>& r_c, const double a_exp, const double c_exp) -> void;

/// @brief Computes R(PC) = P - C distances.
/// @param buffer The SIMD array containing R(PC) distances.
/// @param index_pc The primary row index of R(PC) distances in SIMD array.
/// @param index_p The primary row index of  Cartesian P points coordinates in SIMD array.
/// @param r_c The Cartesian C point coordinates.
auto comp_distances_pc(CSimdArray<double>& buffer, const size_t index_pc, const size_t index_p, const TPoint<double>& r_c) -> void;

/// Computes Boys function arguments.
/// @param bf_data The Boys function data.
/// @param index_args The primary row index of arguments in Boys function data.
/// @param buffer The SIMD array containing R(PC) distances.
/// @param index_pc The primary row index of R(PC) distances in SIMD array.
/// @param a_exp The exponent on A center.
void comp_boys_args(CSimdArray<double>& bf_data, const size_t index_args, CSimdArray<double>& buffer, const size_t index_pc, const double a_exp);

/// Computes Boys function arguments.
/// @param bf_data The Boys function data.
/// @param index_args The primary row index of arguments in Boys function data.
/// @param buffer The SIMD array containing R(PC) distances.
/// @param index_pc The primary row index of R(PC) distances in SIMD array.
/// @param a_exp The exponent on A center.
/// @param omega The range separation factor.
void comp_boys_args(CSimdArray<double>& bf_data,
                    const size_t        index_args,
                    CSimdArray<double>& buffer,
                    const size_t        index_pc,
                    const double        a_exp,
                    const double        omega);

/// Computes Boys function arguments.
/// @param bf_data The Boys function data.
/// @param index_args The primary row index of arguments in Boys function data.
/// @param buffer The SIMD array containing R(AB) distances.
/// @param index_ab The primary row index of R(AB) distances in SIMD array.
/// @param a_exp The exponent on A center.
void comp_boys_args_with_rho(CSimdArray<double>& bf_data,
                             const size_t index_args,
                             CSimdArray<double>& buffer,
                             const size_t index_ab,
                             const double a_exp);

/// @brief Reduces primitive array to contracted array.
/// @param cbuffer The contracted array.
/// @param pbuffer The primitive array.
/// @param position The starting position in primitive array.
/// @param ndims The dimensions of contracted row.
/// @param nblocks The number of blocks in primitive row.
auto reduce(CSimdArray<double>& cbuffer, CSimdArray<double>& pbuffer, const size_t position, const size_t ndims, const size_t nblocks) -> void;

/// @brief Reduces primitive array to contracted array.
/// @param cbuffer The contracted array.
/// @param cposition The contracted array position.
/// @param pbuffer The primitive array.
/// @param pposition The starting position in primitive array.
/// @param nrows The number of rows to contract. 
/// @param ndims The dimensions of contracted row.
/// @param nblocks The number of blocks in primitive row.
auto reduce(CSimdArray<double>& cbuffer,
            const size_t        cposition,
            CSimdArray<double>& pbuffer,
            const size_t pposition,
            const size_t nrows,
            const size_t ndims,
            const size_t nblocks) -> void;

/// @brief Reduces primitive array to contracted array.
/// @param cbuffer The contracted array.
/// @param pbuffer The primitive array.
/// @param position The starting position in primitive array.
/// @param factor The scaling factor of primitive array.
/// @param ndims The dimensions of contracted row.
/// @param nblocks The number of blocks in primitive row.
auto reduce(CSimdArray<double>& cbuffer,
            CSimdArray<double>& pbuffer,
            const size_t        position,
            const double        factor,
            const size_t        ndims,
            const size_t        nblocks) -> void;

/// @brief Reduces primitive array to contracted array.
/// @param cbuffer The contracted array.
/// @param pbuffer The primitive array.
/// @param position The starting position in primitive array.
/// @param facts The vector of scaling factors.
/// @param nfacts The number of scaling factors for single external center.
/// @param index The index of current scalling factor.
/// @param ndims The dimensions of contracted row.
/// @param nblocks The number of blocks in primitive row.
auto reduce(CSimdArray<double>&        cbuffer,
            CSimdArray<double>&        pbuffer,
            const size_t               position,
            const std::vector<double>& facts,
            const size_t               nfacts,
            const size_t               index,
            const size_t               ndims,
            const size_t               nblocks) -> void;

/// @brief Distributes buffer of integrals into given matrix.
/// @param matrix The pointer to matrix for storage of integrals.
/// @param buffer The integrals buffer.
/// @param indices The compressed contracted basis functions indexes.
/// @param bra_comp The angular component of integrals buffer on bra side.
/// @param ket_comp The angular component of integrals buffer on ket side.
/// @param bra_igto The index of basis function on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of basis functions on ket side.
auto distribute(CSubMatrix*                      matrix,
                const double*                    buffer,
                const std::vector<size_t>&       indices,
                const int                        bra_comp,
                const int                        ket_comp,
                const size_t                     bra_igto,
                const std::pair<size_t, size_t>& ket_range) -> void;

/// Distributes buffer of integrals into given matrix.
/// @param matrix The pointer to matrix for storage of integrals.
/// @param buffer The integrals buffer.
/// @param bra_indices The compressed contracted basis functions indexes on bra side.
/// @param ket_indices The compressed contracted basis functions indexes on ket side.
/// @param bra_comp The angular component of integrals buffer on bra side.
/// @param ket_comp The angular component of integrals buffer on ket side.
/// @param bra_igto The index of basis function on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of basis functions on ket side.
/// @param mat_type The matrix type.
auto distribute(CSubMatrix*                      matrix,
                const double*                    buffer,
                const std::vector<size_t>&       bra_indices,
                const std::vector<size_t>&       ket_indices,
                const int                        bra_comp,
                const int                        ket_comp,
                const size_t                     bra_igto,
                const std::pair<size_t, size_t>& ket_range,
                const mat_t                      mat_type) -> void;

/// Distributes buffer of integrals into given matrix.
/// @param matrix The pointer to matrix for storage of integrals.
/// @param buffer The integrals buffer.
/// @param bra_indices The compressed contracted basis functions indexes on bra side.
/// @param ket_indices The compressed contracted basis functions indexes on ket side.
/// @param bra_comp The angular component of integrals buffer on bra side.
/// @param ket_comp The angular component of integrals buffer on ket side.
/// @param bra_igto The index of basis function on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of basis functions on ket side.
/// @param ang_order The angular order of submatrix.
auto distribute(CSubMatrix*                      matrix,
                const double*                    buffer,
                const std::vector<size_t>&       bra_indices,
                const std::vector<size_t>&       ket_indices,
                const int                        bra_comp,
                const int                        ket_comp,
                const size_t                     bra_igto,
                const std::pair<size_t, size_t>& ket_range,
                const bool                       ang_order) -> void;

/// @brief Distributes buffer of integrals into given matrix.
/// @param matrix The pointer to matrix for storage of integrals.
/// @param buffer The integrals buffer.
/// @param offset The offset in integrals buffer.
/// @param indices The compressed contracted basis functions indexes.
/// @param bra_angmom The angular momentum of integrals buffer on bra and ket sides.
/// @param bra_igto The index of basis function on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of basis functions on ket side.
auto distribute(CSubMatrix*                      matrix,
                const CSimdArray<double>&        buffer,
                const size_t                     offset,
                const std::vector<size_t>&       indices,
                const int                        bra_angmom,
                const size_t                     bra_igto,
                const std::pair<size_t, size_t>& ket_range) -> void;

/// @brief Distributes buffer of integrals into given matrix.
/// @param matrix The pointer to matrix for storage of integrals.
/// @param buffer The integrals buffer.
/// @param offset The offset in integrals buffer.
/// @param bra_indices The compressed contracted basis functions indexes on bra side.
/// @param ket_indices The compressed contracted basis functions indexes on ket side.
/// @param bra_angmom The angular momentum of integrals buffer on bra side.
/// @param ket_angmom The angular momentum of integrals buffer on ket side.
/// @param bra_igto The index of basis function on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of basis functions on ket side.
/// @param mat_type The matrix type.
auto distribute(CSubMatrix*                      matrix,
                const CSimdArray<double>&        buffer,
                const size_t                     offset,
                const std::vector<size_t>&       bra_indices,
                const std::vector<size_t>&       ket_indices,
                const int                        bra_angmom,
                const int                        ket_angmom,
                const size_t                     bra_igto,
                const std::pair<size_t, size_t>& ket_range,
                const mat_t                      mat_type) -> void;

/// Distributes buffer of integrals into given matrix.
/// @param matrix The pointer to matrix for storage of integrals.
/// @param buffer The integrals buffer.
/// @param offset The offset in integrals buffer.
/// @param bra_indices The compressed contracted basis functions indexes on bra side.
/// @param ket_indices The compressed contracted basis functions indexes on ket side.
/// @param bra_angmom The angular momentum of integrals buffer on bra side.
/// @param ket_angmom The angular momentum of integrals buffer on ket side.
/// @param bra_igto The index of basis function on bra side.
/// @param ket_range The index of the range [ket_first, ket_last) of basis functions on ket side.
/// @param ang_order The angular order of submatrix.
auto distribute(CSubMatrix*                      matrix,
                const CSimdArray<double>&        buffer,
                const size_t                     offset,
                const std::vector<size_t>&       bra_indices,
                const std::vector<size_t>&       ket_indices,
                const int                        bra_angmom,
                const int                        ket_angmom,
                const size_t                     bra_igto,
                const std::pair<size_t, size_t>& ket_range,
                const bool                       ang_order) -> void;


/// @brief Computes R(PC) = P - C distances on grid.
/// @param buffer The buffer with factors.
/// @param index_pc The primary row index of R(PC) distances in buffer.
/// @param gcoords_x The Cartesian X coordinates of grid points.
/// @param gcoords_y The Cartesian Y coordinates of grid points.
/// @param gcoords_z The Cartesian Z coordinates of grid points.
/// @param p_x The Cartesian X coordinate of P point.
/// @param p_y The Cartesian Y coordinate of P point.
/// @param p_z The Cartesian Z coordinate of P point.
auto comp_distances_pc(      CSubMatrix&          buffer,
                       const size_t               index_pc,
                       const std::vector<double>& gcoords_x,
                       const std::vector<double>& gcoords_y,
                       const std::vector<double>& gcoords_z,
                       const double               p_x,
                       const double               p_y,
                       const double               p_z) -> void;

/// @brief Computes Boys function arguments  on grid.
/// @param buffer The buffer with factors.
/// @param index_args The primary row index of Boys function arguments in buffer.
/// @param index_pc The primary row index of R(PC) distances in buffer.
/// @param factor The Boys function argumnet scaling factor.
auto comp_boys_args(     CSubMatrix& buffer,
                    const size_t     index_args,
                    const size_t     index_pc,
                    const double     factor) -> void;

/// @brief Reduces primitive array to contracted array.
/// @param buffer The buffer with factors.
/// @param index_contr The primary row index  of contracted integrals in buffer.
/// @param index_prim The primary row index  of primitive integrals in buffer.
/// @param ndims The dimensions of contracted rows.
auto reduce(CSubMatrix& buffer, const size_t index_contr, const size_t index_prim, const size_t ndims) -> void;

/// Distributes buffer of integrals into given matrix on grid.
/// @param gmatrix The matrix for storage of contracted integrals on grid.
/// @param buffer The integrals buffer on grid.
/// @param offset The offset in integrals buffer.
/// @param fmatrix The matrix used for contracttion of integrals on grid.
/// @param weights The vector of weights. 
/// @param ao_mask The mask of local indices.
/// @param bra_indices The compressed contracted basis functions indexes on bra side.
/// @param ket_indices The compressed contracted basis functions indexes on ket side.
/// @param bra_angmom The angular momentum of integrals buffer on bra side.
/// @param ket_angmom The angular momentum of integrals buffer on ket side.
/// @param bra_igto The index of basis function on bra side.
/// @param ket_igto The index of basis function on ket side.
/// @param bra_eq_ket The flag for bra and ket basis function blocks being equal.
auto distribute(      CDenseMatrix&             gmatrix,
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
                const bool                      bra_eq_ket) -> void;

}  // namespace t2cfunc

#endif /* T2CUtils_hpp */
