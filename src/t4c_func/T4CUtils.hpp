#ifndef T4CUtils_hpp
#define T4CUtils_hpp

#include <vector>

#include "SubMatrix.hpp"
#include "SimdArray.hpp"
#include "SphericalMomentum.hpp"
#include "TensorComponents.hpp"

namespace t4cfunc {  // t2cfunc namespace

/// Generates vector of masked local indices for given vector of indices.
/// - Parameter q_x: the vector of indices to generate masked indices.
auto masked_indices(const std::vector<int>& indices) -> std::vector<int>;

/// Computes Q center coordinates.
/// - Parameter q_x: the vector of Cartesian X  coordinates of center P.
/// - Parameter q_y: the vector of Cartesian Y  coordinates of center P.
/// - Parameter q_z: the vector of Cartesian Z  coordinates of center P.
/// - Parameter c_x: the vector of Cartesian X  coordinates of center C.
/// - Parameter c_y: the vector of Cartesian Y coordinates of center C.
/// - Parameter c_z: the vector of Cartesian Z  coordinates of center C.
/// - Parameter d_x: the vector of Cartesian X  coordinates of D centers.
/// - Parameter d_y: the vector of Cartesian Y  coordinates of D centers.
/// - Parameter d_z: the vector of Cartesian Z coordinates of BDcenters.
/// - Parameter c_exps: the vector of exponents on C center.
/// - Parameter d_exps: the vector of exponents on D center.
/// - Parameter ndims: the dimensions of vectors.
auto comp_coordinates_q(double*       q_x,
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
                        const int     ndims) -> void;

/// Computes W center coordinates.
/// - Parameter w_x: the vector of Cartesian X  coordinates of center W.
/// - Parameter w_y: the vector of Cartesian Y  coordinates of center W.
/// - Parameter w_z: the vector of Cartesian Z  coordinates of center W.
/// - Parameter p_x: the Cartesian X coordinates of center P.
/// - Parameter p_y: the Cartesian Y coordinates of center P.
/// - Parameter p_z: the  Cartesian Z coordinates of center P.
/// - Parameter q_x: the vector of Cartesian X coordinates of D centers.
/// - Parameter q_y: the vector of Cartesian Y coordinates of D centers.
/// - Parameter q_z: the vector of Cartesian Z coordinates of D centers.
/// - Parameter a_exp: the exponent on A center.
/// - Parameter b_exp: the exponent on B center.
/// - Parameter c_exps: the vector of exponents on C center.
/// - Parameter d_exps: the vector of exponents on D center.
/// - Parameter ndims: the dimensions of vectors.
auto comp_coordinates_w(double*       w_x,
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
                        const int     ndims) -> void;

/// Computes R(PQ) = P - Q distances.
/// - Parameter pq_x: the vector of Cartesian X  distances R(PQ) = P - Q.
/// - Parameter pq_y: the vector of Cartesian Y  distances R(PQ) = P - Q.
/// - Parameter pq_z: the vector of Cartesian Z  distances R(PQ) = P - Q.
/// - Parameter p_x: the Cartesian X  coordinates of center P.
/// - Parameter p_y: the Cartesian Y  coordinates of center P.
/// - Parameter p_z: the Cartesian Z  coordinates of center P.
/// - Parameter q_x: the vector of Cartesian X  coordinate of center Q.
/// - Parameter q_y: the vector of Cartesian Y coordinate of center Q.
/// - Parameter q_z: the vector of Cartesian Z  coordinate of center Q.
/// - Parameter ndims: the dimensions of vectors.
auto comp_distances_pq(double*       pq_x,
                       double*       pq_y,
                       double*       pq_z,
                       const double  p_x,
                       const double  p_y,
                       const double  p_z,
                       const double* q_x,
                       const double* q_y,
                       const double* q_z,
                       const int     ndims) -> void;

/// Computes R(WQ) = W - Q distances.
/// - Parameter wq_x: the vector of Cartesian X  distances R(WQ) = W - Q.
/// - Parameter wq_y: the vector of Cartesian Y  distances R(WQ) = W - Q.
/// - Parameter wq_z: the vector of Cartesian Z  distances R(WQ) = W - Q.
/// - Parameter w_x: the vector of Cartesian X  coordinate of center W.
/// - Parameter w_y: the vector of Cartesian Y  coordinate of center W.
/// - Parameter w_z: the vector of Cartesian Z coordinate of center W.
/// - Parameter q_x: the vector of Cartesian X  coordinate of center Q.
/// - Parameter q_y: the vector of Cartesian Y coordinate of center Q.
/// - Parameter q_z: the vector of Cartesian Z  coordinate of center Q.
/// - Parameter ndims: the dimensions of vectors.
auto comp_distances_wq(double*       wq_x,
                       double*       wq_y,
                       double*       wq_z,
                       const double* w_x,
                       const double* w_y,
                       const double* w_z,
                       const double* q_x,
                       const double* q_y,
                       const double* q_z,
                       const int     ndims) -> void;

/// Computes R(WP) = W - P distances.
/// - Parameter wp_x: the vector of Cartesian X  distances R(WP) = W - P.
/// - Parameter wp_y: the vector of Cartesian Y  distances R(WP) = W - P.
/// - Parameter wp_z: the vector of Cartesian Z  distances R(WP) = W - P.
/// - Parameter w_x: the vector of Cartesian X  coordinate of center W.
/// - Parameter w_y: the vector of Cartesian Y  coordinate of center W.
/// - Parameter w_z: the vector of Cartesian Z coordinate of center W.
/// - Parameter p_x: the Cartesian X  coordinate of center P.
/// - Parameter p_y: the Cartesian Y coordinate of center P.
/// - Parameter p_z: the Cartesian Z  coordinate of center P.
/// - Parameter ndims: the dimensions of vectors.
auto comp_distances_wp(double*       wp_x,
                       double*       wp_y,
                       double*       wp_z,
                       const double* w_x,
                       const double* w_y,
                       const double* w_z,
                       const double  p_x,
                       const double  p_y,
                       const double  p_z,
                       const int     ndims) -> void;

/// Computes R(CD) = C - D distances.
/// - Parameter cd_x: the vector of Cartesian X  distances R(CD) = C - D.
/// - Parameter cd_y: the vector of Cartesian Y  distances R(CD) = C - D.
/// - Parameter cd_z: the vector of Cartesian Z  distances R(CD) = C - D.
/// - Parameter c_x: the vector of Cartesian X  coordinate of center C.
/// - Parameter c_y: the vector of Cartesian Y coordinate of center C.
/// - Parameter c_z: the vector of Cartesian Z  coordinate of center C.
/// - Parameter d_x: the vector of Cartesian X  coordinate of center D.
/// - Parameter d_y: the vector of Cartesian Y coordinate of center D.
/// - Parameter d_z: the vector of Cartesian Z  coordinate of center D.
/// - Parameter ndims: the dimensions of vectors.
auto comp_distances_cd(double*       cd_x,
                       double*       cd_y,
                       double*       cd_z,
                       const double* c_x,
                       const double* c_y,
                       const double* c_z,
                       const double* d_x,
                       const double* d_y,
                       const double* d_z,
                       const int     ndims) -> void;

/// Computes R(QD) = Q - D distances.
/// - Parameter qd_x: the vector of Cartesian X  distances R(QD) = Q - D.
/// - Parameter qd_y: the vector of Cartesian Y  distances R(QD) = Q - D.
/// - Parameter qd_z: the vector of Cartesian Z  distances R(QD) = Q - D.
/// - Parameter q_x: the vector of Cartesian X  coordinate of center Q.
/// - Parameter q_y: the vector of Cartesian Y  coordinate of center Q.
/// - Parameter q_z: the vector of Cartesian Z coordinate of center Q.
/// - Parameter d_x: the Cartesian X  coordinate of center D.
/// - Parameter d_y: the Cartesian Y coordinate of center D.
/// - Parameter d_z: the Cartesian Z  coordinate of center D.
/// - Parameter ndims: the dimensions of vectors.
auto comp_distances_qd(double*       qd_x,
                       double*       qd_y,
                       double*       qd_z,
                       const double* q_x,
                       const double* q_y,
                       const double* q_z,
                       const double* d_x,
                       const double* d_y,
                       const double* d_z,
                       const int     ndims) -> void;

/// Computes R(QC) = Q - C distances.
/// - Parameter qc_x: the vector of Cartesian X  distances R(QC) = Q - C.
/// - Parameter qc_y: the vector of Cartesian Y  distances R(QC) = Q - C.
/// - Parameter qc_z: the vector of Cartesian Z  distances R(QC) = Q - C.
/// - Parameter q_x: the vector of Cartesian X  coordinate of center Q.
/// - Parameter q_y: the vector of Cartesian Y  coordinate of center Q.
/// - Parameter q_z: the vector of Cartesian Z coordinate of center Q.
/// - Parameter c_x: the Cartesian X  coordinate of center C.
/// - Parameter c_y: the Cartesian Y coordinate of center C.
/// - Parameter c_z: the Cartesian Z  coordinate of center C.
/// - Parameter ndims: the dimensions of vectors.
auto comp_distances_qc(double*       qc_x,
                       double*       qc_y,
                       double*       qc_z,
                       const double* q_x,
                       const double* q_y,
                       const double* q_z,
                       const double* c_x,
                       const double* c_y,
                       const double* c_z,
                       const int     ndims) -> void;

/// Computes Boys function arguments.
/// - Parameter bf_args: the vector Boys function arguments.
/// - Parameter pq_x: the vector of Cartesian X  distances R(PQ) = P - Q.
/// - Parameter pq_y: the vector of Cartesian Y  distances R(PQ) = P - Q.
/// - Parameter pq_z: the vector of Cartesian Z  distances R(PQ) = P - Q.
/// - Parameter a_exp: the exponent on A center.
/// - Parameter b_exp: the exponent on B center.
/// - Parameter c_exps: the vector of exponents on C center.
/// - Parameter d_exps: the vector of exponents on D center.
auto comp_boys_args(CSimdArray<double>& bf_args,
                    const double*       pq_x,
                    const double*       pq_y,
                    const double*       pq_z,
                    const double        a_exp,
                    const double        b_exp,
                    const double*       c_exps,
                    const double*       d_exps) -> void;

/// Computes combined overlap factors.
/// - Parameter fss_abcd: the vector of combined overlap factors.
/// - Parameter bra_ovl: the overlap on bra side.
/// - Parameter ket_ovls: the vector of overlaps on ket side.
/// - Parameter bra_norm: the normalization factor on bra side.
/// - Parameter ket_norms: the vector of normalization factors on ket side.
/// - Parameter a_exp: the exponent on A center.
/// - Parameter b_exp: the exponent on B center.
/// - Parameter c_exps: the vector of exponents on C center.
/// - Parameter d_exps: the vector of exponents on D center.
auto comp_ovl_factors(CSimdArray<double>& fss_abcd,
                      const double        bra_ovl,
                      const double*       ket_ovls,
                      const double        bra_norm,
                      const double*       ket_norms,
                      const double        a_exp,
                      const double        b_exp,
                      const double*       c_exps,
                      const double*       d_exps) -> void;

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
    const auto ndims = sbuffer.number_of_columns();

    const auto c_spher_comps = tensor::number_of_spherical_components(N);

    const auto d_spher_comps = tensor::number_of_spherical_components(M);
    
    const auto c_cart_comps = tensor::number_of_cartesian_components(N);
    
    const auto d_cart_comps = tensor::number_of_cartesian_components(M);

    const auto a_cart_comps = tensor::number_of_cartesian_components(a_angmom);

    const auto b_cart_comps = tensor::number_of_cartesian_components(b_angmom);

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
                    auto dst_ptr = sbuffer[spher_off + k * d_spher_comps + l];
                    
                    for (const auto& [c_idx, cfact] : spher_mom::transformation_factors<N>(k))
                    {
                        for (const auto& [d_idx, dfact] : spher_mom::transformation_factors<M>(l))
                        {
                            auto src_ptr = cbuffer[cart_off + c_idx * d_cart_comps + d_idx];
                            
                            const double fact = cfact * dfact;

                            #pragma omp simd aligned(dst_ptr, src_ptr : 64)
                            for (int n = 0; n < ndims; n++)
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
    const auto ndims = sbuffer.number_of_columns();

    const auto c_spher_comps = tensor::number_of_spherical_components(c_angmom);

    const auto d_spher_comps = tensor::number_of_spherical_components(d_angmom);

    const auto a_spher_comps = tensor::number_of_spherical_components(N);

    const auto b_spher_comps = tensor::number_of_spherical_components(M);
    
    const auto b_cart_comps = tensor::number_of_cartesian_components(M);

    for (int i = 0; i < c_spher_comps; i++)
    {
        for (int j = 0; j < d_spher_comps; j++)
        {
            const auto ket_off  = i * d_spher_comps + j;
            
            for (int k = 0; k < a_spher_comps; k++)
            {
                for (int l = 0; l < b_spher_comps; l++)
                {
                    auto dst_ptr = sbuffer[(k * b_spher_comps + l) * c_spher_comps * d_spher_comps + ket_off];
                    
                    for (const auto& [a_idx, afact] : spher_mom::transformation_factors<N>(k))
                    {
                        for (const auto& [b_idx, bfact] : spher_mom::transformation_factors<M>(l))
                        {
                            auto src_ptr = cbuffer[(a_idx * b_cart_comps + b_idx) * c_spher_comps * d_spher_comps + ket_off];
                            
                            const double fact = afact * bfact;

                            #pragma omp simd aligned(dst_ptr, src_ptr : 64)
                            for (int n = 0; n < ndims; n++)
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
