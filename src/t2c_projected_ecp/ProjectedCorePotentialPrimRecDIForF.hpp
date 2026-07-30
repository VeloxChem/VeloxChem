#ifndef ProjectedCorePotentialPrimRecDIForF
#define ProjectedCorePotentialPrimRecDIForF

#include "SimdArray.hpp"
#include "Point.hpp"

namespace t2pecp { // t2pecp namespace

/// @brief Computes primitive [D|U_l|I]_F integrals with F projectors for set of data buffers.
/// @param pbuffer The primitive integrals buffer.
/// @param idx_di_f_0_0_0 The index of integral in primitive integrals buffer.
/// @param idx_si_f_0_0_0 The index of integral in primitive integrals buffer.
/// @param idx_pi_f_0_0_0 The index of integral in primitive integrals buffer.
/// @param idx_ph_d_0_0_1 The index of integral in primitive integrals buffer.
/// @param idx_pi_d_0_0_1 The index of integral in primitive integrals buffer.
/// @param idx_si_f_1_0_0 The index of integral in primitive integrals buffer.
/// @param idx_pi_f_1_0_0 The index of integral in primitive integrals buffer.
/// @param idx_si_p_1_0_1 The index of integral in primitive integrals buffer.
/// @param idx_pi_p_1_0_1 The index of integral in primitive integrals buffer.
/// @param idx_ph_s_1_1_1 The index of integral in primitive integrals buffer.
/// @param idx_pi_s_1_1_1 The index of integral in primitive integrals buffer.
/// @param p The special projector value.
/// @param idx_si_f_0_0_1 The index of integral in primitive integrals buffer.
/// @param idx_pi_f_0_0_1 The index of integral in primitive integrals buffer.
/// @param factors The primitive factors buffer.
/// @param idx_b The vector of Cartesian B points coordinates.
/// @param r_a The Cartesian A point coordinates.
/// @param a_exp The primitive basis function exponent on center A.
/// @param c_exp The core potential exponent on center C.
auto
comp_prim_projected_core_potential_di_f(CSimdArray<double>& pbuffer, 
                                        const size_t idx_di_f_0_0_0,
                                        const size_t idx_si_f_0_0_0,
                                        const size_t idx_pi_f_0_0_0,
                                        const size_t idx_ph_d_0_0_1,
                                        const size_t idx_pi_d_0_0_1,
                                        const size_t idx_si_f_1_0_0,
                                        const size_t idx_pi_f_1_0_0,
                                        const size_t idx_si_p_1_0_1,
                                        const size_t idx_pi_p_1_0_1,
                                        const size_t idx_ph_s_1_1_1,
                                        const size_t idx_pi_s_1_1_1,
                                        const int p,
                                        const size_t idx_si_f_0_0_1,
                                        const size_t idx_pi_f_0_0_1,
                                        const CSimdArray<double>& factors,
                                        const size_t idx_b,
                                        const TPoint<double>& r_a,
                                        const double a_exp,
                                        const double c_exp) -> void;
} // t2pecp namespace

#endif /* ProjectedCorePotentialPrimRecDIForF */
