#ifndef NuclearPotentialGeom010RecDP_hpp
#define NuclearPotentialGeom010RecDP_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"
#include "SimdTypes.hpp"
#include "MatrixType.hpp"

namespace geom_npotrec { // geom_npotrec namespace

/**
 Evaluates <D||P>  integrals for given pair of GTOs blocks.

 @param matrix_xthe pointer to matrix for storage of Cartesian integral component X.
 @param matrix_ythe pointer to matrix for storage of Cartesian integral component Y.
 @param matrix_zthe pointer to matrix for storage of Cartesian integral component Z.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto
compNuclearPotentialGeom010DP(      CSubMatrix* matrix_x,
                                    CSubMatrix* matrix_y,
                                    CSubMatrix* matrix_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                              const CGtoBlock&  bra_gto_block,
                              const CGtoBlock&  ket_gto_block,
                              const bool        ang_order,
                              const int64_t     bra_first,
                              const int64_t     bra_last) -> void;

/**
 Evaluates block of primitive <D_XX||P_X> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DP_XX_X(      TDoubleArray& buffer_x,
                                                  TDoubleArray& buffer_y,
                                                  TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                            const double        bra_exp,
                                            const double        bra_norm,
                                            const TPoint3D&     bra_coord,
                                            const TDoubleArray& ket_exps,
                                            const TDoubleArray& ket_norms,
                                            const TDoubleArray& ket_coords_x,
                                            const TDoubleArray& ket_coords_y,
                                            const TDoubleArray& ket_coords_z,
                                            const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XX||P_Y> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DP_XX_Y(      TDoubleArray& buffer_x,
                                                  TDoubleArray& buffer_y,
                                                  TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                            const double        bra_exp,
                                            const double        bra_norm,
                                            const TPoint3D&     bra_coord,
                                            const TDoubleArray& ket_exps,
                                            const TDoubleArray& ket_norms,
                                            const TDoubleArray& ket_coords_x,
                                            const TDoubleArray& ket_coords_y,
                                            const TDoubleArray& ket_coords_z,
                                            const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XX||P_Z> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DP_XX_Z(      TDoubleArray& buffer_x,
                                                  TDoubleArray& buffer_y,
                                                  TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                            const double        bra_exp,
                                            const double        bra_norm,
                                            const TPoint3D&     bra_coord,
                                            const TDoubleArray& ket_exps,
                                            const TDoubleArray& ket_norms,
                                            const TDoubleArray& ket_coords_x,
                                            const TDoubleArray& ket_coords_y,
                                            const TDoubleArray& ket_coords_z,
                                            const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XY||P_X> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DP_XY_X(      TDoubleArray& buffer_x,
                                                  TDoubleArray& buffer_y,
                                                  TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                            const double        bra_exp,
                                            const double        bra_norm,
                                            const TPoint3D&     bra_coord,
                                            const TDoubleArray& ket_exps,
                                            const TDoubleArray& ket_norms,
                                            const TDoubleArray& ket_coords_x,
                                            const TDoubleArray& ket_coords_y,
                                            const TDoubleArray& ket_coords_z,
                                            const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XY||P_Y> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DP_XY_Y(      TDoubleArray& buffer_x,
                                                  TDoubleArray& buffer_y,
                                                  TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                            const double        bra_exp,
                                            const double        bra_norm,
                                            const TPoint3D&     bra_coord,
                                            const TDoubleArray& ket_exps,
                                            const TDoubleArray& ket_norms,
                                            const TDoubleArray& ket_coords_x,
                                            const TDoubleArray& ket_coords_y,
                                            const TDoubleArray& ket_coords_z,
                                            const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XY||P_Z> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DP_XY_Z(      TDoubleArray& buffer_x,
                                                  TDoubleArray& buffer_y,
                                                  TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                            const double        bra_exp,
                                            const double        bra_norm,
                                            const TPoint3D&     bra_coord,
                                            const TDoubleArray& ket_exps,
                                            const TDoubleArray& ket_norms,
                                            const TDoubleArray& ket_coords_x,
                                            const TDoubleArray& ket_coords_y,
                                            const TDoubleArray& ket_coords_z,
                                            const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XZ||P_X> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DP_XZ_X(      TDoubleArray& buffer_x,
                                                  TDoubleArray& buffer_y,
                                                  TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                            const double        bra_exp,
                                            const double        bra_norm,
                                            const TPoint3D&     bra_coord,
                                            const TDoubleArray& ket_exps,
                                            const TDoubleArray& ket_norms,
                                            const TDoubleArray& ket_coords_x,
                                            const TDoubleArray& ket_coords_y,
                                            const TDoubleArray& ket_coords_z,
                                            const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XZ||P_Y> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DP_XZ_Y(      TDoubleArray& buffer_x,
                                                  TDoubleArray& buffer_y,
                                                  TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                            const double        bra_exp,
                                            const double        bra_norm,
                                            const TPoint3D&     bra_coord,
                                            const TDoubleArray& ket_exps,
                                            const TDoubleArray& ket_norms,
                                            const TDoubleArray& ket_coords_x,
                                            const TDoubleArray& ket_coords_y,
                                            const TDoubleArray& ket_coords_z,
                                            const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_XZ||P_Z> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DP_XZ_Z(      TDoubleArray& buffer_x,
                                                  TDoubleArray& buffer_y,
                                                  TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                            const double        bra_exp,
                                            const double        bra_norm,
                                            const TPoint3D&     bra_coord,
                                            const TDoubleArray& ket_exps,
                                            const TDoubleArray& ket_norms,
                                            const TDoubleArray& ket_coords_x,
                                            const TDoubleArray& ket_coords_y,
                                            const TDoubleArray& ket_coords_z,
                                            const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YY||P_X> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DP_YY_X(      TDoubleArray& buffer_x,
                                                  TDoubleArray& buffer_y,
                                                  TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                            const double        bra_exp,
                                            const double        bra_norm,
                                            const TPoint3D&     bra_coord,
                                            const TDoubleArray& ket_exps,
                                            const TDoubleArray& ket_norms,
                                            const TDoubleArray& ket_coords_x,
                                            const TDoubleArray& ket_coords_y,
                                            const TDoubleArray& ket_coords_z,
                                            const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YY||P_Y> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DP_YY_Y(      TDoubleArray& buffer_x,
                                                  TDoubleArray& buffer_y,
                                                  TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                            const double        bra_exp,
                                            const double        bra_norm,
                                            const TPoint3D&     bra_coord,
                                            const TDoubleArray& ket_exps,
                                            const TDoubleArray& ket_norms,
                                            const TDoubleArray& ket_coords_x,
                                            const TDoubleArray& ket_coords_y,
                                            const TDoubleArray& ket_coords_z,
                                            const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YY||P_Z> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DP_YY_Z(      TDoubleArray& buffer_x,
                                                  TDoubleArray& buffer_y,
                                                  TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                            const double        bra_exp,
                                            const double        bra_norm,
                                            const TPoint3D&     bra_coord,
                                            const TDoubleArray& ket_exps,
                                            const TDoubleArray& ket_norms,
                                            const TDoubleArray& ket_coords_x,
                                            const TDoubleArray& ket_coords_y,
                                            const TDoubleArray& ket_coords_z,
                                            const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YZ||P_X> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DP_YZ_X(      TDoubleArray& buffer_x,
                                                  TDoubleArray& buffer_y,
                                                  TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                            const double        bra_exp,
                                            const double        bra_norm,
                                            const TPoint3D&     bra_coord,
                                            const TDoubleArray& ket_exps,
                                            const TDoubleArray& ket_norms,
                                            const TDoubleArray& ket_coords_x,
                                            const TDoubleArray& ket_coords_y,
                                            const TDoubleArray& ket_coords_z,
                                            const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YZ||P_Y> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DP_YZ_Y(      TDoubleArray& buffer_x,
                                                  TDoubleArray& buffer_y,
                                                  TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                            const double        bra_exp,
                                            const double        bra_norm,
                                            const TPoint3D&     bra_coord,
                                            const TDoubleArray& ket_exps,
                                            const TDoubleArray& ket_norms,
                                            const TDoubleArray& ket_coords_x,
                                            const TDoubleArray& ket_coords_y,
                                            const TDoubleArray& ket_coords_z,
                                            const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_YZ||P_Z> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DP_YZ_Z(      TDoubleArray& buffer_x,
                                                  TDoubleArray& buffer_y,
                                                  TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                            const double        bra_exp,
                                            const double        bra_norm,
                                            const TPoint3D&     bra_coord,
                                            const TDoubleArray& ket_exps,
                                            const TDoubleArray& ket_norms,
                                            const TDoubleArray& ket_coords_x,
                                            const TDoubleArray& ket_coords_y,
                                            const TDoubleArray& ket_coords_z,
                                            const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_ZZ||P_X> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DP_ZZ_X(      TDoubleArray& buffer_x,
                                                  TDoubleArray& buffer_y,
                                                  TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                            const double        bra_exp,
                                            const double        bra_norm,
                                            const TPoint3D&     bra_coord,
                                            const TDoubleArray& ket_exps,
                                            const TDoubleArray& ket_norms,
                                            const TDoubleArray& ket_coords_x,
                                            const TDoubleArray& ket_coords_y,
                                            const TDoubleArray& ket_coords_z,
                                            const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_ZZ||P_Y> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DP_ZZ_Y(      TDoubleArray& buffer_x,
                                                  TDoubleArray& buffer_y,
                                                  TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                            const double        bra_exp,
                                            const double        bra_norm,
                                            const TPoint3D&     bra_coord,
                                            const TDoubleArray& ket_exps,
                                            const TDoubleArray& ket_norms,
                                            const TDoubleArray& ket_coords_x,
                                            const TDoubleArray& ket_coords_y,
                                            const TDoubleArray& ket_coords_z,
                                            const int64_t       ket_dim) -> void;

/**
 Evaluates block of primitive <D_ZZ||P_Z> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_exp the primitive exponent on bra side.
 @param bra_norm the primitive normalization factor on bra side.
 @param bra_coord the 3d coordinate of basis function on bra side.
 @param ket_exps the array of primitive exponents on ket side.
 @param ket_norms the array of primitive normalization factors on ket side.
 @param ket_coords_x the array of Cartesian X coordinates on ket side.
 @param ket_coords_y the array of Cartesian Y coordinates on ket side.
 @param ket_coords_z the array of Cartesian Z coordinates on ket side.
 @param ket_dim the end size of ket arrays.
*/
auto
compPrimitiveNuclearPotentialGeom010DP_ZZ_Z(      TDoubleArray& buffer_x,
                                                  TDoubleArray& buffer_y,
                                                  TDoubleArray& buffer_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                                            const double        bra_exp,
                                            const double        bra_norm,
                                            const TPoint3D&     bra_coord,
                                            const TDoubleArray& ket_exps,
                                            const TDoubleArray& ket_norms,
                                            const TDoubleArray& ket_coords_x,
                                            const TDoubleArray& ket_coords_y,
                                            const TDoubleArray& ket_coords_z,
                                            const int64_t       ket_dim) -> void;

} // geom_npotrec namespace

#endif /* NuclearPotentialGeom010RecDP_hpp */
