#ifndef NuclearPotentialGeom010RecGP_hpp
#define NuclearPotentialGeom010RecGP_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"
#include "SimdTypes.hpp"
#include "MatrixType.hpp"

namespace geom_npotrec { // geom_npotrec namespace

/**
 Evaluates <G||P>  integrals for given pair of GTOs blocks.

 @param matrix_xthe pointer to matrix for storage of Cartesian integral component X.
 @param matrix_ythe pointer to matrix for storage of Cartesian integral component Y.
 @param matrix_zthe pointer to matrix for storage of Cartesian integral component Z.
 @param dipole the charge of external point.
 @param point the coordinates of external point.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto
compNuclearPotentialGeom010GP(      CSubMatrix* matrix_x,
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
 Evaluates block of primitive <G_XXXX||P_X> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XXXX_X(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXX||P_Y> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XXXX_Y(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXX||P_Z> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XXXX_Z(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXY||P_X> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XXXY_X(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXY||P_Y> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XXXY_Y(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXY||P_Z> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XXXY_Z(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXZ||P_X> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XXXZ_X(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXZ||P_Y> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XXXZ_Y(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXXZ||P_Z> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XXXZ_Z(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYY||P_X> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XXYY_X(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYY||P_Y> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XXYY_Y(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYY||P_Z> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XXYY_Z(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYZ||P_X> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XXYZ_X(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYZ||P_Y> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XXYZ_Y(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXYZ||P_Z> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XXYZ_Z(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXZZ||P_X> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XXZZ_X(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXZZ||P_Y> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XXZZ_Y(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XXZZ||P_Z> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XXZZ_Z(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYY||P_X> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XYYY_X(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYY||P_Y> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XYYY_Y(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYY||P_Z> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XYYY_Z(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYZ||P_X> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XYYZ_X(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYZ||P_Y> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XYYZ_Y(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYYZ||P_Z> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XYYZ_Z(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYZZ||P_X> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XYZZ_X(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYZZ||P_Y> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XYZZ_Y(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XYZZ||P_Z> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XYZZ_Z(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XZZZ||P_X> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XZZZ_X(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XZZZ||P_Y> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XZZZ_Y(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_XZZZ||P_Z> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_XZZZ_Z(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYY||P_X> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_YYYY_X(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYY||P_Y> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_YYYY_Y(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYY||P_Z> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_YYYY_Z(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYZ||P_X> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_YYYZ_X(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYZ||P_Y> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_YYYZ_Y(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYYZ||P_Z> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_YYYZ_Z(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYZZ||P_X> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_YYZZ_X(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYZZ||P_Y> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_YYZZ_Y(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YYZZ||P_Z> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_YYZZ_Z(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YZZZ||P_X> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_YZZZ_X(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YZZZ||P_Y> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_YZZZ_Y(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_YZZZ||P_Z> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_YZZZ_Z(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_ZZZZ||P_X> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_ZZZZ_X(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_ZZZZ||P_Y> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_ZZZZ_Y(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <G_ZZZZ||P_Z> integrals.

 @param buffer_x the partial integrals buffer.
 @param buffer_y the partial integrals buffer.
 @param buffer_z the partial integrals buffer.
 @param dipole the charge of external point.
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
compPrimitiveNuclearPotentialGeom010GP_ZZZZ_Z(      TDoubleArray& buffer_x,
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

#endif /* NuclearPotentialGeom010RecGP_hpp */
