#ifndef NuclearPotentialGeom010RecPP_hpp
#define NuclearPotentialGeom010RecPP_hpp

#include <cstdint>

#include "GtoBlock.hpp"
#include "SubMatrix.hpp"
#include "SimdTypes.hpp"
#include "MatrixType.hpp"

namespace geom_npotrec { // geom_npotrec namespace

/**
 Evaluates <P||P>  integrals for given GTOs block.

 @param matrix_xthe pointer to matrix for storage of Cartesian integral component X.
 @param matrix_ythe pointer to matrix for storage of Cartesian integral component Y.
 @param matrix_zthe pointer to matrix for storage of Cartesian integral component Z.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param gto_block the GTOs block.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
*/
auto
compNuclearPotentialGeom010PP(      CSubMatrix* matrix_x,
                                    CSubMatrix* matrix_y,
                                    CSubMatrix* matrix_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                              const CGtoBlock&  gto_block,
                              const int64_t     bra_first,
                              const int64_t     bra_last) -> void;

/**
 Evaluates <P||P>  integrals for given pair of GTOs blocks.

 @param matrix_xthe pointer to matrix for storage of Cartesian integral component X.
 @param matrix_ythe pointer to matrix for storage of Cartesian integral component Y.
 @param matrix_zthe pointer to matrix for storage of Cartesian integral component Z.
 @param dipole the dipole of external point.
 @param point the coordinates of external point.
 @param bra_gto_block the GTOs block on bra side.
 @param ket_gto_block the GTOs block on ket side.
 @param bra_first the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param bra_last the index of the range [bra_first, bra_last) of GTOs on bra side.
 @param mat_type the matrix type.
*/
auto
compNuclearPotentialGeom010PP(      CSubMatrix* matrix_x,
                                    CSubMatrix* matrix_y,
                                    CSubMatrix* matrix_z,
                              const TPoint3D& dipole,
                              const TPoint3D& point,
                              const CGtoBlock&  bra_gto_block,
                              const CGtoBlock&  ket_gto_block,
                              const int64_t     bra_first,
                              const int64_t     bra_last,
                              const mat_t       mat_type) -> void;

/**
 Evaluates block of primitive <P_X||P_X> integrals.

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
compPrimitiveNuclearPotentialGeom010PP_X_X(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <P_X||P_Y> integrals.

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
compPrimitiveNuclearPotentialGeom010PP_X_Y(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <P_X||P_Z> integrals.

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
compPrimitiveNuclearPotentialGeom010PP_X_Z(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <P_Y||P_X> integrals.

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
compPrimitiveNuclearPotentialGeom010PP_Y_X(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <P_Y||P_Y> integrals.

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
compPrimitiveNuclearPotentialGeom010PP_Y_Y(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <P_Y||P_Z> integrals.

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
compPrimitiveNuclearPotentialGeom010PP_Y_Z(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <P_Z||P_X> integrals.

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
compPrimitiveNuclearPotentialGeom010PP_Z_X(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <P_Z||P_Y> integrals.

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
compPrimitiveNuclearPotentialGeom010PP_Z_Y(      TDoubleArray& buffer_x,
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
 Evaluates block of primitive <P_Z||P_Z> integrals.

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
compPrimitiveNuclearPotentialGeom010PP_Z_Z(      TDoubleArray& buffer_x,
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

#endif /* NuclearPotentialGeom010RecPP_hpp */
