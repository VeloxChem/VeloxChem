#ifndef T2CDistributor_hpp
#define T2CDistributor_hpp

#include <cstdint>
#include <vector>

#include "SubMatrix.hpp"
#include "SimdTypes.hpp"
#include "MatrixType.hpp"

namespace t2cfunc { // t2cfunc namespace

/**
 Distributes buffer of integrals into given matrix.
 
 @param matrix the pointer to matrix for storage of integrals.
 @param buffer the integrals buffer.
 @param indexes the compressed contracted GTOs indexes.
 @param bra_comp the angular component of integrals buffer on bra side.
 @param ket_comp the angular component of integrals buffer on ket side.
 @param bra_igto the index of GTO on bra side.
 @param ket_first the index of the range [ket_first, ket_last) of GTOs on bra side.
 @param ket_last the index of the range [ket_first, ket_last) of GTOs on bra side.
 */
auto
distribute(      CSubMatrix*           matrix,
           const TDoubleArray&         buffer,
           const std::vector<int64_t>& indexes,
           const int64_t               bra_comp,
           const int64_t               ket_comp,
           const int64_t               bra_igto,
           const int64_t               ket_first,
           const int64_t               ket_last) -> void;

/**
 Distributes buffer of sclaed integrals into given matrix.
 
 @param matrix the pointer to matrix for storage of integrals.
 @param buffer the integrals buffer.
 @param factor the scaling factor of integrals.
 @param indexes the compressed contracted GTOs indexes.
 @param bra_comp the angular component of integrals buffer on bra side.
 @param ket_comp the angular component of integrals buffer on ket side.
 @param bra_igto the index of GTO on bra side.
 @param ket_first the index of the range [ket_first, ket_last) of GTOs on bra side.
 @param ket_last the index of the range [ket_first, ket_last) of GTOs on bra side.
 */
auto
distribute(      CSubMatrix*           matrix,
           const TDoubleArray&         buffer,
           const double                factor,
           const std::vector<int64_t>& indexes,
           const int64_t               bra_comp,
           const int64_t               ket_comp,
           const int64_t               bra_igto,
           const int64_t               ket_first,
           const int64_t               ket_last) -> void;

/**
 Distributes buffer of integrals into given matrix.
 
 @param matrix the pointer to matrix for storage of integrals.
 @param buffer the integrals buffer.
 @param bra_indexes the compressed contracted GTOs indexes on bra side.
 @param ket_indexes the compressed contracted GTOs indexes on ket side.
 @param bra_comp the angular component of integrals buffer on bra side.
 @param ket_comp the angular component of integrals buffer on ket side.
 @param bra_igto the index of GTO on bra side.
 @param ket_first the index of the range [ket_first, ket_last) of GTOs on bra side.
 @param ket_last the index of the range [ket_first, ket_last) of GTOs on bra side.
 @param mat_type the matrix type.
 */
auto
distribute(      CSubMatrix*           matrix,
           const TDoubleArray&         buffer,
           const std::vector<int64_t>& bra_indexes,
           const std::vector<int64_t>& ket_indexes,
           const int64_t               bra_comp,
           const int64_t               ket_comp,
           const int64_t               bra_igto,
           const int64_t               ket_first,
           const int64_t               ket_last,
           const mat_t                 mat_type) -> void;

/**
 Distributes buffer of scaled integrals into given matrix.
 
 @param matrix the pointer to matrix for storage of integrals.
 @param buffer the integrals buffer.
 @param factor the scaling factor of integrals.
 @param bra_indexes the compressed contracted GTOs indexes on bra side.
 @param ket_indexes the compressed contracted GTOs indexes on ket side.
 @param bra_comp the angular component of integrals buffer on bra side.
 @param ket_comp the angular component of integrals buffer on ket side.
 @param bra_igto the index of GTO on bra side.
 @param ket_first the index of the range [ket_first, ket_last) of GTOs on bra side.
 @param ket_last the index of the range [ket_first, ket_last) of GTOs on bra side.
 @param mat_type the matrix type.
 */
auto
distribute(      CSubMatrix*           matrix,
           const TDoubleArray&         buffer,
           const double                factor,
           const std::vector<int64_t>& bra_indexes,
           const std::vector<int64_t>& ket_indexes,
           const int64_t               bra_comp,
           const int64_t               ket_comp,
           const int64_t               bra_igto,
           const int64_t               ket_first,
           const int64_t               ket_last,
           const mat_t                 mat_type) -> void;

/**
 Distributes buffer of integrals into given matrix.
 
 @param matrix the pointer to matrix for storage of integrals.
 @param buffer the integrals buffer.
 @param bra_indexes the compressed contracted GTOs indexes on bra side.
 @param ket_indexes the compressed contracted GTOs indexes on ket side.
 @param bra_comp the angular component of integrals buffer on bra side.
 @param ket_comp the angular component of integrals buffer on ket side.
 @param bra_igto the index of GTO on bra side.
 @param ket_first the index of the range [ket_first, ket_last) of GTOs on bra side.
 @param ket_last the index of the range [ket_first, ket_last) of GTOs on bra side.
 @param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.
 */
auto
distribute(      CSubMatrix*           matrix,
           const TDoubleArray&         buffer,
           const std::vector<int64_t>& bra_indexes,
           const std::vector<int64_t>& ket_indexes,
           const int64_t               bra_comp,
           const int64_t               ket_comp,
           const int64_t               bra_igto,
           const int64_t               ket_first,
           const int64_t               ket_last,
           const bool                  ang_order) -> void;

/**
 Distributes buffer of sclaed integrals into given matrix.
 
 @param matrix the pointer to matrix for storage of integrals.
 @param buffer the integrals buffer.
 @param factor the scaling factor of integrals.
 @param bra_indexes the compressed contracted GTOs indexes on bra side.
 @param ket_indexes the compressed contracted GTOs indexes on ket side.
 @param bra_comp the angular component of integrals buffer on bra side.
 @param ket_comp the angular component of integrals buffer on ket side.
 @param bra_igto the index of GTO on bra side.
 @param ket_first the index of the range [ket_first, ket_last) of GTOs on bra side.
 @param ket_last the index of the range [ket_first, ket_last) of GTOs on bra side.
 @param ang_order the flag for matching angular order between matrix and pair of GTOs blocks.
 */
auto
distribute(      CSubMatrix*           matrix,
           const TDoubleArray&         buffer,
           const double                factor,
           const std::vector<int64_t>& bra_indexes,
           const std::vector<int64_t>& ket_indexes,
           const int64_t               bra_comp,
           const int64_t               ket_comp,
           const int64_t               bra_igto,
           const int64_t               ket_first,
           const int64_t               ket_last,
           const bool                  ang_order) -> void;

}  // t2cfunc namespace

#endif /* T2CDistributor_hpp */
