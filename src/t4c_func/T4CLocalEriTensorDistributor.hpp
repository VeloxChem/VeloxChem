#ifndef T4CLocalEriTensorDistributor_hpp
#define T4CLocalEriTensorDistributor_hpp

#include <array>
#include <cstddef>

#include "Dense4DTensor.hpp"
#include "Matrices.hpp"
#include "Matrix.hpp"
#include "SimdArray.hpp"

namespace t4cfunc {  // t4cfunc namespace

auto local_distribute_eri_tensor(CDense4DTensor*                  eri_tensor,
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

}  // namespace t4cfunc

#endif /* T4CLocalEriTensorDistributor_hpp */
