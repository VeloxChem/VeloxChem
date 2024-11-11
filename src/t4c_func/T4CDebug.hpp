#ifndef T4CDebug_hpp
#define T4CDebug_hpp

#include <array>
#include <cstddef>
#include <utility>

#include "SimdArray.hpp"
#include "GtoPairBlock.hpp"

namespace t4cfunc {  // t4cfunc namespace

/// @brief Dumps buffer to output stream.
/// @param buffer The buffer to be dumped.
/// @param bra_gto_pair_block The GTOs pair block on bra side.
/// @param ket_gto_pair_block The GTOs pair block on ket side.
/// @param ket_indices The range [ket_first, ket_last) of basis function pairs on ket side.
/// @param bra_index The index of basis function pairs on bra side.
/// @param components The number of geometrical components.
auto
dump_buffer(const CSimdArray<double>&        buffer,
            const CGtoPairBlock&             bra_gto_pair_block,
            const CGtoPairBlock&             ket_gto_pair_block,
            const std::pair<size_t, size_t>& ket_indices,
            const size_t                     bra_index,
            const size_t                     components) -> void;

}  // namespace t4cfunc

#endif /* T4CDebug_hpp */
