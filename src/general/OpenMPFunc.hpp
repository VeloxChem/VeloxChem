#ifndef OpenMPFunc_hpp
#define OpenMPFunc_hpp

#include <array>
#include <cstddef>
#include <vector>

#include "GtoBlock.hpp"
#include "omp.h"

namespace omp {

/// @brief Sets number of OMP threads available.
/// @param nthreads The number of OMP threads.
inline auto
set_number_of_threads(const int nthreads) -> void
{
    omp_set_num_threads(nthreads);
}

/// @brief Gets number of OMP threads available.
/// @return The number of OMP threads.
inline auto
get_number_of_threads() -> int
{
    return omp_get_max_threads();
}

/// @brief Sets static scheduling for parallel region.
inline auto
set_static_scheduler() -> void
{
    omp_set_dynamic(0);
}

/// @brief Generates work groups for OMP tasks manager.
/// @param gto_blocks The vector of basis functions blocks.
/// @return The vector of work tasks.
auto make_work_tasks(const std::vector<CGtoBlock>& gto_blocks) -> std::vector<std::array<size_t, 6>>;

/// @brief Generates work groups for OMP tasks manager.
/// @param bra_gto_blocks The vector of basis functions blocks on bra side.
/// @param ket_gto_blocks The vector of basis functions blocks on ket side.
/// @return The vector of work tasks.
auto make_work_tasks(const std::vector<CGtoBlock>& bra_gto_blocks,
                     const std::vector<CGtoBlock>& ket_gto_blocks) -> std::vector<std::array<size_t, 6>>;

}  // namespace omp

#endif /* OpenMPFunc_hpp */
