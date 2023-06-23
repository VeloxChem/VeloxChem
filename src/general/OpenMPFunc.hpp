#ifndef OpenMPFunc_hpp
#define OpenMPFunc_hpp

#include <cstdint>
#include <vector>

#include "GtoBlock.hpp"
#include "GtoFunc.hpp"
#include "T4Index.hpp"
#include "omp.h"

using TGraph = std::vector<T4Index>;

using TWorkGroup = std::vector<TGraph>;

namespace omp {  // omp namespace

/**
 Sets number of OMP threads available.

 @param nthreads the number of OMP threads.
 */
inline auto
setNumberOfThreads(const int nthreads) -> void
{
    omp_set_num_threads(nthreads);
};

/**
 Gets number of OMP threads available.

 @return the number of OMP threads.
 */
inline auto
getNumberOfThreads() -> int
{
    return omp_get_max_threads();
};

/**
 Sets static scheduling for parallel region.
 */
inline auto
setStaticScheduler() -> void
{
    omp_set_dynamic(0);
};

/**
 Gets thread identifier in parallel region.

 @return the thread identifier.
 */
inline auto
getThreadIdentifier() -> int
{
    return omp_get_thread_num();
};

/**
 Generates work group for OMP tasks manager.

 @param gto_blocks the vector of basis function blocks.
 @return the work group.
 */
auto makeWorkGroup(const std::vector<CGtoBlock>& gto_blocks) -> TWorkGroup;

/**
 Generates work group for OMP tasks manager.

 @param bra_gto_blocks the vector of basis function blocks on bra side.
 @param ket_gto_blocks the vector of basis function blocks on ket side.
 @return the work group.
 */
auto makeWorkGroup(const std::vector<CGtoBlock>& bra_gto_blocks, const std::vector<CGtoBlock>& ket_gto_blocks) -> TWorkGroup;

}  // namespace omp

#endif /* OpenMPFunc_hpp */
