//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//
//  Redistribution and use in source and binary forms, with or without modification,
//  are permitted provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this
//     list of conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//  3. Neither the name of the copyright holder nor the names of its contributors
//     may be used to endorse or promote products derived from this software without
//     specific prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
//  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
//  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
//  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
//  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
//  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
//  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
//  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
//  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
