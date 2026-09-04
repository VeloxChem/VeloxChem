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

#ifndef GpuThreadCheck_hpp
#define GpuThreadCheck_hpp

#include <cstdint>
#include <omp.h>
#include <string>

#include "ErrorHandler.hpp"

namespace gpu {  // gpu namespace

// The GPU integral drivers run one OpenMP thread per GPU: checkNumGpusPerNode
// is called before each parallel region (entered with
// num_threads(num_gpus_per_node)) and checkNumGpuThreads is called inside it;
// a mismatch aborts.

// Pre-region check: reject num_gpus_per_node <= 0 before it reaches the
// num_threads clause.
inline void
checkNumGpusPerNode(const int64_t num_gpus_per_node, const char* func_name)
{
    std::string err_ngpus(std::string("gpu::") + func_name + ": invalid number of GPUs per node (" +
                          std::to_string(num_gpus_per_node) + ")");

    errors::assertMsgCritical(num_gpus_per_node > 0, err_ngpus);
}

// In-region check: abort when the actual team size differs from
// num_gpus_per_node.
inline void
checkNumGpuThreads(const int64_t num_gpus_per_node, const char* func_name)
{
    std::string err_nthreads(std::string("gpu::") + func_name +
                             ": number of OpenMP threads (" + std::to_string(omp_get_num_threads()) +
                             ") does not match the number of GPUs per node (" +
                             std::to_string(num_gpus_per_node) + ")");

    errors::assertMsgCritical(omp_get_num_threads() == num_gpus_per_node, err_nthreads);
}

}  // namespace gpu

#endif /* GpuThreadCheck_hpp */
