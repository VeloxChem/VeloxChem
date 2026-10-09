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

#ifndef SerialVectorMath_hpp
#define SerialVectorMath_hpp

#include <span>

/**
 Serial vector math.

 Elementwise elementary functions over contiguous arrays: the vector
 counterpart of the scalar math functions, for code that is parallelized at a
 higher level with OpenMP or MPI. Each call runs on the calling thread only and
 never starts threads of its own, so a routine may be called from inside an
 OpenMP parallel region. The routines are not dense linear algebra; those live
 in the sdenblas namespace.
 */
namespace svecmath {  // svecmath namespace

/**
 Computes the exponential of the elements of a contiguous array in place.

 values[i] = exp(values[i]) for i < values.size().

 Uses vvexp of Accelerate's vForce on macOS when VLX_USE_MATHLIB is defined, and
 the vectorized exp of Eigen otherwise.

 @param values the values, overwritten by their exponential.
 */
auto exp_in_place(std::span<double> values) -> void;

}  // namespace svecmath

#endif /* SerialVectorMath_hpp */
