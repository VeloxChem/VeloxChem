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

#ifndef GtoValuesMggaGPU_hpp
#define GtoValuesMggaGPU_hpp

#include <cstdint>

namespace gpu {  // gpu namespace

__global__ void gtoValuesMggaRecS(double*        gto_values,
                                  double*        gto_values_x,
                                  double*        gto_values_y,
                                  double*        gto_values_z,
                                  double*        gto_values_xx,
                                  double*        gto_values_xy,
                                  double*        gto_values_xz,
                                  double*        gto_values_yy,
                                  double*        gto_values_yz,
                                  double*        gto_values_zz,
                                  const uint32_t row_offset,
                                  const double*  gto_info,
                                  const double*  grid_x,
                                  const double*  grid_y,
                                  const double*  grid_z,
                                  const uint32_t grid_offset,
                                  const uint32_t nrows,
                                  const uint32_t npgtos,
                                  const uint32_t ncols);

__global__ void gtoValuesMggaRecP(double*        gto_values_p3,
                                  double*        gto_values_p3_x,
                                  double*        gto_values_p3_y,
                                  double*        gto_values_p3_z,
                                  double*        gto_values_p3_xx,
                                  double*        gto_values_p3_xy,
                                  double*        gto_values_p3_xz,
                                  double*        gto_values_p3_yy,
                                  double*        gto_values_p3_yz,
                                  double*        gto_values_p3_zz,
                                  const uint32_t row_offset,
                                  const double*  gto_info,
                                  const double*  grid_x,
                                  const double*  grid_y,
                                  const double*  grid_z,
                                  const uint32_t grid_offset,
                                  const uint32_t nrows,
                                  const uint32_t npgtos,
                                  const uint32_t ncols);

__global__ void gtoValuesMggaRecD(double*        gto_values_d5,
                                  double*        gto_values_d5_x,
                                  double*        gto_values_d5_y,
                                  double*        gto_values_d5_z,
                                  double*        gto_values_d5_xx,
                                  double*        gto_values_d5_xy,
                                  double*        gto_values_d5_xz,
                                  double*        gto_values_d5_yy,
                                  double*        gto_values_d5_yz,
                                  double*        gto_values_d5_zz,
                                  const uint32_t row_offset,
                                  const double   f2_3,
                                  const double*  gto_info,
                                  const double*  grid_x,
                                  const double*  grid_y,
                                  const double*  grid_z,
                                  const uint32_t grid_offset,
                                  const uint32_t nrows,
                                  const uint32_t npgtos,
                                  const uint32_t ncols);

}  // namespace gpu

#endif
