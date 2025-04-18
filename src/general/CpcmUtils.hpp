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

#ifndef CpcmUtils_hpp
#define CpcmUtils_hpp

#include <vector>

namespace cpcm {  // cpcm namespace

auto local_matrix_A_diagonals(const double* ptr_grid_data,
                              const int     row_start,
                              const int     row_end,
                              const int     ncols,
                              const double* ptr_sw_func) -> std::vector<double>;

auto local_matrix_A_dot_vector(const int     npoints,
                               const double* ptr_grid_data,
                               const int     row_start,
                               const int     row_end,
                               const int     ncols,
                               const double* ptr_sw_func,
                               const double* ptr_vector) -> std::vector<double>;

auto comp_grad_Aij(const double* ptr_grid_coords,
                   const double* ptr_zeta,
                   const int*    ptr_atom_indices,
                   const double* ptr_q,
                   const int     row_start,
                   const int     row_end,
                   const int     npoints,
                   const int     natoms) -> std::vector<double>;

auto comp_grad_Aii(const double* ptr_grid_coords,
                   const double* ptr_zeta,
                   const double* ptr_sw_f,
                   const int*    ptr_atom_indices,
                   const double* ptr_q,
                   const double* ptr_atom_coords,
                   const double* ptr_atom_radii,
                   const int     row_start,
                   const int     row_end,
                   const int     npoints,
                   const int     natoms) -> std::vector<double>;

}  // namespace cpcm

#endif /* CpcmUtils_hpp */
