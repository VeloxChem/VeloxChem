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

#ifndef OverlapDriver_hpp
#define OverlapDriver_hpp

#include "Matrix.hpp"
#include "MolecularBasis.hpp"
#include "Molecule.hpp"

/**
 Class COverlapDriver provides methods for computing two-center
 overlap integrals.
 */
class COverlapDriver
{
   public:
    /**
     Creates an overlap integrals driver.
     */
    COverlapDriver() = default;

    /**
     Computes overlap matrix for given molecule and molecular basis.

     @param basis the molecular basis.
     @param molecule the molecule.
     @return the overlap matrix.
     */
    auto compute(const CMolecularBasis& basis, const CMolecule& molecule) const -> CMatrix;
    
    /**
     Computes overlap matrix for given molecule and pair of molecular bases.

     @param bra_basis the molecular basis on bra side.
     @param ket_basis the molecular basis on ket side.
     @param molecule the molecule.
     @return the overlap matrix.
     */
    auto compute(const CMolecularBasis& bra_basis, const CMolecularBasis& ket_basis, const CMolecule& molecule) const -> CMatrix;
};

#endif /* OverlapDriver_hpp */
