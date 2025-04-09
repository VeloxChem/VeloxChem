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

#ifndef OldOneElecIntsDrivers_hpp
#define OldOneElecIntsDrivers_hpp

#include "DenseMatrix.hpp"
#include "Molecule.hpp"
#include "MolecularBasis.hpp"

class COldOneElecIntsMatrix
{
    std::vector<CDenseMatrix> _matrices;

   public:
    COldOneElecIntsMatrix(const std::vector<CDenseMatrix>& matrices);

    COldOneElecIntsMatrix(const COldOneElecIntsMatrix& source);
    COldOneElecIntsMatrix(COldOneElecIntsMatrix&& source) noexcept;

    ~COldOneElecIntsMatrix();

    COldOneElecIntsMatrix& operator=(const COldOneElecIntsMatrix& source);
    COldOneElecIntsMatrix& operator=(COldOneElecIntsMatrix&& source) noexcept;

    CDenseMatrix get_matrix(const int idx) const;
};

class COldElectricDipoleIntegralsDriver
{
   public:
    COldElectricDipoleIntegralsDriver();

    ~COldElectricDipoleIntegralsDriver();

    COldOneElecIntsMatrix compute(const CMolecule& molecule, const CMolecularBasis& basis) const;
};

class COldLinearMomentumIntegralsDriver
{
   public:
    COldLinearMomentumIntegralsDriver();

    ~COldLinearMomentumIntegralsDriver();

    COldOneElecIntsMatrix compute(const CMolecule& molecule, const CMolecularBasis& basis) const;
};

class COldAngularMomentumIntegralsDriver
{
   public:
    COldAngularMomentumIntegralsDriver();

    ~COldAngularMomentumIntegralsDriver();

    COldOneElecIntsMatrix compute(const CMolecule& molecule, const CMolecularBasis& basis) const;
};

#endif /* OldOneElecIntsDrivers_hpp */
