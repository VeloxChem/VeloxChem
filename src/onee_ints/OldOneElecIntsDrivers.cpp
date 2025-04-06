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

#include "OldOneElecIntsDrivers.hpp"

#include <cstring>

#include "AngularMomentumIntegrals.hpp"
#include "ElectricDipoleMomentumDriver.hpp"
#include "LinearMomentumIntegrals.hpp"
#include "Point.hpp"

COldOneElecIntsMatrix::COldOneElecIntsMatrix(const std::vector<CDenseMatrix>& matrices)

    : _matrices(matrices)
{
}

COldOneElecIntsMatrix::COldOneElecIntsMatrix(const COldOneElecIntsMatrix& source)

    : _matrices(source._matrices)
{
}

COldOneElecIntsMatrix::COldOneElecIntsMatrix(COldOneElecIntsMatrix&& source) noexcept

    : _matrices(std::move(source._matrices))
{
}

COldOneElecIntsMatrix::~COldOneElecIntsMatrix()
{
}

COldOneElecIntsMatrix&
COldOneElecIntsMatrix::operator=(const COldOneElecIntsMatrix& source)
{
    if (this == &source) return *this;

    _matrices = source._matrices;

    return *this;
}

COldOneElecIntsMatrix&
COldOneElecIntsMatrix::operator=(COldOneElecIntsMatrix&& source) noexcept
{
    if (this == &source) return *this;

    _matrices = std::move(source._matrices);

    return *this;
}

CDenseMatrix
COldOneElecIntsMatrix::get_matrix(const int idx) const
{
    if ((0 <= idx) && (idx < 3)) return _matrices[idx];

    return CDenseMatrix();
}

COldElectricDipoleIntegralsDriver::COldElectricDipoleIntegralsDriver()
{
}

COldElectricDipoleIntegralsDriver::~COldElectricDipoleIntegralsDriver()
{
}

COldOneElecIntsMatrix
COldElectricDipoleIntegralsDriver::compute(const CMolecule& molecule, const CMolecularBasis& basis) const
{
    CElectricDipoleMomentumDriver dip_drv;

    TPoint<double> origin({0.0, 0.0, 0.0});

    auto dip_mats = dip_drv.compute(basis, molecule, origin);

    std::vector<std::string> keys({"X", "Y", "Z"});

    std::vector<CDenseMatrix> new_dip_mats;

    for (const auto& key : keys)
    {
        const auto matptr = dip_mats.matrix(key);

        const auto mat = matptr->full_matrix();

        const auto nrows = mat.number_of_rows();
        const auto ncols = mat.number_of_columns();

        CDenseMatrix new_mat(nrows, ncols);

        std::memcpy(new_mat.values(), mat.data(), nrows * ncols * sizeof(double));

        new_dip_mats.push_back(new_mat);
    }

    return COldOneElecIntsMatrix(new_dip_mats);
}

COldLinearMomentumIntegralsDriver::COldLinearMomentumIntegralsDriver()
{
}

COldLinearMomentumIntegralsDriver::~COldLinearMomentumIntegralsDriver()
{
}

COldOneElecIntsMatrix
COldLinearMomentumIntegralsDriver::compute(const CMolecule& molecule, const CMolecularBasis& basis) const
{
    const auto lmom_mats = onee::computeLinearMomentumIntegrals(molecule, basis);

    return COldOneElecIntsMatrix(lmom_mats);
}

COldAngularMomentumIntegralsDriver::COldAngularMomentumIntegralsDriver()
{
}

COldAngularMomentumIntegralsDriver::~COldAngularMomentumIntegralsDriver()
{
}

COldOneElecIntsMatrix
COldAngularMomentumIntegralsDriver::compute(const CMolecule& molecule, const CMolecularBasis& basis) const
{
    std::vector<double> origin({0.0, 0.0, 0.0});

    auto amom_mats = onee::computeAngularMomentumIntegrals(molecule, basis, origin);

    for (int i = 0; i < 3; i++)
    {
        amom_mats[i].scale(-1.0);
    }

    return COldOneElecIntsMatrix(amom_mats);
}
