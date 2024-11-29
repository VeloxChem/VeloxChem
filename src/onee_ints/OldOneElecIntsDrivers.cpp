//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

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
