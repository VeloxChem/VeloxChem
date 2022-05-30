//
//                           VELOXCHEM 1.0-RC2
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
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

#ifndef XCMolecularGradient_hpp
#define XCMolecularGradient_hpp

#include <mpi.h>

#include <cstdint>
#include <string>
#include <vector>

#include "AODensityMatrix.hpp"
#include "DenseMatrix.hpp"
#include "DensityGrid.hpp"
#include "GtoContainer.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "XCGradientGrid.hpp"

/**
 Class CXCMolecularGradient implements integration of exchange-correlation contribution to molecular
 gradient.

 @author Z. Rinkevicius
 */
class CXCMolecularGradient
{
    /**
     The rank of associated local MPI process.
     */
    int32_t _locRank;

    /**
     The total number of local MPI processes.
     */
    int32_t _locNodes;

    /**
     The MPI communicator.
     */
    MPI_Comm _locComm;

    /**
     The threshold of density screening.
     */
    double _thresholdOfDensity;

   public:
    /**
     Creates a XC molecular gradient integrator object using MPI info.

     @param comm the MPI communicator.
     */
    CXCMolecularGradient(MPI_Comm comm);

    /**
     Destroys a XC molecular gradient integrator object.
     */
    ~CXCMolecularGradient();

    /**
     Integrates exchnage-correlation functional contribution to molecular gradient.

     @param idsAtomic  the atomic indexes.
     @param aoDensityMatrix the AO density matrix object.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the molecular gradient.
     */
    CDenseMatrix integrate(const std::vector<int32_t>& idsAtomic,
                           const CAODensityMatrix&     aoDensityMatrix,
                           const CMolecule&            molecule,
                           const CMolecularBasis&      basis,
                           const CMolecularGrid&       molecularGrid,
                           const std::string&          xcFuncLabel) const;
};

#endif /* XCMolecularGradient_hpp */
