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

#ifndef XCIntegrator_hpp
#define XCIntegrator_hpp

#include <mpi.h>

#include <cstdint>
#include <string>
#include <tuple>

#include "AODensityMatrix.hpp"
#include "AOFockMatrix.hpp"
#include "AOKohnShamMatrix.hpp"
#include "DensityGrid.hpp"
#include "DensityGridQuad.hpp"
#include "GtoContainer.hpp"
#include "MemBlock.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "OverlapMatrix.hpp"
#include "XCCubicHessianGrid.hpp"
#include "XCGradientGrid.hpp"
#include "XCHessianGrid.hpp"

/**
 Class CXCIntegrator implements exchange-correlation functional and it's derrivatives integraion.

 @author Z. Rinkevicius
 */
class CXCIntegrator
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

    /**
     Computes exchange-correlation energy and number of electrons for given density grid.

     @param xcGradientGrid the exchange-correlation functional gradient grid.
     @param densityGrid the density grid.
     @param molecularGrid the molecular grid.
     @return the tuple (exchange-correlation energy, number of electrons).
     */
    std::tuple<double, double> _compEnergyAndDensityUnrestricted(const CXCGradientGrid& xcGradientGrid,
                                                                 const CDensityGrid&    densityGrid,
                                                                 const CMolecularGrid&  molecularGrid) const;

   public:
    /**
     Creates a XC integrator object using MPI info.

     @param comm the MPI communicator.
     */
    CXCIntegrator(MPI_Comm comm);

    /**
     Destroys a XC integrator object.
     */
    ~CXCIntegrator();

    /**
    Integrates exchange-correlation pair-density functional contribution to energy.

    @param aoDensityMatrix the AO density matrix object.
    @param twoDM the "active" MO two-body density matrix.
    @param activeMOs the active MO coefficients.
    @param nActive the number of active orbitals.
    @param molecule the molecule.
    @param basis the molecular basis.
    @param molecularGrid the molecular grid.
    @param xcFuncLabel the label of exchange-correlation functional.
    @return the XC energy.
    */
    double integratePdft(const CAODensityMatrix& aoDensityMatrix,
                         const double*           twoDM,
                         const double*           activeMOs,
                         const int32_t           nActive,
                         const CMolecule&        molecule,
                         const CMolecularBasis&  basis,
                         const CMolecularGrid&   molecularGrid,
                         const std::string&      xcFuncLabel) const;
};

#endif /* XCIntegrator_hpp */
