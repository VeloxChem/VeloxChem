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

#include "XCIntegrator.hpp"

#include <mpi.h>

#include <cmath>
#include <iostream>

#include "AOKohnShamMatrix.hpp"
#include "AngularMomentum.hpp"
#include "DenseLinearAlgebra.hpp"
#include "DensityGridDriver.hpp"
#include "FunctionalParser.hpp"
#include "GtoFunc.hpp"
#include "MathConst.hpp"
#include "MpiFunc.hpp"
#include "OMPTasks.hpp"

CXCIntegrator::CXCIntegrator(MPI_Comm comm)
{
    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    _locComm = comm;

    _thresholdOfDensity = 1.0e-13;
}

CXCIntegrator::~CXCIntegrator()
{
}

std::tuple<double, double>
CXCIntegrator::_compEnergyAndDensityUnrestricted(const CXCGradientGrid& xcGradientGrid,
                                                 const CDensityGrid&    densityGrid,
                                                 const CMolecularGrid&  molecularGrid) const
{
    // set up pointers to density grid

    auto rhoa = densityGrid.alphaDensity(0);
    auto rhob = densityGrid.betaDensity(0);

    // set up pointer to exchange-correlation energy grid

    auto efunc = xcGradientGrid.xcFunctionalValues();

    // set up pointer to grid weights

    auto gw = molecularGrid.getWeights();

    // set up number of grid points

    auto gpoints = molecularGrid.getNumberOfGridPoints();

    // initialize exchange-correlation energy and gradient

    double xcene = 0.0, nele = 0.0;

    if ((densityGrid.getDensityGridType() == dengrid::limb) && (densityGrid.getNumberOfGridPoints() != 0))
    {
        #pragma omp parallel for reduction(+ : nele) reduction(+ : xcene)
        for (int32_t i = 0; i < gpoints; i++)
        {
            nele += gw[i] * rhoa[i];

            xcene += gw[i] * efunc[i];
        }
    }

    if ((densityGrid.getDensityGridType() == dengrid::lima) && (densityGrid.getNumberOfGridPoints() != 0))
    {
        #pragma omp parallel for reduction(+ : nele) reduction(+ : xcene)
        for (int32_t i = 0; i < gpoints; i++)
        {
            nele += gw[i] * (rhob[i]);

            xcene += gw[i] * efunc[i];
        }
    }

    if ((densityGrid.getDensityGridType() == dengrid::ab) && (densityGrid.getNumberOfGridPoints() != 0))
    {
        #pragma omp parallel for reduction(+ : nele) reduction(+ : xcene)
        for (int32_t i = 0; i < gpoints; i++)
        {
            nele += gw[i] * (rhoa[i] + rhob[i]);

            xcene += gw[i] * efunc[i];
        }
    }

    return std::make_tuple(xcene, nele);
}

double CXCIntegrator::integratePdft(const CAODensityMatrix& aoDensityMatrix,
                                    const double*           twoDM,
                                    const double*           activeMOs,
                                    const int32_t           nActive,
                                    const CMolecule&        molecule,
                                    const CMolecularBasis&  basis,
                                    const CMolecularGrid&   molecularGrid,
                                    const std::string&      xcFuncLabel) const
{
    CAOKohnShamMatrix ksmat;

    // integrator handles single AO density only

    double xc_energy=0.0;

    if (aoDensityMatrix.getNumberOfDensityMatrices() == 1)
    {   
        // parse exchange-correlation functional data

        auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);
        
        // create GTOs container
        
        CGtoContainer* gtovec = new CGtoContainer(molecule, basis);
        
        // generate reference density grid
        
        CDensityGridDriver dgdrv(_locComm);
        
	    auto refdengrid = dgdrv.generatePdftGrid(aoDensityMatrix, twoDM, activeMOs, nActive, molecule, basis, molecularGrid, fvxc.getFunctionalType());

        // generate screened molecular and density grids

        CMolecularGrid mgridab(molecularGrid);

        CMolecularGrid mgrida(molecularGrid);

        CMolecularGrid mgridb(molecularGrid);

        CDensityGrid dgridab;

        CDensityGrid dgrida;

        CDensityGrid dgridb;

        // screen the molecular and density grids 

        refdengrid.getScreenedGridPairUnrestricted(dgridab, dgrida, dgridb, mgridab, mgrida, mgridb,
                                                   0, _thresholdOfDensity, fvxc.getFunctionalType());

        // allocate XC gradient grid

        CXCGradientGrid vxcgridab(mgridab.getNumberOfGridPoints(), dgridab.getDensityGridType(), fvxc.getFunctionalType());

        CXCGradientGrid vxcgrida(mgrida.getNumberOfGridPoints(), dgrida.getDensityGridType(), fvxc.getFunctionalType());

        CXCGradientGrid vxcgridb(mgridb.getNumberOfGridPoints(), dgridb.getDensityGridType(), fvxc.getFunctionalType());

        // compute exchange-correlation functional first derrivatives

        fvxc.compute(vxcgridab, dgridab);

        fvxc.compute(vxcgrida, dgrida);

        fvxc.compute(vxcgridb, dgridb);

        // Compute correlation energy

        auto xcdata = _compEnergyAndDensityUnrestricted(vxcgrida, dgrida, mgrida);

        auto xcdatb = _compEnergyAndDensityUnrestricted(vxcgridb, dgridb, mgridb);

        auto xcdatab = _compEnergyAndDensityUnrestricted(vxcgridab, dgridab, mgridab);

        xc_energy = std::get<0>(xcdata) + std::get<0>(xcdatb) + std::get<0>(xcdatab);

        // delete GTOs container
        
        delete gtovec;
    }

    return xc_energy;
}
