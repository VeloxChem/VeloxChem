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

#include "XCMolecularGradient.hpp"

#include "MpiFunc.hpp"
#include "DensityGridDriver.hpp"
#include "DensityGradientGridDriver.hpp"
#include "FunctionalParser.hpp"
#include "OMPTasks.hpp"
#include "GtoFunc.hpp"

CXCMolecularGradient::CXCMolecularGradient(MPI_Comm comm)
{
    _locRank = mpi::rank(comm);
    
    _locNodes = mpi::nodes(comm);
    
    _locComm = comm;
    
    _thresholdOfDensity = 1.0e-13;
}

CXCMolecularGradient::~CXCMolecularGradient()
{
    
}

CDenseMatrix
CXCMolecularGradient::integrate(const std::vector<int32_t>& idsAtomic,
                                const CAODensityMatrix&     aoDensityMatrix,
                                const CMolecule&            molecule,
                                const CMolecularBasis&      basis,
                                const CMolecularGrid&       molecularGrid,
                                const std::string&          xcFuncLabel) const
{
    // parse exchange-correlation functional data

    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);
    
    // generate reference density grid
    
    CDensityGridDriver dgdrv(_locComm);
    
    auto refdengrid = dgdrv.generate(aoDensityMatrix, molecule, basis, molecularGrid, fvxc.getFunctionalType());
    
    // create molecular gradient
    
    const auto natoms = static_cast<int32_t>(idsAtomic.size());
    
    CDenseMatrix molgrad(3, natoms);
    
    auto mgradx = molgrad.row(0);
    
    auto mgrady = molgrad.row(1);
    
    auto mgradz = molgrad.row(2);
    
    if (aoDensityMatrix.isClosedShell())
    {
        // generate screened molecular and density grids
        
        CMolecularGrid mgrid(molecularGrid);
        
        CDensityGrid dgrid;
        
        refdengrid.getScreenedGridsPair(dgrid, mgrid, 0, _thresholdOfDensity, fvxc.getFunctionalType());

        auto gw = mgrid.getWeights();
        
        const auto gpoints = mgrid.getNumberOfGridPoints();
        
        // allocate XC gradient grid
        
        CXCGradientGrid vxcgrid(mgrid.getNumberOfGridPoints(), dgrid.getDensityGridType(), fvxc.getFunctionalType());
        
        // compute exchange-correlation functional first derrivatives
        
        fvxc.compute(vxcgrid, dgrid);
        
        auto grhoa = vxcgrid.xcGradientValues(xcvars::rhoa);
        
        // set up density gradient grid driver
        
        CDensityGradientGridDriver graddrv(_locComm);
        
        for (int32_t i = 0; i < natoms; i++)
        {
            auto gradgrid = graddrv.generate(aoDensityMatrix, molecule, basis, mgrid,
                                             idsAtomic[i]);
            
            // set up pointers to density gradient grid
            
            const auto gdenx = gradgrid.getComponent(0);
            
            const auto gdeny = gradgrid.getComponent(1);
            
            const auto gdenz = gradgrid.getComponent(2);
            
            // compute molecular gradient for LDA
            
            double gatmx = 0.0;
            
            double gatmy = 0.0;
            
            double gatmz = 0.0;
            
            for (int32_t j = 0; j < gpoints; j++)
            {
                gatmx += gw[j] * grhoa[j] * gdenx[j];
                
                gatmy += gw[j] * grhoa[j] * gdeny[j];
                
                gatmz += gw[j] * grhoa[j] * gdenz[j];
            }

            // factor of 2 from sum of alpha and beta contributions
            
            mgradx[i] = 2.0 * gatmx;
            
            mgrady[i] = 2.0 * gatmy;
            
            mgradz[i] = 2.0 * gatmz;
        }
    }
    else
    {
        // not implemented
    }
    
    // done with molecular gradient
    
    return molgrad;
}
