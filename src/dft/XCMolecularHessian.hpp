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

#ifndef XCMolecularHessian_hpp
#define XCMolecularHessian_hpp

#include <array>
#include <list>
#include <string>

#include "AODensityMatrix.hpp"
#include "AOFockMatrix.hpp"
#include "AOKohnShamMatrix.hpp"
#include "DenseMatrix.hpp"
#include "DensityGridQuad.hpp"
#include "GridBox.hpp"
#include "GtoContainer.hpp"
#include "MemBlock2D.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "MultiTimer.hpp"
#include "XCFunctional.hpp"

/**
 Class CXCMolecularHessian implements XC molecular Hessian.

 @author X. Li
 */
class CXCMolecularHessian
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
     Screening threshold for GTO values on grid points.
     */
    double _screeningThresholdForGTOValues;

    /**
     Integrates LDA exchnage-correlation functional contribution to molecular
     Hessian.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param gsDensityMatrix the ground state AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @return the molecular Hessian.
     */
    CDenseMatrix _integrateVxcHessianForLDA(const CMolecule&        molecule,
                                            const CMolecularBasis&  basis,
                                            const CAODensityMatrix& gsDensityMatrix,
                                            const CMolecularGrid&   molecularGrid,
                                            const CXCFunctional&    xcFunctional) const;

    CDenseMatrix _integrateFxcHessianForLDA(const CMolecule&        molecule,
                                            const CMolecularBasis&  basis,
                                            const CAODensityMatrix& gsDensityMatrix,
                                            const CMolecularGrid&   molecularGrid,
                                            const CXCFunctional&    xcFunctional) const;

    /**
     Integrates GGA exchnage-correlation functional contribution to molecular
     Hessian.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param gsDensityMatrix the ground state AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @return the molecular Hessian.
     */
    CDenseMatrix _integrateVxcHessianForGGA(const CMolecule&        molecule,
                                            const CMolecularBasis&  basis,
                                            const CAODensityMatrix& gsDensityMatrix,
                                            const CMolecularGrid&   molecularGrid,
                                            const CXCFunctional&    xcFunctional) const;

    CDenseMatrix _integrateFxcHessianForGGA(const CMolecule&        molecule,
                                            const CMolecularBasis&  basis,
                                            const CAODensityMatrix& gsDensityMatrix,
                                            const CMolecularGrid&   molecularGrid,
                                            const CXCFunctional&    xcFunctional) const;

    /**
     Integrates LDA exchnage-correlation functional contribution to molecular
     gradient of Vxc matrix element.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param gsDensityMatrix the ground state AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @param atomIdx the index of the atom with respect to which gradient is
     computed.
     @return the molecular Hessian.
     */
    std::vector<CDenseMatrix> _integrateVxcFockGradientForLDA(const CMolecule&        molecule,
                                                              const CMolecularBasis&  basis,
                                                              const CAODensityMatrix& gsDensityMatrix,
                                                              const CMolecularGrid&   molecularGrid,
                                                              const CXCFunctional&    xcFunctional,
                                                              const int32_t           atomIdx) const;

    /**
     Integrates GGA exchnage-correlation functional contribution to molecular
     gradient of Vxc matrix element.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param gsDensityMatrix the ground state AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @param atomIdx the index of the atom with respect to which gradient is
     computed.
     @return the molecular Hessian.
     */
    std::vector<CDenseMatrix> _integrateVxcFockGradientForGGA(const CMolecule&        molecule,
                                                              const CMolecularBasis&  basis,
                                                              const CAODensityMatrix& gsDensityMatrix,
                                                              const CMolecularGrid&   molecularGrid,
                                                              const CXCFunctional&    xcFunctional,
                                                              const int32_t           atomIdx) const;

    /**
     Computes AO-to-atom mapping.

     @param ao_to_atom_ids the vector for storing the mapping.
     @param molecule the molecule.
     @param basis the molecular basis.
     */
    void _computeAOtoAtomMapping(std::vector<int32_t>& ao_to_atom_ids, const CMolecule& molecule, const CMolecularBasis& basis) const;

   public:
    /**
     Creates an XC integrator object using MPI info.

     @param comm the MPI communicator.
     */
    CXCMolecularHessian(MPI_Comm comm);

    /**
     Integrates exchnage-correlation functional contribution to molecular
     Hessian.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param gsDensityMatrix the ground state AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the molecular Hessian.
     */
    CDenseMatrix integrateExcHessian(const CMolecule&        molecule,
                                     const CMolecularBasis&  basis,
                                     const CAODensityMatrix& gsDensityMatrix,
                                     const CMolecularGrid&   molecularGrid,
                                     const std::string&      xcFuncLabel) const;

    /**
     Integrates exchnage-correlation functional contribution to molecular
     gradient of Vxc matrix element.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param gsDensityMatrix the ground state AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @param atomIdx the index of the atom with respect to which gradient is
     computed.
     @return the molecular Hessian.
     */
    std::vector<CDenseMatrix> integrateVxcFockGradient(const CMolecule&        molecule,
                                                       const CMolecularBasis&  basis,
                                                       const CAODensityMatrix& gsDensityMatrix,
                                                       const CMolecularGrid&   molecularGrid,
                                                       const std::string&      xcFuncLabel,
                                                       const int32_t           atomIdx) const;
};

#endif /* XCMolecularHessian_hpp */
