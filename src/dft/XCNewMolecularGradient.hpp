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

#ifndef XCNewMolecularGradient_hpp
#define XCNewMolecularGradient_hpp

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
#include "XCCubicHessianGrid.hpp"
#include "XCFunctional.hpp"
#include "XCGradientGrid.hpp"
#include "XCHessianGrid.hpp"

/**
 Class CXCNewMolecularGradient implements XC molecular gradient.

 @author X. Li
 */
class CXCNewMolecularGradient
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
     Screening threshold for density values on grid points.
     */
    double _screeningThresholdForDensityValues;

    /**
     Integrates first-order LDA exchnage-correlation functional contribution to
     molecular gradient.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the perturbed AO density matrix (to be contracted
            with GTO gradient).
     @param gsDensityMatrix the ground state AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @return the molecular gradient.
     */
    CDenseMatrix _integrateVxcGradientForLDA(const CMolecule&        molecule,
                                             const CMolecularBasis&  basis,
                                             const CAODensityMatrix& rwDensityMatrix,
                                             const CAODensityMatrix& gsDensityMatrix,
                                             const CMolecularGrid&   molecularGrid,
                                             const CXCFunctional&    xcFunctional) const;

    /**
     Integrates first-order GGA exchnage-correlation functional contribution to
     molecular gradient.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the perturbed AO density matrix.
     @param gsDensityMatrix the ground state AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @return the molecular gradient.
     */
    CDenseMatrix _integrateVxcGradientForGGA(const CMolecule&        molecule,
                                             const CMolecularBasis&  basis,
                                             const CAODensityMatrix& rwDensityMatrix,
                                             const CAODensityMatrix& gsDensityMatrix,
                                             const CMolecularGrid&   molecularGrid,
                                             const CXCFunctional&    xcFunctional) const;

    /**
     Integrates second-order LDA exchnage-correlation functional contribution
     to molecular gradient.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrixOne the perturbed AO density matrix.
     @param rwDensityMatrixTwo the perturbed AO density matrix (to be
            contracted with GTO gradient).
     @param gsDensityMatrix the ground state AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @return the molecular gradient.
     */
    CDenseMatrix _integrateFxcGradientForLDA(const CMolecule&        molecule,
                                             const CMolecularBasis&  basis,
                                             const CAODensityMatrix& rwDensityMatrixOne,
                                             const CAODensityMatrix& rwDensityMatrixTwo,
                                             const CAODensityMatrix& gsDensityMatrix,
                                             const CMolecularGrid&   molecularGrid,
                                             const CXCFunctional&    xcFunctional) const;

    /**
     Integrates second-order GGA exchnage-correlation functional contribution
     to molecular gradient.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrixOne the perturbed AO density matrix.
     @param rwDensityMatrixTwo the perturbed AO density matrix (to be
            contracted with GTO gradient).
     @param gsDensityMatrix the ground state AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @return the molecular gradient.
     */
    CDenseMatrix _integrateFxcGradientForGGA(const CMolecule&        molecule,
                                             const CMolecularBasis&  basis,
                                             const CAODensityMatrix& rwDensityMatrixOne,
                                             const CAODensityMatrix& rwDensityMatrixTwo,
                                             const CAODensityMatrix& gsDensityMatrix,
                                             const CMolecularGrid&   molecularGrid,
                                             const CXCFunctional&    xcFunctional) const;

    /**
     Integrates third-order LDA exchnage-correlation functional contribution to
     molecular gradient.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrixOne the perturbed AO density matrix.
     @param rwDensityMatrixTwo the perturbed AO density matrix.
     @param gsDensityMatrix the ground state AO density matrix (to be
            contracted with GTO gradient).
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @return the molecular gradient.
     */
    CDenseMatrix _integrateKxcGradientForLDA(const CMolecule&        molecule,
                                             const CMolecularBasis&  basis,
                                             const CAODensityMatrix& rwDensityMatrixOne,
                                             const CAODensityMatrix& rwDensityMatrixTwo,
                                             const CAODensityMatrix& gsDensityMatrix,
                                             const CMolecularGrid&   molecularGrid,
                                             const CXCFunctional&    xcFunctional) const;

    /**
     Integrates third-order GGA exchnage-correlation functional contribution to
     molecular gradient.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrixOne the perturbed AO density matrix.
     @param rwDensityMatrixTwo the perturbed AO density matrix.
     @param gsDensityMatrix the ground state AO density matrix (to be
            contracted with GTO gradient).
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @return the molecular gradient.
     */
    CDenseMatrix _integrateKxcGradientForGGA(const CMolecule&        molecule,
                                             const CMolecularBasis&  basis,
                                             const CAODensityMatrix& rwDensityMatrixOne,
                                             const CAODensityMatrix& rwDensityMatrixTwo,
                                             const CAODensityMatrix& gsDensityMatrix,
                                             const CMolecularGrid&   molecularGrid,
                                             const CXCFunctional&    xcFunctional) const;

    /**
     Computes AO-to-atom mapping.

     @param ao_to_atom_ids the vector for storing the mapping.
     @param molecule the molecule.
     @param basis the molecular basis.
     */
    void _computeAOtoAtomMapping(std::vector<int32_t>&  ao_to_atom_ids,
                                 const CMolecule&       molecule,
                                 const CMolecularBasis& basis) const;

   public:
    /**
     Creates an XC integrator object using MPI info.

     @param comm the MPI communicator.
     */
    CXCNewMolecularGradient(MPI_Comm comm);

    /**
     Integrates first-order exchnage-correlation functional contribution to
     molecular gradient.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param gsDensityMatrix the ground state AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the molecular gradient.
     */
    CDenseMatrix integrateVxcGradient(const CMolecule&        molecule,
                                      const CMolecularBasis&  basis,
                                      const CAODensityMatrix& gsDensityMatrix,
                                      CMolecularGrid&         molecularGrid,
                                      const std::string&      xcFuncLabel) const;

    /**
     Integrates first-order exchnage-correlation functional contribution to
     molecular gradient.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the perturbed AO density matrix (to be contracted
            with GTO gradient).
     @param gsDensityMatrix the ground state AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the molecular gradient.
     */
    CDenseMatrix integrateVxcGradient(const CMolecule&        molecule,
                                      const CMolecularBasis&  basis,
                                      const CAODensityMatrix& rwDensityMatrix,
                                      const CAODensityMatrix& gsDensityMatrix,
                                      CMolecularGrid&         molecularGrid,
                                      const std::string&      xcFuncLabel) const;

    /**
     Integrates second-order exchnage-correlation functional contribution to
     molecular gradient.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrixOne the perturbed AO density matrix.
     @param rwDensityMatrixTwo the perturbed AO density matrix (to be
            contracted with GTO gradient).
     @param gsDensityMatrix the ground state AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the molecular gradient.
     */
    CDenseMatrix integrateFxcGradient(const CMolecule&        molecule,
                                      const CMolecularBasis&  basis,
                                      const CAODensityMatrix& rwDensityMatrixOne,
                                      const CAODensityMatrix& rwDensityMatrixTwo,
                                      const CAODensityMatrix& gsDensityMatrix,
                                      CMolecularGrid&         molecularGrid,
                                      const std::string&      xcFuncLabel) const;

    /**
     Integrates third-order exchnage-correlation functional contribution to
     molecular gradient.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrixOne the perturbed AO density matrix.
     @param rwDensityMatrixTwo the perturbed AO density matrix.
     @param gsDensityMatrix the ground state AO density matrix (to be
            contracted with GTO gradient).
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the molecular gradient.
     */
    CDenseMatrix integrateKxcGradient(const CMolecule&        molecule,
                                      const CMolecularBasis&  basis,
                                      const CAODensityMatrix& rwDensityMatrixOne,
                                      const CAODensityMatrix& rwDensityMatrixTwo,
                                      const CAODensityMatrix& gsDensityMatrix,
                                      CMolecularGrid&         molecularGrid,
                                      const std::string&      xcFuncLabel) const;
};

#endif /* XCNewMolecularGradient_hpp */
