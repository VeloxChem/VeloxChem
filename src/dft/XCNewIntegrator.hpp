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

#ifndef XCNewIntegrator_hpp
#define XCNewIntegrator_hpp

#include <array>
#include <string>

#include "AODensityMatrix.hpp"
#include "AOFockMatrix.hpp"
#include "AOKohnShamMatrix.hpp"
#include "Dense4DTensor.hpp"
#include "DenseMatrix.hpp"
#include "DensityGridQuad.hpp"
#include "DensityGridCubic.hpp"
#include "GridBox.hpp"
#include "GtoContainer.hpp"
#include "MemBlock2D.hpp"
#include "MolecularBasis.hpp"
#include "MolecularGrid.hpp"
#include "Molecule.hpp"
#include "MultiTimer.hpp"
#include "XCNewFunctional.hpp"
#include "XCPairDensityFunctional.hpp"

/**
 Class CXCNewIntegrator implements XC integrator.

 @author X. Li, K. Ahmadzadeh, M. Delcey
 */
class CXCNewIntegrator
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
     Integrates first-order LDA exchange-correlation functional contribution to
     AO Kohn-Sham matrix.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param densityMatrix the AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @param flag the flag for closed/open shell.
     @return the AO Kohn-Sham matrix.
     */
    CAOKohnShamMatrix _integrateVxcFockForLDA(const CMolecule&        molecule,
                                              const CMolecularBasis&  basis,
                                              const CAODensityMatrix& densityMatrix,
                                              const CMolecularGrid&   molecularGrid,
                                              const CXCNewFunctional& xcFunctional,
                                              const std::string&      flag=std::string("closedshell")) const;

    /**
     Integrates first-order GGA exchange-correlation functional contribution to
     AO Kohn-Sham matrix.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param densityMatrix the AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @param flag the flag for closed/open shell.
     @return the AO Kohn-Sham matrix.
     */
    CAOKohnShamMatrix _integrateVxcFockForGGA(const CMolecule&        molecule,
                                              const CMolecularBasis&  basis,
                                              const CAODensityMatrix& densityMatrix,
                                              const CMolecularGrid&   molecularGrid,
                                              const CXCNewFunctional& xcFunctional,
                                              const std::string&      flag=std::string("closedshell")) const;

    /**
     Integrates first-order meta-GGA exchange-correlation functional
     contribution to AO Kohn-Sham matrix.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param densityMatrix the AO density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @param flag the flag for closed/open shell.
     @return the AO Kohn-Sham matrix.
     */
    CAOKohnShamMatrix _integrateVxcFockForMGGA(const CMolecule&        molecule,
                                               const CMolecularBasis&  basis,
                                               const CAODensityMatrix& densityMatrix,
                                               const CMolecularGrid&   molecularGrid,
                                               const CXCNewFunctional& xcFunctional,
                                               const std::string&      flag=std::string("closedshell")) const;

    /**
     Integrates second-order LDA exchange-correlation functional contribution
     to AO Fock matrix.

     @param aoFockMatrix the AO Fock matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the perturbed density matrix.
     @param gsDensityMatrix the ground state density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     */
    void _integrateFxcFockForLDA(CAOFockMatrix&          aoFockMatrix,
                                 const CMolecule&        molecule,
                                 const CMolecularBasis&  basis,
                                 const CAODensityMatrix& rwDensityMatrix,
                                 const CAODensityMatrix& gsDensityMatrix,
                                 const CMolecularGrid&   molecularGrid,
                                 const CXCNewFunctional& xcFunctional) const;

    /**
     Integrates second-order GGA exchange-correlation functional contribution
     to AO Fock matrix.

     @param aoFockMatrix the AO Fock matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the perturbed density matrix.
     @param gsDensityMatrix the ground state density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     */
    void _integrateFxcFockForGGA(CAOFockMatrix&          aoFockMatrix,
                                 const CMolecule&        molecule,
                                 const CMolecularBasis&  basis,
                                 const CAODensityMatrix& rwDensityMatrix,
                                 const CAODensityMatrix& gsDensityMatrix,
                                 const CMolecularGrid&   molecularGrid,
                                 const CXCNewFunctional& xcFunctional) const;

    /**
     Integrates second-order meta-GGA exchange-correlation functional
     contribution to AO Fock matrix.

     @param aoFockMatrix the AO Fock matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the perturbed density matrix.
     @param gsDensityMatrix the ground state density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     */
    void _integrateFxcFockForMGGA(CAOFockMatrix&         aoFockMatrix,
                                 const CMolecule&        molecule,
                                 const CMolecularBasis&  basis,
                                 const CAODensityMatrix& rwDensityMatrix,
                                 const CAODensityMatrix& gsDensityMatrix,
                                 const CMolecularGrid&   molecularGrid,
                                 const CXCNewFunctional& xcFunctional) const;

    /**
     Integrates third-order LDA exchange-correlation functional contribution
     to AO Fock matrix.

     @param aoFockMatrix the AO Fock matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the perturbed one-time transformed densities.
     @param rw2DensityMatrix the two-time transformed densities.
     @param gsDensityMatrix the ground state density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @param quadMode a string that specifies which densities should be combined.
     */
    void _integrateKxcFockForLDA(CAOFockMatrix&          aoFockMatrix,
                                 const CMolecule&        molecule,
                                 const CMolecularBasis&  basis,
                                 const CAODensityMatrix& rwDensityMatrix,
                                 const CAODensityMatrix& rw2DensityMatrix,
                                 const CAODensityMatrix& gsDensityMatrix,
                                 const CMolecularGrid&   molecularGrid,
                                 const CXCNewFunctional& xcFunctional,
                                 const std::string&      quadMode) const;

    /**
     Integrates third-order GGA exchange-correlation functional contribution
     to AO Fock matrix.

     @param aoFockMatrix the AO Fock matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the perturbed one-time transformed densities.
     @param rw2DensityMatrix the two-time transformed densities.
     @param gsDensityMatrix the ground state density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @param quadMode a string that specifies which densities should be combined.
     */
    void _integrateKxcFockForGGA(CAOFockMatrix&          aoFockMatrix,
                                 const CMolecule&        molecule,
                                 const CMolecularBasis&  basis,
                                 const CAODensityMatrix& rwDensityMatrix,
                                 const CAODensityMatrix& rw2DensityMatrix,
                                 const CAODensityMatrix& gsDensityMatrix,
                                 const CMolecularGrid&   molecularGrid,
                                 const CXCNewFunctional& xcFunctional,
                                 const std::string&      quadMode) const;

    /**
     Integrates third-order meta-GGA exchange-correlation functional
     contribution to AO Fock matrix.

     @param aoFockMatrix the AO Fock matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the perturbed one-time transformed densities.
     @param rw2DensityMatrix the two-time transformed densities.
     @param gsDensityMatrix the ground state density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @param quadMode a string that specifies which densities should be combined.
     */
    void _integrateKxcFockForMGGA(CAOFockMatrix&          aoFockMatrix,
                                  const CMolecule&        molecule,
                                  const CMolecularBasis&  basis,
                                  const CAODensityMatrix& rwDensityMatrix,
                                  const CAODensityMatrix& rw2DensityMatrix,
                                  const CAODensityMatrix& gsDensityMatrix,
                                  const CMolecularGrid&   molecularGrid,
                                  const CXCNewFunctional& xcFunctional,
                                  const std::string&      quadMode) const;

    /**
     Integrates fourth-order LDA exchnage-correlation functional contribution
     to AO Fock matrix.

     @param aoFockMatrix the AO Fock matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the perturbed one-time transformed densities.
     @param rw2DensityMatrix the two-time transformed densities.
     @param rw3DensityMatrix the three-time transformed densities.
     @param gsDensityMatrix the ground state density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @param quadMode a string that specifies which densities should be combined.
     */
    void _integrateLxcFockForLDA(CAOFockMatrix&          aoFockMatrix,
                                 const CMolecule&        molecule,
                                 const CMolecularBasis&  basis,
                                 const CAODensityMatrix& rwDensityMatrix,
                                 const CAODensityMatrix& rw2DensityMatrix,
                                 const CAODensityMatrix& rw3DensityMatrix,
                                 const CAODensityMatrix& gsDensityMatrix,
                                 const CMolecularGrid&   molecularGrid,
                                 const CXCNewFunctional& xcFunctional,
                                 const std::string&      quadMode) const;

    /**
     Integrates fourth-order GGA exchange-correlation functional contribution
     to AO Fock matrix.

     @param aoFockMatrix the AO Fock matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the perturbed one-time transformed densities.
     @param rw2DensityMatrix the two-time transformed densities.
     @param rw3DensityMatrix the three-time transformed densities.
     @param gsDensityMatrix the ground state density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @param cubeMode a string that specifies which densities should be combined.
     */
    void _integrateLxcFockForGGA(CAOFockMatrix&          aoFockMatrix,
                                 const CMolecule&        molecule,
                                 const CMolecularBasis&  basis,
                                 const CAODensityMatrix& rwDensityMatrix,
                                 const CAODensityMatrix& rw2DensityMatrix,
                                 const CAODensityMatrix& rw3DensityMatrix,
                                 const CAODensityMatrix& gsDensityMatrix,
                                 const CMolecularGrid&   molecularGrid,
                                 const CXCNewFunctional& xcFunctional,
                                 const std::string&      cubeMode) const;

    /**
     Integrates fourth-order LDA exchnage-correlation functional contribution
     to AO Fock matrix.

     @param aoFockMatrix the AO Fock matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the perturbed one-time transformed densities.
     @param rw2DensityMatrix the two-time transformed densities.
     @param rw3DensityMatrix the three-time transformed densities.
     @param gsDensityMatrix the ground state density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @param quadMode a string that specifies which densities should be combined.
     */
    void _integrateKxcLxcFockForLDA(CAOFockMatrix&          aoFockMatrix,
                                    const CMolecule&        molecule,
                                    const CMolecularBasis&  basis,
                                    const CAODensityMatrix& rwDensityMatrix,
                                    const CAODensityMatrix& rw2DensityMatrix,
                                    const CAODensityMatrix& rw3DensityMatrix,
                                    const CAODensityMatrix& gsDensityMatrix,
                                    const CMolecularGrid&   molecularGrid,
                                    const CXCNewFunctional& xcFunctional,
                                    const std::string&      quadMode) const;

    /**
     Integrates fourth-order GGA exchange-correlation functional contribution
     to AO Fock matrix.

     @param aoFockMatrix the AO Fock matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the perturbed one-time transformed densities.
     @param rw2DensityMatrix the two-time transformed densities.
     @param rw3DensityMatrix the three-time transformed densities.
     @param gsDensityMatrix the ground state density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @param cubeMode a string that specifies which densities should be combined.
     */
    void _integrateKxcLxcFockForGGA(CAOFockMatrix&          aoFockMatrix,
                                 const CMolecule&        molecule,
                                 const CMolecularBasis&  basis,
                                 const CAODensityMatrix& rwDensityMatrix,
                                 const CAODensityMatrix& rw2DensityMatrix,
                                 const CAODensityMatrix& rw3DensityMatrix,
                                 const CAODensityMatrix& gsDensityMatrix,
                                 const CMolecularGrid&   molecularGrid,
                                 const CXCNewFunctional& xcFunctional,
                                 const std::string&      cubeMode) const;

    /**
     Integrates fourth-order meta-GGA exchange-correlation functional
     contribution to AO Fock matrix.

     @param aoFockMatrix the AO Fock matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the perturbed one-time transformed densities.
     @param rw2DensityMatrix the two-time transformed densities.
     @param rw3DensityMatrix the three-time transformed densities.
     @param gsDensityMatrix the ground state density matrix.
     @param molecularGrid the molecular grid.
     @param xcFunctional the exchange-correlation functional.
     @param cubeMode a string that specifies which densities should be combined.
     */
    void _integrateKxcLxcFockForMGGA(CAOFockMatrix&          aoFockMatrix,
                                     const CMolecule&        molecule,
                                     const CMolecularBasis&  basis,
                                     const CAODensityMatrix& rwDensityMatrix,
                                     const CAODensityMatrix& rw2DensityMatrix,
                                     const CAODensityMatrix& rw3DensityMatrix,
                                     const CAODensityMatrix& gsDensityMatrix,
                                     const CMolecularGrid&   molecularGrid,
                                     const CXCNewFunctional& xcFunctional,
                                     const std::string&      cubeMode) const;

    /**
     Integrates first-order LDA pair-density functional contribution to
     AO Kohn-Sham matrix and MO "Q-matrix".

     @param aoFockMatrix the AO Fock matrix.
     @param moTwoBodyGradient the MO Two-body energy gradient term.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param DensityMatrix the AO density matrix object.
     @param TwoBodyDensityMatrix the MO two-body active density matrix.
     @param ActiveMOs the active molecular orbitals.
     @param molecularGrid the molecular grid.
     @param fvxc the exchange-correlation functional.
     */
    void _integrateVxcPDFTForLDA(CAOKohnShamMatrix&              aoFockMatrix,
                                 CDense4DTensor&                 moTwoBodyGradient,
                                 const CMolecule&                molecule,
                                 const CMolecularBasis&          basis,
                                 const CAODensityMatrix&         DensityMatrix,
                                 const CDense4DTensor&           TwoBodyDensityMatrix,
                                 const CDenseMatrix&             ActiveMOs,
                                 const CMolecularGrid&           molecularGrid,
                                 const CXCPairDensityFunctional& xcFunctional) const;

    /**
     Integrates first-order GGA pair-density functional contribution to
     AO Kohn-Sham matrix and MO "Q-matrix".

     @param aoFockMatrix the AO Fock matrix.
     @param moTwoBodyGradient the MO Two-body energy gradient term.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param DensityMatrix the AO density matrix object.
     @param TwoBodyDensityMatrix the MO two-body active density matrix.
     @param ActiveMOs the active molecular orbitals.
     @param molecularGrid the molecular grid.
     @param fvxc the exchange-correlation functional.
     */
    void _integrateVxcPDFTForGGA(CAOKohnShamMatrix&              aoFockMatrix,
                                 CDense4DTensor&                 moTwoBodyGradient,
                                 const CMolecule&                molecule,
                                 const CMolecularBasis&          basis,
                                 const CAODensityMatrix&         DensityMatrix,
                                 const CDense4DTensor&           TwoBodyDensityMatrix,
                                 const CDenseMatrix&             ActiveMOs,
                                 const CMolecularGrid&           molecularGrid,
                                 const CXCPairDensityFunctional& xcFunctional) const;

    /**
     Integrates LDA contribution to (first-order) Vxc matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param vrho the 1st-order functional derivative wrt density.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialVxcFockForLDA(const int32_t       npoints,
                                                const double*       weights,
                                                const CDenseMatrix& gtoValues,
                                                const double*       vrho,
                                                CMultiTimer&        timer) const;

    /**
     Integrates LDA contribution to (first-order) Vxc matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param vrho the 1st-order functional derivative wrt density.
     @param timer the timer.
     @return the alpha and beta contribution as a list of CDenseMatrix objects.
     */
    std::vector<CDenseMatrix> _integratePartialVxcFockForLDAOpenShell(const int32_t       npoints,
                                                                      const double*       weights,
                                                                      const CDenseMatrix& gtoValues,
                                                                      const double*       vrho,
                                                                      CMultiTimer&        timer) const;

    /**
     Integrates GGA contribution to AO Kohn-Sham matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param gtoValuesX the GTO gradient X values on grid points.
     @param gtoValuesY the GTO gradient Y values on grid points.
     @param gtoValuesZ the GTO gradient Z values on grid points.
     @param rhograd the gradient density.
     @param vrho the 1st-order functional derivative wrt rho.
     @param vsigma the 1st-order functional derivative wrt sigma.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialVxcFockForGGA(const int32_t       npoints,
                                                const double*       weights,
                                                const CDenseMatrix& gtoValues,
                                                const CDenseMatrix& gtoValuesX,
                                                const CDenseMatrix& gtoValuesY,
                                                const CDenseMatrix& gtoValuesZ,
                                                const double*       rhograd,
                                                const double*       vrho,
                                                const double*       vsigma,
                                                CMultiTimer&        timer) const;

    /**
     Integrates GGA contribution to AO Kohn-Sham matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param gtoValuesX the GTO gradient X values on grid points.
     @param gtoValuesY the GTO gradient Y values on grid points.
     @param gtoValuesZ the GTO gradient Z values on grid points.
     @param rhograd the gradient density.
     @param vrho the 1st-order functional derivative wrt rho.
     @param vsigma the 1st-order functional derivative wrt sigma.
     @param timer the timer.
     @return the alpha and beta contribution as a list of CDenseMatrix objects.
     */
    std::vector<CDenseMatrix> _integratePartialVxcFockForGGAOpenShell(const int32_t       npoints,
                                                                      const double*       weights,
                                                                      const CDenseMatrix& gtoValues,
                                                                      const CDenseMatrix& gtoValuesX,
                                                                      const CDenseMatrix& gtoValuesY,
                                                                      const CDenseMatrix& gtoValuesZ,
                                                                      const double*       rhograd,
                                                                      const double*       vrho,
                                                                      const double*       vsigma,
                                                                      CMultiTimer&        timer) const;

    /**
     Integrates meta-GGA contribution to AO Kohn-Sham matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param gtoValuesX the GTO gradient X values on grid points.
     @param gtoValuesY the GTO gradient Y values on grid points.
     @param gtoValuesZ the GTO gradient Z values on grid points.
     @param rhograd the gradient density.
     @param vrho the 1st-order functional derivative wrt rho.
     @param vsigma the 1st-order functional derivative wrt sigma.
     @param vlapl the 1st-order functional derivative wrt laplacian.
     @param vtau the 1st-order functional derivative wrt tau.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialVxcFockForMGGA(const int32_t       npoints,
                                                 const double*       weights,
                                                 const CDenseMatrix& gtoValues,
                                                 const CDenseMatrix& gtoValuesX,
                                                 const CDenseMatrix& gtoValuesY,
                                                 const CDenseMatrix& gtoValuesZ,
                                                 const double*       rhograd,
                                                 const double*       vrho,
                                                 const double*       vsigma,
                                                 const double*       vlapl,
                                                 const double*       vtau,
                                                 CMultiTimer&        timer) const;

    /**
     Integrates meta-GGA contribution to AO Kohn-Sham matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param gtoValuesX the GTO gradient X values on grid points.
     @param gtoValuesY the GTO gradient Y values on grid points.
     @param gtoValuesZ the GTO gradient Z values on grid points.
     @param rhograd the gradient density.
     @param vrho the 1st-order functional derivative wrt rho.
     @param vsigma the 1st-order functional derivative wrt sigma.
     @param vlapl the 1st-order functional derivative wrt laplacian.
     @param vtau the 1st-order functional derivative wrt tau.
     @param timer the timer.
     @return the alpha and beta contribution as a list of CDenseMatrix objects.
     */
    std::vector<CDenseMatrix> _integratePartialVxcFockForMGGAOpenShell(const int32_t       npoints,
                                                                       const double*       weights,
                                                                       const CDenseMatrix& gtoValues,
                                                                       const CDenseMatrix& gtoValuesX,
                                                                       const CDenseMatrix& gtoValuesY,
                                                                       const CDenseMatrix& gtoValuesZ,
                                                                       const double*       rhograd,
                                                                       const double*       vrho,
                                                                       const double*       vsigma,
                                                                       const double*       vlapl,
                                                                       const double*       vtau,
                                                                       CMultiTimer&        timer) const;

    /**
     Integrates LDA contribution to (second-order) Fxc matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param rhow the pointer to perturbed density.
     @param v2rho2 the 2nd-order functional derivative wrt density.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialFxcFockForLDA(const int32_t       npoints,
                                                const double*       weights,
                                                const CDenseMatrix& gtoValues,
                                                const double*       rhow,
                                                const double*       v2rho2,
                                                CMultiTimer&        timer) const;

    /**
     Integrates GGA contribution to (second-order) Fxc matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param gtoValuesX the GTO gradient X values on grid points.
     @param gtoValuesY the GTO gradient Y values on grid points.
     @param gtoValuesZ the GTO gradient Z values on grid points.
     @param rhow the pointer to perturbed density.
     @param rhograd the pointer to density gradient.
     @param rhowgrad the pointer to perturbed density gradient.
     @param v2rho2 the 2nd-order functional derivative wrt density.
     @param v2rhosigma the 2nd-order functional derivative wrt density and
            density gradient.
     @param v2sigma2 the 2nd-order functional derivative wrt density gradient.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialFxcFockForGGA(const int32_t       npoints,
                                                const double*       weights,
                                                const CDenseMatrix& gtoValues,
                                                const CDenseMatrix& gtoValuesX,
                                                const CDenseMatrix& gtoValuesY,
                                                const CDenseMatrix& gtoValuesZ,
                                                const double*       rhow,
                                                const double*       rhograd,
                                                const double*       rhowgrad,
                                                const double*       vsigma,
                                                const double*       v2rho2,
                                                const double*       v2rhosigma,
                                                const double*       v2sigma2,
                                                CMultiTimer&        timer) const;

    /**
     Integrates meta-GGA contribution to (second-order) Fxc matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param gtoValuesX the GTO gradient X values on grid points.
     @param gtoValuesY the GTO gradient Y values on grid points.
     @param gtoValuesZ the GTO gradient Z values on grid points.
     @param rhow the pointer to perturbed density.
     @param rhograd the pointer to density gradient.
     @param rhowgrad the pointer to perturbed density gradient.
     @param tauw , 
     @param laplw ,
     @param vrho , 
     @param vsigma , 
     @param vlapl , 
     @param vtau , 
     @param v2rho2 the 2nd-order functional derivative wrt density.
     @param v2lapl2 , 
     @param v2tau2 , 
     @param v2rholapl , 
     @param v2rhotau ,
     @param v2lapltau , 
     @param v2rhosigma the 2nd-order functional derivative wrt density and
            density gradient.
     @param v2sigmalapl , 
     @param v2sigmatau ,
     @param v2sigma2 the 2nd-order functional derivative wrt density gradient.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialFxcFockForMGGA(const int32_t       npoints, 
                                                 const double*       local_weights, 
                                                 const CDenseMatrix& gtoValues,
                                                 const CDenseMatrix& gtoValuesX, 
                                                 const CDenseMatrix& gtoValuesY, 
                                                 const CDenseMatrix& gtoValuesZ,
                                                 const double*       rhow, 
                                                 const double*       rhograd, 
                                                 const double*       rhowgrad, 
                                                 const double*       tauw, 
                                                 const double*       laplw,
                                                 const double*       vrho, 
                                                 const double*       vsigma, 
                                                 const double*       vlapl, 
                                                 const double*       vtau, 
                                                 const double*       v2rho2,
                                                 const double*       v2lapl2, 
                                                 const double*       v2tau2, 
                                                 const double*       v2rholapl, 
                                                 const double*       v2rhotau,
                                                 const double*       v2lapltau, 
                                                 const double*       v2rhosigma, 
                                                 const double*       v2sigmalapl, 
                                                 const double*       v2sigmatau,
                                                 const double*       v2sigma2, 
                                                 CMultiTimer&        timer) const;

    /**
     Integrates LDA contribution to (third-order) Kxc matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param v2rho2 the 2nd-order functional derivative wrt density.
     @param v3rho3 the 3rd-order functional derivative wrt density.
     @param rwDensityGridQuad the products of one-time transformed densities on grid points.
     @param rw2DensityMatrix the two-time transformed densities on grid points.
     @param iFock the index of the AO Fock matrix.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialKxcFockForLDA(const int32_t              npoints,
                                                const double*              weights,
                                                const CDenseMatrix&        gtoValues,
                                                const double*              v2rho2,
                                                const double*              v3rho3,
                                                const CDensityGridQuad&    rwDensityGridQuad,
                                                const CDensityGrid&        rw2DensityGrid,
                                                const int32_t              iFock,
                                                CMultiTimer&               timer) const;

    /**
     Integrates GGA contribution to (third-order) Kxc matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param gtoValuesX the GTO gradient X values on grid points.
     @param gtoValuesY the GTO gradient Y values on grid points.
     @param gtoValuesZ the GTO gradient Z values on grid points.
     @param rhograd the density gradient.
     @param vsigma the 1st-order functional derivative wrt sigma.
     @param v2rho2 the 2nd-order functional derivative wrt rho.
     @param v2rhosigma the 2nd-order functional derivative wrt rho and sigma.
     @param v2sigma2 the 2nd-order functional derivative wrt sigma.
     @param v3rho3 the 3rd-order functional derivative wrt rho.
     @param v3rho2sigma the 3rd-order functional derivative wrt rho and sigma.
     @param v3rhosigma2 the 3rd-order functional derivative wrt rho and sigma.
     @param v3sigma3 the 3rd-order functional derivative wrt sigma.
     @param rwDensityGridQuad the products of one-time transformed densities on grid points.
     @param rw2DensityMatrix the two-time transformed densities on grid points.
     @param iFock the index of the AO Fock matrix.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialKxcFockForGGA(const int32_t              npoints,
                                                const double*              weights,
                                                const CDenseMatrix&        gtoValues,
                                                const CDenseMatrix&        gtoValuesX,
                                                const CDenseMatrix&        gtoValuesY,
                                                const CDenseMatrix&        gtoValuesZ,
                                                const double*              rhograd,
                                                const double*              vsigma,
                                                const double*              v2rho2,
                                                const double*              v2rhosigma,
                                                const double*              v2sigma2,
                                                const double*              v3rho3,
                                                const double*              v3rho2sigma,
                                                const double*              v3rhosigma2,
                                                const double*              v3sigma3,
                                                const CDensityGridQuad&    rwDensityGridQuad,
                                                const CDensityGrid&        rw2DensityGrid,
                                                const int32_t              iFock,
                                                CMultiTimer&               timer) const;

    /**
     Integrates meta-GGA contribution to (third-order) Kxc matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param gtoValuesX the GTO gradient X values on grid points.
     @param gtoValuesY the GTO gradient Y values on grid points.
     @param gtoValuesZ the GTO gradient Z values on grid points.
     @param rhograd the density gradient.
     @param vsigma the 1st-order functional derivative wrt sigma.
     @param v2rho2 the 2nd-order functional derivative wrt rho.
     @param v2lapl2 ,
     @param v2tau2 ,
     @param v2rholapl ,
     @param v2rhotau ,
     @param v2lapltau ,
     @param v2rhosigma the 2nd-order functional derivative wrt rho and sigma.
     @param v2sigmalapl ,
     @param v2sigmatau ,
     @param v2sigma2 the 2nd-order functional derivative wrt sigma.
     @param v3rho3 the 3rd-order functional derivative wrt rho.
     @param v3rho2sigma the 3rd-order functional derivative wrt rho and sigma.
     @param v3rho2lapl ,
     @param v3rho2tau ,
     @param v3rhosigma2 the 3rd-order functional derivative wrt rho and sigma.
     @param v3rhosigmalapl ,
     @param v3rhosigmatau ,
     @param v3rholapl2 ,
     @param v3rholapltau ,
     @param v3rhotau2 ,
     @param v3sigma3 the 3rd-order functional derivative wrt sigma.
     @param v3sigma2lapl ,
     @param v3sigma2tau ,
     @param v3sigmalapl2 ,
     @param v3sigmalapltau ,
     @param v3sigmatau2 ,
     @param v3lapl3 ,
     @param v3lapl2tau ,
     @param v3lapltau2 ,
     @param v3tau3 ,
     @param rwDensityGridQuad the products of one-time transformed densities on grid points.
     @param rw2DensityMatrix the two-time transformed densities on grid points.
     @param iFock the index of the AO Fock matrix.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialKxcFockForMGGA(const int32_t           npoints,
                                                 const double*           weights,
                                                 const CDenseMatrix&     gtoValues,
                                                 const CDenseMatrix&     gtoValuesX,
                                                 const CDenseMatrix&     gtoValuesY,
                                                 const CDenseMatrix&     gtoValuesZ,
                                                 const double*           rhograd,
                                                 const double*           vsigma,
                                                 const double*           v2rho2,
                                                 const double*           v2lapl2,
                                                 const double*           v2tau2,
                                                 const double*           v2rholapl,
                                                 const double*           v2rhotau,
                                                 const double*           v2lapltau,
                                                 const double*           v2rhosigma,
                                                 const double*           v2sigmalapl,
                                                 const double*           v2sigmatau,
                                                 const double*           v2sigma2,
                                                 const double*           v3rho3,
                                                 const double*           v3rho2sigma,
                                                 const double*           v3rho2lapl,
                                                 const double*           v3rho2tau,
                                                 const double*           v3rhosigma2,
                                                 const double*           v3rhosigmalapl,
                                                 const double*           v3rhosigmatau,
                                                 const double*           v3rholapl2,
                                                 const double*           v3rholapltau,
                                                 const double*           v3rhotau2,
                                                 const double*           v3sigma3,
                                                 const double*           v3sigma2lapl,
                                                 const double*           v3sigma2tau,
                                                 const double*           v3sigmalapl2,
                                                 const double*           v3sigmalapltau,
                                                 const double*           v3sigmatau2,
                                                 const double*           v3lapl3,
                                                 const double*           v3lapl2tau,
                                                 const double*           v3lapltau2,
                                                 const double*           v3tau3,
                                                 const CDensityGridQuad& rwDensityGridQuad,
                                                 const CDensityGrid&     rw2DensityGrid,
                                                 const int32_t           iFock,
                                                 CMultiTimer&            timer) const;

    /**
     Integrates LDA contribution to (third-order) Kxc matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param v2rho2 the 2nd-order functional derivative wrt density.
     @param v3rho3 the 3rd-order functional derivative wrt density.
     @param rwDensityGridCubic the products of one and two-time transformed densities on grid points.
     @param rw2DensityMatrix the two-time transformed densities on grid points.
     @param iFock the index of the AO Fock matrix.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialKxcFockForLDA2(const int32_t              npoints,
                                                 const double*              weights,
                                                 const CDenseMatrix&        gtoValues,
                                                 const double*              v2rho2,
                                                 const double*              v3rho3,
                                                 const CDensityGridCubic&   rwDensityGridCubic,
                                                 const CDensityGrid&        rw2DensityGrid,
                                                 const int32_t              iFock,
                                                 CMultiTimer&               timer) const;

    /**
     Integrates GGA contribution to (third-order) Kxc matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param gtoValuesX the GTO gradient X values on grid points.
     @param gtoValuesY the GTO gradient Y values on grid points.
     @param gtoValuesZ the GTO gradient Z values on grid points.
     @param rhograd the density gradient.
     @param vsigma the 1st-order functional derivative wrt sigma.
     @param v2rho2 the 2nd-order functional derivative wrt rho.
     @param v2rhosigma the 2nd-order functional derivative wrt rho and sigma.
     @param v2sigma2 the 2nd-order functional derivative wrt sigma.
     @param v3rho3 the 3rd-order functional derivative wrt rho.
     @param v3rho2sigma the 3rd-order functional derivative wrt rho and sigma.
     @param v3rhosigma2 the 3rd-order functional derivative wrt rho and sigma.
     @param v3sigma3 the 3rd-order functional derivative wrt sigma.
     @param rwDensityGridCubic the products of one and two-time transformed densities on grid points.
     @param rw2DensityMatrix the two-time transformed densities on grid points.
     @param iFock the index of the AO Fock matrix.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialKxcFockForGGA2(const int32_t            npoints,
                                                 const double*            weights,
                                                 const CDenseMatrix&      gtoValues,
                                                 const CDenseMatrix&      gtoValuesX,
                                                 const CDenseMatrix&      gtoValuesY,
                                                 const CDenseMatrix&      gtoValuesZ,
                                                 const double*            rhograd,
                                                 const double*            vsigma,
                                                 const double*            v2rho2,
                                                 const double*            v2rhosigma,
                                                 const double*            v2sigma2,
                                                 const double*            v3rho3,
                                                 const double*            v3rho2sigma,
                                                 const double*            v3rhosigma2,
                                                 const double*            v3sigma3,
                                                 const CDensityGridCubic& rwDensityGridCubic,
                                                 const CDensityGrid&      rw2DensityGrid,
                                                 const int32_t            iFock,
                                                 CMultiTimer&             timer) const;

    /**
     Integrates meta-GGA contribution to (third-order) Kxc matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param gtoValuesX the GTO gradient X values on grid points.
     @param gtoValuesY the GTO gradient Y values on grid points.
     @param gtoValuesZ the GTO gradient Z values on grid points.
     @param rhograd the density gradient.
     @param vsigma the 1st-order functional derivative wrt sigma.
     @param v2rho2 the 2nd-order functional derivative wrt rho.
     @param v2lapl2 , 
     @param v2tau2 , 
     @param v2rholapl , 
     @param v2rhotau ,
     @param v2lapltau , 
     @param v2rhosigma the 2nd-order functional derivative wrt rho and sigma.
     @param v2sigmalapl , 
     @param v2sigmatau ,
     @param v2sigma2 the 2nd-order functional derivative wrt sigma.
     @param v3rho3 the 3rd-order functional derivative wrt rho.
     @param v3rho2sigma the 3rd-order functional derivative wrt rho and sigma.
     @param v3rho2lapl ,
     @param v3rho2tau ,
     @param v3rhosigma2 the 3rd-order functional derivative wrt rho and sigma.
     @param v3rhosigmalapl ,
     @param v3rhosigmatau ,
     @param v3rholapl2 ,
     @param v3rholapltau ,
     @param v3rhotau2 ,
     @param v3sigma3 the 3rd-order functional derivative wrt sigma.
     @param v3sigma2lapl ,
     @param v3sigma2tau ,
     @param v3sigmalapl2 ,
     @param v3sigmalapltau ,
     @param v3sigmatau2 ,
     @param v3lapl3 ,
     @param v3lapl2tau ,
     @param v3lapltau2 ,
     @param v3tau3 ,
     @param rwDensityGridQuad the products of one-time transformed densities on grid points.
     @param rw2DensityMatrix the two-time transformed densities on grid points.
     @param iFock the index of the AO Fock matrix.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialKxcFockForMGGA2(const int32_t            npoints, 
                                                  const double*            weights, 
                                                  const CDenseMatrix&      gtoValues,
                                                  const CDenseMatrix&      gtoValuesX, 
                                                  const CDenseMatrix&      gtoValuesY, 
                                                  const CDenseMatrix&      gtoValuesZ,
                                                  const double*            rhograd, 
                                                  const double*            vsigma,
                                                  const double*            v2rho2,
                                                  const double*            v2lapl2, 
                                                  const double*            v2tau2, 
                                                  const double*            v2rholapl, 
                                                  const double*            v2rhotau,
                                                  const double*            v2lapltau, 
                                                  const double*            v2rhosigma, 
                                                  const double*            v2sigmalapl, 
                                                  const double*            v2sigmatau,
                                                  const double*            v2sigma2,
                                                  const double*            v3rho3,
                                                  const double*            v3rho2sigma,
                                                  const double*            v3rho2lapl,
                                                  const double*            v3rho2tau,
                                                  const double*            v3rhosigma2,
                                                  const double*            v3rhosigmalapl,
                                                  const double*            v3rhosigmatau,
                                                  const double*            v3rholapl2,
                                                  const double*            v3rholapltau,
                                                  const double*            v3rhotau2,
                                                  const double*            v3sigma3,
                                                  const double*            v3sigma2lapl,
                                                  const double*            v3sigma2tau,
                                                  const double*            v3sigmalapl2,
                                                  const double*            v3sigmalapltau,
                                                  const double*            v3sigmatau2,
                                                  const double*            v3lapl3,
                                                  const double*            v3lapl2tau,
                                                  const double*            v3lapltau2,
                                                  const double*            v3tau3,
                                                  const CDensityGridCubic& rwDensityGridcubic,
                                                  const CDensityGrid&      rw2DensityGrid,
                                                  const int32_t            iFock,
                                                  CMultiTimer&             timer) const;

    /**
     Integrates LDA contribution to (fourth-order) Kxc matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param v2rho2 the 2nd-order functional derivative wrt density.
     @param v3rho3 the 3rd-order functional derivative wrt density.
     @param v4rho4 the 4rd-order functional derivative wrt density.
     @param rwDensityGridCubic the products of one and two-time transformed densities on grid points.
     @param rw3DensityMatrix the three-time transformed densities on grid points.
     @param iFock the index of the AO Fock matrix.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialLxcFockForLDA(const int32_t            npoints,
                                                const double*            weights,
                                                const CDenseMatrix&      gtoValues,
                                                const double*            v2rho2,
                                                const double*            v3rho3,
                                                const double*            v4rho4,
                                                const CDensityGridCubic& rwDensityGridCubic,
                                                const CDensityGrid&      rw3DensityGrid,
                                                const int32_t            iFock,
                                                CMultiTimer&             timer) const;

    /**
     Integrates GGA contribution to (fourth-order) Lxc matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param gtoValuesX the GTO gradient X values on grid points.
     @param gtoValuesY the GTO gradient Y values on grid points.
     @param gtoValuesZ the GTO gradient Z values on grid points.
     @param rhograd the density gradient.
     @param vsigma the 1st-order functional derivative wrt sigma.
     @param v2rho2 the 2nd-order functional derivative wrt rho.
     @param v2rhosigma the 2nd-order functional derivative wrt rho and sigma.
     @param v2sigma2 the 2nd-order functional derivative wrt sigma.
     @param v3rho3 the 3rd-order functional derivative wrt rho.
     @param v3rho2sigma the 3rd-order functional derivative wrt rho and sigma.
     @param v3rhosigma2 the 3rd-order functional derivative wrt rho and sigma.
     @param v3sigma3 the 3rd-order functional derivative wrt sigma.
     @param v4rho4 ,
     @param v4rho3sigma ,
     @param v4rho2sigma2 ,
     @param v4rhosigma3 ,
     @param v4sigma4 ,                   
     @param rwDensityGridCubic the products of one and two-time transformed densities on grid points.
     @param rw3DensityMatrix the three-time transformed densities on grid points.
     @param iFock the index of the AO Fock matrix.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialLxcFockForGGA(const int32_t            npoints,
                                                const double*            weights,
                                                const CDenseMatrix&      gtoValues,
                                                const CDenseMatrix&      gtoValuesX,
                                                const CDenseMatrix&      gtoValuesY,
                                                const CDenseMatrix&      gtoValuesZ,
                                                const double*            rhograd,
                                                const double*            vsigma,
                                                const double*            v2rho2,
                                                const double*            v2rhosigma,
                                                const double*            v2sigma2,
                                                const double*            v3rho3,
                                                const double*            v3rho2sigma,
                                                const double*            v3rhosigma2,
                                                const double*            v3sigma3,
                                                const double*            v4rho4,
                                                const double*            v4rho3sigma,
                                                const double*            v4rho2sigma2,
                                                const double*            v4rhosigma3,
                                                const double*            v4sigma4,                   
                                                const CDensityGridCubic& rwDensityGridCubic,
                                                const CDensityGrid&      rw3DensityGrid,
                                                const int32_t            iFock,
                                                CMultiTimer&             timer) const;

    /**
     Integrates meta-GGA contribution to (fourth-order) Lxc matrix.

     @param npoints the number of grid points.
     @param weights the weights of grid points.
     @param gtoValues the GTO values on grid points.
     @param gtoValuesX the GTO gradient X values on grid points.
     @param gtoValuesY the GTO gradient Y values on grid points.
     @param gtoValuesZ the GTO gradient Z values on grid points.
     @param rhograd the density gradient.
     @param vsigma the 1st-order functional derivative wrt sigma.
     @param v2rho2 ,
     @param v2lapl2 , 
     @param v2tau2 , 
     @param v2rholapl , 
     @param v2rhotau ,
     @param v2lapltau , 
     @param v2rhosigma , 
     @param v2sigmalapl , 
     @param v2sigmatau ,
     @param v2sigma2 ,
     @param v3rho3 ,
     @param v3rho2sigma ,
     @param v3rho2lapl ,
     @param v3rho2tau ,
     @param v3rhosigma2 ,
     @param v3rhosigmalapl ,
     @param v3rhosigmatau ,
     @param v3rholapl2 ,
     @param v3rholapltau ,
     @param v3rhotau2 ,
     @param v3sigma3 ,
     @param v3sigma2lapl ,
     @param v3sigma2tau ,
     @param v3sigmalapl2 ,
     @param v3sigmalapltau ,
     @param v3sigmatau2 ,
     @param v3lapl3 ,
     @param v3lapl2tau ,
     @param v3lapltau2 ,
     @param v3tau3 , 
     @param v4rho4 ,
     @param v4rho3sigma ,
     @param v4rho3lapl ,
     @param v4rho3tau ,
     @param v4rho2sigma2 ,
     @param v4rho2sigmalapl ,
     @param v4rho2sigmatau ,
     @param v4rho2lapl2 ,
     @param v4rho2lapltau ,
     @param v4rho2tau2 ,
     @param v4rhosigma3 ,
     @param v4rhosigma2lapl ,
     @param v4rhosigma2tau ,
     @param v4rhosigmalapl2 ,
     @param v4rhosigmalapltau ,
     @param v4rhosigmatau2 ,
     @param v4rholapl3 ,
     @param v4rholapl2tau ,
     @param v4rholapltau2 ,
     @param v4rhotau3 ,
     @param v4sigma4 ,
     @param v4sigma3lapl ,
     @param v4sigma3tau ,
     @param v4sigma2lapl2 ,
     @param v4sigma2lapltau ,
     @param v4sigma2tau2 ,
     @param v4sigmalapl3 ,
     @param v4sigmalapl2tau ,
     @param v4sigmalapltau2 ,
     @param v4sigmatau3 ,
     @param v4lapl4 ,
     @param v4lapl3tau ,
     @param v4lapl2tau2 ,
     @param v4lapltau3 ,
     @param v4tau4 ,                 
     @param rwDensityGridQuad the products of one-time transformed densities on grid points.
     @param rw3DensityMatrix the three-time transformed densities on grid points.
     @param iFock the index of the AO Fock matrix.
     @param timer the timer.
     @return the contribution as a CDenseMatrix object.
     */
    CDenseMatrix _integratePartialLxcFockForMGGA(const int32_t            npoints, 
                                                 const double*            weights, 
                                                 const CDenseMatrix&      gtoValues,
                                                 const CDenseMatrix&      gtoValuesX, 
                                                 const CDenseMatrix&      gtoValuesY, 
                                                 const CDenseMatrix&      gtoValuesZ,
                                                 const double*            rhograd,
                                                 const double*            vsigma, 
                                                 const double*            v2rho2,
                                                 const double*            v2lapl2, 
                                                 const double*            v2tau2, 
                                                 const double*            v2rholapl, 
                                                 const double*            v2rhotau,
                                                 const double*            v2lapltau, 
                                                 const double*            v2rhosigma, 
                                                 const double*            v2sigmalapl, 
                                                 const double*            v2sigmatau,
                                                 const double*            v2sigma2,
                                                 const double*            v3rho3,
                                                 const double*            v3rho2sigma,
                                                 const double*            v3rho2lapl,
                                                 const double*            v3rho2tau,
                                                 const double*            v3rhosigma2,
                                                 const double*            v3rhosigmalapl,
                                                 const double*            v3rhosigmatau,
                                                 const double*            v3rholapl2,
                                                 const double*            v3rholapltau,
                                                 const double*            v3rhotau2,
                                                 const double*            v3sigma3,
                                                 const double*            v3sigma2lapl,
                                                 const double*            v3sigma2tau,
                                                 const double*            v3sigmalapl2,
                                                 const double*            v3sigmalapltau,
                                                 const double*            v3sigmatau2,
                                                 const double*            v3lapl3,
                                                 const double*            v3lapl2tau,
                                                 const double*            v3lapltau2,
                                                 const double*            v3tau3, 
                                                 const double*            v4rho4,
                                                 const double*            v4rho3sigma,
                                                 const double*            v4rho3lapl,
                                                 const double*            v4rho3tau,
                                                 const double*            v4rho2sigma2,
                                                 const double*            v4rho2sigmalapl,
                                                 const double*            v4rho2sigmatau,
                                                 const double*            v4rho2lapl2,
                                                 const double*            v4rho2lapltau,
                                                 const double*            v4rho2tau2,
                                                 const double*            v4rhosigma3,
                                                 const double*            v4rhosigma2lapl,
                                                 const double*            v4rhosigma2tau,
                                                 const double*            v4rhosigmalapl2,
                                                 const double*            v4rhosigmalapltau,
                                                 const double*            v4rhosigmatau2,
                                                 const double*            v4rholapl3,
                                                 const double*            v4rholapl2tau,
                                                 const double*            v4rholapltau2,
                                                 const double*            v4rhotau3,
                                                 const double*            v4sigma4,
                                                 const double*            v4sigma3lapl,
                                                 const double*            v4sigma3tau,
                                                 const double*            v4sigma2lapl2,
                                                 const double*            v4sigma2lapltau,
                                                 const double*            v4sigma2tau2,
                                                 const double*            v4sigmalapl3,
                                                 const double*            v4sigmalapl2tau,
                                                 const double*            v4sigmalapltau2,
                                                 const double*            v4sigmatau3,
                                                 const double*            v4lapl4,
                                                 const double*            v4lapl3tau,
                                                 const double*            v4lapl2tau2,
                                                 const double*            v4lapltau3,
                                                 const double*            v4tau4,                 
                                                 const CDensityGridCubic& rwDensityGridCubic,
                                                 const CDensityGrid&      rw3DensityGrid,
                                                 const int32_t            iFock,
                                                 CMultiTimer&             timer) const;

   public:
    /**
     Creates an XC integrator object using MPI info.

     @param comm the MPI communicator.
     */
    CXCNewIntegrator(MPI_Comm comm);

    /**
     Integrates first-order exchange-correlation functional contribution to AO
     Kohn-Sham matrix.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param densityMatrix the AO density matrix object.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @return the AO Kohn-Sham matrix.
     */
    CAOKohnShamMatrix integrateVxcFock(const CMolecule&        molecule,
                                       const CMolecularBasis&  basis,
                                       const CAODensityMatrix& densityMatrix,
                                       const CMolecularGrid&   molecularGrid,
                                       const std::string&      xcFuncLabel) const;

    /**
     Integrates second-order exchange-correlation functional contribution to AO
     Fock matrix.

     @param aoFockMatrix the AO Fock matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the perturbed density matrix.
     @param gsDensityMatrix the ground state density matrix.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     */
    void integrateFxcFock(CAOFockMatrix&          aoFockMatrix,
                          const CMolecule&        molecule,
                          const CMolecularBasis&  basis,
                          const CAODensityMatrix& rwDensityMatrix,
                          const CAODensityMatrix& gsDensityMatrix,
                          const CMolecularGrid&   molecularGrid,
                          const std::string&      xcFuncLabel) const;

    /**
     Integrates third-order exchange-correlation functional contribution to AO
     Fock matrix in quadratic response.

     @param aoFockMatrix the AO Fock matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the one-time transformed densities.
     @param rw2DensityMatrix the two-time transformed densities.
     @param gsDensityMatrix the ground state density matrix.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @param quadMode a string that specifies which densities should be combined.
     */
    void integrateKxcFock(CAOFockMatrix&          aoFockMatrix,
                          const CMolecule&        molecule,
                          const CMolecularBasis&  basis,
                          const CAODensityMatrix& rwDensityMatrix,
                          const CAODensityMatrix& rw2DensityMatrix,
                          const CAODensityMatrix& gsDensityMatrix,
                          const CMolecularGrid&   molecularGrid,
                          const std::string&      xcFuncLabel,
                          const std::string&      quadMode) const;

    /**
     Integrates fourth-order exchnage-correlation functional contribution to AO
     Fock matrix in quadratic response.

     @param aoFockMatrix the AO Fock matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the one-time transformed densities.
     @param rw2DensityMatrix the two-time transformed densities.
     @param rw3DensityMatrix the three-time transformed densities.
     @param gsDensityMatrix the ground state density matrix.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @param quadMode a string that specifies which densities should be combined.
     */
    void integrateLxcFock(CAOFockMatrix&          aoFockMatrix,
                          const CMolecule&        molecule,
                          const CMolecularBasis&  basis,
                          const CAODensityMatrix& rwDensityMatrix,
                          const CAODensityMatrix& rw2DensityMatrix,
                          const CAODensityMatrix& rw3DensityMatrix,
                          const CAODensityMatrix& gsDensityMatrix,
                          const CMolecularGrid&   molecularGrid,
                          const std::string&      xcFuncLabel,
                          const std::string&      cubeMode) const;

    /**
     Integrates fourth-order exchnage-correlation functional contribution to AO
     Fock matrix in quadratic response.

     @param aoFockMatrix the AO Fock matrix.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param rwDensityMatrix the one-time transformed densities.
     @param rw2DensityMatrix the two-time transformed densities.
     @param rw3DensityMatrix the three-time transformed densities.
     @param gsDensityMatrix the ground state density matrix.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     @param quadMode a string that specifies which densities should be combined.
     */
    void integrateKxcLxcFock(CAOFockMatrix&          aoFockMatrix,
                             const CMolecule&        molecule,
                             const CMolecularBasis&  basis,
                             const CAODensityMatrix& rwDensityMatrix,
                             const CAODensityMatrix& rw2DensityMatrix,
                             const CAODensityMatrix& rw3DensityMatrix,
                             const CAODensityMatrix& gsDensityMatrix,
                             const CMolecularGrid&   molecularGrid,
                             const std::string&      xcFuncLabel,
                             const std::string&      cubeMode) const;

    /**
     Integrates first-order pair-density functional contribution to AO
     Fock matrix and MO "Q-matrix".

     @param aoFockMatrix the AO Fock matrix.
     @param moTwoBodyGradient the MO Two-body energy gradient term.
     @param molecule the molecule.
     @param basis the molecular basis.
     @param densityMatrix the AO density matrix object.
     @param TwoBodyDensityMatrix the MO two-body active density matrix.
     @param ActiveMOs the active molecular orbitals.
     @param molecularGrid the molecular grid.
     @param xcFuncLabel the label of exchange-correlation functional.
     */
    void integrateVxcPDFT(CAOKohnShamMatrix&      aoFockMatrix,
                          CDense4DTensor&         moTwoBodyGradient,
                          const CMolecule&        molecule,
                          const CMolecularBasis&  basis,
                          const CAODensityMatrix& DensityMatrix,
                          const CDense4DTensor&   TwoBodyDensityMatrix,
                          const CDenseMatrix&     ActiveMOs,
                          const CMolecularGrid&   molecularGrid,
                          const std::string&      xcFuncLabel) const;

    /**
     Computes GTOs values on grid points.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @return the GTO values on grid points.
     */
    CDenseMatrix computeGtoValuesOnGridPoints(const CMolecule& molecule, const CMolecularBasis& basis, const CMolecularGrid& molecularGrid) const;

    /**
     Computes GTOs values and derivatives on grid points.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param molecularGrid the molecular grid.
     @return the GTO values and derivatives on grid points.
     */
    std::vector<CDenseMatrix> computeGtoValuesAndDerivativesOnGridPoints(const CMolecule&       molecule,
                                                                         const CMolecularBasis& basis,
                                                                         const CMolecularGrid&  molecularGrid) const;

    /**
     Computes GTOs values and derivatives on grid points.

     @param molecule the molecule.
     @param basis the molecular basis.
     @param npoints the number of grid points.
     @param xcoords the X coordinates of grid points.
     @param ycoords the Y coordinates of grid points.
     @param zcoords the Z coordinates of grid points.
     @return the GTO values and derivatives on grid points.
     */
    std::vector<CDenseMatrix> computeGtoValuesAndDerivativesOnGridPoints(const CMolecule&       molecule,
                                                                         const CMolecularBasis& basis,
                                                                         const int32_t          npoints,
                                                                         const double*          xcoords,
                                                                         const double*          ycoords,
                                                                         const double*          zcoords) const;
};

#endif /* XCNewIntegrator_hpp */
