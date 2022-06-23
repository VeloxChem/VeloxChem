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

#include <string>

#include "AODensityMatrix.hpp"
#include "DenseMatrix.hpp"
#include "DensityGradientGridDriver.hpp"
#include "DensityGridDriver.hpp"
#include "DensityGridQuad.hpp"
#include "DensityMatrixType.hpp"
#include "ErrorHandler.hpp"
#include "FunctionalParser.hpp"
#include "GtoFunc.hpp"
#include "MpiFunc.hpp"
#include "OMPTasks.hpp"

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
CXCMolecularGradient::integrateVxcGradient(const std::vector<int32_t>& idsAtomic,
                                           const CAODensityMatrix&     aoDensityMatrix,
                                           const CMolecule&            molecule,
                                           const CMolecularBasis&      basis,
                                           const CMolecularGrid&       molecularGrid,
                                           const std::string&          xcFuncLabel) const
{
    return integrateVxcGradient(idsAtomic, aoDensityMatrix, aoDensityMatrix, molecule, basis, molecularGrid, xcFuncLabel);
}

CDenseMatrix
CXCMolecularGradient::integrateVxcGradient(const std::vector<int32_t>& idsAtomic,
                                           const CAODensityMatrix&     rwDensityMatrix,
                                           const CAODensityMatrix&     gsDensityMatrix,
                                           const CMolecule&            molecule,
                                           const CMolecularBasis&      basis,
                                           const CMolecularGrid&       molecularGrid,
                                           const std::string&          xcFuncLabel) const
{
    // parse exchange-correlation functional data

    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    // generate reference density grid

    CDensityGridDriver dgdrv(_locComm);

    auto refdengrid = dgdrv.generate(gsDensityMatrix, molecule, basis, molecularGrid, fvxc.getFunctionalType());

    // create molecular gradient

    const auto natoms = static_cast<int32_t>(idsAtomic.size());

    CDenseMatrix molgrad(3, natoms);

    auto mgradx = molgrad.row(0);

    auto mgrady = molgrad.row(1);

    auto mgradz = molgrad.row(2);

    if (rwDensityMatrix.isClosedShell())
    {
        // generate screened molecular and density grids

        CMolecularGrid mgrid(molecularGrid);

        CDensityGrid dgrid;

        refdengrid.getScreenedGridsPair(dgrid, mgrid, 0, _thresholdOfDensity, fvxc.getFunctionalType());

        auto gw = mgrid.getWeights();

        const auto gpoints = mgrid.getNumberOfGridPoints();

        // allocate XC gradient grid

        CXCGradientGrid vxcgrid(mgrid.getNumberOfGridPoints(), dgrid.getDensityGridType(), fvxc.getFunctionalType());

        // compute exchange-correlation functional derivative

        fvxc.compute(vxcgrid, dgrid);

        // set up pointers to exchange-correlation functional derivatives

        auto grhoa = vxcgrid.xcGradientValues(xcvars::rhoa);

        auto ggrada = vxcgrid.xcGradientValues(xcvars::grada);

        auto ggradab = vxcgrid.xcGradientValues(xcvars::gradab);

        // set up pointers to density gradient norms

        auto ngrada = dgrid.alphaDensityGradient(0);

        auto gradax = dgrid.alphaDensityGradientX(0);

        auto graday = dgrid.alphaDensityGradientY(0);

        auto gradaz = dgrid.alphaDensityGradientZ(0);

        // set up density gradient grid driver

        CDensityGradientGridDriver graddrv(_locComm);

        for (int32_t i = 0; i < natoms; i++)
        {
            auto gradgrid = graddrv.generate(rwDensityMatrix, molecule, basis, mgrid, fvxc.getFunctionalType(), idsAtomic[i]);

            double gatmx = 0.0;

            double gatmy = 0.0;

            double gatmz = 0.0;

            // compute LDA contribution to molecular gradient

            if (fvxc.getFunctionalType() == xcfun::lda)
            {
                // set up pointers to density gradient grid

                const auto gdenx = gradgrid.getComponent(0);

                const auto gdeny = gradgrid.getComponent(1);

                const auto gdenz = gradgrid.getComponent(2);

                for (int32_t j = 0; j < gpoints; j++)
                {
                    gatmx += gw[j] * grhoa[j] * gdenx[j];

                    gatmy += gw[j] * grhoa[j] * gdeny[j];

                    gatmz += gw[j] * grhoa[j] * gdenz[j];
                }
            }

            // compute GGA contribution to molecular gradient

            if (fvxc.getFunctionalType() == xcfun::gga)
            {
                // set up pointers to density gradient grid

                const auto gdenx = gradgrid.getComponent(0);

                const auto gdeny = gradgrid.getComponent(1);

                const auto gdenz = gradgrid.getComponent(2);

                const auto gdenxx = gradgrid.getComponent(3);

                const auto gdenxy = gradgrid.getComponent(4);

                const auto gdenxz = gradgrid.getComponent(5);

                const auto gdenyx = gradgrid.getComponent(6);

                const auto gdenyy = gradgrid.getComponent(7);

                const auto gdenyz = gradgrid.getComponent(8);

                const auto gdenzx = gradgrid.getComponent(9);

                const auto gdenzy = gradgrid.getComponent(10);

                const auto gdenzz = gradgrid.getComponent(11);

                for (int32_t j = 0; j < gpoints; j++)
                {
                    // contribution from \nabla_A (\phi_mu \phi_nu)

                    double prefac = gw[j] * grhoa[j];

                    gatmx += prefac * gdenx[j];

                    gatmy += prefac * gdeny[j];

                    gatmz += prefac * gdenz[j];

                    // contribution from \nabla_A (\nabla (\phi_mu \phi_nu))

                    prefac = gw[j] * (ggrada[j] / ngrada[j] + ggradab[j]);

                    gatmx += prefac * (gradax[j] * gdenxx[j] + graday[j] * gdenxy[j] + gradaz[j] * gdenxz[j]);

                    gatmy += prefac * (gradax[j] * gdenyx[j] + graday[j] * gdenyy[j] + gradaz[j] * gdenyz[j]);

                    gatmz += prefac * (gradax[j] * gdenzx[j] + graday[j] * gdenzy[j] + gradaz[j] * gdenzz[j]);
                }
            }

            if (fvxc.getFunctionalType() == xcfun::mgga)
            {
                // not implemented

                std::string errmgga("XCMolecularGradient.integrateVxcGradient: Not implemented for meta-GGA");

                errors::assertMsgCritical(false, errmgga);
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

        std::string erropenshell("XCMolecularGradient.integrateVxcGradient: Not implemented for open-shell");

        errors::assertMsgCritical(false, erropenshell);
    }

    // done with molecular gradient

    CDenseMatrix molgrad_T(natoms, 3);

    for (int32_t i = 0; i < natoms; i++)
    {
        for (int32_t d = 0; d < 3; d++)
        {
            molgrad_T.row(i)[d] = molgrad.row(d)[i];
        }
    }

    return molgrad_T;
}

CDenseMatrix
CXCMolecularGradient::integrateFxcGradient(const CAODensityMatrix& rwDensityMatrixOne,
                                           const CAODensityMatrix& rwDensityMatrixTwo,
                                           const CAODensityMatrix& gsDensityMatrix,
                                           const CMolecule&        molecule,
                                           const CMolecularBasis&  basis,
                                           const CMolecularGrid&   molecularGrid,
                                           const std::string&      xcFuncLabel) const
{
    // parse exchange-correlation functional data

    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    // generate reference density grid

    CDensityGridDriver dgdrv(_locComm);

    auto refdengrid = dgdrv.generate(gsDensityMatrix, molecule, basis, molecularGrid, fvxc.getFunctionalType());

    // create molecular gradient

    const auto natoms = molecule.getNumberOfAtoms();

    CDenseMatrix molgrad(3, natoms);

    molgrad.zero();

    if (rwDensityMatrixOne.isClosedShell())
    {
        // generate screened molecular and density grids

        CMolecularGrid mgrid(molecularGrid);

        CDensityGrid gsdengrid;

        refdengrid.getScreenedGridsPair(gsdengrid, mgrid, 0, _thresholdOfDensity, fvxc.getFunctionalType());

        // allocate XC gradient/hessian grids

        CXCGradientGrid vxcgrid(mgrid.getNumberOfGridPoints(), gsdengrid.getDensityGridType(), fvxc.getFunctionalType());

        CXCHessianGrid vxc2grid(mgrid.getNumberOfGridPoints(), gsdengrid.getDensityGridType(), fvxc.getFunctionalType());

        // compute exchange-correlation functional derivatives

        fvxc.compute(vxcgrid, gsdengrid);

        fvxc.compute(vxc2grid, gsdengrid);

        // create perturbed density grid

        auto rwdengrid = dgdrv.generate(rwDensityMatrixOne, molecule, basis, mgrid, fvxc.getFunctionalType());

        // compute Fxc contribution to molecular gradient

        _compFxcContrib(molgrad, molecule, basis, fvxc.getFunctionalType(), rwDensityMatrixTwo, mgrid, gsdengrid, rwdengrid, vxcgrid, vxc2grid);
    }
    else
    {
        // not implemented

        std::string erropenshell("XCMolecularGradient.integrateFxcGradient: Not implemented for open-shell");

        errors::assertMsgCritical(false, erropenshell);
    }

    // done with molecular gradient

    return molgrad.transpose();
}

CDenseMatrix
CXCMolecularGradient::integrateGxcGradient(const CAODensityMatrix& rwDensityMatrixOne,
                                           const CAODensityMatrix& rwDensityMatrixTwo,
                                           const CAODensityMatrix& gsDensityMatrix,
                                           const CMolecule&        molecule,
                                           const CMolecularBasis&  basis,
                                           const CMolecularGrid&   molecularGrid,
                                           const std::string&      xcFuncLabel) const
{
    // parse exchange-correlation functional data

    auto fvxc = vxcfuncs::getExchangeCorrelationFunctional(xcFuncLabel);

    // generate reference density grid

    CDensityGridDriver dgdrv(_locComm);

    auto refdengrid = dgdrv.generate(gsDensityMatrix, molecule, basis, molecularGrid, fvxc.getFunctionalType());

    // create molecular gradient

    const auto natoms = molecule.getNumberOfAtoms();

    CDenseMatrix molgrad(3, natoms);

    molgrad.zero();

    if (rwDensityMatrixOne.isClosedShell())
    {
        // generate screened molecular and density grids

        CMolecularGrid mgrid(molecularGrid);

        CDensityGrid gsdengrid;

        refdengrid.getScreenedGridsPair(gsdengrid, mgrid, 0, _thresholdOfDensity, fvxc.getFunctionalType());

        // allocate XC gradient/hessian/cubic hessian grids

        CXCGradientGrid vxcgrid(mgrid.getNumberOfGridPoints(), gsdengrid.getDensityGridType(), fvxc.getFunctionalType());

        CXCHessianGrid vxc2grid(mgrid.getNumberOfGridPoints(), gsdengrid.getDensityGridType(), fvxc.getFunctionalType());

        CXCCubicHessianGrid vxc3grid(mgrid.getNumberOfGridPoints(), gsdengrid.getDensityGridType(), fvxc.getFunctionalType());

        // compute exchange-correlation functional derivatives

        fvxc.compute(vxcgrid, gsdengrid);

        fvxc.compute(vxc2grid, gsdengrid);

        fvxc.compute(vxc3grid, gsdengrid);

        // prepare rwDensityMatrix for quadratic response

        auto rwdenmat1 = rwDensityMatrixOne.getReferenceToDensity(0);

        auto rwdenmat2 = rwDensityMatrixTwo.getReferenceToDensity(0);

        CDenseMatrix zerodenmat1(rwdenmat1);

        CDenseMatrix zerodenmat2(rwdenmat2);

        zerodenmat1.zero();

        zerodenmat2.zero();

        CAODensityMatrix rwDensityMatrix(std::vector<CDenseMatrix>({rwdenmat1, zerodenmat1, rwdenmat2, zerodenmat2}), denmat::rest);

        // Compute all and store all products of first-order transformed denisites

        // Note: We use quadratic response (quadMode == "QRF") to calculate
        // third-order functional derivative contribution. The rw2DensityMatrix
        // contains zero matrices and is therefore removed from the following code.
        // Same for rw2dengrid.

        // For "QRF" we have rwDensityMatrix.getNumberOfDensityMatrices() ==
        // 2 * rw2DensityMatrix.getNumberOfDensityMatrices()

        std::string quadMode("QRF");

        int32_t rw2NumberOfDensityMatrices = rwDensityMatrix.getNumberOfDensityMatrices() / 2;

        auto rwdengrid = dgdrv.generate(rwDensityMatrix, molecule, basis, mgrid, fvxc.getFunctionalType());

        auto rwdengridc = CDensityGridQuad(mgrid.getNumberOfGridPoints(), rw2NumberOfDensityMatrices, fvxc.getFunctionalType(), dengrid::ab);

        rwdengridc.DensityProd(rwdengridc, mgrid, rwdengrid, fvxc.getFunctionalType(), rw2NumberOfDensityMatrices, quadMode);

        // compute Gxc contribution to molecular gradient

        _compGxcContrib(
            molgrad, molecule, basis, fvxc.getFunctionalType(), gsDensityMatrix, mgrid, gsdengrid, rwdengridc, vxcgrid, vxc2grid, vxc3grid);
    }
    else
    {
        // not implemented

        std::string erropenshell("XCMolecularGradient.integrateGxcGradient: Not implemented for open-shell");

        errors::assertMsgCritical(false, erropenshell);
    }

    // done with molecular gradient

    return molgrad.transpose();
}

void
CXCMolecularGradient::_compFxcContrib(CDenseMatrix&           molecularGradient,
                                      const CMolecule&        molecule,
                                      const CMolecularBasis&  basis,
                                      const xcfun             xcFuncType,
                                      const CAODensityMatrix& densityMatrix,
                                      const CMolecularGrid&   molecularGrid,
                                      const CDensityGrid&     gsDensityGrid,
                                      const CDensityGrid&     rwDensityGrid,
                                      const CXCGradientGrid&  xcGradientGrid,
                                      const CXCHessianGrid&   xcHessianGrid) const
{
    auto gridWeights = molecularGrid.getWeights();

    const auto nGridPoints = molecularGrid.getNumberOfGridPoints();

    const auto nAtoms = molecule.getNumberOfAtoms();

    // set up pointers to Cartesian components of molecular gradient

    auto mgradx = molecularGradient.row(0);

    auto mgrady = molecularGradient.row(1);

    auto mgradz = molecularGradient.row(2);

    // set up pointers to exchange-correlation functional derivatives

    auto ggrad_a = xcGradientGrid.xcGradientValues(xcvars::grada);

    auto ggrad_c = xcGradientGrid.xcGradientValues(xcvars::gradab);

    auto grho_aa = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::rhoa);

    auto grho_ab = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::rhob);

    auto gmix_aa = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::grada);

    auto gmix_ab = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::gradb);

    auto gmix_ac = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::gradab);

    auto gmix_bc = xcHessianGrid.xcHessianValues(xcvars::rhob, xcvars::gradab);

    auto ggrad_aa = xcHessianGrid.xcHessianValues(xcvars::grada, xcvars::grada);

    auto ggrad_ab = xcHessianGrid.xcHessianValues(xcvars::grada, xcvars::gradb);

    auto ggrad_ac = xcHessianGrid.xcHessianValues(xcvars::grada, xcvars::gradab);

    auto ggrad_bc = xcHessianGrid.xcHessianValues(xcvars::gradb, xcvars::gradab);

    auto ggrad_cc = xcHessianGrid.xcHessianValues(xcvars::gradab, xcvars::gradab);

    // set up pointers to ground state density gradient norms

    auto ngrada = gsDensityGrid.alphaDensityGradient(0);

    auto ngradb = gsDensityGrid.betaDensityGradient(0);

    auto grada_x = gsDensityGrid.alphaDensityGradientX(0);

    auto grada_y = gsDensityGrid.alphaDensityGradientY(0);

    auto grada_z = gsDensityGrid.alphaDensityGradientZ(0);

    auto gradb_x = gsDensityGrid.betaDensityGradientX(0);

    auto gradb_y = gsDensityGrid.betaDensityGradientY(0);

    auto gradb_z = gsDensityGrid.betaDensityGradientZ(0);

    // set up pointers to perturbed density gradient norms

    auto rhowa = rwDensityGrid.alphaDensity(0);

    auto rhowb = rwDensityGrid.betaDensity(0);

    auto gradwa_x = rwDensityGrid.alphaDensityGradientX(0);

    auto gradwa_y = rwDensityGrid.alphaDensityGradientY(0);

    auto gradwa_z = rwDensityGrid.alphaDensityGradientZ(0);

    auto gradwb_x = rwDensityGrid.betaDensityGradientX(0);

    auto gradwb_y = rwDensityGrid.betaDensityGradientY(0);

    auto gradwb_z = rwDensityGrid.betaDensityGradientZ(0);

    // set up density gradient grid driver

    CDensityGradientGridDriver graddrv(_locComm);

    for (int32_t i = 0; i < nAtoms; i++)
    {
        auto gradgrid = graddrv.generate(densityMatrix, molecule, basis, molecularGrid, xcFuncType, i);

        double gatmx = 0.0;

        double gatmy = 0.0;

        double gatmz = 0.0;

        // compute LDA contribution to molecular gradient

        if (xcFuncType == xcfun::lda)
        {
            // set up pointers to density gradient grid

            const auto gdenx = gradgrid.getComponent(0);

            const auto gdeny = gradgrid.getComponent(1);

            const auto gdenz = gradgrid.getComponent(2);

            for (int32_t j = 0; j < nGridPoints; j++)
            {
                double prefac = gridWeights[j] * (grho_aa[j] * rhowa[j] + grho_ab[j] * rhowb[j]);

                gatmx += prefac * gdenx[j];

                gatmy += prefac * gdeny[j];

                gatmz += prefac * gdenz[j];
            }
        }

        // compute GGA contribution to molecular gradient

        if (xcFuncType == xcfun::gga)
        {
            // set up pointers to density gradient grid

            const auto gdenx = gradgrid.getComponent(0);

            const auto gdeny = gradgrid.getComponent(1);

            const auto gdenz = gradgrid.getComponent(2);

            const auto gdenxx = gradgrid.getComponent(3);

            const auto gdenxy = gradgrid.getComponent(4);

            const auto gdenxz = gradgrid.getComponent(5);

            const auto gdenyx = gradgrid.getComponent(6);

            const auto gdenyy = gradgrid.getComponent(7);

            const auto gdenyz = gradgrid.getComponent(8);

            const auto gdenzx = gradgrid.getComponent(9);

            const auto gdenzy = gradgrid.getComponent(10);

            const auto gdenzz = gradgrid.getComponent(11);

            for (int32_t j = 0; j < nGridPoints; j++)
            {
                double w = gridWeights[j];

                double znva = 1.0 / ngrada[j];

                double znvb = 1.0 / ngradb[j];

                double rxa = znva * grada_x[j];

                double rya = znva * grada_y[j];

                double rza = znva * grada_z[j];

                double rxb = znvb * gradb_x[j];

                double ryb = znvb * gradb_y[j];

                double rzb = znvb * gradb_z[j];

                double rxwa = gradwa_x[j];

                double rywa = gradwa_y[j];

                double rzwa = gradwa_z[j];

                double rxwb = gradwb_x[j];

                double rywb = gradwb_y[j];

                double rzwb = gradwb_z[j];

                // GTOs values
                // a0 = bgaos[m] * kgaos[m];
                // ax = bgaox[m] * kgaos[m] + bgaos[m] * kgaox[m];
                // ay = bgaoy[m] * kgaos[m] + bgaos[m] * kgaoy[m];
                // az = bgaoz[m] * kgaos[m] + bgaos[m] * kgaoz[m];

                //  variations of functionals variables

                double zetaa = rxwa * rxa + rywa * rya + rzwa * rza;

                double zetab = rxwb * rxb + rywb * ryb + rzwb * rzb;

                double zetac = grada_x[j] * rxwb + grada_y[j] * rywb

                               + grada_z[j] * rzwb + gradb_x[j] * rxwa

                               + gradb_y[j] * rywa + gradb_z[j] * rzwa;

                // contribution from \nabla_A (\phi_mu \phi_nu)

                // first contribution

                double fac0 = gmix_aa[j] * zetaa + gmix_ab[j] * zetab

                              + gmix_ac[j] * zetac + grho_aa[j] * rhowa[j]

                              + grho_ab[j] * rhowb[j];

                // w * a0 * fac0;   a0 = bgaos[m] * kgaos[m];

                double prefac = w * fac0;

                gatmx += prefac * gdenx[j];

                gatmy += prefac * gdeny[j];

                gatmz += prefac * gdenz[j];

                // contribution from \nabla_A (\nabla (\phi_mu \phi_nu))

                double xcomp = 0.0;

                double ycomp = 0.0;

                double zcomp = 0.0;

                // second contribution

                double facr = gmix_aa[j] * rhowa[j] + gmix_ab[j] * rhowb[j]

                              + ggrad_aa[j] * zetaa + ggrad_ab[j] * zetab + ggrad_ac[j] * zetac;

                // w * facr * ar;   ar = ax * rxa + ay * rya + az * rza;

                prefac = w * facr;

                xcomp += prefac * rxa;

                ycomp += prefac * rya;

                zcomp += prefac * rza;

                // third contribution

                double facz = gmix_ac[j] * rhowa[j] + gmix_bc[j] * rhowb[j]

                              + ggrad_ac[j] * zetaa + ggrad_bc[j] * zetab + ggrad_cc[j] * zetac;

                // w * facz * arb;   arb = ax * grada_x[j] + ay * grada_y[j] + az * grada_z[j];

                prefac = w * facz;

                xcomp += prefac * grada_x[j];

                ycomp += prefac * grada_y[j];

                zcomp += prefac * grada_z[j];

                // fourth contribution

                // w * znva * ggrad_a[j] * ab;
                // ab = ax * rxwa + ay * rywa + az * rzwa - ar * zetaa;
                // ar = ax * rxa + ay * rya + az * rza;

                prefac = w * znva * ggrad_a[j];

                xcomp += prefac * (rxwa - rxa * zetaa);

                ycomp += prefac * (rywa - rya * zetaa);

                zcomp += prefac * (rzwa - rza * zetaa);

                // fifth contribution

                // w * ggrad_c[j] * abw;
                // abw = ax * rxwa + ay * rywa + az * rzwa;

                prefac = w * ggrad_c[j];

                xcomp += prefac * rxwa;

                ycomp += prefac * rywa;

                zcomp += prefac * rzwa;

                gatmx += (xcomp * gdenxx[j] + ycomp * gdenxy[j] + zcomp * gdenxz[j]);

                gatmy += (xcomp * gdenyx[j] + ycomp * gdenyy[j] + zcomp * gdenyz[j]);

                gatmz += (xcomp * gdenzx[j] + ycomp * gdenzy[j] + zcomp * gdenzz[j]);
            }
        }

        if (xcFuncType == xcfun::mgga)
        {
            // not implemented

            std::string errmgga("XCMolecularGradient._compFxcContrib: Not implemented for meta-GGA");

            errors::assertMsgCritical(false, errmgga);
        }

        // factor of 2 from sum of alpha and beta contributions

        mgradx[i] = 2.0 * gatmx;

        mgrady[i] = 2.0 * gatmy;

        mgradz[i] = 2.0 * gatmz;
    }
}

void
CXCMolecularGradient::_compGxcContrib(CDenseMatrix&              molecularGradient,
                                      const CMolecule&           molecule,
                                      const CMolecularBasis&     basis,
                                      const xcfun                xcFuncType,
                                      const CAODensityMatrix&    densityMatrix,
                                      const CMolecularGrid&      molecularGrid,
                                      const CDensityGrid&        gsDensityGrid,
                                      const CDensityGridQuad&    rwDensityGridQuad,
                                      const CXCGradientGrid&     xcGradientGrid,
                                      const CXCHessianGrid&      xcHessianGrid,
                                      const CXCCubicHessianGrid& xcCubicHessianGrid) const
{
    auto gridWeights = molecularGrid.getWeights();

    const auto nGridPoints = molecularGrid.getNumberOfGridPoints();

    const auto nAtoms = molecule.getNumberOfAtoms();

    // set up pointers to Cartesian components of molecular gradient

    auto mgradx = molecularGradient.row(0);

    auto mgrady = molecularGradient.row(1);

    auto mgradz = molecularGradient.row(2);

    // set up pointers to exchange-correlation functional derivatives

    auto df0010 = xcGradientGrid.xcGradientValues(xcvars::grada);

    auto df1010 = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::grada);

    auto df1001 = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::gradb);

    auto df10001 = xcHessianGrid.xcHessianValues(xcvars::rhoa, xcvars::gradab);

    auto df0020 = xcHessianGrid.xcHessianValues(xcvars::grada, xcvars::grada);

    auto df0011 = xcHessianGrid.xcHessianValues(xcvars::grada, xcvars::gradb);

    auto df00101 = xcHessianGrid.xcHessianValues(xcvars::grada, xcvars::gradab);

    auto df00002 = xcHessianGrid.xcHessianValues(xcvars::gradab, xcvars::gradab);

    auto df00011 = xcHessianGrid.xcHessianValues(xcvars::gradb, xcvars::gradab);

    auto df01001 = xcHessianGrid.xcHessianValues(xcvars::rhob, xcvars::gradab);

    auto df0110 = xcHessianGrid.xcHessianValues(xcvars::rhob, xcvars::grada);

    auto df3000 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::rhoa);

    auto df2100 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::rhob);

    auto df1200 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhob, xcvars::rhob);

    auto df2010 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::grada);

    auto df0030 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::grada, xcvars::grada, xcvars::grada);

    auto df0021 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::grada, xcvars::grada, xcvars::gradb);

    auto df0012 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::grada, xcvars::gradb, xcvars::gradb);

    auto df00201 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::grada, xcvars::grada, xcvars::gradab);

    auto df00111 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::grada, xcvars::gradb, xcvars::gradab);

    auto df00102 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::grada, xcvars::gradab, xcvars::gradab);

    auto df00003 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::gradab, xcvars::gradab, xcvars::gradab);

    auto df2001 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::gradb);

    auto df1110 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhob, xcvars::grada);

    auto df1101 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhob, xcvars::gradb);

    auto df20001 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhoa, xcvars::gradab);

    auto df11001 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::rhob, xcvars::gradab);

    auto df1020 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::grada, xcvars::grada);

    auto df1011 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::grada, xcvars::gradb);

    auto df1002 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::gradb, xcvars::gradb);

    auto df10101 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::grada, xcvars::gradab);

    auto df10002 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::gradab, xcvars::gradab);

    auto df01002 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhob, xcvars::gradab, xcvars::gradab);

    auto df0120 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhob, xcvars::grada, xcvars::grada);

    auto df0111 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhob, xcvars::grada, xcvars::gradb);

    auto df01101 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhob, xcvars::grada, xcvars::gradab);

    auto df10011 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhoa, xcvars::gradb, xcvars::gradab);

    auto df01011 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhob, xcvars::gradb, xcvars::gradab);

    auto df0210 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhob, xcvars::rhob, xcvars::grada);

    auto df02001 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::rhob, xcvars::rhob, xcvars::gradab);

    auto df00021 = xcCubicHessianGrid.xcCubicHessianValues(xcvars::gradb, xcvars::gradb, xcvars::gradab);

    // set up pointers to ground state density gradient norms

    auto ngrada = gsDensityGrid.alphaDensityGradient(0);

    auto grada_x = gsDensityGrid.alphaDensityGradientX(0);

    auto grada_y = gsDensityGrid.alphaDensityGradientY(0);

    auto grada_z = gsDensityGrid.alphaDensityGradientZ(0);

    // set up pointers to perturbed density gradient norms

    auto rhow1a = rwDensityGridQuad.rhow1rhow2(0);

    auto rhow1rhow2 = rwDensityGridQuad.rhow1rhow2(0);

    auto rxw1rhow2 = rwDensityGridQuad.rxw1rhow2(0);

    auto ryw1rhow2 = rwDensityGridQuad.ryw1rhow2(0);

    auto rzw1rhow2 = rwDensityGridQuad.rzw1rhow2(0);

    auto rxw1rxw2 = rwDensityGridQuad.rxw1rxw2(0);

    auto rxw1ryw2 = rwDensityGridQuad.rxw1ryw2(0);

    auto rxw1rzw2 = rwDensityGridQuad.rxw1rzw2(0);

    auto ryw1rxw2 = rwDensityGridQuad.ryw1rxw2(0);

    auto ryw1ryw2 = rwDensityGridQuad.ryw1ryw2(0);

    auto ryw1rzw2 = rwDensityGridQuad.ryw1rzw2(0);

    auto rzw1rxw2 = rwDensityGridQuad.rzw1rxw2(0);

    auto rzw1ryw2 = rwDensityGridQuad.rzw1ryw2(0);

    auto rzw1rzw2 = rwDensityGridQuad.rzw1rzw2(0);

    // set up density gradient grid driver

    CDensityGradientGridDriver graddrv(_locComm);

    for (int32_t i = 0; i < nAtoms; i++)
    {
        auto gradgrid = graddrv.generate(densityMatrix, molecule, basis, molecularGrid, xcFuncType, i);

        double gatmx = 0.0;

        double gatmy = 0.0;

        double gatmz = 0.0;

        // compute LDA contribution to molecular gradient

        if (xcFuncType == xcfun::lda)
        {
            // set up pointers to density gradient grid

            const auto gdenx = gradgrid.getComponent(0);

            const auto gdeny = gradgrid.getComponent(1);

            const auto gdenz = gradgrid.getComponent(2);

            for (int32_t j = 0; j < nGridPoints; j++)
            {
                double prefac = gridWeights[j] * (df3000[j] + df2100[j] + df2100[j] + df1200[j]) * rhow1a[j];

                gatmx += prefac * gdenx[j];

                gatmy += prefac * gdeny[j];

                gatmz += prefac * gdenz[j];
            }
        }

        // compute GGA contribution to molecular gradient

        if (xcFuncType == xcfun::gga)
        {
            // set up pointers to density gradient grid

            const auto gdenx = gradgrid.getComponent(0);

            const auto gdeny = gradgrid.getComponent(1);

            const auto gdenz = gradgrid.getComponent(2);

            const auto gdenxx = gradgrid.getComponent(3);

            const auto gdenxy = gradgrid.getComponent(4);

            const auto gdenxz = gradgrid.getComponent(5);

            const auto gdenyx = gradgrid.getComponent(6);

            const auto gdenyy = gradgrid.getComponent(7);

            const auto gdenyz = gradgrid.getComponent(8);

            const auto gdenzx = gradgrid.getComponent(9);

            const auto gdenzy = gradgrid.getComponent(10);

            const auto gdenzz = gradgrid.getComponent(11);

            for (int32_t j = 0; j < nGridPoints; j++)
            {
                // auto omega = bgaos[m] * kgaos[m];
                // auto xomega = bgaox[m] * kgaos[m] + bgaos[m] * kgaox[m];
                // auto yomega = bgaoy[m] * kgaos[m] + bgaos[m] * kgaoy[m];
                // auto zomega = bgaoz[m] * kgaos[m] + bgaos[m] * kgaoz[m];

                double w = gridWeights[j];

                double znva = 1.0 / ngrada[j];

                double znva3 = 1.0 / std::pow(ngrada[j], 3.0);

                double znva5 = 1.0 / std::pow(ngrada[j], 5.0);

                double xigrad_x = znva * grada_x[j];

                double xigrad_y = znva * grada_y[j];

                double xigrad_z = znva * grada_z[j];

                double xigrad_xx = (znva - grada_x[j] * grada_x[j] * znva3);

                double xigrad_yy = (znva - grada_y[j] * grada_y[j] * znva3);

                double xigrad_zz = (znva - grada_z[j] * grada_z[j] * znva3);

                double xigrad_xy = -grada_x[j] * grada_y[j] * znva3;

                double xigrad_xz = -grada_x[j] * grada_z[j] * znva3;

                double xigrad_yz = -grada_y[j] * grada_z[j] * znva3;

                double xigrad_xxy = 3.0 * grada_x[j] * grada_x[j] * grada_y[j] * znva5 - grada_y[j] * znva3;

                double xigrad_xxz = 3.0 * grada_x[j] * grada_x[j] * grada_z[j] * znva5 - grada_z[j] * znva3;

                double xigrad_xyy = 3.0 * grada_x[j] * grada_y[j] * grada_y[j] * znva5 - grada_x[j] * znva3;

                double xigrad_xzz = 3.0 * grada_x[j] * grada_z[j] * grada_z[j] * znva5 - grada_x[j] * znva3;

                double xigrad_yzz = 3.0 * grada_y[j] * grada_z[j] * grada_z[j] * znva5 - grada_y[j] * znva3;

                double xigrad_yyz = 3.0 * grada_y[j] * grada_y[j] * grada_z[j] * znva5 - grada_z[j] * znva3;

                double xigrad_xyz = 3.0 * grada_x[j] * grada_y[j] * grada_z[j] * znva5;

                double xigrad_xxx = 3.0 * grada_x[j] * grada_x[j] * grada_x[j] * znva5 - 3.0 * grada_x[j] * znva3;

                double xigrad_yyy = 3.0 * grada_y[j] * grada_y[j] * grada_y[j] * znva5 - 3.0 * grada_y[j] * znva3;

                double xigrad_zzz = 3.0 * grada_z[j] * grada_z[j] * grada_z[j] * znva5 - 3.0 * grada_z[j] * znva3;

                // Various required quantities

                // xigrad_dot_omega = (xigrad_x * xomega + xigrad_y * yomega + xigrad_z * zomega);

                double xigrad_dot_rw1rw2 = xigrad_x * rxw1rhow2[j] + xigrad_y * ryw1rhow2[j] + xigrad_z * rzw1rhow2[j];

                double rw1_dot_rw2 = rxw1rxw2[j] + ryw1ryw2[j] + rzw1rzw2[j];

                double xigrad_dot_rw1rhow2 = xigrad_x * rxw1rhow2[j] + xigrad_y * ryw1rhow2[j] + xigrad_z * rzw1rhow2[j];

                // grad_dot_omega = grada_x[j] * xomega + grada_y[j] * yomega + grada_z[j] * zomega;

                double grad_dot_rw1rw2 = grada_x[j] * rxw1rhow2[j] + grada_y[j] * ryw1rhow2[j] + grada_z[j] * rzw1rhow2[j];

                // omega_dot_rw1rhow2 = xomega * rxw1rhow2[j] + yomega * ryw1rhow2[j] + zomega * rzw1rhow2[j];

                double grad_dot_rw1rhow2 = grada_x[j] * rxw1rhow2[j] + grada_y[j] * ryw1rhow2[j] + grada_z[j] * rzw1rhow2[j];

                double xigrad_dot_rw1_xigrad_dot_rw2 =
                    xigrad_x * xigrad_x * rxw1rxw2[j] + xigrad_x * xigrad_y * rxw1ryw2[j] + xigrad_x * xigrad_z * rxw1rzw2[j] +
                    xigrad_y * xigrad_x * ryw1rxw2[j] + xigrad_y * xigrad_y * ryw1ryw2[j] + xigrad_y * xigrad_z * ryw1rzw2[j] +
                    xigrad_z * xigrad_x * rzw1rxw2[j] + xigrad_z * xigrad_y * rzw1ryw2[j] + xigrad_z * xigrad_z * rzw1rzw2[j];

                // twelthfifth_gam = (xigrad_x * grada_x[j] + grada_x[j] * xigrad_x) * rxw1rxw2[j] +
                //                   (xigrad_x * grada_y[j] + grada_x[j] * xigrad_y) * rxw1ryw2[j] +
                //                   (xigrad_x * grada_z[j] + grada_x[j] * xigrad_z) * rxw1rzw2[j] +
                //                   (xigrad_y * grada_x[j] + grada_y[j] * xigrad_x) * rxw1rxw2[j] +
                //                   (xigrad_y * grada_y[j] + grada_y[j] * xigrad_y) * rxw1ryw2[j] +
                //                   (xigrad_y * grada_z[j] + grada_y[j] * xigrad_z) * rxw1rzw2[j] +
                //                   (xigrad_z * grada_x[j] + grada_z[j] * xigrad_x) * rxw1rxw2[j] +
                //                   (xigrad_z * grada_y[j] + grada_z[j] * xigrad_y) * rxw1ryw2[j] +
                //                   (xigrad_z * grada_z[j] + grada_z[j] * xigrad_z) * rxw1rzw2[j];

                double twelthfifth_gam =
                    (xigrad_x * grada_x[j] + grada_x[j] * xigrad_x) * rxw1rxw2[j] + (xigrad_x * grada_y[j] + grada_x[j] * xigrad_y) * rxw1ryw2[j] +
                    (xigrad_x * grada_z[j] + grada_x[j] * xigrad_z) * rxw1rzw2[j] + (xigrad_y * grada_x[j] + grada_y[j] * xigrad_x) * rxw1rxw2[j] +
                    (xigrad_y * grada_y[j] + grada_y[j] * xigrad_y) * rxw1ryw2[j] + (xigrad_y * grada_z[j] + grada_y[j] * xigrad_z) * rxw1rzw2[j] +
                    (xigrad_z * grada_x[j] + grada_z[j] * xigrad_x) * rxw1rxw2[j] + (xigrad_z * grada_y[j] + grada_z[j] * xigrad_y) * rxw1ryw2[j] +
                    (xigrad_z * grada_z[j] + grada_z[j] * xigrad_z) * rxw1rzw2[j];

                // contribution from \nabla_A (\phi_mu \phi_nu)

                double prefac;

                // fifth = w * (df3000[j] + 2.0 * df2100[j] + df1200[j]) * rhow1rhow2[j] * omega;
                // seventh = w * (df2010[j] + df2001[j]) * xigrad_dot_rw1rw2 * omega;
                // seventh += w * (df1110[j] + df1101[j]) * xigrad_dot_rw1rw2 * omega;
                // seventh += w * 2.0 * (df20001[j] + df11001[j]) * grad_dot_rw1rw2 * omega;

                // eighth = w * (df1020[j] + 2.0 * df1011[j] + df1002[j]) * xigrad_dot_rw1_xigrad_dot_rw2 * omega;
                // eighth += w * (df1010[j] + df1001[j]) *
                //           (xigrad_xx * rxw1rxw2[j] + xigrad_xy * rxw1ryw2[j] + xigrad_xz * rxw1rzw2[j]
                //          + xigrad_xy * ryw1rxw2[j] + xigrad_yy * ryw1ryw2[j] + xigrad_yz * ryw1rzw2[j]
                //          + xigrad_xz * rzw1rxw2[j] + xigrad_yz * rzw1ryw2[j] + xigrad_zz * rzw1rzw2[j]) * omega;
                // eighth += w * 2.0 * (df10101[j] + df10101[j]) * ngrada[j] * xigrad_dot_rw1_xigrad_dot_rw2 * omega;
                // eighth += w * 4.0 * df10002[j] * ngrada[j] * ngrada[j] * xigrad_dot_rw1_xigrad_dot_rw2 * omega;
                // eighth += w * 2.0 * df10001[j] * rw1_dot_rw2 * omega;

                prefac = w * (df3000[j] + 2.0 * df2100[j] + df1200[j]) * rhow1rhow2[j]

                         + w * (df2010[j] + df2001[j]) * xigrad_dot_rw1rw2

                         + w * (df1110[j] + df1101[j]) * xigrad_dot_rw1rw2

                         + w * 2.0 * (df20001[j] + df11001[j]) * grad_dot_rw1rw2

                         + w * (df1020[j] + 2.0 * df1011[j] + df1002[j]) * xigrad_dot_rw1_xigrad_dot_rw2

                         + w * (df1010[j] + df1001[j]) *

                               (xigrad_xx * rxw1rxw2[j] + xigrad_xy * rxw1ryw2[j] + xigrad_xz * rxw1rzw2[j]

                                + xigrad_xy * ryw1rxw2[j] + xigrad_yy * ryw1ryw2[j] + xigrad_yz * ryw1rzw2[j]

                                + xigrad_xz * rzw1rxw2[j] + xigrad_yz * rzw1ryw2[j] + xigrad_zz * rzw1rzw2[j])

                         + w * 2.0 * (df10101[j] + df10101[j]) * ngrada[j] * xigrad_dot_rw1_xigrad_dot_rw2

                         + w * 4.0 * df10002[j] * ngrada[j] * ngrada[j] * xigrad_dot_rw1_xigrad_dot_rw2

                         + w * 2.0 * df10001[j] * rw1_dot_rw2;

                gatmx += prefac * gdenx[j];

                gatmy += prefac * gdeny[j];

                gatmz += prefac * gdenz[j];

                // contribution from \nabla_A (\nabla (\phi_mu \phi_nu))

                double xcomp = 0.0;

                double ycomp = 0.0;

                double zcomp = 0.0;

                // ninth = w * (df2010[j] + 2.0 * df1110[j] + df0210[j]) * rhow1rhow2[j] * xigrad_dot_omega;

                // tenth += w * (df1020[j] + df1011[j] + df0120[j] + df0111[j]) * xigrad_dot_rw1rhow2 * xigrad_dot_omega;
                // tenth += w * (df10101[j] + df10011[j] + df01101[j] + df0111[j]) * grad_dot_rw1rhow2 * xigrad_dot_omega;

                // twelfth += w * (df0030[j] + 2.0 * df0021[j] + df0012[j]) * xigrad_dot_rw1_xigrad_dot_rw2 * xigrad_dot_omega;
                // twelfth += w * (df00101[j] + df00011[j]) * xigrad_dot_omega * rw1_dot_rw2;
                // twelfth += w * (df00201[j] + df00111[j]) * twelthfifth_gam * xigrad_dot_omega;
                // twelfth += w * df00102[j] * xigrad_dot_rw1_xigrad_dot_rw2 * xigrad_dot_omega;

                // xigrad_dot_omega = (xigrad_x * xomega + xigrad_y * yomega + xigrad_z * zomega);

                prefac = w * (df2010[j] + 2.0 * df1110[j] + df0210[j]) * rhow1rhow2[j]

                         + w * (df1020[j] + df1011[j] + df0120[j] + df0111[j]) * xigrad_dot_rw1rhow2

                         + w * (df10101[j] + df10011[j] + df01101[j] + df0111[j]) * grad_dot_rw1rhow2

                         + w * (df0030[j] + 2.0 * df0021[j] + df0012[j]) * xigrad_dot_rw1_xigrad_dot_rw2

                         + w * (df00101[j] + df00011[j]) * rw1_dot_rw2

                         + w * (df00201[j] + df00111[j]) * twelthfifth_gam

                         + w * df00102[j] * xigrad_dot_rw1_xigrad_dot_rw2;

                xcomp += prefac * xigrad_x;

                ycomp += prefac * xigrad_y;

                zcomp += prefac * xigrad_z;

                // ninth += w * (df20001[j] + 2.0 * df11001[j] + df02001[j]) * grad_dot_omega * rhow1rhow2[j];

                // tenth += w * (df10101[j] + df10011[j] + df01101[j] + df0111[j] + df01011[j]) * xigrad_dot_rw1rhow2 * grad_dot_omega;
                // tenth += w * (df10002[j] + df01002[j]) * grad_dot_rw1rhow2 * grad_dot_omega;

                // twelfth += w * df00002[j] * grad_dot_omega * rw1_dot_rw2;
                // twelfth += w * (df00201[j] + 2 * df00111[j] + df00021[j]) * xigrad_dot_rw1_xigrad_dot_rw2 * grad_dot_omega;
                // twelfth += w * (df00102[j] + df00011[j]) * ngrada[j] * xigrad_dot_rw1_xigrad_dot_rw2 * grad_dot_omega;
                // twelfth += w * df00003[j] * ngrada[j] * ngrada[j] * grad_dot_omega;

                // grad_dot_omega = grada_x[j] * xomega + grada_y[j] * yomega + grada_z[j] * zomega;

                prefac = w * (df20001[j] + 2.0 * df11001[j] + df02001[j]) * rhow1rhow2[j]

                         + w * (df10101[j] + df10011[j] + df01101[j] + df0111[j] + df01011[j]) * xigrad_dot_rw1rhow2

                         + w * (df10002[j] + df01002[j]) * grad_dot_rw1rhow2

                         + w * df00002[j] * rw1_dot_rw2

                         + w * (df00201[j] + 2 * df00111[j] + df00021[j]) * xigrad_dot_rw1_xigrad_dot_rw2

                         + w * (df00102[j] + df00011[j]) * ngrada[j] * xigrad_dot_rw1_xigrad_dot_rw2

                         + w * df00003[j] * ngrada[j] * ngrada[j];

                xcomp += prefac * grada_x[j];

                ycomp += prefac * grada_y[j];

                zcomp += prefac * grada_z[j];

                // tenth += w * (df10001[j] + df01001[j]) * omega_dot_rw1rhow2;

                // omega_dot_rw1rhow2 = xomega * rxw1rhow2[j] + yomega * ryw1rhow2[j] + zomega * rzw1rhow2[j];

                prefac = w * (df10001[j] + df01001[j]);

                xcomp += prefac * rxw1rhow2[j];

                ycomp += prefac * ryw1rhow2[j];

                zcomp += prefac * rzw1rhow2[j];

                // tenth = w * (df1010[j] + df0110[j]) *
                //          ((xigrad_xx * rxw1rhow2[j] + xigrad_xy * ryw1rhow2[j] + xigrad_xz * rzw1rhow2[j]) * xomega
                //         + (xigrad_xy * rxw1rhow2[j] + xigrad_yy * ryw1rhow2[j] + xigrad_yz * rzw1rhow2[j]) * yomega
                //         + (xigrad_xz * rxw1rhow2[j] + xigrad_yz * ryw1rhow2[j] + xigrad_zz * rzw1rhow2[j]) * zomega);

                prefac = w * (df1010[j] + df0110[j]);

                xcomp += prefac * (xigrad_xx * rxw1rhow2[j] + xigrad_xy * ryw1rhow2[j] + xigrad_xz * rzw1rhow2[j]);

                ycomp += prefac * (xigrad_xy * rxw1rhow2[j] + xigrad_yy * ryw1rhow2[j] + xigrad_yz * rzw1rhow2[j]);

                zcomp += prefac * (xigrad_xz * rxw1rhow2[j] + xigrad_yz * ryw1rhow2[j] + xigrad_zz * rzw1rhow2[j]);

                // twelfth = w * df0010[j] * twelthfirst;
                // twelthfirst = xigrad_xxx * xomega * rxw1rxw2[j] + xigrad_xxy * xomega * rxw1ryw2[j] +
                //               xigrad_xxz * xomega * rxw1rzw2[j] + xigrad_xxy * xomega * ryw1rxw2[j] +
                //               xigrad_xyy * xomega * ryw1ryw2[j] + xigrad_xyz * xomega * ryw1rzw2[j] +
                //               xigrad_xxz * xomega * rzw1rxw2[j] + xigrad_xyz * xomega * rzw1ryw2[j] +
                //               xigrad_xzz * xomega * rzw1rzw2[j] + xigrad_xxy * yomega * rxw1rxw2[j] +
                //               xigrad_xyy * yomega * rxw1ryw2[j] + xigrad_xyz * yomega * rxw1rzw2[j] +
                //               xigrad_xyy * yomega * ryw1rxw2[j] + xigrad_yyy * yomega * ryw1ryw2[j] +
                //               xigrad_yyz * yomega * ryw1rzw2[j] + xigrad_xyz * yomega * rzw1rxw2[j] +
                //               xigrad_yyz * yomega * rzw1ryw2[j] + xigrad_yzz * yomega * rzw1rzw2[j] +
                //               xigrad_xxz * zomega * rxw1rxw2[j] + xigrad_xyz * zomega * rxw1ryw2[j] +
                //               xigrad_xzz * zomega * rxw1rzw2[j] + xigrad_xyz * zomega * ryw1rxw2[j] +
                //               xigrad_yyz * zomega * ryw1ryw2[j] + xigrad_yzz * zomega * ryw1rzw2[j] +
                //               xigrad_xzz * zomega * rzw1rxw2[j] + xigrad_yzz * zomega * rzw1ryw2[j] +
                //               xigrad_zzz * zomega * rzw1rzw2[j];

                prefac = w * df0010[j];

                xcomp += prefac * (xigrad_xxx * rxw1rxw2[j] + xigrad_xxy * rxw1ryw2[j] + xigrad_xxz * rxw1rzw2[j]

                                   + xigrad_xxy * ryw1rxw2[j] + xigrad_xyy * ryw1ryw2[j] + xigrad_xyz * ryw1rzw2[j]

                                   + xigrad_xxz * rzw1rxw2[j] + xigrad_xyz * rzw1ryw2[j] + xigrad_xzz * rzw1rzw2[j]);

                ycomp += prefac * (xigrad_xxy * rxw1rxw2[j] + xigrad_xyy * rxw1ryw2[j] + xigrad_xyz * rxw1rzw2[j]

                                   + xigrad_xyy * ryw1rxw2[j] + xigrad_yyy * ryw1ryw2[j] + xigrad_yyz * ryw1rzw2[j]

                                   + xigrad_xyz * rzw1rxw2[j] + xigrad_yyz * rzw1ryw2[j] + xigrad_yzz * rzw1rzw2[j]);

                zcomp += prefac * (xigrad_xxz * rxw1rxw2[j] + xigrad_xyz * rxw1ryw2[j] + xigrad_xzz * rxw1rzw2[j]

                                   + xigrad_xyz * ryw1rxw2[j] + xigrad_yyz * ryw1ryw2[j] + xigrad_yzz * ryw1rzw2[j]

                                   + xigrad_xzz * rzw1rxw2[j] + xigrad_yzz * rzw1ryw2[j] + xigrad_zzz * rzw1rzw2[j]);

                // twelfth += w * (df0020[j] + df0011[j]) * twelthsecond;
                // twelfth += w * (df00101[j] + df00011[j]) * ngrada[j] * twelthsecond;

                // twelthsecond = xigrad_xx * xigrad_x * xomega * rxw1rxw2[j]
                //              + xigrad_xx * xigrad_y * yomega * rxw1rxw2[j]
                //              + xigrad_xx * xigrad_z * zomega * rxw1rxw2[j]
                //              + xigrad_xy * xigrad_x * xomega * rxw1ryw2[j]
                //              + xigrad_xy * xigrad_y * yomega * rxw1ryw2[j]
                //              + xigrad_xy * xigrad_z * zomega * rxw1ryw2[j]
                //              + xigrad_xz * xigrad_x * xomega * rxw1rzw2[j]
                //              + xigrad_xz * xigrad_y * yomega * rxw1rzw2[j]
                //              + xigrad_xz * xigrad_z * zomega * rxw1rzw2[j]
                //              + xigrad_xy * xigrad_x * xomega * ryw1rxw2[j]
                //              + xigrad_xy * xigrad_y * yomega * ryw1rxw2[j]
                //              + xigrad_xy * xigrad_z * zomega * ryw1rxw2[j]
                //              + xigrad_yy * xigrad_x * xomega * ryw1ryw2[j]
                //              + xigrad_yy * xigrad_y * yomega * ryw1ryw2[j]
                //              + xigrad_yy * xigrad_z * zomega * ryw1ryw2[j]
                //              + xigrad_yz * xigrad_x * xomega * ryw1rzw2[j]
                //              + xigrad_yz * xigrad_y * yomega * ryw1rzw2[j]
                //              + xigrad_yz * xigrad_z * zomega * ryw1rzw2[j]
                //              + xigrad_xz * xigrad_x * xomega * rzw1rxw2[j]
                //              + xigrad_xz * xigrad_y * yomega * rzw1rxw2[j]
                //              + xigrad_xz * xigrad_z * zomega * rzw1rxw2[j]
                //              + xigrad_yz * xigrad_x * xomega * rzw1ryw2[j]
                //              + xigrad_yz * xigrad_y * yomega * rzw1ryw2[j]
                //              + xigrad_yz * xigrad_z * zomega * rzw1ryw2[j]
                //              + xigrad_zz * xigrad_x * xomega * rzw1rzw2[j]
                //              + xigrad_zz * xigrad_y * yomega * rzw1rzw2[j]
                //              + xigrad_zz * xigrad_z * zomega * rzw1rzw2[j];

                prefac = w * (df0020[j] + df0011[j])

                         + w * (df00101[j] + df00011[j]) * ngrada[j];

                xcomp += prefac * (xigrad_xx * xigrad_x * rxw1rxw2[j]

                                   + xigrad_xy * xigrad_x * rxw1ryw2[j]

                                   + xigrad_xz * xigrad_x * rxw1rzw2[j]

                                   + xigrad_xy * xigrad_x * ryw1rxw2[j]

                                   + xigrad_yy * xigrad_x * ryw1ryw2[j]

                                   + xigrad_yz * xigrad_x * ryw1rzw2[j]

                                   + xigrad_xz * xigrad_x * rzw1rxw2[j]

                                   + xigrad_yz * xigrad_x * rzw1ryw2[j]

                                   + xigrad_zz * xigrad_x * rzw1rzw2[j]);

                ycomp += prefac * (xigrad_xx * xigrad_y * rxw1rxw2[j]

                                   + xigrad_xy * xigrad_y * rxw1ryw2[j]

                                   + xigrad_xz * xigrad_y * rxw1rzw2[j]

                                   + xigrad_xy * xigrad_y * ryw1rxw2[j]

                                   + xigrad_yy * xigrad_y * ryw1ryw2[j]

                                   + xigrad_yz * xigrad_y * ryw1rzw2[j]

                                   + xigrad_xz * xigrad_y * rzw1rxw2[j]

                                   + xigrad_yz * xigrad_y * rzw1ryw2[j]

                                   + xigrad_zz * xigrad_y * rzw1rzw2[j]);

                zcomp += prefac * (xigrad_xx * xigrad_z * rxw1rxw2[j]

                                   + xigrad_xy * xigrad_z * rxw1ryw2[j]

                                   + xigrad_xz * xigrad_z * rxw1rzw2[j]

                                   + xigrad_xy * xigrad_z * ryw1rxw2[j]

                                   + xigrad_yy * xigrad_z * ryw1ryw2[j]

                                   + xigrad_yz * xigrad_z * ryw1rzw2[j]

                                   + xigrad_xz * xigrad_z * rzw1rxw2[j]

                                   + xigrad_yz * xigrad_z * rzw1ryw2[j]

                                   + xigrad_zz * xigrad_z * rzw1rzw2[j]);

                // twelfth += w * (df0020[j] + df0011[j]) * twelththird;
                // twelfth += w * df00101[j] * ngrada[j] * twelththird;

                // twelththird = xigrad_xx * xigrad_x * xomega * (rxw1rxw2[j] + rxw1rxw2[j]) +
                //               xigrad_xx * xigrad_y * xomega * (ryw1rxw2[j] + rxw1ryw2[j]) +
                //               xigrad_xx * xigrad_z * xomega * (rzw1rxw2[j] + rxw1rzw2[j]) +
                //               xigrad_xy * xigrad_x * xomega * (rxw1ryw2[j] + ryw1rxw2[j]) +
                //               xigrad_xy * xigrad_y * xomega * (ryw1ryw2[j] + ryw1ryw2[j]) +
                //               xigrad_xy * xigrad_z * xomega * (rzw1ryw2[j] + ryw1rzw2[j]) +
                //               xigrad_xz * xigrad_x * xomega * (rxw1rzw2[j] + rzw1rxw2[j]) +
                //               xigrad_xz * xigrad_y * xomega * (ryw1rzw2[j] + rzw1ryw2[j]) +
                //               xigrad_xz * xigrad_z * xomega * (rzw1rzw2[j] + rzw1rzw2[j]) +
                //               xigrad_xy * xigrad_x * yomega * (rxw1rxw2[j] + rxw1rxw2[j]) +
                //               xigrad_xy * xigrad_y * yomega * (ryw1rxw2[j] + rxw1ryw2[j]) +
                //               xigrad_xy * xigrad_z * yomega * (rzw1rxw2[j] + rxw1rzw2[j]) +
                //               xigrad_yy * xigrad_x * yomega * (rxw1ryw2[j] + ryw1rxw2[j]) +
                //               xigrad_yy * xigrad_y * yomega * (ryw1ryw2[j] + ryw1ryw2[j]) +
                //               xigrad_yy * xigrad_z * yomega * (rzw1ryw2[j] + ryw1rzw2[j]) +
                //               xigrad_yz * xigrad_x * yomega * (rxw1rzw2[j] + rzw1rxw2[j]) +
                //               xigrad_yz * xigrad_y * yomega * (ryw1rzw2[j] + rzw1ryw2[j]) +
                //               xigrad_yz * xigrad_z * yomega * (rzw1rzw2[j] + rzw1rzw2[j]) +
                //               xigrad_xz * xigrad_x * zomega * (rxw1rxw2[j] + rxw1rxw2[j]) +
                //               xigrad_xz * xigrad_y * zomega * (ryw1rxw2[j] + rxw1ryw2[j]) +
                //               xigrad_xz * xigrad_z * zomega * (rzw1rxw2[j] + rxw1rzw2[j]) +
                //               xigrad_yz * xigrad_x * zomega * (rxw1ryw2[j] + ryw1rxw2[j]) +
                //               xigrad_yz * xigrad_y * zomega * (ryw1ryw2[j] + ryw1ryw2[j]) +
                //               xigrad_yz * xigrad_z * zomega * (rzw1ryw2[j] + ryw1rzw2[j]) +
                //               xigrad_zz * xigrad_x * zomega * (rxw1rzw2[j] + rzw1rxw2[j]) +
                //               xigrad_zz * xigrad_y * zomega * (ryw1rzw2[j] + rzw1ryw2[j]) +
                //               xigrad_zz * xigrad_z * zomega * (rzw1rzw2[j] + rzw1rzw2[j]);

                prefac = w * (df0020[j] + df0011[j])

                         + w * df00101[j] * ngrada[j];

                xcomp += prefac * (xigrad_xx * xigrad_x * (rxw1rxw2[j] + rxw1rxw2[j])

                                   + xigrad_xx * xigrad_y * (ryw1rxw2[j] + rxw1ryw2[j])

                                   + xigrad_xx * xigrad_z * (rzw1rxw2[j] + rxw1rzw2[j])

                                   + xigrad_xy * xigrad_x * (rxw1ryw2[j] + ryw1rxw2[j])

                                   + xigrad_xy * xigrad_y * (ryw1ryw2[j] + ryw1ryw2[j])

                                   + xigrad_xy * xigrad_z * (rzw1ryw2[j] + ryw1rzw2[j])

                                   + xigrad_xz * xigrad_x * (rxw1rzw2[j] + rzw1rxw2[j])

                                   + xigrad_xz * xigrad_y * (ryw1rzw2[j] + rzw1ryw2[j])

                                   + xigrad_xz * xigrad_z * (rzw1rzw2[j] + rzw1rzw2[j]));

                ycomp += prefac * (xigrad_xy * xigrad_x * (rxw1rxw2[j] + rxw1rxw2[j])

                                   + xigrad_xy * xigrad_y * (ryw1rxw2[j] + rxw1ryw2[j])

                                   + xigrad_xy * xigrad_z * (rzw1rxw2[j] + rxw1rzw2[j])

                                   + xigrad_yy * xigrad_x * (rxw1ryw2[j] + ryw1rxw2[j])

                                   + xigrad_yy * xigrad_y * (ryw1ryw2[j] + ryw1ryw2[j])

                                   + xigrad_yy * xigrad_z * (rzw1ryw2[j] + ryw1rzw2[j])

                                   + xigrad_yz * xigrad_x * (rxw1rzw2[j] + rzw1rxw2[j])

                                   + xigrad_yz * xigrad_y * (ryw1rzw2[j] + rzw1ryw2[j])

                                   + xigrad_yz * xigrad_z * (rzw1rzw2[j] + rzw1rzw2[j]));

                zcomp += prefac * (xigrad_xz * xigrad_x * (rxw1rxw2[j] + rxw1rxw2[j])

                                   + xigrad_xz * xigrad_y * (ryw1rxw2[j] + rxw1ryw2[j])

                                   + xigrad_xz * xigrad_z * (rzw1rxw2[j] + rxw1rzw2[j])

                                   + xigrad_yz * xigrad_x * (rxw1ryw2[j] + ryw1rxw2[j])

                                   + xigrad_yz * xigrad_y * (ryw1ryw2[j] + ryw1ryw2[j])

                                   + xigrad_yz * xigrad_z * (rzw1ryw2[j] + ryw1rzw2[j])

                                   + xigrad_zz * xigrad_x * (rxw1rzw2[j] + rzw1rxw2[j])

                                   + xigrad_zz * xigrad_y * (ryw1rzw2[j] + rzw1ryw2[j])

                                   + xigrad_zz * xigrad_z * (rzw1rzw2[j] + rzw1rzw2[j]));

                // twelfth += w * (df00101[j] + df00011[j]) * twelthfourth_gam;
                // twelfth += w * df00002[j] * ngrada[j] * twelthfourth_gam;
                // twelthfourth_gam =
                //   xigrad_x * xomega * rxw1rxw2[j] + xigrad_x * yomega * rxw1ryw2[j] + xigrad_x * zomega * rxw1rzw2[j] +
                //   xigrad_y * xomega * ryw1rxw2[j] + xigrad_y * yomega * ryw1ryw2[j] + xigrad_y * zomega * ryw1rzw2[j] +
                //   xigrad_z * xomega * rzw1rxw2[j] + xigrad_z * yomega * rzw1ryw2[j] + xigrad_z * zomega * rzw1rzw2[j];

                prefac = w * (df00101[j] + df00011[j])

                         + w * df00002[j] * ngrada[j];

                xcomp += prefac * (xigrad_x * rxw1rxw2[j]

                                   + xigrad_y * ryw1rxw2[j]

                                   + xigrad_z * rzw1rxw2[j]);

                ycomp += prefac * (xigrad_x * rxw1ryw2[j]

                                   + xigrad_y * ryw1ryw2[j]

                                   + xigrad_z * rzw1ryw2[j]);

                zcomp += prefac * (xigrad_x * rxw1rzw2[j]

                                   + xigrad_y * ryw1rzw2[j]

                                   + xigrad_z * rzw1rzw2[j]);

                gatmx += (xcomp * gdenxx[j] + ycomp * gdenxy[j] + zcomp * gdenxz[j]);

                gatmy += (xcomp * gdenyx[j] + ycomp * gdenyy[j] + zcomp * gdenyz[j]);

                gatmz += (xcomp * gdenzx[j] + ycomp * gdenzy[j] + zcomp * gdenzz[j]);
            }
        }

        if (xcFuncType == xcfun::mgga)
        {
            // not implemented

            std::string errmgga("XCMolecularGradient._compGxcContrib: Not implemented for meta-GGA");

            errors::assertMsgCritical(false, errmgga);
        }

        // factor of 2 from sum of alpha and beta contributions

        mgradx[i] = 2.0 * gatmx;

        mgrady[i] = 2.0 * gatmy;

        mgradz[i] = 2.0 * gatmz;
    }
}
