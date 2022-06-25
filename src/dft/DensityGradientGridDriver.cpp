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

#include "DensityGradientGridDriver.hpp"

#include "GtoFunc.hpp"
#include "MpiFunc.hpp"
#include "OMPTasks.hpp"

CDensityGradientGridDriver::CDensityGradientGridDriver(MPI_Comm comm)
{
    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    _locComm = comm;
}

CDensityGradientGridDriver::~CDensityGradientGridDriver()
{
}

CDensityGrid
CDensityGradientGridDriver::generate(const CAODensityMatrix& aoDensityMatrix,
                                     const CMolecule&        molecule,
                                     const CMolecularBasis&  basis,
                                     const CMolecularGrid&   molecularGrid,
                                     const xcfun             xcFunctional,
                                     const int32_t           iAtom)
{
    if (xcFunctional == xcfun::lda)
    {
        // initialize density grid

        CDensityGrid dgrid(molecularGrid.getNumberOfGridPoints(), aoDensityMatrix.getNumberOfDensityMatrices(), 3, dengrid::ab);

        dgrid.zero();

        // LDA

        _genRestrictedDensityForLda(dgrid, aoDensityMatrix, molecule, basis, molecularGrid, iAtom);

        return dgrid;
    }

    if (xcFunctional == xcfun::gga)
    {
        // initialize density grid

        CDensityGrid dgrid(molecularGrid.getNumberOfGridPoints(), aoDensityMatrix.getNumberOfDensityMatrices(), 12, dengrid::ab);

        dgrid.zero();

        // GGA

        _genRestrictedDensityForGga(dgrid, aoDensityMatrix, molecule, basis, molecularGrid, iAtom);

        return dgrid;
    }

    if (xcFunctional == xcfun::mgga)
    {
        // not implemented

        std::string errmgga("DensityGradientGridDriver.integrate: Not implemented for meta-GGA");

        errors::assertMsgCritical(false, errmgga);
    }

    return CDensityGrid();
}

void
CDensityGradientGridDriver::_genRestrictedDensityForLda(CDensityGrid&           densityGrid,
                                                        const CAODensityMatrix& aoDensityMatrix,
                                                        const CMolecule&        molecule,
                                                        const CMolecularBasis&  basis,
                                                        const CMolecularGrid&   molecularGrid,
                                                        const int32_t           iAtom) const
{
    // set up OMP tasks

    COMPTasks omptaks(5);

    omptaks.set(molecularGrid.getNumberOfGridPoints());

    auto ntasks = omptaks.getNumberOfTasks();

    auto tbsizes = omptaks.getTaskSizes();

    auto tbpositions = omptaks.getTaskPositions();

    // create GTOs container

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    CGtoContainer* xgtovec = new CGtoContainer(molecule, basis, iAtom, 1);

    // set up molecular grid data

    auto mgx = molecularGrid.getCoordinatesX();

    auto mgy = molecularGrid.getCoordinatesY();

    auto mgz = molecularGrid.getCoordinatesZ();

    // set up pointer to density matrix

    auto denptr = &aoDensityMatrix;

    // set up poinet to density grid

    auto dgridptr = &densityGrid;

    // generate density on grid points

    #pragma omp parallel shared(tbsizes, tbpositions, ntasks, mgx, mgy, mgz, gtovec, xgtovec, denptr, dgridptr)
    {
        #pragma omp single nowait
        {
            for (int32_t i = 0; i < ntasks; i++)
            {
                // set up task parameters

                auto tbsize = tbsizes[i];

                auto tbposition = tbpositions[i];

                // generate task

                #pragma omp task firstprivate(tbsize, tbposition)
                {
                    _genBatchOfRestrictedDensityGridPointsForLda(dgridptr, denptr, gtovec, xgtovec, mgx, mgy, mgz, tbposition, tbsize);
                }
            }
        }
    }

    // finalize density grid

    densityGrid.updateBetaDensities();

    // destroy GTOs container

    delete gtovec;

    delete xgtovec;
}

void
CDensityGradientGridDriver::_genBatchOfRestrictedDensityGridPointsForLda(CDensityGrid*           densityGrid,
                                                                         const CAODensityMatrix* aoDensityMatrix,
                                                                         const CGtoContainer*    gtoContainer,
                                                                         const CGtoContainer*    atmGtoContainer,
                                                                         const double*           gridCoordinatesX,
                                                                         const double*           gridCoordinatesY,
                                                                         const double*           gridCoordinatesZ,
                                                                         const int32_t           gridOffset,
                                                                         const int32_t           nGridPoints) const
{
    // set up number of AOs

    auto naos = gtoContainer->getNumberOfAtomicOrbitals();

    auto atmnaos = atmGtoContainer->getNumberOfAtomicOrbitals();

    // determine number of grid blocks

    auto blockdim = _getSizeOfBlock();

    auto nblocks = nGridPoints / blockdim;

    // set up current grid point

    int32_t igpnt = 0;

    // loop over grid points blocks

    if (nblocks > 0)
    {
        CMemBlock2D<double> gaos(blockdim, naos);

        CMemBlock2D<double> xgaos(blockdim, atmnaos);

        CMemBlock2D<double> xgaox(blockdim, atmnaos);

        CMemBlock2D<double> xgaoy(blockdim, atmnaos);

        CMemBlock2D<double> xgaoz(blockdim, atmnaos);

        CMemBlock<int32_t> aoidx(atmnaos);

        for (int32_t i = 0; i < nblocks; i++)
        {
            gaos.zero();

            xgaos.zero();

            xgaox.zero();

            xgaoy.zero();

            xgaoz.zero();

            gtorec::computeGtosValuesForLDA(gaos, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, igpnt, blockdim);

            gtorec::computeGtosValuesForLDA2(
                aoidx, xgaos, xgaox, xgaoy, xgaoz, atmGtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, igpnt, blockdim);

            _distRestrictedDensityValuesForLda(densityGrid, aoDensityMatrix, gaos, aoidx, xgaox, xgaoy, xgaoz, gridOffset, igpnt, blockdim);

            igpnt += blockdim;
        }
    }

    // comopute remaining grid points block

    blockdim = nGridPoints % blockdim;

    if (blockdim > 0)
    {
        CMemBlock2D<double> gaos(blockdim, naos);

        CMemBlock2D<double> xgaos(blockdim, atmnaos);

        CMemBlock2D<double> xgaox(blockdim, atmnaos);

        CMemBlock2D<double> xgaoy(blockdim, atmnaos);

        CMemBlock2D<double> xgaoz(blockdim, atmnaos);

        CMemBlock<int32_t> aoidx(atmnaos);

        gaos.zero();

        xgaos.zero();

        xgaox.zero();

        xgaoy.zero();

        xgaoz.zero();

        gtorec::computeGtosValuesForLDA(gaos, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, igpnt, blockdim);

        gtorec::computeGtosValuesForLDA2(
            aoidx, xgaos, xgaox, xgaoy, xgaoz, atmGtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, igpnt, blockdim);

        _distRestrictedDensityValuesForLda(densityGrid, aoDensityMatrix, gaos, aoidx, xgaox, xgaoy, xgaoz, gridOffset, igpnt, blockdim);
    }
}

void
CDensityGradientGridDriver::_distRestrictedDensityValuesForLda(CDensityGrid*              densityGrid,
                                                               const CAODensityMatrix*    aoDensityMatrix,
                                                               const CMemBlock2D<double>& braGtoValues,
                                                               const CMemBlock<int32_t>&  aoIdentifiers,
                                                               const CMemBlock2D<double>& ketGtoValuesX,
                                                               const CMemBlock2D<double>& ketGtoValuesY,
                                                               const CMemBlock2D<double>& ketGtoValuesZ,
                                                               const int32_t              gridOffset,
                                                               const int32_t              gridBlockPosition,
                                                               const int32_t              nGridPoints) const
{
    if (aoDensityMatrix->getNumberOfDensityMatrices() == 1)
    {
        // set up pointer to density grid data

        auto rhoax = densityGrid->getComponent(0);

        auto rhoay = densityGrid->getComponent(1);

        auto rhoaz = densityGrid->getComponent(2);

        // set up pointer to density matrix data

        auto denmat = aoDensityMatrix->alphaDensity(0);

        auto naos = aoDensityMatrix->getNumberOfRows(0);

        auto atmnaos = aoIdentifiers.size();

        // loop over density matrix

        for (int32_t i = 0; i < naos; i++)
        {
            auto bgaos = braGtoValues.data(i);

            for (int32_t j = 0; j < atmnaos; j++)
            {
                auto kgaox = ketGtoValuesX.data(j);

                auto kgaoy = ketGtoValuesY.data(j);

                auto kgaoz = ketGtoValuesZ.data(j);

                const auto dval = 2.0 * denmat[i * naos + aoIdentifiers.at(j)];

                #pragma omp simd
                for (int32_t k = 0; k < nGridPoints; k++)
                {
                    rhoax[gridOffset + gridBlockPosition + k] -= dval * bgaos[k] * kgaox[k];

                    rhoay[gridOffset + gridBlockPosition + k] -= dval * bgaos[k] * kgaoy[k];

                    rhoaz[gridOffset + gridBlockPosition + k] -= dval * bgaos[k] * kgaoz[k];
                }
            }
        }
    }
}

void
CDensityGradientGridDriver::_genRestrictedDensityForGga(CDensityGrid&           densityGrid,
                                                        const CAODensityMatrix& aoDensityMatrix,
                                                        const CMolecule&        molecule,
                                                        const CMolecularBasis&  basis,
                                                        const CMolecularGrid&   molecularGrid,
                                                        const int32_t           iAtom) const
{
    // set up OMP tasks

    COMPTasks omptaks(5);

    omptaks.set(molecularGrid.getNumberOfGridPoints());

    auto ntasks = omptaks.getNumberOfTasks();

    auto tbsizes = omptaks.getTaskSizes();

    auto tbpositions = omptaks.getTaskPositions();

    // create GTOs container

    CGtoContainer* gtovec = new CGtoContainer(molecule, basis);

    CGtoContainer* xgtovec = new CGtoContainer(molecule, basis, iAtom, 1);

    // set up molecular grid data

    auto mgx = molecularGrid.getCoordinatesX();

    auto mgy = molecularGrid.getCoordinatesY();

    auto mgz = molecularGrid.getCoordinatesZ();

    // set up pointer to density matrix

    auto denptr = &aoDensityMatrix;

    // set up poinet to density grid

    auto dgridptr = &densityGrid;

    // generate density on grid points

    #pragma omp parallel shared(tbsizes, tbpositions, ntasks, mgx, mgy, mgz, gtovec, xgtovec, denptr, dgridptr)
    {
    #pragma omp single nowait
        {
            for (int32_t i = 0; i < ntasks; i++)
            {
                // set up task parameters

                auto tbsize = tbsizes[i];

                auto tbposition = tbpositions[i];

                // generate task

                #pragma omp task firstprivate(tbsize, tbposition)
                {
                    _genBatchOfRestrictedDensityGridPointsForGga(dgridptr, denptr, gtovec, xgtovec, mgx, mgy, mgz, tbposition, tbsize);
                }
            }
        }
    }

    // finalize density grid

    densityGrid.updateBetaDensities();

    // destroy GTOs container

    delete gtovec;

    delete xgtovec;
}

void
CDensityGradientGridDriver::_genBatchOfRestrictedDensityGridPointsForGga(CDensityGrid*           densityGrid,
                                                                         const CAODensityMatrix* aoDensityMatrix,
                                                                         const CGtoContainer*    gtoContainer,
                                                                         const CGtoContainer*    atmGtoContainer,
                                                                         const double*           gridCoordinatesX,
                                                                         const double*           gridCoordinatesY,
                                                                         const double*           gridCoordinatesZ,
                                                                         const int32_t           gridOffset,
                                                                         const int32_t           nGridPoints) const
{
    // set up number of AOs

    auto naos = gtoContainer->getNumberOfAtomicOrbitals();

    auto atmnaos = atmGtoContainer->getNumberOfAtomicOrbitals();

    // determine number of grid blocks

    auto blockdim = _getSizeOfBlock();

    auto nblocks = nGridPoints / blockdim;

    // set up current grid point

    int32_t igpnt = 0;

    // loop over grid points blocks

    if (nblocks > 0)
    {
        CMemBlock2D<double> bgaos(blockdim, naos);

        CMemBlock2D<double> bgaox(blockdim, naos);

        CMemBlock2D<double> bgaoy(blockdim, naos);

        CMemBlock2D<double> bgaoz(blockdim, naos);

        CMemBlock2D<double> kgaos(blockdim, naos);

        CMemBlock2D<double> kgaox(blockdim, naos);

        CMemBlock2D<double> kgaoy(blockdim, naos);

        CMemBlock2D<double> kgaoz(blockdim, naos);

        CMemBlock2D<double> kgaoxx(blockdim, naos);

        CMemBlock2D<double> kgaoxy(blockdim, naos);

        CMemBlock2D<double> kgaoxz(blockdim, naos);

        CMemBlock2D<double> kgaoyy(blockdim, naos);

        CMemBlock2D<double> kgaoyz(blockdim, naos);

        CMemBlock2D<double> kgaozz(blockdim, naos);

        CMemBlock<int32_t> aoidx(atmnaos);

        for (int32_t i = 0; i < nblocks; i++)
        {
            bgaos.zero();

            bgaox.zero();

            bgaoy.zero();

            bgaoz.zero();

            kgaos.zero();

            kgaox.zero();

            kgaoy.zero();

            kgaoz.zero();

            kgaoxx.zero();

            kgaoxy.zero();

            kgaoxz.zero();

            kgaoyy.zero();

            kgaoyz.zero();

            kgaozz.zero();

            gtorec::computeGtosValuesForGGA(
                bgaos, bgaox, bgaoy, bgaoz, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, igpnt, blockdim);

            gtorec::computeGtosValuesForGGA2(aoidx,
                                             kgaos,
                                             kgaox,
                                             kgaoy,
                                             kgaoz,
                                             kgaoxx,
                                             kgaoxy,
                                             kgaoxz,
                                             kgaoyy,
                                             kgaoyz,
                                             kgaozz,
                                             atmGtoContainer,
                                             gridCoordinatesX,
                                             gridCoordinatesY,
                                             gridCoordinatesZ,
                                             gridOffset,
                                             igpnt,
                                             blockdim);

            _distRestrictedDensityValuesForGga(
                densityGrid, aoDensityMatrix, bgaos, bgaox, bgaoy, bgaoz, aoidx, kgaox, kgaoy, kgaoz, kgaoxx, kgaoxy, kgaoxz, kgaoyy, kgaoyz, kgaozz, gridOffset, igpnt, blockdim);

            igpnt += blockdim;
        }
    }

    // comopute remaining grid points block

    blockdim = nGridPoints % blockdim;

    if (blockdim > 0)
    {
        CMemBlock2D<double> bgaos(blockdim, naos);

        CMemBlock2D<double> bgaox(blockdim, naos);

        CMemBlock2D<double> bgaoy(blockdim, naos);

        CMemBlock2D<double> bgaoz(blockdim, naos);

        CMemBlock2D<double> kgaos(blockdim, naos);

        CMemBlock2D<double> kgaox(blockdim, naos);

        CMemBlock2D<double> kgaoy(blockdim, naos);

        CMemBlock2D<double> kgaoz(blockdim, naos);

        CMemBlock2D<double> kgaoxx(blockdim, naos);

        CMemBlock2D<double> kgaoxy(blockdim, naos);

        CMemBlock2D<double> kgaoxz(blockdim, naos);

        CMemBlock2D<double> kgaoyy(blockdim, naos);

        CMemBlock2D<double> kgaoyz(blockdim, naos);

        CMemBlock2D<double> kgaozz(blockdim, naos);

        CMemBlock<int32_t> aoidx(atmnaos);

        bgaos.zero();

        bgaox.zero();

        bgaoy.zero();

        bgaoz.zero();

        kgaos.zero();

        kgaox.zero();

        kgaoy.zero();

        kgaoz.zero();

        kgaoxx.zero();

        kgaoxy.zero();

        kgaoxz.zero();

        kgaoyy.zero();

        kgaoyz.zero();

        kgaozz.zero();

        gtorec::computeGtosValuesForGGA(
            bgaos, bgaox, bgaoy, bgaoz, gtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, igpnt, blockdim);

        gtorec::computeGtosValuesForGGA2(
            aoidx, kgaos, kgaox, kgaoy, kgaoz, kgaoxx, kgaoxy, kgaoxz, kgaoyy, kgaoyz, kgaozz,
            atmGtoContainer, gridCoordinatesX, gridCoordinatesY, gridCoordinatesZ, gridOffset, igpnt, blockdim);

        _distRestrictedDensityValuesForGga(
            densityGrid, aoDensityMatrix, bgaos, bgaox, bgaoy, bgaoz, aoidx, kgaox, kgaoy, kgaoz, kgaoxx, kgaoxy, kgaoxz, kgaoyy, kgaoyz, kgaozz, gridOffset, igpnt, blockdim);
    }
}

void
CDensityGradientGridDriver::_distRestrictedDensityValuesForGga(CDensityGrid*              densityGrid,
                                                               const CAODensityMatrix*    aoDensityMatrix,
                                                               const CMemBlock2D<double>& braGtoValues,
                                                               const CMemBlock2D<double>& braGtoValuesX,
                                                               const CMemBlock2D<double>& braGtoValuesY,
                                                               const CMemBlock2D<double>& braGtoValuesZ,
                                                               const CMemBlock<int32_t>&  aoIdentifiers,
                                                               const CMemBlock2D<double>& ketGtoValuesX,
                                                               const CMemBlock2D<double>& ketGtoValuesY,
                                                               const CMemBlock2D<double>& ketGtoValuesZ,
                                                               const CMemBlock2D<double>& ketGtoValuesXX,
                                                               const CMemBlock2D<double>& ketGtoValuesXY,
                                                               const CMemBlock2D<double>& ketGtoValuesXZ,
                                                               const CMemBlock2D<double>& ketGtoValuesYY,
                                                               const CMemBlock2D<double>& ketGtoValuesYZ,
                                                               const CMemBlock2D<double>& ketGtoValuesZZ,
                                                               const int32_t              gridOffset,
                                                               const int32_t              gridBlockPosition,
                                                               const int32_t              nGridPoints) const
{
    if (aoDensityMatrix->getNumberOfDensityMatrices() == 1)
    {
        // set up pointer to density grid data

        auto rhoax = densityGrid->getComponent(0);

        auto rhoay = densityGrid->getComponent(1);

        auto rhoaz = densityGrid->getComponent(2);

        auto rhoaxx = densityGrid->getComponent(3);

        auto rhoaxy = densityGrid->getComponent(4);

        auto rhoaxz = densityGrid->getComponent(5);

        auto rhoayx = densityGrid->getComponent(6);

        auto rhoayy = densityGrid->getComponent(7);

        auto rhoayz = densityGrid->getComponent(8);

        auto rhoazx = densityGrid->getComponent(9);

        auto rhoazy = densityGrid->getComponent(10);

        auto rhoazz = densityGrid->getComponent(11);

        // set up pointer to density matrix data

        auto denmat = aoDensityMatrix->alphaDensity(0);

        auto naos = aoDensityMatrix->getNumberOfRows(0);

        auto atmnaos = aoIdentifiers.size();

        // loop over density matrix

        for (int32_t i = 0; i < naos; i++)
        {
            auto bgaos = braGtoValues.data(i);

            auto bgaox = braGtoValuesX.data(i);

            auto bgaoy = braGtoValuesY.data(i);

            auto bgaoz = braGtoValuesZ.data(i);

            for (int32_t j = 0; j < atmnaos; j++)
            {
                auto kgaox = ketGtoValuesX.data(j);

                auto kgaoy = ketGtoValuesY.data(j);

                auto kgaoz = ketGtoValuesZ.data(j);

                auto kgaoxx = ketGtoValuesXX.data(j);

                auto kgaoxy = ketGtoValuesXY.data(j);

                auto kgaoxz = ketGtoValuesXZ.data(j);

                auto kgaoyy = ketGtoValuesYY.data(j);

                auto kgaoyz = ketGtoValuesYZ.data(j);

                auto kgaozz = ketGtoValuesZZ.data(j);

                const auto dval = 2.0 * denmat[i * naos + aoIdentifiers.at(j)];

                #pragma omp simd
                for (int32_t k = 0; k < nGridPoints; k++)
                {
                    rhoax[gridOffset + gridBlockPosition + k] -= dval * bgaos[k] * kgaox[k];

                    rhoay[gridOffset + gridBlockPosition + k] -= dval * bgaos[k] * kgaoy[k];

                    rhoaz[gridOffset + gridBlockPosition + k] -= dval * bgaos[k] * kgaoz[k];

                    rhoaxx[gridOffset + gridBlockPosition + k] -= dval * (bgaos[k] * kgaoxx[k] + bgaox[k] * kgaox[k]);

                    rhoaxy[gridOffset + gridBlockPosition + k] -= dval * (bgaos[k] * kgaoxy[k] + bgaoy[k] * kgaox[k]);

                    rhoaxz[gridOffset + gridBlockPosition + k] -= dval * (bgaos[k] * kgaoxz[k] + bgaoz[k] * kgaox[k]);

                    rhoayx[gridOffset + gridBlockPosition + k] -= dval * (bgaos[k] * kgaoxy[k] + bgaox[k] * kgaoy[k]);

                    rhoayy[gridOffset + gridBlockPosition + k] -= dval * (bgaos[k] * kgaoyy[k] + bgaoy[k] * kgaoy[k]);

                    rhoayz[gridOffset + gridBlockPosition + k] -= dval * (bgaos[k] * kgaoyz[k] + bgaoz[k] * kgaoy[k]);

                    rhoazx[gridOffset + gridBlockPosition + k] -= dval * (bgaos[k] * kgaoxz[k] + bgaox[k] * kgaoz[k]);

                    rhoazy[gridOffset + gridBlockPosition + k] -= dval * (bgaos[k] * kgaoyz[k] + bgaoy[k] * kgaoz[k]);

                    rhoazz[gridOffset + gridBlockPosition + k] -= dval * (bgaos[k] * kgaozz[k] + bgaoz[k] * kgaoz[k]);
                }
            }
        }
    }
}

int32_t
CDensityGradientGridDriver::_getSizeOfBlock() const
{
    return 500;
}
