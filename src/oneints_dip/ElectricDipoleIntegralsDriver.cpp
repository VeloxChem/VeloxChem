//
//                           VELOXCHEM 1.0-RC
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

#include "ElectricDipoleIntegralsDriver.hpp"

#include <mpi.h>

#include "AngularMomentum.hpp"
#include "GenFunc.hpp"
#include "GenIntsFunc.hpp"
#include "MemBlock.hpp"
#include "MpiFunc.hpp"
#include "OneIntsDistributor.hpp"
#include "OneIntsFunc.hpp"
#include "RecursionFunctionsList.hpp"
#include "TwoCentersRecursionFunctions.hpp"

#include "ElectricDipoleRecFuncForDX.hpp"
#include "ElectricDipoleRecFuncForFF.hpp"
#include "ElectricDipoleRecFuncForFG.hpp"
#include "ElectricDipoleRecFuncForGF.hpp"
#include "ElectricDipoleRecFuncForGG.hpp"
#include "ElectricDipoleRecFuncForPX.hpp"
#include "ElectricDipoleRecFuncForSX.hpp"
#include "OverlapRecFuncForDX.hpp"
#include "OverlapRecFuncForFF.hpp"
#include "OverlapRecFuncForFG.hpp"
#include "OverlapRecFuncForPX.hpp"
#include "OverlapRecFuncForSX.hpp"

CElectricDipoleIntegralsDriver::CElectricDipoleIntegralsDriver(MPI_Comm comm)

    : _xOrigin(0.0)

    , _yOrigin(0.0)

    , _zOrigin(0.0)
{
    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    _locComm = comm;
}

CElectricDipoleIntegralsDriver::~CElectricDipoleIntegralsDriver()
{
}

void
CElectricDipoleIntegralsDriver::setElectricDipoleOrigin(const double xOrigin, const double yOrigin, const double zOrigin)
{
    _xOrigin = xOrigin;

    _yOrigin = yOrigin;

    _zOrigin = zOrigin;
}

std::array<double, 3>
CElectricDipoleIntegralsDriver::getElectricDipoleOrigin() const noexcept
{
    return std::array<double, 3>{{_xOrigin, _yOrigin, _zOrigin}};
}

CElectricDipoleMatrix
CElectricDipoleIntegralsDriver::compute(const CMolecule& molecule, const CMolecularBasis& basis) const
{
    CElectricDipoleMatrix dipmat;

    if (_locRank == mpi::master())
    {
        // set up GTOs container

        CGtoContainer bracontr(molecule, basis);

        // compute electric dipole integrals

        dipmat = _compElectricDipoleIntegrals(&bracontr, &bracontr);
    }

    return dipmat;
}

CElectricDipoleMatrix
CElectricDipoleIntegralsDriver::compute(const CMolecule& molecule, const CMolecularBasis& braBasis, const CMolecularBasis& ketBasis) const
{
    CElectricDipoleMatrix dipmat;

    if (_locRank == mpi::master())
    {
        // set up GTOs containers

        CGtoContainer bracontr(molecule, braBasis);

        CGtoContainer ketcontr(molecule, ketBasis);

        // compute electric dipole integrals

        dipmat = _compElectricDipoleIntegrals(&bracontr, &ketcontr);
    }

    return dipmat;
}

CElectricDipoleMatrix
CElectricDipoleIntegralsDriver::compute(const CMolecule& braMolecule, const CMolecule& ketMolecule, const CMolecularBasis& basis) const
{
    CElectricDipoleMatrix dipmat;

    if (_locRank == mpi::master())
    {
        // set up GTOs containers

        CGtoContainer bracontr(braMolecule, basis);

        CGtoContainer ketcontr(ketMolecule, basis);

        // compute electric dipole integrals

        dipmat = _compElectricDipoleIntegrals(&bracontr, &ketcontr);
    }

    return dipmat;
}

CElectricDipoleMatrix
CElectricDipoleIntegralsDriver::compute(const CMolecule&       braMolecule,
                                        const CMolecule&       ketMolecule,
                                        const CMolecularBasis& braBasis,
                                        const CMolecularBasis& ketBasis) const
{
    CElectricDipoleMatrix dipmat;

    if (_locRank == mpi::master())
    {
        // set up GTOs containers

        CGtoContainer bracontr(braMolecule, braBasis);

        CGtoContainer ketcontr(ketMolecule, ketBasis);

        // compute electric dipole integrals

        dipmat = _compElectricDipoleIntegrals(&bracontr, &ketcontr);
    }

    return dipmat;
}

void
CElectricDipoleIntegralsDriver::compute(double*          intsBatchX,
                                        double*          intsBatchY,
                                        double*          intsBatchZ,
                                        const CGtoBlock& braGtoBlock,
                                        const CGtoBlock& ketGtoBlock) const
{
    // determine dimensions of integrals batch

    auto nrow = angmom::to_SphericalComponents(braGtoBlock.getAngularMomentum())

                * braGtoBlock.getNumberOfContrGtos();

    auto ncol = angmom::to_SphericalComponents(ketGtoBlock.getAngularMomentum())

                * ketGtoBlock.getNumberOfContrGtos();

    // set up distribution pattern

    COneIntsDistribution distx(intsBatchX, nrow, ncol, dist1e::batch);

    COneIntsDistribution disty(intsBatchY, nrow, ncol, dist1e::batch);

    COneIntsDistribution distz(intsBatchZ, nrow, ncol, dist1e::batch);

    // compute electric dipole integrals

    _compElectricDipoleForGtoBlocks(&distx, &disty, &distz, _xOrigin, _yOrigin, _zOrigin, braGtoBlock, ketGtoBlock);
}

CElectricDipoleMatrix
CElectricDipoleIntegralsDriver::_compElectricDipoleIntegrals(const CGtoContainer* braGtoContainer, const CGtoContainer* ketGtoContainer) const
{
    // check if GTOs containers are same on bra and ket sides

    auto symbk = ((*braGtoContainer) == (*ketGtoContainer));

    // determine dimensions of overlap matrix

    auto nrow = braGtoContainer->getNumberOfAtomicOrbitals();

    auto ncol = ketGtoContainer->getNumberOfAtomicOrbitals();

    // allocate dense matrix for electric dipole integrals

    CDenseMatrix dipxmat(nrow, ncol);

    CDenseMatrix dipymat(nrow, ncol);

    CDenseMatrix dipzmat(nrow, ncol);

    // set up distributio pattern

    dist1e dstyp = (symbk) ? dist1e::symsq : dist1e::rect;

    COneIntsDistribution* distpatx = new COneIntsDistribution(dipxmat.values(), nrow, ncol, dstyp);

    COneIntsDistribution* distpaty = new COneIntsDistribution(dipymat.values(), nrow, ncol, dstyp);

    COneIntsDistribution* distpatz = new COneIntsDistribution(dipzmat.values(), nrow, ncol, dstyp);

    // copy origin coordinates

    auto origx = _xOrigin;

    auto origy = _yOrigin;

    auto origz = _zOrigin;

    // compute electric dipole integral blocks

    #pragma omp parallel shared(braGtoContainer, ketGtoContainer, distpatx,\
                                distpaty, distpatz, origx, origy, origz, symbk)
    {
        #pragma omp single nowait
        {
            // determine number of GTOs blocks in bra/ket sides

            auto nbra = braGtoContainer->getNumberOfGtoBlocks();

            auto nket = ketGtoContainer->getNumberOfGtoBlocks();

            // loop over pairs of GTOs blocks

            for (int32_t i = 0; i < nbra; i++)
            {
                auto bgtos = braGtoContainer->getGtoBlock(i);

                auto joff = (symbk) ? i : 0;

                for (int32_t j = joff; j < nket; j++)
                {
                    #pragma omp task firstprivate(j)
                    {
                        auto kgtos = ketGtoContainer->getGtoBlock(j);

                        _compElectricDipoleForGtoBlocks(distpatx, distpaty, distpatz, origx, origy, origz, bgtos, kgtos);
                    }
                }
            }
        }
    }

    // deallocate distribution pattern

    delete distpatx;

    delete distpaty;

    delete distpatz;

    return CElectricDipoleMatrix(dipxmat, dipymat, dipzmat, _xOrigin, _yOrigin, _zOrigin);
}

void
CElectricDipoleIntegralsDriver::_compElectricDipoleForGtoBlocks(COneIntsDistribution* distPatternX,
                                                                COneIntsDistribution* distPatternY,
                                                                COneIntsDistribution* distPatternZ,
                                                                const double          xOrigin,
                                                                const double          yOrigin,
                                                                const double          zOrigin,
                                                                const CGtoBlock&      braGtoBlock,
                                                                const CGtoBlock&      ketGtoBlock) const
{
    // copy GTOs blocks for bra and ket sides

    auto bragtos = braGtoBlock;

    auto ketgtos = ketGtoBlock;

    // copy distribution pattern

    auto distpatx = *distPatternX;

    auto distpaty = *distPatternY;

    auto distpatz = *distPatternZ;

    // set up angular momentum data

    auto bang = bragtos.getAngularMomentum();

    auto kang = ketgtos.getAngularMomentum();

    // set up spherical angular momentum for bra and ket sides

    CSphericalMomentum bmom(bang);

    CSphericalMomentum kmom(kang);

    // allocate prefactors used in Obara-Saika recursion

    auto pdim = ketgtos.getNumberOfPrimGtos();

    CMemBlock2D<double> rab(pdim, 3);

    auto pmax = bragtos.getMaxContractionDepth();

    CMemBlock2D<double> rfacts(pdim, 2 * pmax);

    // allocate P center coordinates

    CMemBlock2D<double> rp(pdim, 3 * pmax);

    // set up PA and PB distances

    auto rpa = (bang > 0) ? CMemBlock2D<double>(pdim, 3 * pmax) : CMemBlock2D<double>();

    auto rpb = (kang > 0) ? CMemBlock2D<double>(pdim, 3 * pmax) : CMemBlock2D<double>();

    // set up PC distances

    auto rpc = CMemBlock2D<double>(pdim, 3 * pmax);

    // allocate primitives and auxilary integrals buffer

    auto recmap = _setRecursionMap(bang, kang, pmax);

    auto nblock = recmap.getNumberOfComponents();

    CMemBlock2D<double> primbuffer(pdim, nblock);

    // set up contracted GTOs dimensions

    auto kdim = ketgtos.getNumberOfContrGtos();

    // allocate contracted Cartesian integrals buffer

    auto ncart = angmom::to_CartesianComponents(bang, kang);

    CMemBlock2D<double> cartbufferx(kdim, ncart);

    CMemBlock2D<double> cartbuffery(kdim, ncart);

    CMemBlock2D<double> cartbufferz(kdim, ncart);

    // allocate contracted spherical integrals buffer

    auto nspher = angmom::to_SphericalComponents(bang, kang);

    CMemBlock2D<double> spherbufferx(kdim, nspher);

    CMemBlock2D<double> spherbuffery(kdim, nspher);

    CMemBlock2D<double> spherbufferz(kdim, nspher);

    // set up indexes for contraction

    auto pidx = recmap.getIndexOfTerm(CRecursionTerm({"Electric Dipole"}, 1, true, {bang, -1, -1, -1}, {kang, -1, -1, -1}, 1, 1, 0));

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // determine bra and ket sides symmetry

    bool symbk = (bragtos == ketgtos);

    for (int32_t i = 0; i < bragtos.getNumberOfContrGtos(); i++)
    {
        // compute bra dimensions and shift

        auto bdim = epos[i] - spos[i];

        auto poff = bdim * ncart;

        // compute distances: R(AB) = A - B

        intsfunc::compDistancesAB(rab, bragtos, ketgtos, i);

        // compute Obara-Saika recursion factors

        intsfunc::compFactorsForOverlap(rfacts, bragtos, ketgtos, i);

        // compute coordinates of center P

        intsfunc::compCoordinatesForP(rp, rfacts, 2, bragtos, ketgtos, i);

        // compute distances: R(PA) = P - A

        intsfunc::compDistancesPA(rpa, rp, bragtos, ketgtos, i);

        // compute distances: R(PB) = P - B

        intsfunc::compDistancesPB(rpb, rp, bragtos, ketgtos, i);

        // compute distances: R(PC) = P - C

        intsfunc::compDistancesPC(rpc, rp, xOrigin, yOrigin, zOrigin, bragtos, ketgtos, i);

        // compute primitive electric dipole integrals

        _compPrimElectricDipoleInts(primbuffer, recmap, rfacts, rab, rpa, rpb, rpc, bragtos, ketgtos, i);

        // contract primitive linear momentum integrals

        genfunc::contract(cartbufferx, primbuffer, pidx, bragtos, ketgtos, i);

        genfunc::contract(cartbuffery, primbuffer, pidx + poff, bragtos, ketgtos, i);

        genfunc::contract(cartbufferz, primbuffer, pidx + 2 * poff, bragtos, ketgtos, i);

        // transform Cartesian to spherical integrals

        genfunc::transform(spherbufferx, cartbufferx, bmom, kmom, 0, 0, kdim);

        genfunc::transform(spherbuffery, cartbuffery, bmom, kmom, 0, 0, kdim);

        genfunc::transform(spherbufferz, cartbufferz, bmom, kmom, 0, 0, kdim);

        // add batch of integrals to integrals matrix

        distpatx.distribute(spherbufferx, bragtos, ketgtos, symbk, i);

        distpaty.distribute(spherbuffery, bragtos, ketgtos, symbk, i);

        distpatz.distribute(spherbufferz, bragtos, ketgtos, symbk, i);
    }
}

void
CElectricDipoleIntegralsDriver::_compPrimElectricDipoleInts(CMemBlock2D<double>&       primBuffer,
                                                            const CRecursionMap&       recursionMap,
                                                            const CMemBlock2D<double>& osFactors,
                                                            const CMemBlock2D<double>& abDistances,
                                                            const CMemBlock2D<double>& paDistances,
                                                            const CMemBlock2D<double>& pbDistances,
                                                            const CMemBlock2D<double>& pcDistances,
                                                            const CGtoBlock&           braGtoBlock,
                                                            const CGtoBlock&           ketGtoBlock,
                                                            const int32_t              iContrGto) const
{
    ovlrecfunc::compOverlapForSS(primBuffer, recursionMap, osFactors, 2, abDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForSS(primBuffer, recursionMap, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForSP(primBuffer, recursionMap, osFactors, 2, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForPS(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForSP(primBuffer, recursionMap, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForPS(primBuffer, recursionMap, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForSD(primBuffer, recursionMap, osFactors, 2, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForDS(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForSD(primBuffer, recursionMap, osFactors, 2, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForDS(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForSF(primBuffer, recursionMap, osFactors, 2, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForFS(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForSF(primBuffer, recursionMap, osFactors, 2, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForFS(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForSG(primBuffer, recursionMap, osFactors, 2, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForGS(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForPP(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForPP(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForPD(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForDP(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForPD(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForDP(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForPF(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForFP(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForPF(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForFP(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForSG(primBuffer, recursionMap, osFactors, 2, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForPG(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForGP(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForDD(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForDD(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForDF(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForFD(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForDF(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForFD(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForPG(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForDG(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForGD(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForFF(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForFF(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForDG(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForFG(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForGF(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForFG(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ediprecfunc::compElectricDipoleForGG(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    // NOTE: add l > 4 recursion here
}

CRecursionMap
CElectricDipoleIntegralsDriver::_setRecursionMap(const int32_t braAngularMomentum,
                                                 const int32_t ketAngularMomentum,
                                                 const int32_t maxNumberOfPrimitives) const
{
    CRecursionFunctionsList recfuncs;

    recfuncs.add(CRecursionFunction({"Overlap"}, &t2crecfunc::obRecursionForOverlap));

    recfuncs.add(CRecursionFunction({"Electric Dipole"}, &t2crecfunc::obRecursionForElectricDipole));

    auto rterm = gintsfunc::genIntegral({"Electric Dipole"}, braAngularMomentum, ketAngularMomentum, 0);

    return gintsfunc::genRecursionMap(rterm, recblock::cc, maxNumberOfPrimitives, recfuncs);
}
