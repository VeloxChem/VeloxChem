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

#include "OverlapIntegralsDriver.hpp"

#include <mpi.h>

#include "AngularMomentum.hpp"
#include "GenFunc.hpp"
#include "GenIntsFunc.hpp"
#include "MemBlock.hpp"
#include "MpiFunc.hpp"
#include "OneIntsFunc.hpp"
#include "RecursionFunctionsList.hpp"
#include "StringFormat.hpp"
#include "TwoCentersRecursionFunctions.hpp"

#include "OverlapRecFuncForDX.hpp"
#include "OverlapRecFuncForFF.hpp"
#include "OverlapRecFuncForFG.hpp"
#include "OverlapRecFuncForGF.hpp"
#include "OverlapRecFuncForGG.hpp"
#include "OverlapRecFuncForPX.hpp"
#include "OverlapRecFuncForSX.hpp"
#include "GeomRecFunc.hpp"

COverlapIntegralsDriver::COverlapIntegralsDriver(MPI_Comm comm)
{
    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    _locComm = comm;
}

COverlapIntegralsDriver::~COverlapIntegralsDriver()
{
}

COverlapMatrix
COverlapIntegralsDriver::compute(const CMolecule&       molecule,
                                 const CMolecularBasis& basis) const
{
    COverlapMatrix ovlmat;

    if (_locRank == mpi::master())
    {
        // set up GTOs container

        CGtoContainer bracontr(molecule, basis);

        // compute overlap integrals

        ovlmat = _compOverlapIntegrals(&bracontr, &bracontr);
    }

    return ovlmat;
}

COverlapMatrix
COverlapIntegralsDriver::compute(const CMolecule&       molecule,
                                 const CMolecularBasis& basis,
                                 const int32_t          iAtom,
                                 const char             axis) const
{
    COverlapMatrix ovlmat;

    if (_locRank == mpi::master())
    {
        // set up GTOs container

        CGtoContainer bracontr(molecule, basis, iAtom, 1);
        
        CGtoContainer ketcontr(molecule, basis);

        // compute overlap integrals

        ovlmat = _compGeomOverlapIntegrals(&bracontr, &ketcontr, axis);
    }

    return ovlmat;
}

COverlapMatrix
COverlapIntegralsDriver::compute(const CMolecule&       molecule,
                                 const CMolecularBasis& braBasis,
                                 const CMolecularBasis& ketBasis) const
{
    COverlapMatrix ovlmat;

    if (_locRank == mpi::master())
    {
        // set up GTOs containers

        CGtoContainer bracontr(molecule, braBasis);

        CGtoContainer ketcontr(molecule, ketBasis);

        // compute overlap integrals

        ovlmat = _compOverlapIntegrals(&bracontr, &ketcontr);
    }

    return ovlmat;
}

COverlapMatrix
COverlapIntegralsDriver::compute(const CMolecule&       braMolecule,
                                 const CMolecule&       ketMolecule,
                                 const CMolecularBasis& basis) const
{
    COverlapMatrix ovlmat;

    if (_locRank == mpi::master())
    {
        // set up GTOs containers

        CGtoContainer bracontr(braMolecule, basis);

        CGtoContainer ketcontr(ketMolecule, basis);

        // compute overlap integrals

        ovlmat = _compOverlapIntegrals(&bracontr, &ketcontr);
    }

    return ovlmat;
}

COverlapMatrix
COverlapIntegralsDriver::compute(const CMolecule&       braMolecule,
                                 const CMolecule&       ketMolecule,
                                 const CMolecularBasis& braBasis,
                                 const CMolecularBasis& ketBasis) const
{
    COverlapMatrix ovlmat;

    if (_locRank == mpi::master())
    {
        // set up GTOs containers

        CGtoContainer bracontr(braMolecule, braBasis);

        CGtoContainer ketcontr(ketMolecule, ketBasis);

        // compute overlap integrals

        ovlmat = _compOverlapIntegrals(&bracontr, &ketcontr);
    }

    return ovlmat;
}

void
COverlapIntegralsDriver::compute(      double*    intsBatch,
                                 const CGtoBlock& braGtoBlock,
                                 const CGtoBlock& ketGtoBlock) const
{
    // determine dimensions of integrals batch

    auto nrow = angmom::to_SphericalComponents(braGtoBlock.getAngularMomentum())

                * braGtoBlock.getNumberOfContrGtos();

    auto ncol = angmom::to_SphericalComponents(ketGtoBlock.getAngularMomentum())

                * ketGtoBlock.getNumberOfContrGtos();

    // set up distribution pattern

    COneIntsDistribution dist(intsBatch, nrow, ncol, dist1e::batch);

    // compute overlap integrals

    _compOverlapForGtoBlocks(&dist, braGtoBlock, ketGtoBlock);
}

COverlapMatrix
COverlapIntegralsDriver::_compOverlapIntegrals(const CGtoContainer* braGtoContainer,
                                               const CGtoContainer* ketGtoContainer) const
{
    // check if GTOs containers are same on bra and ket sides

    auto symbk = ((*braGtoContainer) == (*ketGtoContainer));

    // determine dimensions of overlap matrix

    auto nrow = braGtoContainer->getNumberOfAtomicOrbitals();

    auto ncol = ketGtoContainer->getNumberOfAtomicOrbitals();

    // allocate dense matrix for overlap integrals

    CDenseMatrix ovlmat(nrow, ncol);

    // set up distributio pattern

    dist1e dstyp = (symbk) ? dist1e::symsq : dist1e::rect;

    COneIntsDistribution* distpat = new COneIntsDistribution(ovlmat.values(), nrow, ncol, dstyp);

    // compute overlap integral blocks

    #pragma omp parallel shared(braGtoContainer, ketGtoContainer, distpat, symbk)
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

                        _compOverlapForGtoBlocks(distpat, bgtos, kgtos);
                    }
                }
            }
        }
    }

    // deallocate distribution pattern

    delete distpat;

    return COverlapMatrix(ovlmat);
}

COverlapMatrix
COverlapIntegralsDriver::_compGeomOverlapIntegrals(const CGtoContainer* braGtoContainer,
                                                   const CGtoContainer* ketGtoContainer,
                                                   const char           axis) const
{
    // determine dimensions of overlap matrix

    auto nrow = ketGtoContainer->getNumberOfAtomicOrbitals();
    
    auto ncol = nrow; 

    // allocate dense matrix for overlap integrals

    CDenseMatrix ovlmat(nrow, ncol);
    
    ovlmat.zero(); 

    // set up distributio pattern

    COneIntsDistribution* distpat = new COneIntsDistribution(ovlmat.values(), nrow, ncol, dist1e::symsq);

    // compute overlap integral blocks

    #pragma omp parallel shared(braGtoContainer, ketGtoContainer, distpat)
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

                for (int32_t j = 0; j < nket; j++)
                {
                    #pragma omp task firstprivate(j)
                    {
                        auto kgtos = ketGtoContainer->getGtoBlock(j);

                        _compGeomOverlapForGtoBlocks(distpat, bgtos, kgtos, axis);
                    }
                }
            }
        }
    }

    // deallocate distribution pattern

    delete distpat;

    return COverlapMatrix(ovlmat);
}

void
COverlapIntegralsDriver::_compOverlapForGtoBlocks(      COneIntsDistribution* distPattern,
                                                  const CGtoBlock&            braGtoBlock,
                                                  const CGtoBlock&            ketGtoBlock) const
{
    // copy GTOs blocks for bra and ket sides

    auto bragtos = braGtoBlock;

    auto ketgtos = ketGtoBlock;

    // copy distribution pattern

    auto distpat = *distPattern;

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

    // set up tensors of PA and PB distances

    auto rpa = (bang > 0) ? CMemBlock2D<double>(pdim, 3 * pmax) : CMemBlock2D<double>();

    auto rpb = (kang > 0) ? CMemBlock2D<double>(pdim, 3 * pmax) : CMemBlock2D<double>();

    // allocate primitives and auxilary integrals buffer

    auto recmap = _setRecursionMap(bang, kang, pmax);

    auto nblock = recmap.getNumberOfComponents();

    CMemBlock2D<double> primbuffer(pdim, nblock);

    // allocate contracted Cartesian integrals buffer

    auto kdim = ketgtos.getNumberOfContrGtos();

    auto ncart = angmom::to_CartesianComponents(bang, kang);

    CMemBlock2D<double> cartbuffer(kdim, ncart);

    // allocate contracted spherical integrals buffer

    auto nspher = angmom::to_SphericalComponents(bang, kang);

    CMemBlock2D<double> spherbuffer(kdim, nspher);

    // determine bra and ket sides symmetry

    bool symbk = (bragtos == ketgtos);

    for (int32_t i = 0; i < bragtos.getNumberOfContrGtos(); i++)
    {
        // compute distances: R(AB) = A - B

        intsfunc::compDistancesAB(rab, bragtos, ketgtos, i);

        // compute Obara-Saika recursion factors

        intsfunc::compFactorsForOverlap(rfacts, bragtos, ketgtos, i);

        // compute tensors of distances: R(PA) = P - A

        intsfunc::compDistancesPA(rpa, rab, rfacts, 2, bragtos, ketgtos, i);

        // compute tensors of distances: R(PB) = P - B

        intsfunc::compDistancesPB(rpb, rab, rfacts, 2, bragtos, ketgtos, i);

        // compite primitive overlap integrals

        _compPrimOverlapInts(primbuffer, recmap, rfacts, 2, rab, rpa, rpb, bragtos, ketgtos, i);

        // contract primitive overlap integrals

        genfunc::contract(cartbuffer, primbuffer, 0, bragtos, ketgtos, i);

        // transform Cartesian to spherical integrals

        genfunc::transform(spherbuffer, cartbuffer, bmom, kmom, 0, 0, kdim);

        // add batch of integrals to integrals matrix

        distpat.distribute(spherbuffer, bragtos, ketgtos, symbk, i);
    }
}

void
COverlapIntegralsDriver::_compGeomOverlapForGtoBlocks(      COneIntsDistribution* distPattern,
                                                      const CGtoBlock&            braGtoBlock,
                                                      const CGtoBlock&            ketGtoBlock,
                                                      const char                  axis) const
{
    // copy GTOs blocks for bra and ket sides

    auto bragtos = braGtoBlock;

    auto ketgtos = ketGtoBlock;

    // copy distribution pattern

    auto distpat = *distPattern;

    // set up angular momentum data

    auto bang = bragtos.getAngularMomentum();

    auto kang = ketgtos.getAngularMomentum();
    
    auto mang = bang + 1;

    // set up spherical angular momentum for bra and ket sides

    CSphericalMomentum bmom(bang);

    CSphericalMomentum kmom(kang);

    // allocate prefactors used in Obara-Saika recursion

    auto pdim = ketgtos.getNumberOfPrimGtos();

    CMemBlock2D<double> rab(pdim, 3);

    auto pmax = bragtos.getMaxContractionDepth();

    CMemBlock2D<double> rfacts(pdim, 3 * pmax);

    // set up tensors of PA and PB distances

    auto rpa = CMemBlock2D<double>(pdim, 3 * pmax);

    auto rpb = (kang > 0) ? CMemBlock2D<double>(pdim, 3 * pmax) : CMemBlock2D<double>();

    // allocate primitives and auxilary integrals buffer

    auto recmap = _setRecursionMap(mang, kang, pmax);

    auto nblock = recmap.getNumberOfComponents();
    
    const auto ridx = recmap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                           {bang, -1, -1, -1},
                                                           {kang, -1, -1, -1},
                                                           1, 1, 0));
    
//    std::cout << "(" << bang << "," << kang << "): " << pmax << " -> " << nblock  << " :: " << ridx << std::endl;

    CMemBlock2D<double> primbuffer(pdim, nblock);

    // allocate contracted Cartesian integrals buffer

    auto kdim = ketgtos.getNumberOfContrGtos();

    auto ncart = angmom::to_CartesianComponents(bang, kang);

    CMemBlock2D<double> cartbuffer(kdim, ncart);

    // allocate contracted spherical integrals buffer

    auto nspher = angmom::to_SphericalComponents(bang, kang);

    CMemBlock2D<double> spherbuffer(kdim, nspher);
    
    for (int32_t i = 0; i < bragtos.getNumberOfContrGtos(); i++)
    {
        // compute distances: R(AB) = A - B

        intsfunc::compDistancesAB(rab, bragtos, ketgtos, i);

        // compute Obara-Saika recursion factors

        intsfunc::compFactorsForGeomOverlap(rfacts, bragtos, ketgtos, i);

        // compute tensors of distances: R(PA) = P - A

        intsfunc::compDistancesPA(rpa, rab, rfacts, 3, bragtos, ketgtos, i);

        // compute tensors of distances: R(PB) = P - B

        intsfunc::compDistancesPB(rpb, rab, rfacts, 3, bragtos, ketgtos, i);

        // compite primitive overlap integrals

        _compPrimOverlapInts(primbuffer, recmap, rfacts, 3, rab, rpa, rpb, bragtos, ketgtos, i);
        
        // compite primitive overlap integrals

        _compPrimGeomOverlapInts(primbuffer, recmap, rfacts, 3, bragtos, ketgtos, i, axis);

        // contract primitive overlap integrals

        genfunc::contract(cartbuffer, primbuffer, ridx, bragtos, ketgtos, i);

        // transform Cartesian to spherical integrals

        genfunc::transform(spherbuffer, cartbuffer, bmom, kmom, 0, 0, kdim);

        // add batch of integrals to integrals matrix

        distpat.distribute(spherbuffer, bragtos, ketgtos, false, i);
    }
}

void
COverlapIntegralsDriver::_compPrimOverlapInts(CMemBlock2D<double>&       primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const int32_t              nFactors, 
                                              const CMemBlock2D<double>& abDistances,
                                              const CMemBlock2D<double>& paDistances,
                                              const CMemBlock2D<double>& pbDistances,
                                              const CGtoBlock&           braGtoBlock,
                                              const CGtoBlock&           ketGtoBlock,
                                              const int32_t              iContrGto) const
{
    ovlrecfunc::compOverlapForSS(primBuffer, recursionMap, osFactors, nFactors, abDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForSP(primBuffer, recursionMap, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForPS(primBuffer, recursionMap, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForSD(primBuffer, recursionMap, osFactors, nFactors, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForDS(primBuffer, recursionMap, osFactors, nFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForSF(primBuffer, recursionMap, osFactors, nFactors, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForFS(primBuffer, recursionMap, osFactors, nFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForSG(primBuffer, recursionMap, osFactors, nFactors, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForGS(primBuffer, recursionMap, osFactors, nFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForPP(primBuffer, recursionMap, osFactors, nFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForPD(primBuffer, recursionMap, osFactors, nFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForDP(primBuffer, recursionMap, osFactors, nFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForPF(primBuffer, recursionMap, osFactors, nFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForFP(primBuffer, recursionMap, osFactors, nFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForPG(primBuffer, recursionMap, osFactors, nFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForGP(primBuffer, recursionMap, osFactors, nFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForDD(primBuffer, recursionMap, osFactors, nFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForDF(primBuffer, recursionMap, osFactors, nFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForFD(primBuffer, recursionMap, osFactors, nFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForDG(primBuffer, recursionMap, osFactors, nFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForGD(primBuffer, recursionMap, osFactors, nFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForFF(primBuffer, recursionMap, osFactors, nFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForFG(primBuffer, recursionMap, osFactors, nFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForGF(primBuffer, recursionMap, osFactors, nFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForGG(primBuffer, recursionMap, osFactors, nFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    // NOTE: Add l > 4 terms if needed
}

void
COverlapIntegralsDriver::_compPrimGeomOverlapInts(CMemBlock2D<double>&       primBuffer,
                                                  const CRecursionMap&       recursionMap,
                                                  const CMemBlock2D<double>& osFactors,
                                                  const int32_t              nFactors,
                                                  const CGtoBlock&           braGtoBlock,
                                                  const CGtoBlock&           ketGtoBlock,
                                                  const int32_t              iContrGto,
                                                  const char                 axis) const
{
    // get angular momentum data
    
    const auto bang = braGtoBlock.getAngularMomentum();
    
    const auto kang = ketGtoBlock.getAngularMomentum();
    
    // get indexes of integrals
    
    const auto ridx = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                                 {bang, -1, -1, -1},
                                                                 {kang, -1, -1, -1},
                                                                 1, 1, 0));
    
    const auto uidx = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                                 {bang + 1, -1, -1, -1},
                                                                 {kang, -1, -1, -1},
                                                                 1, 1, 0));
    
    const auto lidx = recursionMap.getIndexOfTerm(CRecursionTerm({"Overlap"}, 0, true,
                                                                 {bang - 1, -1, -1, -1},
                                                                 {kang, -1, -1, -1},
                                                                 1, 1, 0));
    
//    std::cout << "bra: " << bang << " ket: " << kang << " indexes: "  << ridx << " : " << uidx << " : " << lidx << std::endl;
    
    // compute geometrical derivatives
    
    if (bang == 0)
    {
        geomrecfunc::compGeomForSX(primBuffer, osFactors, nFactors, ridx, uidx, braGtoBlock, ketGtoBlock, iContrGto, axis);
    }
    
    if (bang == 1)
    {
        geomrecfunc::compGeomForPX(primBuffer, osFactors, nFactors, ridx, uidx, lidx, braGtoBlock, ketGtoBlock, iContrGto, axis);
    }
    
    if (bang == 2)
    {
        geomrecfunc::compGeomForDX(primBuffer, osFactors, nFactors, ridx, uidx, lidx, braGtoBlock, ketGtoBlock, iContrGto, axis);
    }
    
    // need l > 2 implemented....
}

CRecursionMap
COverlapIntegralsDriver::_setRecursionMap(const int32_t braAngularMomentum,
                                          const int32_t ketAngularMomentum,
                                          const int32_t maxNumberOfPrimitives) const
{
    CRecursionFunctionsList recfuncs;

    recfuncs.add(CRecursionFunction({"Overlap"}, &t2crecfunc::obRecursionForOverlap));

    auto rterm = gintsfunc::genIntegral({"Overlap"}, braAngularMomentum, ketAngularMomentum, 0);

    return gintsfunc::genRecursionMap(rterm, recblock::cc, maxNumberOfPrimitives, recfuncs);
}
