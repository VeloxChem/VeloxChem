//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "OverlapIntegralsDriver.hpp"

#include "AngularMomentum.hpp"
#include "GenFunc.hpp"
#include "GenIntsFunc.hpp"
#include "MemBlock.hpp"
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

COverlapIntegralsDriver::COverlapIntegralsDriver(MPI_Comm comm)
{
    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    mpi::duplicate(comm, &_locComm);
}

COverlapIntegralsDriver::~COverlapIntegralsDriver()
{
    mpi::destroy(&_locComm);
}

COverlapMatrix
COverlapIntegralsDriver::compute(const CMolecule& molecule, const CMolecularBasis& basis) const
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
COverlapIntegralsDriver::compute(const CMolecule& molecule, const CMolecularBasis& braBasis, const CMolecularBasis& ketBasis) const
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
COverlapIntegralsDriver::compute(const CMolecule& braMolecule, const CMolecule& ketMolecule, const CMolecularBasis& basis) const
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
COverlapIntegralsDriver::compute(double* intsBatch, const CGtoBlock& braGtoBlock, const CGtoBlock& ketGtoBlock) const
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
COverlapIntegralsDriver::_compOverlapIntegrals(const CGtoContainer* braGtoContainer, const CGtoContainer* ketGtoContainer) const
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

void
COverlapIntegralsDriver::_compOverlapForGtoBlocks(COneIntsDistribution* distPattern, const CGtoBlock& braGtoBlock, const CGtoBlock& ketGtoBlock) const
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

        _compPrimOverlapInts(primbuffer, recmap, rfacts, rab, rpa, rpb, bragtos, ketgtos, i);

        // contract primitive overlap integrals

        genfunc::contract(cartbuffer, primbuffer, 0, bragtos, ketgtos, i);

        // transform Cartesian to spherical integrals

        genfunc::transform(spherbuffer, cartbuffer, bmom, kmom, 0, 0, kdim);

        // add batch of integrals to integrals matrix

        distpat.distribute(spherbuffer, bragtos, ketgtos, symbk, i);
    }
}

void
COverlapIntegralsDriver::_compPrimOverlapInts(CMemBlock2D<double>&       primBuffer,
                                              const CRecursionMap&       recursionMap,
                                              const CMemBlock2D<double>& osFactors,
                                              const CMemBlock2D<double>& abDistances,
                                              const CMemBlock2D<double>& paDistances,
                                              const CMemBlock2D<double>& pbDistances,
                                              const CGtoBlock&           braGtoBlock,
                                              const CGtoBlock&           ketGtoBlock,
                                              const int32_t              iContrGto) const
{
    ovlrecfunc::compOverlapForSS(primBuffer, recursionMap, osFactors, 2, abDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForSP(primBuffer, recursionMap, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForPS(primBuffer, recursionMap, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForSD(primBuffer, recursionMap, osFactors, 2, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForDS(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForSF(primBuffer, recursionMap, osFactors, 2, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForFS(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForSG(primBuffer, recursionMap, osFactors, 2, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForGS(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForPP(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForPD(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForDP(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForPF(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForFP(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForPG(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForGP(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForDD(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForDF(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForFD(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForDG(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForGD(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForFF(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForFG(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForGF(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForGG(primBuffer, recursionMap, osFactors, 2, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    // NOTE: Add l > 4 terms if needed
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
