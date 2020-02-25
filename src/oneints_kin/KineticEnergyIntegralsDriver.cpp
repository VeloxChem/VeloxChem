//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "KineticEnergyIntegralsDriver.hpp"

#include "AngularMomentum.hpp"
#include "GenFunc.hpp"
#include "GenIntsFunc.hpp"
#include "OneIntsFunc.hpp"
#include "RecursionFunctionsList.hpp"
#include "StringFormat.hpp"
#include "TwoCentersRecursionFunctions.hpp"

#include "KineticEnergyRecFuncForDX.hpp"
#include "KineticEnergyRecFuncForFF.hpp"
#include "KineticEnergyRecFuncForFG.hpp"
#include "KineticEnergyRecFuncForGF.hpp"
#include "KineticEnergyRecFuncForGG.hpp"
#include "KineticEnergyRecFuncForPX.hpp"
#include "KineticEnergyRecFuncForSX.hpp"
#include "OverlapRecFuncForDX.hpp"
#include "OverlapRecFuncForFF.hpp"
#include "OverlapRecFuncForFG.hpp"
#include "OverlapRecFuncForGF.hpp"
#include "OverlapRecFuncForGG.hpp"
#include "OverlapRecFuncForPX.hpp"
#include "OverlapRecFuncForSX.hpp"

CKineticEnergyIntegralsDriver::CKineticEnergyIntegralsDriver(MPI_Comm comm)
{
    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    mpi::duplicate(comm, &_locComm);
}

CKineticEnergyIntegralsDriver::~CKineticEnergyIntegralsDriver()
{
    mpi::destroy(&_locComm);
}

CKineticEnergyMatrix
CKineticEnergyIntegralsDriver::compute(const CMolecule& molecule, const CMolecularBasis& basis) const
{
    CKineticEnergyMatrix kinmat;

    if (_locRank == mpi::master())
    {
        // set up GTOs container

        CGtoContainer bracontr(molecule, basis);

        // compute kinetic energy integrals

        kinmat = _compKineticEnergyIntegrals(&bracontr, &bracontr);
    }

    return kinmat;
}

CKineticEnergyMatrix
CKineticEnergyIntegralsDriver::compute(const CMolecule& molecule, const CMolecularBasis& braBasis, const CMolecularBasis& ketBasis) const
{
    CKineticEnergyMatrix kinmat;

    if (_locRank == mpi::master())
    {
        // set up GTOs containers

        CGtoContainer bracontr(molecule, braBasis);

        CGtoContainer ketcontr(molecule, ketBasis);

        // compute kinetic energy integrals

        kinmat = _compKineticEnergyIntegrals(&bracontr, &ketcontr);
    }

    return kinmat;
}

CKineticEnergyMatrix
CKineticEnergyIntegralsDriver::compute(const CMolecule& braMolecule, const CMolecule& ketMolecule, const CMolecularBasis& basis) const
{
    CKineticEnergyMatrix kinmat;

    if (_locRank == mpi::master())
    {
        // set up GTOs containers

        CGtoContainer bracontr(braMolecule, basis);

        CGtoContainer ketcontr(ketMolecule, basis);

        // compute kinetic energy integrals

        kinmat = _compKineticEnergyIntegrals(&bracontr, &ketcontr);
    }

    return kinmat;
}

CKineticEnergyMatrix
CKineticEnergyIntegralsDriver::compute(const CMolecule&       braMolecule,
                                       const CMolecule&       ketMolecule,
                                       const CMolecularBasis& braBasis,
                                       const CMolecularBasis& ketBasis) const
{
    CKineticEnergyMatrix kinmat;

    if (_locRank == mpi::master())
    {
        // set up GTOs containers

        CGtoContainer bracontr(braMolecule, braBasis);

        CGtoContainer ketcontr(ketMolecule, ketBasis);

        // compute kinetic energy integrals

        kinmat = _compKineticEnergyIntegrals(&bracontr, &ketcontr);
    }

    return kinmat;
}

void
CKineticEnergyIntegralsDriver::compute(double* intsBatch, const CGtoBlock& braGtoBlock, const CGtoBlock& ketGtoBlock) const
{
    // determine dimensions of integrals batch

    auto nrow = angmom::to_SphericalComponents(braGtoBlock.getAngularMomentum())

                * braGtoBlock.getNumberOfContrGtos();

    auto ncol = angmom::to_SphericalComponents(ketGtoBlock.getAngularMomentum())

                * ketGtoBlock.getNumberOfContrGtos();

    // set up distribution pattern

    COneIntsDistribution dist(intsBatch, nrow, ncol, dist1e::batch);

    // compute kinetic energy integrals

    _compKineticEnergyForGtoBlocks(&dist, braGtoBlock, ketGtoBlock);
}

CKineticEnergyMatrix
CKineticEnergyIntegralsDriver::_compKineticEnergyIntegrals(const CGtoContainer* braGtoContainer, const CGtoContainer* ketGtoContainer) const
{
    // check if GTOs containers are same on bra and ket sides

    auto symbk = ((*braGtoContainer) == (*ketGtoContainer));

    // determine dimensions of overlap matrix

    auto nrow = braGtoContainer->getNumberOfAtomicOrbitals();

    auto ncol = ketGtoContainer->getNumberOfAtomicOrbitals();

    // allocate dense matrix for kinetic energy integrals

    CDenseMatrix kinmat(nrow, ncol);

    // set up distributio pattern

    dist1e dstyp = (symbk) ? dist1e::symsq : dist1e::rect;

    COneIntsDistribution* distpat = new COneIntsDistribution(kinmat.values(), nrow, ncol, dstyp);

    // compute kinetic energy integral blocks

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

                        _compKineticEnergyForGtoBlocks(distpat, bgtos, kgtos);
                    }
                }
            }
        }
    }

    // deallocate distribution pattern

    delete distpat;

    return CKineticEnergyMatrix(kinmat);
}

void
CKineticEnergyIntegralsDriver::_compKineticEnergyForGtoBlocks(COneIntsDistribution* distPattern,
                                                              const CGtoBlock&      braGtoBlock,
                                                              const CGtoBlock&      ketGtoBlock) const
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

    CMemBlock2D<double> rfacts(pdim, 4 * pmax);

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

        intsfunc::compFactorsForKineticEnergy(rfacts, bragtos, ketgtos, i);

        // compute distances: R(PA) = P - A

        intsfunc::compDistancesPA(rpa, rab, rfacts, 4, bragtos, ketgtos, i);

        // compute distances: R(PB) = P - B

        intsfunc::compDistancesPB(rpb, rab, rfacts, 4, bragtos, ketgtos, i);

        // compite primitive kinetic energy integrals

        _compPrimKineticEnergyInts(primbuffer, recmap, rfacts, rab, rpa, rpb, bragtos, ketgtos, i);

        // contract primitive kinetic energy integrals

        genfunc::contract(cartbuffer, primbuffer, 0, bragtos, ketgtos, i);

        // transform Cartesian to spherical integrals

        genfunc::transform(spherbuffer, cartbuffer, bmom, kmom, 0, 0, kdim);

        // add batch of integrals to integrals matrix

        distpat.distribute(spherbuffer, bragtos, ketgtos, symbk, i);
    }
}

void
CKineticEnergyIntegralsDriver::_compPrimKineticEnergyInts(CMemBlock2D<double>&       primBuffer,
                                                          const CRecursionMap&       recursionMap,
                                                          const CMemBlock2D<double>& osFactors,
                                                          const CMemBlock2D<double>& abDistances,
                                                          const CMemBlock2D<double>& paDistances,
                                                          const CMemBlock2D<double>& pbDistances,
                                                          const CGtoBlock&           braGtoBlock,
                                                          const CGtoBlock&           ketGtoBlock,
                                                          const int32_t              iContrGto) const
{
    ovlrecfunc::compOverlapForSS(primBuffer, recursionMap, osFactors, 4, abDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForSS(primBuffer, recursionMap, osFactors, abDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForSP(primBuffer, recursionMap, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForSP(primBuffer, recursionMap, osFactors, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForPS(primBuffer, recursionMap, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForPS(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForSD(primBuffer, recursionMap, osFactors, 4, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForSD(primBuffer, recursionMap, osFactors, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForDS(primBuffer, recursionMap, osFactors, 4, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForDS(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForSF(primBuffer, recursionMap, osFactors, 4, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForSF(primBuffer, recursionMap, osFactors, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForFS(primBuffer, recursionMap, osFactors, 4, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForFS(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForSG(primBuffer, recursionMap, osFactors, 4, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForSG(primBuffer, recursionMap, osFactors, pbDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForGS(primBuffer, recursionMap, osFactors, 4, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForGS(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForPP(primBuffer, recursionMap, osFactors, 4, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForPP(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForPD(primBuffer, recursionMap, osFactors, 4, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForPD(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForDP(primBuffer, recursionMap, osFactors, 4, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForDP(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForPF(primBuffer, recursionMap, osFactors, 4, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForPF(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForFP(primBuffer, recursionMap, osFactors, 4, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForFP(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForPG(primBuffer, recursionMap, osFactors, 4, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForPG(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForGP(primBuffer, recursionMap, osFactors, 4, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForGP(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForDD(primBuffer, recursionMap, osFactors, 4, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForDD(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForDF(primBuffer, recursionMap, osFactors, 4, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForDF(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForFD(primBuffer, recursionMap, osFactors, 4, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForFD(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForDG(primBuffer, recursionMap, osFactors, 4, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForDG(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForGD(primBuffer, recursionMap, osFactors, 4, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForGD(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForFF(primBuffer, recursionMap, osFactors, 4, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForFF(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForFG(primBuffer, recursionMap, osFactors, 4, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForFG(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForGF(primBuffer, recursionMap, osFactors, 4, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForGF(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    ovlrecfunc::compOverlapForGG(primBuffer, recursionMap, osFactors, 4, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    kinrecfunc::compKineticEnergyForGG(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);

    // NOTE: Add l > 4 terms if needed
}

CRecursionMap
CKineticEnergyIntegralsDriver::_setRecursionMap(const int32_t braAngularMomentum,
                                                const int32_t ketAngularMomentum,
                                                const int32_t maxNumberOfPrimitives) const
{
    CRecursionFunctionsList recfuncs;

    recfuncs.add(CRecursionFunction({"Overlap"}, &t2crecfunc::obRecursionForOverlap));

    recfuncs.add(CRecursionFunction({"Kinetic Energy"}, &t2crecfunc::obRecursionForKineticEnergy));

    auto rterm = gintsfunc::genIntegral({"Kinetic Energy"}, braAngularMomentum, ketAngularMomentum, 0);

    return gintsfunc::genRecursionMap(rterm, recblock::cc, maxNumberOfPrimitives, recfuncs);
}
