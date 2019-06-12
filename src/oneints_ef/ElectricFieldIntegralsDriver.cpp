//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "ElectricFieldIntegralsDriver.hpp"

#include "AngularMomentum.hpp"
#include "GenFunc.hpp"
#include "GenIntsFunc.hpp"
#include "OneIntsFunc.hpp"
#include "RecursionFunctionsList.hpp"
#include "TwoCentersRecursionFunctions.hpp"

#include "ElectricFieldRecFuncForDX.hpp"
#include "ElectricFieldRecFuncForFF.hpp"
#include "ElectricFieldRecFuncForFG.hpp"
#include "ElectricFieldRecFuncForGF.hpp"
#include "ElectricFieldRecFuncForGG.hpp"
#include "ElectricFieldRecFuncForPX.hpp"
#include "ElectricFieldRecFuncForSX.hpp"
#include "NuclearPotentialRecFuncForDX.hpp"
#include "NuclearPotentialRecFuncForFF.hpp"
#include "NuclearPotentialRecFuncForFG.hpp"
#include "NuclearPotentialRecFuncForPX.hpp"
#include "NuclearPotentialRecFuncForSX.hpp"
#include "OverlapRecFuncForSX.hpp"

CElectricFieldIntegralsDriver::CElectricFieldIntegralsDriver(MPI_Comm comm)
{
    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    mpi::duplicate(comm, &_locComm);
}

CElectricFieldIntegralsDriver::~CElectricFieldIntegralsDriver()
{
    mpi::destroy(&_locComm);
}

CElectricFieldMatrix
CElectricFieldIntegralsDriver::compute(const CMolecule&       molecule,
                                       const CMolecularBasis& basis,
                                       const double           coordinateX,
                                       const double           coordinateY,
                                       const double           coordinateZ) const
{
    CElectricFieldMatrix efieldmat;

    if (_locRank == mpi::master())
    {
        // set up GTOs container

        CGtoContainer bracontr(molecule, basis);

        // set up point dipoles data

        auto dipoles = CMemBlock2D<double>({1.0, 1.0, 1.0}, 1, 3);

        auto coords = CMemBlock2D<double>({coordinateX, coordinateY, coordinateZ}, 1, 3);

        // compute electric field integrals

        efieldmat = _compElectricFieldIntegrals(&dipoles, &coords, &bracontr, &bracontr);
    }

    return efieldmat;
}

CElectricFieldMatrix
CElectricFieldIntegralsDriver::compute(const CMolecule&           molecule,
                                       const CMolecularBasis&     basis,
                                       const CMemBlock2D<double>* dipoles,
                                       const CMemBlock2D<double>* coordinates) const
{
    CElectricFieldMatrix efieldmat;

    if (_locRank == mpi::master())
    {
        // set up GTOs container

        CGtoContainer bracontr(molecule, basis);

        // compute electric field integrals

        efieldmat = _compElectricFieldIntegrals(dipoles, coordinates, &bracontr, &bracontr);
    }

    return efieldmat;
}

CElectricFieldMatrix
CElectricFieldIntegralsDriver::compute(const CMolecule&       molecule,
                                       const CMolecularBasis& braBasis,
                                       const CMolecularBasis& ketBasis,
                                       const double           coordinateX,
                                       const double           coordinateY,
                                       const double           coordinateZ) const
{
    CElectricFieldMatrix efieldmat;

    if (_locRank == mpi::master())
    {
        // set up GTOs container

        CGtoContainer bracontr(molecule, braBasis);

        CGtoContainer ketcontr(molecule, ketBasis);

        // set up point dipoles data

        auto dipoles = CMemBlock2D<double>({1.0, 1.0, 1.0}, 1, 3);

        auto coords = CMemBlock2D<double>({coordinateX, coordinateY, coordinateZ}, 1, 3);

        // compute electric field integrals

        efieldmat = _compElectricFieldIntegrals(&dipoles, &coords, &bracontr, &ketcontr);
    }

    return efieldmat;
}

CElectricFieldMatrix
CElectricFieldIntegralsDriver::compute(const CMolecule&           molecule,
                                       const CMolecularBasis&     braBasis,
                                       const CMolecularBasis&     ketBasis,
                                       const CMemBlock2D<double>* dipoles,
                                       const CMemBlock2D<double>* coordinates) const
{
    CElectricFieldMatrix efieldmat;

    if (_locRank == mpi::master())
    {
        // set up GTOs container

        CGtoContainer bracontr(molecule, braBasis);

        CGtoContainer ketcontr(molecule, ketBasis);

        // compute electric field integrals

        efieldmat = _compElectricFieldIntegrals(dipoles, coordinates, &bracontr, &ketcontr);
    }

    return efieldmat;
}

CElectricFieldMatrix
CElectricFieldIntegralsDriver::compute(const CMolecule&       braMolecule,
                                       const CMolecule&       ketMolecule,
                                       const CMolecularBasis& basis,
                                       const double           coordinateX,
                                       const double           coordinateY,
                                       const double           coordinateZ) const
{
    CElectricFieldMatrix efieldmat;

    if (_locRank == mpi::master())
    {
        // set up GTOs container

        CGtoContainer bracontr(braMolecule, basis);

        CGtoContainer ketcontr(ketMolecule, basis);

        // set up point dipoles data

        auto dipoles = CMemBlock2D<double>({1.0, 1.0, 1.0}, 1, 3);

        auto coords = CMemBlock2D<double>({coordinateX, coordinateY, coordinateZ}, 1, 3);

        // compute electric field integrals

        efieldmat = _compElectricFieldIntegrals(&dipoles, &coords, &bracontr, &ketcontr);
    }

    return efieldmat;
}

CElectricFieldMatrix
CElectricFieldIntegralsDriver::compute(const CMolecule&           braMolecule,
                                       const CMolecule&           ketMolecule,
                                       const CMolecularBasis&     basis,
                                       const CMemBlock2D<double>* dipoles,
                                       const CMemBlock2D<double>* coordinates) const
{
    CElectricFieldMatrix efieldmat;

    if (_locRank == mpi::master())
    {
        // set up GTOs container

        CGtoContainer bracontr(braMolecule, basis);

        CGtoContainer ketcontr(ketMolecule, basis);

        // compute electric field integrals

        efieldmat = _compElectricFieldIntegrals(dipoles, coordinates, &bracontr, &ketcontr);
    }

    return efieldmat;
}

CElectricFieldMatrix
CElectricFieldIntegralsDriver::compute(const CMolecule&       braMolecule,
                                       const CMolecule&       ketMolecule,
                                       const CMolecularBasis& braBasis,
                                       const CMolecularBasis& ketBasis,
                                       const double           coordinateX,
                                       const double           coordinateY,
                                       const double           coordinateZ) const
{
    CElectricFieldMatrix efieldmat;

    if (_locRank == mpi::master())
    {
        // set up GTOs container

        CGtoContainer bracontr(braMolecule, braBasis);

        CGtoContainer ketcontr(ketMolecule, ketBasis);

        // set up point dipoles data

        auto dipoles = CMemBlock2D<double>({1.0, 1.0, 1.0}, 1, 3);

        auto coords = CMemBlock2D<double>({coordinateX, coordinateY, coordinateZ}, 1, 3);

        // compute electric field integrals

        efieldmat = _compElectricFieldIntegrals(&dipoles, &coords, &bracontr, &ketcontr);
    }

    return efieldmat;
}

CElectricFieldMatrix
CElectricFieldIntegralsDriver::compute(const CMolecule&           braMolecule,
                                       const CMolecule&           ketMolecule,
                                       const CMolecularBasis&     braBasis,
                                       const CMolecularBasis&     ketBasis,
                                       const CMemBlock2D<double>* dipoles,
                                       const CMemBlock2D<double>* coordinates) const
{
    CElectricFieldMatrix efieldmat;

    if (_locRank == mpi::master())
    {
        // set up GTOs container

        CGtoContainer bracontr(braMolecule, braBasis);

        CGtoContainer ketcontr(ketMolecule, ketBasis);

        // compute electric field integrals

        efieldmat = _compElectricFieldIntegrals(dipoles, coordinates, &bracontr, &ketcontr);
    }

    return efieldmat;
}

void
CElectricFieldIntegralsDriver::compute(double*                    intsBatchX,
                                       double*                    intsBatchY,
                                       double*                    intsBatchZ,
                                       const CMemBlock2D<double>* dipoles,
                                       const CMemBlock2D<double>* coordinates,
                                       const CGtoBlock&           braGtoBlock,
                                       const CGtoBlock&           ketGtoBlock) const
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

    // compute electric field integrals

    _compElectricFieldForGtoBlocks(&distx, &disty, &distz, dipoles, coordinates, braGtoBlock, ketGtoBlock);
}

CElectricFieldMatrix
CElectricFieldIntegralsDriver::_compElectricFieldIntegrals(const CMemBlock2D<double>* dipoles,
                                                           const CMemBlock2D<double>* coordinates,
                                                           const CGtoContainer*       braGtoContainer,
                                                           const CGtoContainer*       ketGtoContainer) const
{
    // check if GTOs containers are same on bra and ket sides

    auto symbk = ((*braGtoContainer) == (*ketGtoContainer));

    // determine dimensions of overlap matrix

    auto nrow = braGtoContainer->getNumberOfAtomicOrbitals();

    auto ncol = ketGtoContainer->getNumberOfAtomicOrbitals();

    // allocate dense matrix for electric field integrals

    CDenseMatrix efxmat(nrow, ncol);

    CDenseMatrix efymat(nrow, ncol);

    CDenseMatrix efzmat(nrow, ncol);

    // set up distributio pattern

    dist1e dstyp = (symbk) ? dist1e::symsq : dist1e::rect;

    COneIntsDistribution* distpatx = new COneIntsDistribution(efxmat.values(), nrow, ncol, dstyp);

    COneIntsDistribution* distpaty = new COneIntsDistribution(efymat.values(), nrow, ncol, dstyp);

    COneIntsDistribution* distpatz = new COneIntsDistribution(efzmat.values(), nrow, ncol, dstyp);

    // compute electric field integral blocks

    #pragma omp parallel shared(braGtoContainer, ketGtoContainer, dipoles,\
                                coordinates, distpatx, distpaty, distpatz, symbk)
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

                        _compElectricFieldForGtoBlocks(distpatx, distpaty, distpatz, dipoles, coordinates, bgtos, kgtos);
                    }
                }
            }
        }
    }

    // deallocate distribution pattern

    delete distpatx;

    delete distpaty;

    delete distpatz;

    return CElectricFieldMatrix(efxmat, efymat, efzmat);
}

void
CElectricFieldIntegralsDriver::_compElectricFieldForGtoBlocks(COneIntsDistribution*      distPatternX,
                                                              COneIntsDistribution*      distPatternY,
                                                              COneIntsDistribution*      distPatternZ,
                                                              const CMemBlock2D<double>* dipoles,
                                                              const CMemBlock2D<double>* coordinates,
                                                              const CGtoBlock&           braGtoBlock,
                                                              const CGtoBlock&           ketGtoBlock) const
{
    // copy GTOs blocks for bra and ket sides

    auto bragtos = braGtoBlock;

    auto ketgtos = ketGtoBlock;

    // copy distribution pattern

    auto distpatx = *distPatternX;

    auto distpaty = *distPatternY;

    auto distpatz = *distPatternZ;

    // copy dipoles and their coordinates

    auto dipvalues = *dipoles;

    auto dipcoords = *coordinates;

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

    CMemBlock2D<double> rfacts(pdim, 3 * pmax);

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

    // allocate primitive integrals accumulation buffer

    auto ncart = angmom::to_CartesianComponents(bang, kang);

    CMemBlock2D<double> accbuffer(pdim, 3 * ncart * pmax);

    // set up contracted GTOs dimensions

    auto kdim = ketgtos.getNumberOfContrGtos();

    // allocate contracted Cartesian integrals buffer

    CMemBlock2D<double> cartbufferx(kdim, ncart);

    CMemBlock2D<double> cartbuffery(kdim, ncart);

    CMemBlock2D<double> cartbufferz(kdim, ncart);

    // allocate contracted spherical integrals buffer

    auto nspher = angmom::to_SphericalComponents(bang, kang);

    CMemBlock2D<double> spherbufferx(kdim, nspher);

    CMemBlock2D<double> spherbuffery(kdim, nspher);

    CMemBlock2D<double> spherbufferz(kdim, nspher);

    // set up indexes for contraction

    auto pidx = recmap.getIndexOfTerm(CRecursionTerm({"Electric Field"}, 1, true, {bang, -1, -1, -1}, {kang, -1, -1, -1}, 1, 1, 0));

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    // set uo Boys function data

    CBoysFunction bftab(bang + kang + 1);

    CMemBlock<double> bargs(pdim);

    CMemBlock2D<double> bvals(pdim, bang + kang + 2);

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

        intsfunc::compFactorsForNuclearPotential(rfacts, bragtos, ketgtos, i);

        // compute coordinates of center P

        intsfunc::compCoordinatesForP(rp, rfacts, 3, bragtos, ketgtos, i);

        // compute distances: R(PA) = P - A

        intsfunc::compDistancesPA(rpa, rp, bragtos, ketgtos, i);

        // compute distances: R(PB) = P - B

        intsfunc::compDistancesPB(rpb, rp, bragtos, ketgtos, i);

        // reset accumulation buffer

        accbuffer.zero();

        // loop over charges

        for (int32_t j = 0; j < dipvalues.size(0); j++)
        {
            // compute distances: R(PC) = P - C

            intsfunc::compDistancesPC(rpc, rp, dipcoords, bragtos, ketgtos, i, j);

            // compute primitive integrals

            _compPrimElectricFieldInts(primbuffer, recmap, bftab, bargs, bvals, bang + kang + 1, rfacts, rab, rpa, rpb, rpc, bragtos, ketgtos, i);

            // add scaled contribution to accumulation buffer

            _addPointDipoleContribution(accbuffer, primbuffer, pidx, dipvalues, bragtos, ketgtos, i, j);
        }

        // contract primitive overlap integrals

        genfunc::contract(cartbufferx, accbuffer, 0, bragtos, ketgtos, i);

        genfunc::contract(cartbuffery, accbuffer, poff, bragtos, ketgtos, i);

        genfunc::contract(cartbufferz, accbuffer, 2 * poff, bragtos, ketgtos, i);

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
CElectricFieldIntegralsDriver::_addPointDipoleContribution(CMemBlock2D<double>&       accBuffer,
                                                           const CMemBlock2D<double>& primBuffer,
                                                           const int32_t              primIndex,
                                                           const CMemBlock2D<double>& dipoles,
                                                           const CGtoBlock&           braGtoBlock,
                                                           const CGtoBlock&           ketGtoBlock,
                                                           const int32_t              iContrGto,
                                                           const int32_t              iPointDipole) const
{
    // set up angular momentum for bra and ket sides

    auto bang = braGtoBlock.getAngularMomentum();

    auto kang = ketGtoBlock.getAngularMomentum();

    auto ncart = angmom::to_CartesianComponents(bang, kang);

    // set up pointers to primitives data on bra side

    auto spos = braGtoBlock.getStartPositions();

    auto epos = braGtoBlock.getEndPositions();

    auto bdim = epos[iContrGto] - spos[iContrGto];

    // set up pointers to primitives data on ket side

    auto nprim = ketGtoBlock.getNumberOfPrimGtos();

    for (int32_t i = 0; i < 3; i++)
    {
        // set up point dipole factor

        auto fact = (dipoles.data(i))[iPointDipole];

        // set up electric field component offset

        auto poff = i * bdim * ncart;

        for (int32_t j = 0; j < bdim; j++)
        {
            for (int32_t k = 0; k < ncart; k++)
            {
                auto abuf = accBuffer.data(poff + j * ncart + k);

                auto pbuf = primBuffer.data(primIndex + poff + j * ncart + k);

                #pragma omp simd aligned(abuf, pbuf: VLX_ALIGN)
                for (int32_t l = 0; l < nprim; l++)
                {
                    abuf[l] += fact * pbuf[l];
                }
            }
        }
    }
}

void
CElectricFieldIntegralsDriver::_compPrimElectricFieldInts(CMemBlock2D<double>&       primBuffer,
                                                          const CRecursionMap&       recursionMap,
                                                          const CBoysFunction&       bfTable,
                                                          CMemBlock<double>&         bfArguments,
                                                          CMemBlock2D<double>&       bfValues,
                                                          const int32_t              bfOrder,
                                                          const CMemBlock2D<double>& osFactors,
                                                          const CMemBlock2D<double>& abDistances,
                                                          const CMemBlock2D<double>& paDistances,
                                                          const CMemBlock2D<double>& pbDistances,
                                                          const CMemBlock2D<double>& pcDistances,
                                                          const CGtoBlock&           braGtoBlock,
                                                          const CGtoBlock&           ketGtoBlock,
                                                          const int32_t              iContrGto) const
{
    ovlrecfunc::compOverlapForSS(primBuffer, recursionMap, osFactors, 3, abDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForSS(
        primBuffer, recursionMap, bfTable, bfArguments, bfValues, bfOrder, osFactors, abDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForSS(primBuffer, recursionMap, osFactors, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForSP(primBuffer, recursionMap, pbDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForPS(primBuffer, recursionMap, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForSP(primBuffer, recursionMap, pbDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForSD(primBuffer, recursionMap, osFactors, pbDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForPS(primBuffer, recursionMap, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForDS(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForSD(primBuffer, recursionMap, osFactors, pbDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForSF(primBuffer, recursionMap, osFactors, pbDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForDS(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForFS(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForSF(primBuffer, recursionMap, osFactors, pbDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForSG(primBuffer, recursionMap, osFactors, pbDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForFS(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGS(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForPP(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForPD(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForPP(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForDP(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForPF(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForDP(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForFP(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForSG(primBuffer, recursionMap, osFactors, pbDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForPG(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForFP(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGP(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForPD(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForDD(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForPF(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForDF(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForDD(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForFD(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForPG(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForDG(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForFD(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGD(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForDF(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForFF(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForFF(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGF(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForDG(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForFG(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    npotrecfunc::compNuclearPotentialForFG(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);

    efieldrecfunc::compElectricFieldForGG(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
}

CRecursionMap
CElectricFieldIntegralsDriver::_setRecursionMap(const int32_t braAngularMomentum,
                                                const int32_t ketAngularMomentum,
                                                const int32_t maxNumberOfPrimitives) const
{
    CRecursionFunctionsList recfuncs;

    recfuncs.add(CRecursionFunction({"Nuclear Potential"}, &t2crecfunc::obRecursionForNuclearPotential));

    recfuncs.add(CRecursionFunction({"Electric Field"}, &t2crecfunc::obRecursionForElectricField));

    auto rterm = gintsfunc::genIntegral({"Electric Field"}, braAngularMomentum, ketAngularMomentum, 0);

    return gintsfunc::genRecursionMap(rterm, recblock::cc, maxNumberOfPrimitives, recfuncs);
}
