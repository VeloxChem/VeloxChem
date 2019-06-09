//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "AngularMomentumIntegralsDriver.hpp"

#include "AngularMomentum.hpp"
#include "GenFunc.hpp"
#include "MemBlock.hpp"
#include "OneIntsFunc.hpp"
#include "OneIntsDistributor.hpp"
#include "RecursionFunctionsList.hpp"
#include "GenIntsFunc.hpp"
#include "TwoCentersRecursionFunctions.hpp"

#include "OverlapRecFuncForSX.hpp"
#include "LinearMomentumRecFuncForSX.hpp"
#include "ElectricDipoleRecFuncForSX.hpp"
#include "AngularMomentumRecFuncForSX.hpp"

CAngularMomentumIntegralsDriver::CAngularMomentumIntegralsDriver(MPI_Comm comm)

    : _xOrigin(0.0)

    , _yOrigin(0.0)

    , _zOrigin(0.0)
{
    _locRank  = mpi::rank(comm);
    
    _locNodes = mpi::nodes(comm);

    mpi::duplicate(comm, &_locComm);
}

CAngularMomentumIntegralsDriver::~CAngularMomentumIntegralsDriver()
{
    mpi::destroy(&_locComm);
}

void
CAngularMomentumIntegralsDriver::setAngularMomentumOrigin(const double xOrigin,
                                                          const double yOrigin,
                                                          const double zOrigin)
{
    _xOrigin = xOrigin;
    
    _yOrigin = yOrigin;
    
    _zOrigin = zOrigin;
}

CAngularMomentumMatrix
CAngularMomentumIntegralsDriver::compute(const CMolecule&       molecule,
                                         const CMolecularBasis& basis) const
{
    CAngularMomentumMatrix dipmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs container
        
        CGtoContainer bracontr(molecule, basis);
        
        // compute angular momentum integrals
        
        dipmat = _compAngularMomentumIntegrals(&bracontr, &bracontr);
    }
    
    return dipmat;
}

CAngularMomentumMatrix
CAngularMomentumIntegralsDriver::compute(const CMolecule&       molecule,
                                         const CMolecularBasis& braBasis,
                                         const CMolecularBasis& ketBasis) const
{
    CAngularMomentumMatrix dipmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs containers
        
        CGtoContainer bracontr(molecule, braBasis);
        
        CGtoContainer ketcontr(molecule, ketBasis);
        
        // compute angular momentum integrals
        
        dipmat = _compAngularMomentumIntegrals(&bracontr, &ketcontr);
    }
    
    return dipmat;
}

CAngularMomentumMatrix
CAngularMomentumIntegralsDriver::compute(const CMolecule&       braMolecule,
                                         const CMolecule&       ketMolecule,
                                         const CMolecularBasis& basis) const
{
    CAngularMomentumMatrix dipmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs containers
        
        CGtoContainer bracontr(braMolecule, basis);
        
        CGtoContainer ketcontr(ketMolecule, basis);
        
        // compute angular momentum integrals
        
        dipmat = _compAngularMomentumIntegrals(&bracontr, &ketcontr);
    }
    
    return dipmat;
}

CAngularMomentumMatrix
CAngularMomentumIntegralsDriver::compute(const CMolecule&       braMolecule,
                                         const CMolecule&       ketMolecule,
                                         const CMolecularBasis& braBasis,
                                         const CMolecularBasis& ketBasis) const
{
    CAngularMomentumMatrix dipmat;
    
    if (_locRank == mpi::master())
    {
        // set up GTOs containers
        
        CGtoContainer bracontr(braMolecule, braBasis);
        
        CGtoContainer ketcontr(ketMolecule, ketBasis);
        
        // compute angular momentum integrals
        
        dipmat = _compAngularMomentumIntegrals(&bracontr, &ketcontr);
    }
    
    return dipmat;
}

void
CAngularMomentumIntegralsDriver::compute(   double*    intsBatchX,
                                            double*    intsBatchY,
                                            double*    intsBatchZ,
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
    
    // compute angular momentum integrals
    
    _compAngularMomentumForGtoBlocks(&distx, &disty, &distz, _xOrigin, _yOrigin,
                                    _zOrigin, braGtoBlock, ketGtoBlock);
}

CAngularMomentumMatrix
CAngularMomentumIntegralsDriver::_compAngularMomentumIntegrals(const CGtoContainer* braGtoContainer,
                                                               const CGtoContainer* ketGtoContainer) const
{
    // check if GTOs containers are same on bra and ket sides
    
    auto symbk = ((*braGtoContainer) == (*ketGtoContainer));
    
    // determine dimensions of overlap matrix
    
    auto nrow = braGtoContainer->getNumberOfAtomicOrbitals();
    
    auto ncol = ketGtoContainer->getNumberOfAtomicOrbitals();
    
    // allocate dense matrix for angular momentum integrals
    
    CDenseMatrix dipxmat(nrow, ncol);
    
    CDenseMatrix dipymat(nrow, ncol);
    
    CDenseMatrix dipzmat(nrow, ncol);
    
    // set up distributio pattern
    
    dist1e dstyp = (symbk) ? dist1e::antisq : dist1e::rect;
    
    COneIntsDistribution* distpatx = new COneIntsDistribution(dipxmat.values(),
                                                              nrow, ncol, dstyp);
    
    COneIntsDistribution* distpaty = new COneIntsDistribution(dipymat.values(),
                                                              nrow, ncol, dstyp);
    
    COneIntsDistribution* distpatz = new COneIntsDistribution(dipzmat.values(),
                                                              nrow, ncol, dstyp);
    
    // copy origin coordinates
    
    auto origx = _xOrigin;
    
    auto origy = _yOrigin;
    
    auto origz = _zOrigin;
    
    // compute angular momentum integral blocks
    
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
                        
                        _compAngularMomentumForGtoBlocks(distpatx, distpaty,
                                                         distpatz, origx, origy,
                                                         origz, bgtos, kgtos);
                    }
                }
            }
        }
    }
    
    // deallocate distribution pattern
    
    delete distpatx;
    
    delete distpaty;
    
    delete distpatz;
    
    return CAngularMomentumMatrix(dipxmat, dipymat, dipzmat, _xOrigin, _yOrigin,
                                 _zOrigin);
}

void
CAngularMomentumIntegralsDriver::_compAngularMomentumForGtoBlocks(      COneIntsDistribution* distPatternX,
                                                                        COneIntsDistribution* distPatternY,
                                                                        COneIntsDistribution* distPatternZ,
                                                                  const double                xOrigin,
                                                                  const double                yOrigin,
                                                                  const double                zOrigin,
                                                                  const CGtoBlock&            braGtoBlock,
                                                                  const CGtoBlock&            ketGtoBlock) const
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
    
    CMemBlock2D<double> rfacts(pdim, 4 * pmax);
    
    // allocate P center coordinates
    
    CMemBlock2D<double> rp(pdim, 3 * pmax);
    
    // set up PA and PB distances
    
    auto rpa = CMemBlock2D<double>(pdim, 3 * pmax);
    
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
    
    auto pidx = recmap.getIndexOfTerm(CRecursionTerm({"Angular Momentum"}, 1, true,
                                                     {bang, -1, -1, -1}, {kang, -1, -1, -1}, 1, 1, 0));
    
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
        
        intsfunc::compFactorsForLinearMomentum(rfacts, bragtos, ketgtos, i);
        
        // compute coordinates of center P
        
        intsfunc::compCoordinatesForP(rp, rfacts, 4, bragtos, ketgtos, i);
        
        // compute distances: R(PA) = P - A
        
        intsfunc::compDistancesPA(rpa, rp, bragtos, ketgtos, i);
        
        // compute distances: R(PB) = P - B
        
        intsfunc::compDistancesPB(rpb, rp, bragtos, ketgtos, i);
        
        // compute distances: R(PC) = P - C
        
        intsfunc::compDistancesPC(rpc, rp, xOrigin, yOrigin, zOrigin,
                                  bragtos, ketgtos, i);
        
        // compute primitive angular momentum integrals
        
        _compPrimAngularMomentumInts(primbuffer, recmap, rfacts, rab, rpa, rpb, rpc,
                                     bragtos, ketgtos, i);
        
        // contract primitive linear momentum integrals
        
        genfunc::contract(cartbufferx, primbuffer, pidx, bragtos, ketgtos, i);
        
        genfunc::contract(cartbuffery, primbuffer, pidx + poff, bragtos,
                          ketgtos, i);
        
        genfunc::contract(cartbufferz, primbuffer, pidx + 2 * poff, bragtos,
                          ketgtos, i);
        
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
CAngularMomentumIntegralsDriver::_compPrimAngularMomentumInts(      CMemBlock2D<double>&  primBuffer,
                                                              const CRecursionMap&        recursionMap,
                                                              const CMemBlock2D<double>&  osFactors,
                                                              const CMemBlock2D<double>&  abDistances,
                                                              const CMemBlock2D<double>&  paDistances,
                                                              const CMemBlock2D<double>&  pbDistances,
                                                              const CMemBlock2D<double>&  pcDistances,
                                                              const CGtoBlock&            braGtoBlock,
                                                              const CGtoBlock&            ketGtoBlock,
                                                              const int32_t               iContrGto) const
{
    ovlrecfunc::compOverlapForSS(primBuffer, recursionMap, osFactors, 4, abDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    amomrecfunc::compAngularMomentumForSS(primBuffer, recursionMap, osFactors, paDistances, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    ediprecfunc::compElectricDipoleForSS(primBuffer, recursionMap, pcDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    lmomrecfunc::compLinearMomentumForSS(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    amomrecfunc::compAngularMomentumForPS(primBuffer, recursionMap, osFactors, paDistances, braGtoBlock, ketGtoBlock, iContrGto);
    
    // NOTE: add l > 4 recursion here
}

CRecursionMap
CAngularMomentumIntegralsDriver::_setRecursionMap(const int32_t braAngularMomentum,
                                                 const int32_t ketAngularMomentum,
                                                 const int32_t maxNumberOfPrimitives) const
{
    CRecursionFunctionsList recfuncs;
    
    recfuncs.add(CRecursionFunction({"Overlap"}, &t2crecfunc::obRecursionForOverlap));
    
    recfuncs.add(CRecursionFunction({"Electric Dipole"}, &t2crecfunc::obRecursionForElectricDipole));
    
    recfuncs.add(CRecursionFunction({"Linear Momentum"}, &t2crecfunc::obRecursionForLinearMomentum));
    
    recfuncs.add(CRecursionFunction({"Angular Momentum"}, &t2crecfunc::obRecursionForAngularMomentum));
    
    auto rterm = gintsfunc::genIntegral({"Angular Momentum"}, braAngularMomentum,
                                        ketAngularMomentum, 0);
    
    return gintsfunc::genRecursionMap(rterm, recblock::cc, maxNumberOfPrimitives,
                                      recfuncs);
}
