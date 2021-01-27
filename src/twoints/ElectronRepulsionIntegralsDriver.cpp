//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "ElectronRepulsionIntegralsDriver.hpp"

#include <sstream>

#include <mpi.h>

#include "GtoPairsContainer.hpp"
#include "SphericalMomentum.hpp"
#include "GenFunc.hpp"
#include "TwoIntsFunc.hpp"
#include "AngularMomentum.hpp"
#include "MpiFunc.hpp"
#include "StringFormat.hpp"
#include "TwoCentersRecursionFunctions.hpp"
#include "GenIntsFunc.hpp"

#include "ElectronRepulsionRecFuncForSX.hpp"
#include "ElectronRepulsionRecFuncForPX.hpp"
#include "ElectronRepulsionRecFuncForDX.hpp"
#include "ElectronRepulsionRecFuncForFF.hpp"
#include "ElectronRepulsionRecFuncForFG.hpp"
#include "ElectronRepulsionRecFuncForGF.hpp"
#include "ElectronRepulsionRecFuncForFH.hpp"
#include "ElectronRepulsionRecFuncForHF.hpp"
#include "ElectronRepulsionRecFuncForFI.hpp"
#include "ElectronRepulsionRecFuncForIF.hpp"
#include "ElectronRepulsionRecFuncForFK.hpp"
#include "ElectronRepulsionRecFuncForKF.hpp"
#include "ElectronRepulsionRecFuncForFL.hpp"
#include "ElectronRepulsionRecFuncForLF.hpp"
#include "ElectronRepulsionRecFuncForGG.hpp"
#include "ElectronRepulsionRecFuncForGH.hpp"
#include "ElectronRepulsionRecFuncForHG.hpp"
#include "ElectronRepulsionRecFuncForGI.hpp"
#include "ElectronRepulsionRecFuncForIG.hpp"
#include "ElectronRepulsionRecFuncForGK.hpp"
#include "ElectronRepulsionRecFuncForKG.hpp"
#include "ElectronRepulsionRecFuncForGL.hpp"
#include "ElectronRepulsionRecFuncForLG.hpp"
#include "ElectronRepulsionRecFuncForHH.hpp"
#include "ElectronRepulsionRecFuncForHI.hpp"
#include "ElectronRepulsionRecFuncForIH.hpp"
#include "ElectronRepulsionRecFuncForHK.hpp"
#include "ElectronRepulsionRecFuncForKH.hpp"
#include "ElectronRepulsionRecFuncForHL.hpp"
#include "ElectronRepulsionRecFuncForLH.hpp"
#include "ElectronRepulsionRecFuncForII.hpp"
#include "ElectronRepulsionRecFuncForIK.hpp"
#include "ElectronRepulsionRecFuncForKI.hpp"
#include "ElectronRepulsionRecFuncForIL.hpp"
#include "ElectronRepulsionRecFuncForLI.hpp"
#include "ElectronRepulsionRecFuncForKK.hpp"
#include "ElectronRepulsionRecFuncForKL.hpp"
#include "ElectronRepulsionRecFuncForLK.hpp"
#include "ElectronRepulsionRecFuncForLL.hpp"

#include "ElectronRepulsionKRRRecFuncForSXPY.hpp"
#include "ElectronRepulsionKRRRecFuncForSXDY.hpp"
#include "ElectronRepulsionKRRRecFuncForSXFY.hpp"
#include "ElectronRepulsionKRRRecFuncForSXGG.hpp"

#include "ElectronRepulsionBRRRecFuncForPXYY.hpp"
#include "ElectronRepulsionBRRRecFuncForDXYY.hpp"
#include "ElectronRepulsionBRRRecFuncForFXYY.hpp"
#include "ElectronRepulsionBRRRecFuncForGGYY.hpp"

CElectronRepulsionIntegralsDriver::CElectronRepulsionIntegralsDriver(MPI_Comm comm)
{
    _locRank  = mpi::rank(comm);
    
    _locNodes = mpi::nodes(comm);

    _locComm = comm;
}

CElectronRepulsionIntegralsDriver::~CElectronRepulsionIntegralsDriver()
{
}

void
CElectronRepulsionIntegralsDriver::compute(      CAOFockMatrix&       aoFockMatrix,
                                           const CAODensityMatrix&    aoDensityMatrix,
                                           const CMolecule&           molecule,
                                           const CMolecularBasis&     aoBasis,
                                           const CScreeningContainer& screeningContainer) const
{
    // generate GTOs pairs blocks for AO basis on bra side
    
    CGtoPairsContainer bgtopairs(molecule, aoBasis, 1.0e-15);
    
    // compute repulsion integrals
    
    aoFockMatrix.zero(); 
    
    // compute integrals on CPUs
    
    auto bbpairs = bgtopairs.split(_locNodes);
    
    _compElectronRepulsionIntegrals(aoFockMatrix, aoDensityMatrix, &bbpairs,
                                    &bbpairs, &screeningContainer);
    
    aoFockMatrix.symmetrize(); 
}

CScreeningContainer
CElectronRepulsionIntegralsDriver::compute(const ericut           screeningScheme,
                                           const double           threshold,
                                           const CMolecule&       molecule,
                                           const CMolecularBasis& aoBasis) const
{
    // generate GTOs pairs blocks for AO basis on bra side
    
    CGtoPairsContainer bgtopairs(molecule, aoBasis, 1.0e-15);
    
    // split GTOs pairs into batches on bra side
    
    auto bbpairs = bgtopairs.split(_locNodes);
    
    // allocate temporary buffer for Q values on bra side
    
    auto bqbuff = _getQValuesBuffer(bbpairs);
    
    CVecMemBlock<double> kqbuff;
    
    // compute Q values on bra side
    
    computeMaxQValues(&bqbuff, &kqbuff, &bbpairs, &bbpairs);
    
    // copy Q values from bra to ket side
    
    kqbuff = bqbuff;
    
    // initialize screening container
    
    CScreeningContainer qcont(bqbuff, kqbuff, bbpairs, bbpairs, screeningScheme,
                              threshold);
    
    return qcont;
}

void
CElectronRepulsionIntegralsDriver::compute(      double*         intsBatch,
                                           const CGtoPairsBlock& braGtoPairsBlock,
                                           const CGtoPairsBlock& ketGtoPairsBlock) const
{
    // set up dimensions of integrals batch
    
    auto nrow = braGtoPairsBlock.getNumberOfScreenedContrPairs();
    
    auto ncol = ketGtoPairsBlock.getNumberOfScreenedContrPairs();
    
    // initialize two electron distributor
    
    CTwoIntsDistribution distpat(intsBatch, nrow, ncol, dist2e::batch);
    
    // set up empty screener
    
    CCauchySchwarzScreener qqdat;
    
    // compute batch of two electron integrals
    
    _compElectronRepulsionForGtoPairsBlocks(distpat, qqdat, braGtoPairsBlock,
                                            ketGtoPairsBlock);
}

void
CElectronRepulsionIntegralsDriver::computeMaxQValues(      CVecMemBlock<double>* braQValuesBuffer,
                                                           CVecMemBlock<double>* ketQValuesBuffer,
                                                     const CGtoPairsContainer*   braGtoPairsContainer,
                                                     const CGtoPairsContainer*   ketGtoPairsContainer) const
{
    // determine symmetry of GTOs pairs containers on bra and ket sides
    
    auto symbk = (*braGtoPairsContainer == *ketGtoPairsContainer);
    
    #pragma omp parallel shared(braGtoPairsContainer, ketGtoPairsContainer)
    {
        #pragma omp single nowait
        {
            // Q values for bra side
            
            auto nbra = braGtoPairsContainer->getNumberOfGtoPairsBlocks();
            
            for (int32_t i = 0; i < nbra; i++)
            {
                #pragma omp task firstprivate(i)
                {
                    auto bqprt = (*braQValuesBuffer)[i].data();
                    
                    auto bpairs = braGtoPairsContainer->getGtoPairsBlock(i);
                    
                    _compMaxQValuesForGtoPairsBlock(bqprt, bpairs);
                }
            }
            
            // Q values for ket side if needed
            
            if (!symbk)
            {
                auto nket = ketGtoPairsContainer->getNumberOfGtoPairsBlocks();
                
                for (int32_t i = 0; i < nket; i++)
                {
                    #pragma omp task firstprivate(i)
                    {
                        auto kqprt = (*ketQValuesBuffer)[i].data();
                        
                        auto kpairs = ketGtoPairsContainer->getGtoPairsBlock(i);
                        
                        _compMaxQValuesForGtoPairsBlock(kqprt, kpairs);
                    }
                }
            }
        }
    }
}

void
CElectronRepulsionIntegralsDriver::computeInMemory(const CMolecule&       molecule,
                                                   const CMolecularBasis& aoBasis,
                                                         double*          eriData) const
{
    auto natoms = molecule.getNumberOfAtoms();

    auto max_angl = aoBasis.getMolecularMaxAngularMomentum(molecule);

    int32_t nao = 0;

    for (int32_t angl = 0; angl <= max_angl; angl++)
    {
        for (int32_t s = -angl; s <= angl; s++)
        {
            for (int32_t atomidx = 0; atomidx < natoms; atomidx++)
            {
                int32_t idelem = molecule.getIdsElemental()[atomidx];

                nao += aoBasis.getNumberOfBasisFunctions(idelem, angl);
            }
        }
    }

    int64_t nao2 = nao * nao;

    CGtoPairsContainer bgtopairs(molecule, aoBasis, 1.0e-15);

    CGtoPairsContainer kgtopairs(molecule, aoBasis, 1.0e-15);

    auto nbra = bgtopairs.getNumberOfGtoPairsBlocks();

    auto nket = kgtopairs.getNumberOfGtoPairsBlocks();

    // loop over pairs of GTOs blocks

    for (int32_t ibra = (nbra - 1); ibra >= 0; ibra--)
    {
        auto bpairs = bgtopairs.getGtoPairsBlock(ibra);

        auto nrow = bpairs.getNumberOfScreenedContrPairs();

        auto acomp = angmom::to_SphericalComponents(bpairs.getBraAngularMomentum());

        auto bcomp = angmom::to_SphericalComponents(bpairs.getKetAngularMomentum());

        auto ketoff = ibra; // make use of symmetry (bgtopairs == kgtopairs)

        for (int32_t iket = (nket - 1); iket >= ketoff; iket--)
        {
            auto kpairs = kgtopairs.getGtoPairsBlock(iket);

            bool symbk = (bpairs == kpairs);

            auto ncol = kpairs.getNumberOfScreenedContrPairs();

            auto ccomp = angmom::to_SphericalComponents(kpairs.getBraAngularMomentum());

            auto dcomp = angmom::to_SphericalComponents(kpairs.getKetAngularMomentum());

            auto ncomp = acomp * bcomp * ccomp * dcomp;

            int32_t npairs = 0;

            for (int32_t irow = 0; irow < nrow; irow++)
            {
                int32_t numcols = symbk ? (irow + 1) : ncol;

                npairs += numcols;
            }

            CMemBlock<double> intsBatch(npairs * ncomp);

            CTwoIntsDistribution distpat(intsBatch.data(), nrow, ncol, dist2e::batch);

            CCauchySchwarzScreener qqdat;

            _compElectronRepulsionForGtoPairsBlocks(distpat, qqdat, bpairs, kpairs);

            for (int32_t irow = 0, intsId = 0; irow < nrow; irow++)
            {
                int32_t numcols = symbk ? (irow + 1) : ncol;

                for (int32_t icol = 0; icol < numcols; icol++)
                {
                    for (int32_t i = 0; i < acomp; i++)
                    {
                        int32_t idp = (bpairs.getBraIdentifiers(i))[irow];

                        for (int32_t j = 0; j < bcomp; j++)
                        {
                            int32_t idq = (bpairs.getKetIdentifiers(j))[irow];

                            int64_t pq = idp * nao + idq;

                            int64_t qp = idq * nao + idp;

                            for (int32_t k = 0; k < ccomp; k++)
                            {
                                int32_t idr = (kpairs.getBraIdentifiers(k))[icol];

                                for (int32_t l = 0; l < dcomp; l++, intsId++)
                                {
                                    int32_t ids = (kpairs.getKetIdentifiers(l))[icol];

                                    int64_t rs = idr * nao + ids;

                                    int64_t sr = ids * nao + idr;

                                    auto fval = intsBatch.data()[intsId];

                                    eriData[pq * nao2 + rs] = fval;

                                    eriData[pq * nao2 + sr] = fval;

                                    eriData[qp * nao2 + rs] = fval;

                                    eriData[qp * nao2 + sr] = fval;

                                    eriData[rs * nao2 + pq] = fval;

                                    eriData[rs * nao2 + qp] = fval;

                                    eriData[sr * nao2 + pq] = fval;

                                    eriData[sr * nao2 + qp] = fval;
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void
CElectronRepulsionIntegralsDriver::_compElectronRepulsionForGtoPairsBlocks(      CTwoIntsDistribution&   distPattern,
                                                                           const CCauchySchwarzScreener& intsScreener,
                                                                           const CGtoPairsBlock&         braGtoPairsBlock,
                                                                           const CGtoPairsBlock&         ketGtoPairsBlock) const
{
    // copy GTOs pairs blocks for bra and ket sides
    
    auto brapairs = braGtoPairsBlock;
    
    auto ketpairs = ketGtoPairsBlock;
    
    // determine symmetry of bra and ket sides
    
    bool symbk = (brapairs == ketpairs);
    
    // set up angular momentum for four centers
    
    auto anga = brapairs.getBraAngularMomentum();
    
    auto angb = brapairs.getKetAngularMomentum();
    
    auto angc = ketpairs.getBraAngularMomentum();
    
    auto angd = ketpairs.getKetAngularMomentum();
    
    // set up spherical angular momentum for bra and ket sides
    
    CSphericalMomentum amom(anga);
    
    CSphericalMomentum bmom(angb);
    
    CSphericalMomentum cmom(angc);
    
    CSphericalMomentum dmom(angd);
    
    // allocate prefactors used in Obara-Saika recursion
    
    auto pdim = ketpairs.getNumberOfScreenedPrimPairs();
    
    auto pmax = brapairs.getMaxContractionDepth();
    
    CMemBlock2D<double> rpq(pdim, 3 * pmax);
    
    CMemBlock2D<double> rfacts(pdim, 4 * pmax);
    
    CMemBlock2D<double> rw(pdim, 3 * pmax);
    
    CMemBlock2D<double> rwp(pdim, 3 * pmax);
    
    CMemBlock2D<double> rwq(pdim, 3 * pmax);
    
    // generate horizontal recursion patterna for bra and ket
    
    auto bhrrmap = _setBraHorizontalRecursionPattern(anga, angb, angc, angd);
    
    auto khrrmap = _setKetHorizontalRecursionPattern(anga, angb, angc, angd);
    
    // generate vertical recursion map
    
    auto vrrmap = _setVerticalRecursionMap(khrrmap, pmax);

    auto nblk = vrrmap.getNumberOfComponents();
    
    // allocate primitives integrals buffer
    
    auto pbuffer = vrrmap.createBuffer(pdim);
    
    // set up horizontal recursion buffer for ket side
    
    auto cdim = ketpairs.getNumberOfScreenedContrPairs();
    
    nblk = khrrmap.getNumberOfComponents();
    
    CMemBlock2D<double> khrrbuffer(cdim, nblk);
    
    // initialize R(CD) = C - D distance for horizontal recursion
    
    CMemBlock2D<double> rcd(cdim, 3);
    
    // set up horizontal recursion buffer for bra side
    
    nblk = bhrrmap.getNumberOfComponents();
    
    CMemBlock2D<double> bhrrbuffer(cdim, nblk);
    
    // initialize R(AB) = A - B distance for horizontal recursion
    
    auto rab = brapairs.getDistancesAB();
    
    // allocate spherical integrals buffer
    
    nblk = angmom::to_SphericalComponents(anga, angb)
         * angmom::to_SphericalComponents(angc, angd);
    
    CMemBlock2D<double> spherbuffer(cdim, nblk);
    
    // initialize Boys function evaluator
    
    auto bord = anga + angb + angc + angd;
    
    CBoysFunction bftab(bord);
    
    CMemBlock<double> bargs(pdim);
    
    CMemBlock2D<double> bvals(pdim, bord + 1);
    
    // set up integrals screening
    
    bool useqq = !intsScreener.isEmpty();
    
    auto qqpairs = ketpairs;
    
    auto ddpairs = ketpairs;
    
    CMemBlock<int32_t> qqvec(cdim);
    
    CMemBlock<int32_t> qqidx(cdim);
    
    CMemBlock<double> distpq(cdim);
    
    CMemBlock<double> qqden(cdim);
    
    // loop over contracted GTOs ob bra side
    
    for (int32_t i = 0; i < brapairs.getNumberOfScreenedContrPairs(); i++)
    {
        // determine GTOs pairs  effective dimensions on ket side
        
        auto nqpdim = (symbk) ? ketpairs.getNumberOfPrimPairs(i) : pdim;
        
        auto nqcdim = (symbk) ? i + 1 : cdim;
        
        // integrals screening: QQ or QQR scheme 
        
        if (useqq)
        {
            // compute effective distances between GTOs pairs on bra and ket sides
            
            if (intsScreener.getScreeningScheme() == ericut::qqr)
            {
                twointsfunc::compEffectiveDistancesPQ(distpq, brapairs, ketpairs,
                                                      symbk, i);
            }
            
            intsScreener.setScreeningVector(qqvec, distpq, symbk, i);
            
            mathfunc::ordering(qqidx.data(), qqvec.data(), nqcdim);
            
            nqcdim = qqpairs.compress(ketpairs, qqvec, nqcdim);
            
            if (nqcdim > 0) nqpdim = qqpairs.getNumberOfPrimPairs(nqcdim - 1);
        }
        
        // all integrals are vanishing in batch, skip computations
        
        if (nqcdim == 0) continue;
        
        // density based screeing of integrals batches
        
        if (useqq)
        {
            // determine max. density element for GTO pairs on ket side
            
            if ((intsScreener.getScreeningScheme() == ericut::qqden) ||
                (intsScreener.getScreeningScheme() == ericut::qqrden))
            {
                distPattern.getMaxDensityElements(qqden, brapairs, qqpairs,
                                                  symbk, nqcdim, i);
                
                intsScreener.setScreeningVector(qqvec, qqidx, qqden, distpq,
                                                nqcdim, i);
                
                nqcdim = ddpairs.compress(qqpairs, qqvec, nqcdim);
                
                if (nqcdim > 0) nqpdim = ddpairs.getNumberOfPrimPairs(nqcdim - 1);
            }
        }
        
        // all integrals are vanishing in batch, skip computations
        
        if (nqcdim == 0) continue;
        
        // density screened scheme
        
        if ((intsScreener.getScreeningScheme() == ericut::qqden) ||
            (intsScreener.getScreeningScheme() == ericut::qqrden))
        {
            // compute distances: R(PQ) = P - Q
            
            twointsfunc::compDistancesPQ(rpq, brapairs, ddpairs, nqpdim, i);
            
            // compute Obara-Saika recursion factors
            
            twointsfunc::compFactorsForElectronRepulsion(rfacts, brapairs,
                                                         ddpairs, nqpdim, i);
            
            // compute coordinates of center W
            
            twointsfunc::compCoordinatesForW(rw, rfacts, 4, brapairs, ddpairs,
                                             nqpdim, i);
            
            // compute distances: R(WP) = W - P
            
            twointsfunc::compDistancesWP(rwp, rw, brapairs, ddpairs, nqpdim, i);
            
            // compute distances: R(WQ) = W - Q;
            
            twointsfunc::compDistancesWQ(rwq, rw, brapairs, ddpairs, nqpdim, i);
            
            // compute primitive electron repulsion integrals
            
            _compPrimElectronRepulsionInts(pbuffer, vrrmap, bftab,
                                           bargs, bvals, bord, rfacts, rpq, rwp,
                                           rwq, brapairs, ddpairs, nqpdim, i);
            
            // contract primitive electron repulsion integrals
            
            genfunc::contract(khrrbuffer, pbuffer, khrrmap, vrrmap, brapairs, ddpairs, nqpdim, nqcdim, i);
            
            // apply horizontal recursion on ket side
            
            ddpairs.getDistancesAB(rcd, nqcdim);
            
            _applyHRRonKet(khrrbuffer, khrrmap, rcd, brapairs, ddpairs, nqcdim, i);
            
            // transform ket side to spherical form
            
            genfunc::transform_ket(bhrrbuffer, khrrbuffer, cmom, dmom, bhrrmap,
                                   khrrmap, ddpairs, nqcdim, i);
            
            // apply horizontal recursion on bra side
            
            _applyHRRonBra(bhrrbuffer, bhrrmap, rab, brapairs, ddpairs, nqcdim, i);
            
            // transform bra side to spherical form
            
            genfunc::transform_bra(spherbuffer, bhrrbuffer, amom, bmom, bhrrmap,
                                   ddpairs, nqcdim, i);
            
            // distribute integrals: add distribution or Fock formation code
            
            distPattern.distribute(spherbuffer, brapairs, ddpairs, symbk,
                                   nqcdim, i);
        }
        else
        {
            // compute distances: R(PQ) = P - Q
        
            twointsfunc::compDistancesPQ(rpq, brapairs, qqpairs, nqpdim, i);
        
            // compute Obara-Saika recursion factors
        
            twointsfunc::compFactorsForElectronRepulsion(rfacts, brapairs,
                                                         qqpairs, nqpdim, i);
        
            // compute coordinates of center W
        
            twointsfunc::compCoordinatesForW(rw, rfacts, 4, brapairs, qqpairs,
                                             nqpdim, i);
        
            // compute distances: R(WP) = W - P
        
            twointsfunc::compDistancesWP(rwp, rw, brapairs, qqpairs, nqpdim, i);
        
            // compute distances: R(WQ) = W - Q;
        
            twointsfunc::compDistancesWQ(rwq, rw, brapairs, qqpairs, nqpdim, i);
        
            // compute primitive electron repulsion integrals
        
            _compPrimElectronRepulsionInts(pbuffer, vrrmap, bftab,
                                           bargs, bvals, bord, rfacts, rpq, rwp,
                                           rwq, brapairs, qqpairs, nqpdim, i);
        
            // contract primitive electron repulsion integrals
        
            genfunc::contract(khrrbuffer, pbuffer, khrrmap, vrrmap, brapairs, qqpairs, nqpdim, nqcdim, i);
        
            // apply horizontal recursion on ket side
        
            qqpairs.getDistancesAB(rcd, nqcdim);
        
            _applyHRRonKet(khrrbuffer, khrrmap, rcd, brapairs, qqpairs, nqcdim, i);
        
            // transform ket side to spherical form
        
            genfunc::transform_ket(bhrrbuffer, khrrbuffer, cmom, dmom, bhrrmap,
                                   khrrmap, qqpairs, nqcdim, i);
        
            // apply horizontal recursion on bra side
        
            _applyHRRonBra(bhrrbuffer, bhrrmap, rab, brapairs, qqpairs, nqcdim, i);
        
            // transform bra side to spherical form
        
            genfunc::transform_bra(spherbuffer, bhrrbuffer, amom, bmom, bhrrmap,
                                   qqpairs, nqcdim, i);
        
            // distribute integrals: add distribution or Fock formation code
        
            distPattern.distribute(spherbuffer, brapairs, qqpairs, symbk,
                                   nqcdim, i);
        }
    }
    
    // deallocate recursion buffers
    
    vrrmap.destroyBuffer(pbuffer);
}

CRecursionMap
CElectronRepulsionIntegralsDriver::_setBraHorizontalRecursionPattern(const int32_t angularMomentumA,
                                                                     const int32_t angularMomentumB,
                                                                     const int32_t angularMomentumC,
                                                                     const int32_t angularMomentumD) const
{
    CRecursionMap recmap(recblock::cs, 1);
    
    for (int32_t i = 0; i <= angularMomentumA; i++)
    {
        for (int32_t j = angularMomentumB; j <= (angularMomentumB + angularMomentumA - i); j++)
        {
            recmap.add(gintsfunc::genElectronRepulsionIntegral(i, j, angularMomentumC, angularMomentumD));
        }
    }
    
    return recmap;
}

CRecursionMap
CElectronRepulsionIntegralsDriver::_setKetHorizontalRecursionPattern(const int32_t angularMomentumA,
                                                                     const int32_t angularMomentumB,
                                                                     const int32_t angularMomentumC,
                                                                     const int32_t angularMomentumD) const
{
    CRecursionMap recmap(recblock::cc, 1);
    
    for (int32_t i = 0; i <= angularMomentumA + angularMomentumB; i++)
    {
        for (int32_t j = 0; j <= angularMomentumC; j++)
        {
            for (int32_t k = angularMomentumD; k <= (angularMomentumC + angularMomentumD - j); k++)
            {
                recmap.add(gintsfunc::genElectronRepulsionIntegral(i, j, k));
            }
        }
    }
    
    return recmap;
}

CRecursionMap
CElectronRepulsionIntegralsDriver::_setVerticalRecursionMap(const CRecursionMap& leadTerms,
                                                            const int32_t        maxNumberOfPrimPairs) const
{
    CRecursionFunctionsList recfuncs;
    
    recfuncs.add(CRecursionFunction({"Electron Repulsion"}, &t2crecfunc::obRecursionForElectronRepulsion));
    
    CRecursionMap recmap(recblock::cc, maxNumberOfPrimPairs);
    
    for (int32_t i = 0; i < leadTerms.getNumberOfTerms(); i++)
    {
        auto rterm = leadTerms.getTerm(i);
        
        if (rterm.getKetAngularMomentum(0) == 0)
        {
            recmap.append(gintsfunc::genRecursionMap(CRecursionTerm(std::string("Electron Repulsion"), 0, true,
                                                                    {rterm.getBraAngularMomentum(0), -1, -1, -1},
                                                                    {rterm.getKetAngularMomentum(1), -1, -1, -1},
                                                                    1, 1, 0),
                                                     recblock::cc, maxNumberOfPrimPairs, recfuncs));
        }
    }
    
    return recmap;
}

int32_t
CElectronRepulsionIntegralsDriver::_getBufferDimensionsForVerticalRecursion(const int32_t braAngularMomentum,
                                                                            const int32_t ketAngularMomentum,
                                                                            const int32_t maxNumberOfPrimPairs) const
{
    int32_t vrrdim = braAngularMomentum +  ketAngularMomentum + 1;
    
    if ((braAngularMomentum ==  0) && (ketAngularMomentum == 1)) vrrdim = 5;
    
    if ((braAngularMomentum ==  1) && (ketAngularMomentum == 0)) vrrdim = 5;
    
    // FIX ME: Add other angular momentum terms
    
    return maxNumberOfPrimPairs * vrrdim;
}

void
CElectronRepulsionIntegralsDriver::_compPrimElectronRepulsionInts(      CMemBlock2D<double>*  primBuffer,
                                                                  const CRecursionMap&        recursionMap, 
                                                                  const CBoysFunction&        bfTable,
                                                                        CMemBlock<double>&    bfArguments,
                                                                        CMemBlock2D<double>&  bfValues,
                                                                  const int32_t               bfOrder,
                                                                  const CMemBlock2D<double>&  osFactors,
                                                                  const CMemBlock2D<double>&  pqDistances,
                                                                  const CMemBlock2D<double>&  wpDistances,
                                                                  const CMemBlock2D<double>&  wqDistances,
                                                                  const CGtoPairsBlock&       braGtoPairsBlock,
                                                                  const CGtoPairsBlock&       ketGtoPairsBlock,
                                                                  const int32_t               nKetPrimPairs,
                                                                  const int32_t               iContrPair) const
{
    // set up angular momentum
    
    auto bang = braGtoPairsBlock.getBraAngularMomentum() + braGtoPairsBlock.getKetAngularMomentum();
    
    auto kang = ketGtoPairsBlock.getBraAngularMomentum() + ketGtoPairsBlock.getKetAngularMomentum();
 
    // primitive integrals
    
    erirecfunc::compElectronRepulsionForSSSS(primBuffer, recursionMap, bfTable, bfArguments, bfValues, bfOrder, osFactors, pqDistances,
                                             braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 0) && (kang == 0)) return;
    
    erirecfunc::compElectronRepulsionForSSSP(primBuffer, recursionMap, wqDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 0) && (kang == 1)) return;
    
    erirecfunc::compElectronRepulsionForSPSS(primBuffer, recursionMap, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 1) && (kang == 0)) return;
    
    erirecfunc::compElectronRepulsionForSPSP(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 1) && (kang == 1)) return;
    
    erirecfunc::compElectronRepulsionForSDSS(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 2) && (kang == 0)) return;
    
    erirecfunc::compElectronRepulsionForSDSP(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 2) && (kang == 1)) return;
    
    erirecfunc::compElectronRepulsionForSFSS(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 3) && (kang == 0)) return;
    
    erirecfunc::compElectronRepulsionForSFSP(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 3) && (kang == 1)) return;
    
    erirecfunc::compElectronRepulsionForSGSS(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 4) && (kang == 0)) return;
    
    erirecfunc::compElectronRepulsionForSGSP(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 4) && (kang == 1)) return;
    
    erirecfunc::compElectronRepulsionForSSSD(primBuffer, recursionMap, osFactors, wqDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 0) && (kang == 2)) return;
    
    erirecfunc::compElectronRepulsionForSPSD(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 1) && (kang == 2)) return;
    
    erirecfunc::compElectronRepulsionForSDSD(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 2) && (kang == 2)) return;
   
    erirecfunc::compElectronRepulsionForSSSF(primBuffer, recursionMap, osFactors, wqDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 0) && (kang == 3)) return;
    
    erirecfunc::compElectronRepulsionForSPSF(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 1) && (kang == 3)) return;
    
    erirecfunc::compElectronRepulsionForSDSF(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 2) && (kang == 3)) return;
    
    erirecfunc::compElectronRepulsionForSFSD(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 3) && (kang == 2)) return;
    
    erirecfunc::compElectronRepulsionForSFSF(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 3) && (kang == 3)) return;
    
    erirecfunc::compElectronRepulsionForSSSG(primBuffer, recursionMap, osFactors, wqDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 0) && (kang == 4)) return;
    
    erirecfunc::compElectronRepulsionForSPSG(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 1) && (kang == 4)) return;
    
    erirecfunc::compElectronRepulsionForSDSG(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 2) && (kang == 4)) return;
    
    erirecfunc::compElectronRepulsionForSGSD(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 4) && (kang == 2)) return;
    
    erirecfunc::compElectronRepulsionForSFSG(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 3) && (kang == 4)) return;
    
    erirecfunc::compElectronRepulsionForSGSF(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 4) && (kang == 3)) return;
    
    erirecfunc::compElectronRepulsionForSGSG(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 4) && (kang == 4)) return;
    
    erirecfunc::compElectronRepulsionForSSSH(primBuffer, recursionMap, osFactors, wqDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 0) && (kang == 5)) return;
    
    erirecfunc::compElectronRepulsionForSPSH(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 1) && (kang == 5)) return;
    
    erirecfunc::compElectronRepulsionForSDSH(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 2) && (kang == 5)) return;
   
    erirecfunc::compElectronRepulsionForSFSH(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 3) && (kang == 5)) return;
    
    erirecfunc::compElectronRepulsionForSGSH(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 4) && (kang == 5)) return;
    
    erirecfunc::compElectronRepulsionForSHSH(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 5) && (kang == 5)) return;
    
    erirecfunc::compElectronRepulsionForSHSS(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 5) && (kang == 0)) return;
    
    erirecfunc::compElectronRepulsionForSHSP(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 5) && (kang == 1)) return;
    
    erirecfunc::compElectronRepulsionForSHSD(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 5) && (kang == 2)) return;
    
    erirecfunc::compElectronRepulsionForSHSF(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 5) && (kang == 3)) return;
    
    erirecfunc::compElectronRepulsionForSHSG(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 5) && (kang == 4)) return;
    
    erirecfunc::compElectronRepulsionForSSSI(primBuffer, recursionMap, osFactors, wqDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 0) && (kang == 6)) return;
    
    erirecfunc::compElectronRepulsionForSPSI(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 1) && (kang == 6)) return;
    
    erirecfunc::compElectronRepulsionForSDSI(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 2) && (kang == 6)) return;
    
    erirecfunc::compElectronRepulsionForSFSI(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 3) && (kang == 6)) return;
    
    erirecfunc::compElectronRepulsionForSGSI(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 4) && (kang == 6)) return;
    
    erirecfunc::compElectronRepulsionForSHSI(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 5) && (kang == 6)) return;
    
    erirecfunc::compElectronRepulsionForSISI(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 6) && (kang == 6)) return;
    
    erirecfunc::compElectronRepulsionForSISS(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 6) && (kang == 0)) return;
    
    erirecfunc::compElectronRepulsionForSISP(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 6) && (kang == 1)) return;
    
    erirecfunc::compElectronRepulsionForSISD(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 6) && (kang == 2)) return;
    
    erirecfunc::compElectronRepulsionForSISF(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 6) && (kang == 3)) return;
    
    erirecfunc::compElectronRepulsionForSISG(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 6) && (kang == 4)) return;
    
    erirecfunc::compElectronRepulsionForSISH(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 6) && (kang == 5)) return;
    
    erirecfunc::compElectronRepulsionForSSSK(primBuffer, recursionMap, osFactors, wqDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 0) && (kang == 7)) return;
    
    erirecfunc::compElectronRepulsionForSPSK(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 1) && (kang == 7)) return;
    
    erirecfunc::compElectronRepulsionForSDSK(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 2) && (kang == 7)) return;
    
    erirecfunc::compElectronRepulsionForSFSK(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 3) && (kang == 7)) return;
    
    erirecfunc::compElectronRepulsionForSGSK(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 4) && (kang == 7)) return;
    
    erirecfunc::compElectronRepulsionForSHSK(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 5) && (kang == 7)) return;
    
    erirecfunc::compElectronRepulsionForSISK(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 6) && (kang == 7)) return;
    
    erirecfunc::compElectronRepulsionForSKSK(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 7) && (kang == 7)) return;
    
    erirecfunc::compElectronRepulsionForSKSS(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 7) && (kang == 0)) return;
    
    erirecfunc::compElectronRepulsionForSKSP(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 7) && (kang == 1)) return;
    
    erirecfunc::compElectronRepulsionForSKSD(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 7) && (kang == 2)) return;
    
    erirecfunc::compElectronRepulsionForSKSF(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 7) && (kang == 3)) return;
    
    erirecfunc::compElectronRepulsionForSKSG(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 7) && (kang == 4)) return;
    
    erirecfunc::compElectronRepulsionForSKSH(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 7) && (kang == 5)) return;
    
    erirecfunc::compElectronRepulsionForSKSI(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 7) && (kang == 6)) return;
    
    erirecfunc::compElectronRepulsionForSSSL(primBuffer, recursionMap, osFactors, wqDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 0) && (kang == 8)) return;
    
    erirecfunc::compElectronRepulsionForSPSL(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 1) && (kang == 8)) return;
    
    erirecfunc::compElectronRepulsionForSDSL(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 2) && (kang == 8)) return;
    
    erirecfunc::compElectronRepulsionForSFSL(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 3) && (kang == 8)) return;
    
    erirecfunc::compElectronRepulsionForSGSL(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 4) && (kang == 8)) return;
    
    erirecfunc::compElectronRepulsionForSHSL(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 5) && (kang == 8)) return;
    
    erirecfunc::compElectronRepulsionForSISL(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 6) && (kang == 8)) return;
    
    erirecfunc::compElectronRepulsionForSKSL(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 7) && (kang == 8)) return;
    
    erirecfunc::compElectronRepulsionForSLSL(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 8) && (kang == 8)) return;
    
    erirecfunc::compElectronRepulsionForSLSS(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 8) && (kang == 0)) return;
    
    erirecfunc::compElectronRepulsionForSLSP(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 8) && (kang == 1)) return;
    
    erirecfunc::compElectronRepulsionForSLSD(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 8) && (kang == 2)) return;
    
    erirecfunc::compElectronRepulsionForSLSF(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 8) && (kang == 3)) return;
    
    erirecfunc::compElectronRepulsionForSLSG(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 8) && (kang == 4)) return;
    
    erirecfunc::compElectronRepulsionForSLSH(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 8) && (kang == 5)) return;
    
    erirecfunc::compElectronRepulsionForSLSI(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
    
    if ((bang == 8) && (kang == 6)) return;
    
    erirecfunc::compElectronRepulsionForSLSK(primBuffer, recursionMap, osFactors, wpDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetPrimPairs, iContrPair);
}

void
CElectronRepulsionIntegralsDriver::_applyHRRonKet(      CMemBlock2D<double>&  ketBuffer,
                                                  const CRecursionMap&        recursionMap,
                                                  const CMemBlock2D<double>&  cdDistances,
                                                  const CGtoPairsBlock&       braGtoPairsBlock,
                                                  const CGtoPairsBlock&       ketGtoPairsBlock,
                                                  const int32_t               nKetContrPairs,
                                                  const int32_t               iContrPair) const
{
    erikrrfunc::compElectronRepulsionForSXPP(ketBuffer, recursionMap, cdDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    erikrrfunc::compElectronRepulsionForSXPD(ketBuffer, recursionMap, cdDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    erikrrfunc::compElectronRepulsionForSXPF(ketBuffer, recursionMap, cdDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    erikrrfunc::compElectronRepulsionForSXPG(ketBuffer, recursionMap, cdDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    erikrrfunc::compElectronRepulsionForSXPH(ketBuffer, recursionMap, cdDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    erikrrfunc::compElectronRepulsionForSXPI(ketBuffer, recursionMap, cdDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    erikrrfunc::compElectronRepulsionForSXPK(ketBuffer, recursionMap, cdDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    erikrrfunc::compElectronRepulsionForSXDD(ketBuffer, recursionMap, cdDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    erikrrfunc::compElectronRepulsionForSXDF(ketBuffer, recursionMap, cdDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    erikrrfunc::compElectronRepulsionForSXDG(ketBuffer, recursionMap, cdDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    erikrrfunc::compElectronRepulsionForSXDH(ketBuffer, recursionMap, cdDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    erikrrfunc::compElectronRepulsionForSXDI(ketBuffer, recursionMap, cdDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    erikrrfunc::compElectronRepulsionForSXFF(ketBuffer, recursionMap, cdDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    erikrrfunc::compElectronRepulsionForSXFG(ketBuffer, recursionMap, cdDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    erikrrfunc::compElectronRepulsionForSXFH(ketBuffer, recursionMap, cdDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    erikrrfunc::compElectronRepulsionForSXGG(ketBuffer, recursionMap, cdDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
}

void
CElectronRepulsionIntegralsDriver::_applyHRRonBra(      CMemBlock2D<double>& braBuffer,
                                                  const CRecursionMap&       recursionMap,
                                                  const CMemBlock2D<double>& abDistances,
                                                  const CGtoPairsBlock&      braGtoPairsBlock,
                                                  const CGtoPairsBlock&      ketGtoPairsBlock,
                                                  const int32_t              nKetContrPairs,
                                                  const int32_t              iContrPair) const
{
    eribrrfunc::compElectronRepulsionForPPXY(braBuffer, recursionMap, abDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    eribrrfunc::compElectronRepulsionForPDXY(braBuffer, recursionMap, abDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    eribrrfunc::compElectronRepulsionForPFXY(braBuffer, recursionMap, abDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    eribrrfunc::compElectronRepulsionForPGXY(braBuffer, recursionMap, abDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    eribrrfunc::compElectronRepulsionForPHXY(braBuffer, recursionMap, abDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    eribrrfunc::compElectronRepulsionForPIXY(braBuffer, recursionMap, abDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    eribrrfunc::compElectronRepulsionForPKXY(braBuffer, recursionMap, abDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    eribrrfunc::compElectronRepulsionForDDXY(braBuffer, recursionMap, abDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    eribrrfunc::compElectronRepulsionForDFXY(braBuffer, recursionMap, abDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    eribrrfunc::compElectronRepulsionForDGXY(braBuffer, recursionMap, abDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    eribrrfunc::compElectronRepulsionForDHXY(braBuffer, recursionMap, abDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    eribrrfunc::compElectronRepulsionForDIXY(braBuffer, recursionMap, abDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    eribrrfunc::compElectronRepulsionForFFXY(braBuffer, recursionMap, abDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    eribrrfunc::compElectronRepulsionForFGXY(braBuffer, recursionMap, abDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    eribrrfunc::compElectronRepulsionForFHXY(braBuffer, recursionMap, abDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
    
    eribrrfunc::compElectronRepulsionForGGXY(braBuffer, recursionMap, abDistances, braGtoPairsBlock, ketGtoPairsBlock, nKetContrPairs, iContrPair);
}

void
CElectronRepulsionIntegralsDriver::_compElectronRepulsionIntegrals(      CAOFockMatrix&       aoFockMatrix,
                                                                   const CAODensityMatrix&    aoDensityMatrix,
                                                                   const CGtoPairsContainer*  braGtoPairsContainer,
                                                                   const CGtoPairsContainer*  ketGtoPairsContainer,
                                                                   const CScreeningContainer* screeningContainer) const
{
    // set up pointers to AO Fock and AO density matrices
    
    auto pfock = &aoFockMatrix;
    
    auto pden = &aoDensityMatrix;
    
    // determine number of GTOs pairs blocks in bra/ket sides
    
    auto nbra = braGtoPairsContainer->getNumberOfGtoPairsBlocks();
    
    auto nket = ketGtoPairsContainer->getNumberOfGtoPairsBlocks();
    
    // determine symmetry of bra/ket GTOs pairs containers
    
    auto symbk = ((*braGtoPairsContainer) == (*ketGtoPairsContainer));
    
    // initialize tasks grid for each MPI process
    
    auto nodpatt = _setTasksGrid(nbra, nket, symbk);
    
    auto ptgrid = nodpatt.data();
    
    #pragma omp parallel shared(braGtoPairsContainer, ketGtoPairsContainer, pfock, pden, nbra, nket, symbk, ptgrid)
    {
        #pragma omp single nowait
        {
            // screeners counter
            
            int32_t idx = 0;
            
            // loop over pairs of GTOs blocks
            
            for (int32_t i = (nbra - 1); i >= 0; i--)
            {
                auto bpairs = braGtoPairsContainer->getGtoPairsBlock(i);
                
                auto joff = (symbk) ? i : 0;
                
                for (int32_t j = (nket - 1); j >= joff; j--)
                {
                    if (ptgrid[idx] == 1)
                    {
                        #pragma omp task firstprivate(j, idx)
                        {
                            auto kpairs = ketGtoPairsContainer->getGtoPairsBlock(j);
                        
                            auto qqdat = screeningContainer->getScreener(idx);
                        
                            CTwoIntsDistribution distpat(pfock, pden);
                        
                            distpat.setFockContainer(bpairs, kpairs);
                        
                            _compElectronRepulsionForGtoPairsBlocks(distpat, qqdat,
                                                                    bpairs, kpairs);
                        
                            // accumulate AO Fock matrix
                        
                            #pragma omp critical (fockacc)
                            distpat.accumulate();
                        }
                    }
                    
                    // update screeners counter
                    
                    idx++; 
                }
            }
        }
    }
}


void
CElectronRepulsionIntegralsDriver::_compMaxQValuesForGtoPairsBlock(      double*         qValuesBuffer,
                                                                   const CGtoPairsBlock& gtoPairsBlock) const
{
    // determine number of contracted GTOs pairs
    
    auto ngto = gtoPairsBlock.getNumberOfScreenedContrPairs();
    
    // set up empty integrals screener
    
    CCauchySchwarzScreener qqdat;
    
    // loop over contracted GTOs pairs
    
    for (int32_t i = 0; i < ngto; i++)
    {
        auto cpair = gtoPairsBlock.pick(i);
        
        CTwoIntsDistribution cdist(qValuesBuffer, 1, 1, i, dist2e::qvalues);
        
        _compElectronRepulsionForGtoPairsBlocks(cdist, qqdat, cpair, cpair);
    }
}

CVecMemBlock<double>
CElectronRepulsionIntegralsDriver::_getQValuesBuffer(const CGtoPairsContainer& gtoPairsContainer) const
{
    CVecMemBlock<double> buffvec;
    
    auto nppblk = gtoPairsContainer.getNumberOfGtoPairsBlocks();
    
    for (int32_t i = 0; i < nppblk; i++)
    {
        auto cpairs = gtoPairsContainer.getGtoPairsBlock(i);
        
        buffvec.push_back(CMemBlock<double>(cpairs.getNumberOfScreenedContrPairs()));
    }
    
    return buffvec;
}

CMemBlock<int32_t>
CElectronRepulsionIntegralsDriver::_setTasksGrid(const int32_t nBraGtoPairsBlocks,
                                                 const int32_t nKetGtoPairsBlocks,
                                                 const bool    isBraEqualKet) const
{
    // determine number of tasks
    
    int32_t nelem = 0;
    
    if (isBraEqualKet)
    {
        nelem = nBraGtoPairsBlocks * (nBraGtoPairsBlocks + 1) / 2;
    }
    else
    {
        nelem = nBraGtoPairsBlocks * nKetGtoPairsBlocks;
    }
    
    // initialize tasks grid
    
    CMemBlock<int32_t> tgrid(nelem);
    
    if (_locNodes == 1)
    {
        mathfunc::set_to(tgrid.data(), 1, nelem);
    }
    else
    {
        mathfunc::set_to(tgrid.data(), 0, nelem);
        
        auto ndim = nelem / _locNodes;
        
        for (int32_t i = 0; i < ndim; i++)
        {
            tgrid.at(i * _locNodes + _locRank) = 1;
        }
        
        auto nrem = nelem % _locNodes;
        
        if (_locRank < nrem)
        {
            tgrid.at(ndim * _locNodes + _locRank) = 1;
        }
    }
    
    return tgrid;
}
