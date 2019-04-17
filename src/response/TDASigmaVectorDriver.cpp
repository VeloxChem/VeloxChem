//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "TDASigmaVectorDriver.hpp"

#include "DensityMatrixType.hpp"
#include "AOFockMatrix.hpp"
#include "ElectronRepulsionIntegralsDriver.hpp"

CTDASigmaVectorDriver::CTDASigmaVectorDriver(const int32_t  globRank,
                                             const int32_t  globNodes,
                                                   MPI_Comm comm)

    : _globRank(globRank)

    , _globNodes(globNodes)

    , _isLocalMode(false)
{
    _locRank  = mpi::rank(comm);
    
    _locNodes = mpi::nodes(comm);
    
    _isLocalMode = !mpi::compare(comm, MPI_COMM_WORLD);
}

CTDASigmaVectorDriver::~CTDASigmaVectorDriver()
{
    
}

std::vector<CDenseMatrix>
CTDASigmaVectorDriver::compute(const std::vector<CExcitationVector>& zVectors,
                               const CScreeningContainer&            screeningContainer,
                               const CMolecularOrbitals&             molecularOrbitals,
                               const CMolecule&                      molecule,
                               const CMolecularBasis&                basis,
                                     MPI_Comm                        comm) const
{
    auto sig_vecs = _allocSigmaVectors(zVectors);
 
    _addCanonicalFockContribution(sig_vecs, zVectors, molecularOrbitals);
    
    _addFirstOrderFockContribution(sig_vecs, zVectors, screeningContainer,
                                   molecularOrbitals, molecule, basis, comm);
    
    return sig_vecs;
}

std::vector<CDenseMatrix>
CTDASigmaVectorDriver::_allocSigmaVectors(const std::vector<CExcitationVector>& zVectors) const
{
    std::vector<CDenseMatrix> sig_vecs;
    
    // determine number of sigma vectors
    
    auto nvecs = static_cast<int32_t>(zVectors.size());
    
    if (nvecs > 0)
    {
        // allocate sigma vectors
        
        for (int32_t i = 0; i < nvecs; i++)
        {
            auto ndim = zVectors[i].getNumberOfExcitations();
            
            sig_vecs.push_back(CDenseMatrix(ndim, 1));
        }
        
        // zero sigma vectors
        
        for (int32_t i = 0; i < nvecs; i++)
        {
            sig_vecs[i].zero();
        }
    }
    
    return sig_vecs;
}

void
CTDASigmaVectorDriver::_addCanonicalFockContribution(      std::vector<CDenseMatrix>&      sigmaVectors,
                                                     const std::vector<CExcitationVector>& zVectors,
                                                     const CMolecularOrbitals&             molecularOrbitals) const
{
    for (size_t i = 0; i < sigmaVectors.size(); i++)
    {
        // compute approximate diagonal of A matrix
        
        auto diagmat = zVectors[i].getApproximateDiagonal(molecularOrbitals);
        
        auto devals = diagmat.data();
        
        // set up pointer to sigma and Z vector values
        
        auto sigdat = sigmaVectors[i].values();
        
        auto zdat = zVectors[i].getCoefficientsZ();
        
        // add diagonal (e_a - e_i) z_ia contribution
        
        auto ndim = sigmaVectors[i].getNumberOfRows();
        
        #pragma omp simd aligned(sigdat, devals, zdat)
        for (int32_t j = 0; j < ndim; j++)
        {
            sigdat[j] += devals[j] * zdat[j];
        }
    }
}

void
CTDASigmaVectorDriver::_addFirstOrderFockContribution(      std::vector<CDenseMatrix>&      sigmaVectors,
                                                      const std::vector<CExcitationVector>& zVectors,
                                                      const CScreeningContainer&            screeningContainer,
                                                      const CMolecularOrbitals&             molecularOrbitals,
                                                      const CMolecule&                      molecule,
                                                      const CMolecularBasis&                basis,
                                                            MPI_Comm                        comm) const
{
    auto nvecs = static_cast<int32_t>(sigmaVectors.size());
    
    // create first order transformed density
    
    CAODensityMatrix dmat;
    
    dmat.setDensityType(denmat::rgen);
    
    for (int32_t i = 0; i < nvecs; i++)
    {
        dmat.append(zVectors[i].getDensityZ(molecularOrbitals));
    }
    
    dmat.broadcast(_locRank, comm);
    
    // compute AO Fock matrices
    
    CAOFockMatrix faomat(dmat);
    
    CElectronRepulsionIntegralsDriver eri_drv(comm);
    
    eri_drv.compute(faomat, dmat, molecule, basis, screeningContainer);
    
    faomat.reduce_sum(_locRank, _locNodes, comm);
    
    // add contributions to sigma vectors on master node
    
    if (_locRank == mpi::master())
    {
        for (int32_t i = 0; i < nvecs; i++)
        {
            // compute MO Fock matrix
            
            auto fmomat = molecularOrbitals.transform(faomat.getReferenceToFock(i),
                                                      szblock::aa);
            
            // set up pointers to creation/anihilation operator indexes
            
            auto bidx = zVectors[i].getBraIndexes();
            
            auto kidx = zVectors[i].getKetIndexes();
            
            // set up pointer to sigma vector values
            
            auto sigdat = sigmaVectors[i].values();
            
            // add first order Fock contribution
            
            auto fdat = fmomat.values();
            
            auto ncol = fmomat.getNumberOfColumns();
            
            auto ndim = sigmaVectors[i].getNumberOfRows();
            
            for (int32_t j = 0; j < ndim; j++)
            {
                sigdat[j] += fdat[bidx[j] * ncol + kidx[j]];
            }
        }
    }
}
