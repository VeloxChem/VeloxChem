//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "TDASigmaVectorDriver.hpp"

#include "AOFockMatrix.hpp"
#include "DensityMatrixType.hpp"
#include "ElectronRepulsionIntegralsDriver.hpp"
#include "XCIntegrator.hpp"

CTDASigmaVectorDriver::CTDASigmaVectorDriver(MPI_Comm comm)
{
    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    mpi::duplicate(comm, &_locComm);
}

CTDASigmaVectorDriver::~CTDASigmaVectorDriver()
{
    mpi::destroy(&_locComm);
}

std::vector<CDenseMatrix>
CTDASigmaVectorDriver::compute(const std::vector<CExcitationVector>& zVectors,
                               const bool                            isTripletStates,
                               const CScreeningContainer&            screeningContainer,
                               const CMolecularGrid&                 molecularGrid,
                               const CXCFunctional&                  xcFunctional, 
                               const CMolecularOrbitals&             molecularOrbitals,
                               const CMolecule&                      molecule,
                               const CMolecularBasis&                basis) const
{
    auto sig_vecs = _allocSigmaVectors(zVectors);

    _addCanonicalFockContribution(sig_vecs, zVectors, molecularOrbitals);

    _addFirstOrderFockContribution(sig_vecs, zVectors, isTripletStates, screeningContainer,
                                   molecularGrid, xcFunctional, molecularOrbitals, molecule, basis);

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
CTDASigmaVectorDriver::_addCanonicalFockContribution(std::vector<CDenseMatrix>&            sigmaVectors,
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
CTDASigmaVectorDriver::_addFirstOrderFockContribution(std::vector<CDenseMatrix>&            sigmaVectors,
                                                      const std::vector<CExcitationVector>& zVectors,
                                                      const bool                            isTripletStates,
                                                      const CScreeningContainer&            screeningContainer,
                                                      const CMolecularGrid&                 molecularGrid,
                                                      const CXCFunctional&                  xcFunctional,
                                                      const CMolecularOrbitals&             molecularOrbitals,
                                                      const CMolecule&                      molecule,
                                                      const CMolecularBasis&                basis) const
{
    auto nvecs = static_cast<int32_t>(sigmaVectors.size());

    // create first order transformed density

    CAODensityMatrix dmat;

    dmat.setDensityType(denmat::rgen);

    for (int32_t i = 0; i < nvecs; i++)
    {
        dmat.append(zVectors[i].getDensityZ(molecularOrbitals));
    }

    dmat.broadcast(_locRank, _locComm);

    // compute AO Fock matrices

    CAOFockMatrix faomat(dmat);

    double fock_prefactor = 1.0;

    if (isTripletStates)
    {
        fock_prefactor = -1.0;

        for (int32_t i = 0; i < nvecs; i++)
        {
            faomat.setFockType(fockmat::rgenk, i);
        }
    }
    
    // 2e-contribution to Fock matrix

    CElectronRepulsionIntegralsDriver eri_drv(_locComm);

    eri_drv.compute(faomat, dmat, molecule, basis, screeningContainer);
    
    faomat.reduce_sum(_locRank, _locNodes, _locComm);
    
    // exchange-correlation contribution to Fock matrix
    
    if (!xcFunctional.isUndefined())
    {
        CXCIntegrator xcdrv(_locComm);
        
        xcdrv.integrate(faomat, dmat, molecule, basis, molecularGrid, xcFunctional.getLabel()); 
    }

    // add contributions to sigma vectors on master node

    if (_locRank == mpi::master())
    {
        for (int32_t i = 0; i < nvecs; i++)
        {
            // compute MO Fock matrix

            auto fmomat = molecularOrbitals.transform(faomat.getReferenceToFock(i), szblock::aa);

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
                sigdat[j] += fdat[bidx[j] * ncol + kidx[j]] * fock_prefactor;
            }
        }
    }
}
