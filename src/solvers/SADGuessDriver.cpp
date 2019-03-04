//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "SADGuessDriver.hpp"

#include "SystemClock.hpp"
#include "StringFormat.hpp"
#include "DenseLinearAlgebra.hpp"
#include "DenseDiagonalizer.hpp"
#include "AODensityMatrix.hpp"
#include "DensityMatrixType.hpp"
#include "ErrorHandler.hpp"

CSADGuessDriver::CSADGuessDriver(const int32_t  globRank,
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

CSADGuessDriver::~CSADGuessDriver()
{
    
}

std::vector<double>
CSADGuessDriver::_getOcc1s(double occ) const
{
    //                           1s
    return std::vector<double>({ occ });
}

std::vector<double>
CSADGuessDriver::_getOcc2s(double occ) const
{
    //                           1s   2s
    return std::vector<double>({ 1.0, occ });
}

std::vector<double>
CSADGuessDriver::_getOcc2s2p(double occ) const
{
    //                           1s   2s   2p-1 2p0  2p+1
    return std::vector<double>({ 1.0, occ, occ, occ, occ });
}

std::vector<double>
CSADGuessDriver::_getOcc3s(double occ) const
{
    //                           1s   2s   3s   2p-1 2p0  2p+1
    return std::vector<double>({ 1.0, 1.0, occ, 1.0, 1.0, 1.0 });
}

std::vector<double>
CSADGuessDriver::_getOcc3s3p(double occ) const
{
    //                           1s   2s   3s   2p-1 3p-1 2p0  3p0  2p+1 3p+1
    return std::vector<double>({ 1.0, 1.0, occ, 1.0, occ, 1.0, occ, 1.0, occ });
}

std::vector<double>
CSADGuessDriver::_getOcc4s(double occ) const
{
    //                           1s   2s   3s   4s   2p-1 3p-1 2p0  3p0  2p+1 3p+1
    return std::vector<double>({ 1.0, 1.0, 1.0, occ, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0 });
}

std::vector<double>
CSADGuessDriver::_getOcc3d(double occ) const
{
    //                           1s   2s   3s   4s   2p-1 3p-1 2p0  3p0  2p+1 3p+1
    return std::vector<double>({ 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,

    //                           3d-2 3d-1 3d0  3d+1 3d+2
                                 occ, occ, occ, occ, occ });
}

std::vector<double>
CSADGuessDriver::_getOcc4s4p(double occ) const
{
    //                           1s   2s   3s   4s   2p-1 3p-1 4p-1 2p0  3p0  4p0  2p+1 3p+1 4p+1
    return std::vector<double>({ 1.0, 1.0, 1.0, occ, 1.0, 1.0, occ, 1.0, 1.0, occ, 1.0, 1.0, occ,

    //                           3d-2 3d-1 3d0  3d+1 3d+2
                                 1.0, 1.0, 1.0, 1.0, 1.0 });
}

std::vector< std::vector<double> >
CSADGuessDriver::_buildQocc() const
{
    std::vector< std::vector<double> > qocc;

    // dummy atom

    qocc.push_back(std::vector<double> ());

    // H,He

    qocc.push_back(_getOcc1s(0.5));

    qocc.push_back(_getOcc1s(1.0));

    // Li,Be

    qocc.push_back(_getOcc2s(0.5));

    qocc.push_back(_getOcc2s(1.0));

    // B,C,N,O,F,Ne

    qocc.push_back(_getOcc2s2p(0.375));

    qocc.push_back(_getOcc2s2p(0.500));

    qocc.push_back(_getOcc2s2p(0.625));

    qocc.push_back(_getOcc2s2p(0.750));

    qocc.push_back(_getOcc2s2p(0.875));

    qocc.push_back(_getOcc2s2p(1.000));

    // Na,Mg

    qocc.push_back(_getOcc3s(0.5));

    qocc.push_back(_getOcc3s(1.0));

    // Al,Si,P,S,Cl,Ar

    qocc.push_back(_getOcc3s3p(0.375));

    qocc.push_back(_getOcc3s3p(0.500));

    qocc.push_back(_getOcc3s3p(0.625));

    qocc.push_back(_getOcc3s3p(0.750));

    qocc.push_back(_getOcc3s3p(0.875));

    qocc.push_back(_getOcc3s3p(1.000));

    // K,Ca

    qocc.push_back(_getOcc4s(0.5));

    qocc.push_back(_getOcc4s(1.0));

    // Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn

    qocc.push_back(_getOcc3d(0.1));

    qocc.push_back(_getOcc3d(0.2));

    qocc.push_back(_getOcc3d(0.3));

    qocc.push_back(_getOcc3d(0.4));

    qocc.push_back(_getOcc3d(0.5));

    qocc.push_back(_getOcc3d(0.6));

    qocc.push_back(_getOcc3d(0.7));

    qocc.push_back(_getOcc3d(0.8));

    qocc.push_back(_getOcc3d(0.9));

    qocc.push_back(_getOcc3d(1.0));

    // Ga,Ge,As,Se,Br,Kr

    qocc.push_back(_getOcc4s4p(0.375));

    qocc.push_back(_getOcc4s4p(0.500));

    qocc.push_back(_getOcc4s4p(0.625));

    qocc.push_back(_getOcc4s4p(0.750));

    qocc.push_back(_getOcc4s4p(0.875));

    qocc.push_back(_getOcc4s4p(1.000));

    return qocc;
}

std::vector< std::vector<int32_t> >
CSADGuessDriver::getAOIndicesOfAtoms(const CMolecule&       molecule,
                                     const CMolecularBasis& basis) const
{
    std::vector< std::vector<int32_t> > aoinds_atoms;

    int32_t natoms = molecule.getNumberOfAtoms();

    for (int32_t atomidx = 0; atomidx < natoms; atomidx++)
    {
        aoinds_atoms.push_back(std::vector<int32_t>());
    }

    int32_t max_angl = basis.getMolecularMaxAngularMomentum(molecule);

    for (int32_t aoidx = 0, angl = 0; angl <= max_angl; angl++)
    {
        for (int32_t s = -angl; s <= angl; s++)
        {
            for (int32_t atomidx = 0; atomidx < natoms; atomidx++)
            {
                int32_t idelem = molecule.getIdsElemental()[atomidx];

                int32_t nao = basis.getNumberOfBasisFunctions(idelem, angl);

                for (int32_t i = 0; i < nao; i++, aoidx++)
                {
                    aoinds_atoms[atomidx].push_back(aoidx);
                }
            }
        }
    }

    return aoinds_atoms;
}

CAODensityMatrix
CSADGuessDriver::compute(const CMolecule&       molecule,
                         const CMolecularBasis& basis_1,
                         const CMolecularBasis& basis_2,
                         const COverlapMatrix&  S12,
                         const COverlapMatrix&  S22,
                               MPI_Comm         comm) const 
{
    CSystemClock timer;
    
    CAODensityMatrix dsad;
    
    if (_locRank == mpi::master())
    {
        // generate SAD guess
        
        dsad = _compSADGuess(molecule, basis_1, basis_2, S12, S22);
    }
    
    return dsad;
}

CAODensityMatrix
CSADGuessDriver::_compSADGuess(const CMolecule&       molecule,
                               const CMolecularBasis& basis_1,
                               const CMolecularBasis& basis_2,
                               const COverlapMatrix&  S12,
                               const COverlapMatrix&  S22) const
{
    const bool closed_shell = (molecule.getMultiplicity() == 1);

    auto natoms = molecule.getNumberOfAtoms();

    auto nao_1 = S12.getNumberOfRows();

    auto nao_2 = S12.getNumberOfColumns();

    // sanity checks

    std::string err_ovl_size("SADGuessDriver - Mismatch between overlap matrices");

    errors::assertMsgCritical(nao_2 == S22.getNumberOfRows(), err_ovl_size);

    errors::assertMsgCritical(nao_2 == S22.getNumberOfColumns(), err_ovl_size);

    // AO indices for atoms

    auto aoinds_atoms_1 = getAOIndicesOfAtoms(molecule, basis_1);

    auto aoinds_atoms_2 = getAOIndicesOfAtoms(molecule, basis_2);

    // more sanity checks

    int32_t count_ao_1 = 0;

    int32_t count_ao_2 = 0;

    for (int32_t atomidx = 0; atomidx < natoms; atomidx++)
    {
        count_ao_1 += static_cast<int32_t>(aoinds_atoms_1[atomidx].size());

        count_ao_2 += static_cast<int32_t>(aoinds_atoms_2[atomidx].size());
    }

    std::string err_bas_size("SADGuessDriver - Mismatch between basis set & overlap matrix");

    errors::assertMsgCritical(count_ao_1 == nao_1 && count_ao_2 == nao_2, err_bas_size);

    errors::assertMsgCritical(count_ao_2 == S22.getNumberOfRows(), err_bas_size);

    // occupation numbers

    auto qocc = _buildQocc();

    // take care of electrons from net charge and/or multiplicity

    double charge = molecule.getCharge();

    double mult_1 = static_cast<double>(molecule.getMultiplicity() - 1);

    double alpha_elec = 0.5 * (-charge + mult_1);

    double beta_elec  = 0.5 * (-charge - mult_1);

    alpha_elec /= static_cast<double>(nao_1);

    beta_elec  /= static_cast<double>(nao_1);

    // C_SAD matrix

    CDenseMatrix csad_alpha (nao_2, nao_1);

    CDenseMatrix csad_beta  (nao_2, nao_1);

    csad_alpha.zero();

    csad_beta.zero();

    #pragma omp parallel for schedule(dynamic)
    for (int32_t atomidx = 0; atomidx < natoms; atomidx++)
    {
        // elemental index (nuclear charge) for this atom

        const int32_t idelem = molecule.getIdsElemental()[atomidx];

        // AO indices for this atom

        const std::vector<int32_t>& aoinds_1 = aoinds_atoms_1[atomidx];

        const std::vector<int32_t>& aoinds_2 = aoinds_atoms_2[atomidx];
        
        // set up AO indices dimensions
        
        auto naodim_1 = static_cast<int32_t>(aoinds_1.size());
        
        auto naodim_2 = static_cast<int32_t>(aoinds_2.size());

        // size checking
        
        std::string err_ao_size("SADGuessDriver - Mismatch between basis set & occupation number");

        errors::assertMsgCritical(qocc[idelem].size() == aoinds_1.size(), err_ao_size);
        
        // atomic block of AOs

        CDenseMatrix block_12 (naodim_1, naodim_2);

        CDenseMatrix block_22 (naodim_2, naodim_2);

        block_12.zero();

        block_22.zero();

        for (int32_t i = 0; i < naodim_1; i++)
        {
            for (int32_t j = 0; j < naodim_2; j++)
            {
                block_12.values()[i * naodim_2 + j] = 
                    
                    S12.values()[aoinds_1[i] * nao_2 + aoinds_2[j]];
            }
        }

        for (int32_t i = 0; i < naodim_2; i++)
        {
            for (int32_t j = 0; j < naodim_2; j++)
            {
                block_22.values()[i * naodim_2 + j] = 
                    
                    S22.values()[aoinds_2[i] * nao_2 + aoinds_2[j]];
            }
        }

        // A = S12' C1(identity)

        CDenseMatrix mat_c1 (naodim_1, naodim_1);

        mat_c1.zero();

        for (int32_t i = 0; i < naodim_1; i++)
        {
            mat_c1.values()[i * naodim_1 + i] = 1.0;
        }

        auto mat_a = denblas::multAtB(block_12, mat_c1);

        // S22^-1

        CDenseDiagonalizer diagdrv;

        diagdrv.diagonalize(block_22);

        std::string err_diag("SADGuessDriver - Matrix diagonalization failed");

        errors::assertMsgCritical(diagdrv.getState(), err_diag);

        auto block_22_inv = diagdrv.getInvertedMatrix();

        // M = A' S22^-1 A

        auto prod = denblas::multAB(block_22_inv, mat_a);

        auto mat_m = denblas::multAtB(mat_a, prod);

        // M^-1/2

        diagdrv.diagonalize(mat_m);

        errors::assertMsgCritical(diagdrv.getState(), err_diag);

        auto mat_m_invsqrt = diagdrv.getInvertedSqrtMatrix();

        // C2 = S22^-1 A M^-1/2

        prod = denblas::multAB(mat_a, mat_m_invsqrt);

        auto mat_c2 = denblas::multAB(block_22_inv, prod);

        // update csad_alpha

        for (int32_t j = 0; j < naodim_2; j++)
        {
            for (int32_t i = 0; i < naodim_1; i++)
            {
                csad_alpha.values()[aoinds_2[j] * nao_1 + aoinds_1[i]] = 

                    mat_c2.values()[j * naodim_1 + i] * sqrt(qocc[idelem][i] + alpha_elec);
            }
        }

        // update csad_beta

        if (! closed_shell)
        {
            for (int32_t j = 0; j < naodim_2; j++)
            {
                for (int32_t i = 0; i < naodim_1; i++)
                {
                    csad_beta.values()[aoinds_2[j] * nao_1 + aoinds_1[i]] = 

                        mat_c2.values()[j * naodim_1 + i] * sqrt(qocc[idelem][i] + beta_elec);
                }
            }
        }
    }

    // D_SAD density matrix

    std::vector<CDenseMatrix> dsad;

    if (closed_shell)
    {
        dsad.push_back(denblas::multABt(csad_alpha, csad_alpha));

        return CAODensityMatrix(dsad, denmat::rest);
    }
    else
    {
        dsad.push_back(denblas::multABt(csad_alpha, csad_alpha));

        dsad.push_back(denblas::multABt(csad_beta,  csad_beta));

        return CAODensityMatrix(dsad, denmat::unrest);
    }
}
