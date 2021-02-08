//
//                           VELOXCHEM 1.0-RC
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2020 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact

#include "SADGuessDriver.hpp"

#include <mpi.h>

#include "AODensityMatrix.hpp"
#include "DenseDiagonalizer.hpp"
#include "DenseLinearAlgebra.hpp"
#include "DensityMatrixType.hpp"
#include "ErrorHandler.hpp"
#include "MpiFunc.hpp"
#include "PartialCharges.hpp"
#include "StringFormat.hpp"

CSADGuessDriver::CSADGuessDriver(MPI_Comm comm)
{
    _locRank = mpi::rank(comm);

    _locNodes = mpi::nodes(comm);

    _locComm = comm;
}

CSADGuessDriver::~CSADGuessDriver()
{
}

std::vector<double>
CSADGuessDriver::_getOcc1s(double occ) const
{
    //                           1s
    return std::vector<double>({occ});
}

std::vector<double>
CSADGuessDriver::_getOcc2s(double occ) const
{
    //                           1s   2s
    return std::vector<double>({1.0, occ});
}

std::vector<double>
CSADGuessDriver::_getOcc2s2p(double occ) const
{
    //                           1s   2s   2p-1 2p0  2p+1
    return std::vector<double>({1.0, occ, occ, occ, occ});
}

std::vector<double>
CSADGuessDriver::_getOcc3s(double occ) const
{
    //                           1s   2s   3s   2p-1 2p0  2p+1
    return std::vector<double>({1.0, 1.0, occ, 1.0, 1.0, 1.0});
}

std::vector<double>
CSADGuessDriver::_getOcc3s3p(double occ) const
{
    //                           1s   2s   3s   2p-1 3p-1 2p0  3p0  2p+1 3p+1
    return std::vector<double>({1.0, 1.0, occ, 1.0, occ, 1.0, occ, 1.0, occ});
}

std::vector<double>
CSADGuessDriver::_getOcc4s(double occ) const
{
    //                           1s   2s   3s   4s   2p-1 3p-1 2p0  3p0  2p+1 3p+1
    return std::vector<double>({1.0, 1.0, 1.0, occ, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0});
}

std::vector<double>
CSADGuessDriver::_getOcc3d(double occ) const
{
    //                           1s   2s   3s   4s  2p-1 3p-1 2p0  3p0  2p+1 3p+1 3d-2 3d-1 3d0  3d+1 3d+2
    return std::vector<double>({1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, occ, occ, occ, occ, occ});
}

std::vector<double>
CSADGuessDriver::_getOcc4s4p(double occ) const
{
    //                           1s   2s   3s   4s  2p-1 3p-1 4p-1 2p0  3p0  4p0  2p+1 3p+1 4p+1 3d-2 3d-1 3d0  3d+1 3d+2
    return std::vector<double>({1.0, 1.0, 1.0, occ, 1.0, 1.0, occ, 1.0, 1.0, occ, 1.0, 1.0, occ, 1.0, 1.0, 1.0, 1.0, 1.0});
}

std::vector<double>
CSADGuessDriver::getOccupationNumbersForElement(const int32_t elem_id, const double num_electrons) const
{
    double nelec = std::min(std::max(-0.5, num_electrons), 0.5);

    switch (elem_id)
    {
            // dummy atom

        case 0:
            return std::vector<double>();

            // H,He

        case 1:
            return _getOcc1s(0.5 + nelec);

        case 2:
            return _getOcc1s(1.0 + nelec);

            // Li,Be

        case 3:
            return _getOcc2s(0.5 + nelec);

        case 4:
            return _getOcc2s(1.0 + nelec);

            // B,C,N,O,F,Ne

        case 5:
            return _getOcc2s2p(0.375 + nelec / 4.0);

        case 6:
            return _getOcc2s2p(0.500 + nelec / 4.0);

        case 7:
            return _getOcc2s2p(0.625 + nelec / 4.0);

        case 8:
            return _getOcc2s2p(0.750 + nelec / 4.0);

        case 9:
            return _getOcc2s2p(0.875 + nelec / 4.0);

        case 10:
            return _getOcc2s2p(1.000 + nelec / 4.0);

            // Na,Mg

        case 11:
            return _getOcc3s(0.5 + nelec);

        case 12:
            return _getOcc3s(1.0 + nelec);

            // Al,Si,P,S,Cl,Ar

        case 13:
            return _getOcc3s3p(0.375 + nelec / 4.0);

        case 14:
            return _getOcc3s3p(0.500 + nelec / 4.0);

        case 15:
            return _getOcc3s3p(0.625 + nelec / 4.0);

        case 16:
            return _getOcc3s3p(0.750 + nelec / 4.0);

        case 17:
            return _getOcc3s3p(0.875 + nelec / 4.0);

        case 18:
            return _getOcc3s3p(1.000 + nelec / 4.0);

            // K,Ca

        case 19:
            return _getOcc4s(0.5 + nelec);

        case 20:
            return _getOcc4s(1.0 + nelec);

            // Sc,Ti,V,Cr,Mn,Fe,Co,Ni,Cu,Zn

        case 21:
            return _getOcc3d(0.1 + nelec / 5.0);

        case 22:
            return _getOcc3d(0.2 + nelec / 5.0);

        case 23:
            return _getOcc3d(0.3 + nelec / 5.0);

        case 24:
            return _getOcc3d(0.4 + nelec / 5.0);

        case 25:
            return _getOcc3d(0.5 + nelec / 5.0);

        case 26:
            return _getOcc3d(0.6 + nelec / 5.0);

        case 27:
            return _getOcc3d(0.7 + nelec / 5.0);

        case 28:
            return _getOcc3d(0.8 + nelec / 5.0);

        case 29:
            return _getOcc3d(0.9 + nelec / 5.0);

        case 30:
            return _getOcc3d(1.0 + nelec / 5.0);

            // Ga,Ge,As,Se,Br,Kr

        case 31:
            return _getOcc4s4p(0.375 + nelec / 4.0);

        case 32:
            return _getOcc4s4p(0.500 + nelec / 4.0);

        case 33:
            return _getOcc4s4p(0.625 + nelec / 4.0);

        case 34:
            return _getOcc4s4p(0.750 + nelec / 4.0);

        case 35:
            return _getOcc4s4p(0.875 + nelec / 4.0);

        case 36:
            return _getOcc4s4p(1.000 + nelec / 4.0);

            // other elements: not implemented

        default:
            return std::vector<double>();
    }
}

std::vector<std::vector<double>>
CSADGuessDriver::getOccupationNumbersForMolecule(const CMolecule& molecule, const double nelec) const
{
    std::vector<std::vector<double>> occnumbers;

    auto natoms = molecule.getNumberOfAtoms();

    auto partialcharges = parchg::getPartialCharges(molecule, -nelec * 2.0);

    auto idselem = molecule.getIdsElemental();

    for (int32_t i = 0; i < natoms; i++)
    {
        auto occ = getOccupationNumbersForElement(idselem[i], -partialcharges[i] / 2.0);

        occnumbers.push_back(occ);
    }

    return occnumbers;
}

std::vector<std::vector<int32_t>>
CSADGuessDriver::getAOIndicesOfAtoms(const CMolecule& molecule, const CMolecularBasis& basis) const
{
    std::vector<std::vector<int32_t>> aoinds_atoms;

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
                         const bool             restricted) const
{
    CAODensityMatrix dsad;

    if (_locRank == mpi::master())
    {
        // generate SAD guess

        dsad = _compSADGuess(molecule, basis_1, basis_2, S12, S22, restricted);
    }

    return dsad;
}

CAODensityMatrix
CSADGuessDriver::_compSADGuess(const CMolecule&       molecule,
                               const CMolecularBasis& basis_1,
                               const CMolecularBasis& basis_2,
                               const COverlapMatrix&  S12,
                               const COverlapMatrix&  S22,
                               const bool             restricted) const
{
    auto natoms = molecule.getNumberOfAtoms();

    auto nao_1 = S12.getNumberOfRows();

    auto nao_2 = S12.getNumberOfColumns();

    // sanity checks

    std::string err_ovl_size("SADGuessDriver.compute: Mismatch between overlap matrices");

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

    std::string err_bas_size("SADGuessDriver.compute: Mismatch between basis set and overlap matrix");

    errors::assertMsgCritical(count_ao_1 == nao_1 && count_ao_2 == nao_2, err_bas_size);

    errors::assertMsgCritical(count_ao_2 == S22.getNumberOfRows(), err_bas_size);

    // number of excessive electrons

    double charge = molecule.getCharge();

    double mult_1 = static_cast<double>(molecule.getMultiplicity() - 1);

    double alpha_elec = 0.5 * (-charge + mult_1);

    double beta_elec = 0.5 * (-charge - mult_1);

    // occupation numbers

    auto alpha_occ = getOccupationNumbersForMolecule(molecule, alpha_elec);

    auto beta_occ = getOccupationNumbersForMolecule(molecule, beta_elec);

    // C_SAD matrix

    CDenseMatrix csad_alpha(nao_2, nao_1);

    CDenseMatrix csad_beta(nao_2, nao_1);

    csad_alpha.zero();

    csad_beta.zero();

    #pragma omp parallel for schedule(dynamic)
    for (int32_t atomidx = 0; atomidx < natoms; atomidx++)
    {
        // AO indices for this atom

        const std::vector<int32_t>& aoinds_1 = aoinds_atoms_1[atomidx];

        const std::vector<int32_t>& aoinds_2 = aoinds_atoms_2[atomidx];

        // set up AO indices dimensions

        auto naodim_1 = static_cast<int32_t>(aoinds_1.size());

        auto naodim_2 = static_cast<int32_t>(aoinds_2.size());

        // size checking

        std::string err_ao_size("SADGuessDriver.compute: Mismatch between basis set and occupation number");

        errors::assertMsgCritical(alpha_occ[atomidx].size() == aoinds_1.size(), err_ao_size);

        errors::assertMsgCritical(beta_occ[atomidx].size() == aoinds_1.size(), err_ao_size);

        // atomic block of AOs

        CDenseMatrix block_12(naodim_1, naodim_2);

        CDenseMatrix block_22(naodim_2, naodim_2);

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

        CDenseMatrix mat_c1(naodim_1, naodim_1);

        mat_c1.zero();

        for (int32_t i = 0; i < naodim_1; i++)
        {
            mat_c1.values()[i * naodim_1 + i] = 1.0;
        }

        auto mat_a = denblas::multAtB(block_12, mat_c1);

        // S22^-1

        CDenseDiagonalizer diagdrv;

        diagdrv.diagonalize(block_22);

        std::string err_diag("SADGuessDriver.compute: Matrix diagonalization failed");

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

                    mat_c2.values()[j * naodim_1 + i] * sqrt(alpha_occ[atomidx][i]);
            }
        }

        // update csad_beta

        if (!restricted)
        {
            for (int32_t j = 0; j < naodim_2; j++)
            {
                for (int32_t i = 0; i < naodim_1; i++)
                {
                    csad_beta.values()[aoinds_2[j] * nao_1 + aoinds_1[i]] =

                        mat_c2.values()[j * naodim_1 + i] * sqrt(beta_occ[atomidx][i]);
                }
            }
        }
    }

    // D_SAD density matrix

    std::vector<CDenseMatrix> dsad;

    if (restricted)
    {
        dsad.push_back(denblas::multABt(csad_alpha, csad_alpha));

        return CAODensityMatrix(dsad, denmat::rest);
    }
    else
    {
        dsad.push_back(denblas::multABt(csad_alpha, csad_alpha));

        dsad.push_back(denblas::multABt(csad_beta, csad_beta));

        return CAODensityMatrix(dsad, denmat::unrest);
    }
}
