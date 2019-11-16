//
//                             VELOXCHEM
//      ---------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2019 by VeloxChem developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "PartialCharges.hpp"
#include "Molecule.hpp"
#include "DenseMatrix.hpp"
#include "DenseLinearAlgebra.hpp"
#include "DenseDiagonalizer.hpp"
#include "ErrorHandler.hpp"

namespace parchg {  // parchg namespace

std::vector<double>
getPartialCharges(const CMolecule& molecule, double netcharge)
{
    auto electronegativity = getElectronegativity();

    auto hardness = getHardness();

    auto natoms = molecule.getNumberOfAtoms();

    auto idselem = molecule.getIdsElemental();

    // form left-hand side matrix

    CDenseMatrix matrix(natoms + 1, natoms + 1);

    for (int32_t i = 0; i < natoms; i++)
    {
        matrix.values()[i * (natoms + 1) + i] = 2.0 * hardness[idselem[i]];

        matrix.values()[i * (natoms + 1) + natoms] = 1.0;

        matrix.values()[natoms * (natoms + 1) + i] = 1.0;
    }

    matrix.values()[natoms * (natoms + 1) + natoms] = 0.0;

    auto rx = molecule.getCoordinatesX();

    auto ry = molecule.getCoordinatesY();

    auto rz = molecule.getCoordinatesZ();

    for (int32_t i = 0; i < natoms; i++)
    {
        for (int32_t j = i + 1; j < natoms; j++)
        {
            double rij = std::sqrt((rx[i] - rx[j]) * (rx[i] - rx[j]) +
                                   (ry[i] - ry[j]) * (ry[i] - ry[j]) +
                                   (rz[i] - rz[j]) * (rz[i] - rz[j]));

            matrix.values()[i * (natoms + 1) + j] = 1.0 / rij;

            matrix.values()[j * (natoms + 1) + i] = 1.0 / rij;
        }
    }

    // form right-hand side vector

    CDenseMatrix vector(natoms + 1, 1);

    for (int32_t i = 0; i < natoms; i++)
    {
        vector.values()[i] = electronegativity[idselem[i]];
    }

    vector.values()[natoms] = netcharge;

    // get partial charges

    CDenseDiagonalizer diagdrv;

    diagdrv.diagonalize(matrix);

    std::string err_diag("PartialCharges - Matrix diagonalization failed");

    errors::assertMsgCritical(diagdrv.getState(), err_diag);

    auto matrix_inv = diagdrv.getInvertedMatrix();

    auto solution = denblas::multAB(matrix_inv, vector);

    std::vector<double> partialcharges(natoms);

    for (int32_t i = 0; i < natoms; i++)
    {
        partialcharges[i] = solution.values()[i];
    }

    return partialcharges;
}

std::vector<double>
getElectronegativity()
{
    return std::vector<double>({

        // dummy
        0.00,

        // H-B
        -0.193,
        -0.242,
         0.003,
         0.125,
        -0.141,

        // C-Ne
        -0.188,
        -0.232,
        -0.293,
        -0.349,
        -0.398,

        // Na-P
         0.005,
         0.099,
        -0.086,
        -0.096,
        -0.099,

        // S-Ca
        -0.106,
        -0.105,
        -0.092,
         0.010,
         0.058,

        // Sc-Mn
        -0.029,
        -0.062,
         0.018,
         0.063,
        -0.106,

        // Fe-Zn
        -0.020,
        -0.003,
         0.015,
         0.065,
         0.054,

        // Ga-Br
        -0.093,
        -0.088,
        -0.074,
        -0.062,
        -0.040,

        // Kr
        -0.006});
}

std::vector<double>
getHardness()
{
    return std::vector<double>({

        // dummy
        0.000,

        // H-B
        0.443,
        1.156,
        0.094,
        0.180,
        0.193,

        // C-Ne
        0.334,
        0.517,
        0.712,
        0.952,
        1.237,

        // Na-P
        0.085,
        0.153,
        0.121,
        0.195,
        0.293,

        // S-Ca
        0.393,
        0.518,
        0.669,
        0.063,
        0.137,

        // Sc-Mn
        0.064,
        0.102,
        0.131,
        0.171,
        0.013,

        // Fe-Zn
        0.159,
        0.177,
        0.196,
        0.187,
        0.234,

        // Ga-Br
        0.127,
        0.183,
        0.258,
        0.328,
        0.417,

        // Kr
        0.524});
}

}  // namespace parchg
