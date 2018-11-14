//
//                     V.E.L.O.X. C.H.E.M. MP
//      ---------------------------------------------------
//           An Electronic Structure Code for Nanoscale
//
//  Copyright © 2018 by Velox Chem MP developers. All rights reserved.
//  Contact: Zilvinas Rinkevicius (rinkevic@kth.se), KTH, Sweden.

#include "CubeGenerator.hpp"

#include <cmath>

#include "SphericalMomentum.hpp"
#include "BasisFunction.hpp"
#include "MolecularBasis.hpp"

const std::vector<std::vector<int32_t>>
cubes::buildCartesianAngularMomentum(int32_t angl)
{
    std::vector<std::vector<int32_t>> lmn;

    // lexical order of Cartesian angular momentum
    // 1S: 0
    // 3P: x,y,z
    // 6D: xx,xy,xz,yy,yz,zz
    // 10F: ...

    for (int32_t index = 0, i = 0; i <= angl; i++)
    {
        int32_t lx = angl - i;

        for (int32_t j = 0; j <= i; j++, index++)
        {
            int32_t ly = i - j;

            int32_t lz = j;

            lmn.push_back(std::vector<int32_t>({lx, ly, lz}));
        }
    }

    return lmn;
}

double
cubes::getPsiMolecularOrbital(const CMolecule&          molecule,
                              const CMolecularBasis&    basis,
                              const CMolecularOrbitals& molorb,
                              const int32_t             moidx,
                              const double              xp,
                              const double              yp,
                              const double              zp)
{
    auto natoms = molecule.getNumberOfAtoms();

    auto max_angl = basis.getMolecularMaxAngularMomentum(molecule);

    //auto morows = molorb.getNumberOfRows();

    auto mocols = molorb.getNumberOfColumns();

    auto mocoefs = molorb.alphaOrbitals();

    double psi = 0.0;

    // azimuthal quantum number: s,p,d,f,...

    for (int32_t aoidx = 0, angl = 0; angl <= max_angl; angl++)
    {
        CSphericalMomentum sphmom (angl);

        auto nsph = sphmom.getNumberOfComponents();

        auto lmn = buildCartesianAngularMomentum(angl);

        // magnetic quantum number: s,p-1,p0,p+1,d-2,d-1,d0,d+1,d+2,...

        for (int32_t isph = 0; isph < nsph; isph++)
        {
            // prepare Cartesian components

            std::vector<double> fcarts, lx, ly, lz;

            int32_t ncomp = sphmom.getNumberOfFactors(isph);
            
            for (int32_t icomp = 0; icomp < ncomp; icomp++)
            {
                fcarts.push_back(sphmom.getFactors(isph)[icomp]);

                int32_t cartind = sphmom.getIndexes(isph)[icomp];

                lx.push_back(lmn[cartind][0]);

                ly.push_back(lmn[cartind][1]);

                lz.push_back(lmn[cartind][2]);
            }

            // go through atoms

            for (int32_t atomidx = 0; atomidx < natoms; atomidx++)
            {
                // process coordinates

                double rx = xp - molecule.getCoordinatesX()[atomidx];

                double ry = yp - molecule.getCoordinatesY()[atomidx];

                double rz = zp - molecule.getCoordinatesZ()[atomidx];

                double r2 = rx * rx + ry * ry + rz * rz;

                // process atomic orbitals

                int32_t idelem = molecule.getIdsElemental()[atomidx];

                int32_t nao = basis.getNumberOfBasisFunctions(idelem, angl);

                auto basisfunc = basis.getBasisFunctions(idelem, angl);

                for (int32_t i = 0; i < nao; i++, aoidx++)
                {
                    // get MO coefficient

                    double mocoef = mocoefs[aoidx * mocols + moidx];

                    // process primitives

                    auto nprims = basisfunc[i].getNumberOfPrimitiveFunctions();

                    auto exponents = basisfunc[i].getExponents();

                    auto normcoefs = basisfunc[i].getNormalizationFactors();

                    for (int32_t iprim = 0; iprim < nprims; iprim++)
                    {
                        double alpha = exponents[iprim];

                        double coef1 = mocoef * normcoefs[iprim];

                        // transform from Cartesian to spherical harmonics

                        for (int32_t icomp = 0; icomp < ncomp; icomp++)
                        {
                            double coef2 = coef1 * fcarts[icomp];

                            double powxyz = pow(rx, lx[icomp])
                                          * pow(ry, ly[icomp])
                                          * pow(rz, lz[icomp]);

                            psi += coef2 * powxyz * exp(-alpha*r2);
                        }
                    }
                }
            }
        }
    }

    return psi;
}
