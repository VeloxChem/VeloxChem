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

    for (int32_t i = 0; i <= angl; i++)
    {
        int32_t lx = angl - i;

        for (int32_t j = 0; j <= i; j++)
        {
            int32_t ly = i - j;

            int32_t lz = j;

            lmn.push_back(std::vector<int32_t>({lx, ly, lz}));
        }
    }

    return lmn;
}

const std::vector<double>
cubes::getPhiAO(const CMolecule&          molecule,
                const CMolecularBasis&    basis,
                const double              xp,
                const double              yp,
                const double              zp)
{
    auto natoms = molecule.getNumberOfAtoms();

    auto max_angl = basis.getMolecularMaxAngularMomentum(molecule);

    std::vector<double> phi;

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

            auto ncomp = sphmom.getNumberOfFactors(isph);
            
            for (int32_t icomp = 0; icomp < ncomp; icomp++)
            {
                fcarts.push_back(sphmom.getFactors(isph)[icomp]);

                auto cartind = sphmom.getIndexes(isph)[icomp];

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

                auto idelem = molecule.getIdsElemental()[atomidx];

                auto nao = basis.getNumberOfBasisFunctions(idelem, angl);

                auto basisfunc = basis.getBasisFunctions(idelem, angl);

                for (int32_t i = 0; i < nao; i++, aoidx++)
                {
                    double phiao = 0.0;

                    // process primitives

                    auto nprims = basisfunc[i].getNumberOfPrimitiveFunctions();

                    auto exponents = basisfunc[i].getExponents();

                    auto normcoefs = basisfunc[i].getNormalizationFactors();

                    for (int32_t iprim = 0; iprim < nprims; iprim++)
                    {
                        double expon = std::exp(-exponents[iprim] * r2);

                        double coef1 = normcoefs[iprim];

                        // transform from Cartesian to spherical harmonics

                        for (int32_t icomp = 0; icomp < ncomp; icomp++)
                        {
                            double coef2 = coef1 * fcarts[icomp];

                            double powxyz = std::pow(rx, lx[icomp])
                                          * std::pow(ry, ly[icomp])
                                          * std::pow(rz, lz[icomp]);

                            phiao += coef2 * powxyz * expon;
                        }
                    }

                    phi.push_back(phiao);
                }
            }
        }
    }

    return phi;
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
    auto phi = getPhiAO(molecule, basis, xp, yp, zp);

    double psi = 0.0;

    auto mocoefs = molorb.alphaOrbitals();

    auto mocols = molorb.getNumberOfColumns();

    int32_t nao = (int32_t)(phi.size());

    for (int32_t aoidx = 0; aoidx < nao; aoidx++)
    {
        double mocoef = mocoefs[aoidx * mocols + moidx];

        psi += mocoef * phi[aoidx];
    }

    return psi;
}

double
cubes::getPsiDensity(const CMolecule&        molecule,
                     const CMolecularBasis&  basis,
                     const CAODensityMatrix& density,
                     const int32_t           densityIndex,
                     const double            xp,
                     const double            yp,
                     const double            zp)
{
    auto phi = getPhiAO(molecule, basis, xp, yp, zp);

    double psi = 0.0;

    auto rho = density.totalDensity(densityIndex);

    int32_t nao = (int32_t)(phi.size());

    for (int32_t iao = 0; iao < nao; iao++)
    {
        for (int32_t jao = 0; jao < nao; jao++)
        {
            psi += phi[iao] * rho[iao * nao + jao] * phi[jao];
        }
    }

    return psi;
}
