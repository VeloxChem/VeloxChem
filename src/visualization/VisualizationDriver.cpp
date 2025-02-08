//
//                              VELOXCHEM
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2024 by VeloxChem developers. All rights reserved.
//
//  SPDX-License-Identifier: LGPL-3.0-or-later
//
//  This file is part of VeloxChem.
//
//  VeloxChem is free software: you can redistribute it and/or modify it under
//  the terms of the GNU Lesser General Public License as published by the Free
//  Software Foundation, either version 3 of the License, or (at your option)
//  any later version.
//
//  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
//  License for more details.
//
//  You should have received a copy of the GNU Lesser General Public License
//  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

#include "VisualizationDriver.hpp"

#include <array>
#include <cmath>
#include <cstring>

#include "BasisFunction.hpp"
#include "CubicGrid.hpp"
#include "ErrorHandler.hpp"
#include "MolecularBasis.hpp"
#include "SphericalMomentum.hpp"
#include "StringFormat.hpp"

CVisualizationDriver::CVisualizationDriver()
{
}

std::vector<std::vector<int>>
CVisualizationDriver::_buildCartesianAngularMomentum(int angl) const
{
    std::vector<std::vector<int>> lmn;

    // lexical order of Cartesian angular momentum
    // 1S: 0
    // 3P: x,y,z
    // 6D: xx,xy,xz,yy,yz,zz
    // 10F: ...

    for (int i = 0; i <= angl; i++)
    {
        int lx = angl - i;

        for (int j = 0; j <= i; j++)
        {
            int ly = i - j;

            int lz = j;

            lmn.push_back(std::vector<int>({lx, ly, lz}));
        }
    }

    return lmn;
}

std::vector<double>
CVisualizationDriver::_compPhiAtomicOrbitals(const CMolecule&       molecule,
                                             const CMolecularBasis& basis,
                                             const double           xp,
                                             const double           yp,
                                             const double           zp) const
{
    auto natoms = molecule.number_of_atoms();

    auto max_angl = basis.max_angular_momentum();

    std::vector<double> phi;

    // azimuthal quantum number: s,p,d,f,...

    for (int aoidx = 0, angl = 0; angl <= max_angl; angl++)
    {
        auto nsph = angl * 2 + 1;

        auto lmn = _buildCartesianAngularMomentum(angl);

        // magnetic quantum number: s,p-1,p0,p+1,d-2,d-1,d0,d+1,d+2,...

        for (int isph = 0; isph < nsph; isph++)
        {
            // prepare Cartesian components

            std::vector<double> fcarts, lx, ly, lz;

            std::vector<std::pair<int, double>> sphmom;

            switch (angl)
            {
                case 0:
                    sphmom = spher_mom::transformation_factors<0>(isph);
                    break;
                case 1:
                    sphmom = spher_mom::transformation_factors<1>(isph);
                    break;
                case 2:
                    sphmom = spher_mom::transformation_factors<2>(isph);
                    break;
                case 3:
                    sphmom = spher_mom::transformation_factors<3>(isph);
                    break;
                case 4:
                    sphmom = spher_mom::transformation_factors<4>(isph);
                    break;
                default:
                    sphmom = std::vector<std::pair<int, double>>();
            }

            auto ncomp = static_cast<int>(sphmom.size());

            for (int icomp = 0; icomp < ncomp; icomp++)
            {
                fcarts.push_back(sphmom[icomp].second);

                auto cartind = sphmom[icomp].first;

                lx.push_back(lmn[cartind][0]);

                ly.push_back(lmn[cartind][1]);

                lz.push_back(lmn[cartind][2]);
            }

            // go through atoms

            for (int atomidx = 0; atomidx < natoms; atomidx++)
            {
                // process coordinates

                auto rxyz = molecule.coordinates()[atomidx].coordinates();

                double rx = xp - rxyz[0];

                double ry = yp - rxyz[1];

                double rz = zp - rxyz[2];

                double r2 = rx * rx + ry * ry + rz * rz;

                // process atomic orbitals

                auto nao = basis.number_of_basis_functions(std::vector<int>({atomidx}), angl);

                auto basisfunc = basis.basis_functions(std::vector<int>({atomidx}), angl);

                for (int i = 0; i < nao; i++, aoidx++)
                {
                    double phiao = 0.0;

                    // process primitives

                    auto nprims = basisfunc[i].number_of_primitive_functions();

                    auto exponents = basisfunc[i].get_exponents();

                    auto normcoefs = basisfunc[i].get_normalization_factors();

                    for (int iprim = 0; iprim < nprims; iprim++)
                    {
                        double expon = std::exp(-exponents[iprim] * r2);

                        double coef1 = normcoefs[iprim];

                        // transform from Cartesian to spherical harmonics

                        for (int icomp = 0; icomp < ncomp; icomp++)
                        {
                            double coef2 = coef1 * fcarts[icomp];

                            double powxyz = std::pow(rx, lx[icomp]) * std::pow(ry, ly[icomp]) * std::pow(rz, lz[icomp]);

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

std::vector<std::vector<int>>
CVisualizationDriver::getAtomicOrbitalInformation(const CMolecule& molecule, const CMolecularBasis& basis) const
{
    std::vector<std::vector<int>> aoinfo;

    auto natoms = molecule.number_of_atoms();

    auto max_angl = basis.max_angular_momentum();

    // azimuthal quantum number: s,p,d,f,...

    for (int angl = 0; angl <= max_angl; angl++)
    {
        auto nsph = angl * 2 + 1;

        // magnetic quantum number: s,p-1,p0,p+1,d-2,d-1,d0,d+1,d+2,...

        for (int isph = 0; isph < nsph; isph++)
        {
            // atoms

            for (int atomidx = 0; atomidx < natoms; atomidx++)
            {
                auto idelem = molecule.identifiers()[atomidx];

                auto nao = basis.number_of_basis_functions(std::vector<int>({atomidx}), angl);

                // atomic orbitals

                for (int iao = 0; iao < nao; iao++)
                {
                    aoinfo.push_back(std::vector<int>({idelem, angl, isph, iao, atomidx}));
                }
            }
        }
    }

    return aoinfo;
}

std::vector<std::vector<int>>
CVisualizationDriver::mapAtomToAtomicOrbitals(const CMolecule& molecule, const CMolecularBasis& basis) const
{
    auto natoms = molecule.number_of_atoms();

    std::vector<std::vector<int>> atomToAO(natoms, std::vector<int>());

    auto max_angl = basis.max_angular_momentum();

    // azimuthal quantum number: s,p,d,f,...

    for (int aoidx = 0, angl = 0; angl <= max_angl; angl++)
    {
        auto nsph = angl * 2 + 1;

        // magnetic quantum number: s,p-1,p0,p+1,d-2,d-1,d0,d+1,d+2,...

        for (int isph = 0; isph < nsph; isph++)
        {
            // atoms

            for (int atomidx = 0; atomidx < natoms; atomidx++)
            {
                auto nao = basis.number_of_basis_functions(std::vector<int>({atomidx}), angl);

                // atomic orbitals

                for (int iao = 0; iao < nao; iao++, aoidx++)
                {
                    atomToAO[atomidx].push_back(aoidx);
                }
            }
        }
    }

    return atomToAO;
}

void
CVisualizationDriver::computeAtomicOrbitalForGrid(CCubicGrid& grid, const CMolecularBasis& basis, const std::vector<int>& aoinfo) const
{
    // atomic orbital information

    auto idelem = aoinfo[0];

    auto angl = aoinfo[1];

    auto isph = aoinfo[2];

    auto iao = aoinfo[3];

    auto atomidx = aoinfo[4];

    // prepare Cartesian components

    std::vector<double> fcarts, lx, ly, lz;

    std::vector<std::pair<int, double>> sphmom;

    switch (angl)
    {
        case 0:
            sphmom = spher_mom::transformation_factors<0>(isph);
            break;
        case 1:
            sphmom = spher_mom::transformation_factors<1>(isph);
            break;
        case 2:
            sphmom = spher_mom::transformation_factors<2>(isph);
            break;
        case 3:
            sphmom = spher_mom::transformation_factors<3>(isph);
            break;
        case 4:
            sphmom = spher_mom::transformation_factors<4>(isph);
            break;
        default:
            sphmom = std::vector<std::pair<int, double>>();
    }

    auto ncomp = static_cast<int>(sphmom.size());

    auto lmn = _buildCartesianAngularMomentum(angl);

    for (int icomp = 0; icomp < ncomp; icomp++)
    {
        fcarts.push_back(sphmom[icomp].second);

        auto cartind = sphmom[icomp].first;

        lx.push_back(lmn[cartind][0]);

        ly.push_back(lmn[cartind][1]);

        lz.push_back(lmn[cartind][2]);
    }

    // calculate atomic orbital on grid points

    auto origin = grid.getOrigin();

    auto stepsize = grid.getStepSize();

    auto numpoints = grid.getNumPoints();

    #pragma omp parallel for schedule(dynamic)
    for (int ix = 0; ix < numpoints[0]; ix++)
    {
        double rx = origin[0] + stepsize[0] * ix;

        int xstride = ix * numpoints[1] * numpoints[2];

        for (int iy = 0; iy < numpoints[1]; iy++)
        {
            double ry = origin[1] + stepsize[1] * iy;

            int ystride = iy * numpoints[2];

            for (int iz = 0; iz < numpoints[2]; iz++)
            {
                double rz = origin[2] + stepsize[2] * iz;

                int zstride = iz;

                // note that the AO is centered at origin

                double r2 = rx * rx + ry * ry + rz * rz;

                // process primitives in atomic orbital

                double phiao = 0.0;

                auto basisfuncs = basis.basis_functions(std::vector<int>({atomidx}), angl);

                auto nprims = basisfuncs[iao].number_of_primitive_functions();

                auto exponents = basisfuncs[iao].get_exponents();

                auto normcoefs = basisfuncs[iao].get_normalization_factors();

                for (int iprim = 0; iprim < nprims; iprim++)
                {
                    double expon = std::exp(-exponents[iprim] * r2);

                    double coef1 = normcoefs[iprim];

                    // transform from Cartesian to spherical harmonics

                    for (int icomp = 0; icomp < ncomp; icomp++)
                    {
                        double coef2 = coef1 * fcarts[icomp];

                        double powxyz = std::pow(rx, lx[icomp]) * std::pow(ry, ly[icomp]) * std::pow(rz, lz[icomp]);

                        phiao += coef2 * powxyz * expon;
                    }
                }

                grid.values()[xstride + ystride + zstride] = phiao;
            }
        }
    }
}

std::vector<std::vector<int>>
CVisualizationDriver::getCountsAndDisplacements(int nx, int nnodes) const
{
    int ave = nx / nnodes;

    int res = nx % nnodes;

    std::vector<int> counts;

    for (int i = 0; i < nnodes; i++)
    {
        if (i < res)
        {
            counts.push_back(ave + 1);
        }
        else
        {
            counts.push_back(ave);
        }
    }

    std::vector<int> displs;

    for (int i = 0, disp = 0; i < nnodes; i++)
    {
        displs.push_back(disp);

        disp += counts[i];
    }

    return std::vector<std::vector<int>>({counts, displs});
}

CCubicGrid
CVisualizationDriver::create_local_cubic_grid(const CCubicGrid& grid,
                                              const int         rank,
                                              const int         nnodes) const
{
    // grid information

    auto origin = grid.getOrigin();

    auto stepsize = grid.getStepSize();

    auto numpoints = grid.getNumPoints();

    // compute local grid on this MPI process

    auto xcntdsp = getCountsAndDisplacements(numpoints[0], nnodes);

    auto xcounts = xcntdsp[0];

    auto xdispls = xcntdsp[1];

    std::array localorigin{origin[0] + stepsize[0] * xdispls[rank], origin[1], origin[2]};

    std::array localnumpoints{xcounts[rank], numpoints[1], numpoints[2]};

    return CCubicGrid(localorigin, stepsize, localnumpoints);
}

void
CVisualizationDriver::compute_local_grid(CCubicGrid&               grid,
                                         const CMolecule&          molecule,
                                         const CMolecularBasis&    basis,
                                         const int                 nao,
                                         const int                 nmo,
                                         const double*             mocoefs,
                                         const int                 moidx) const
{
    // grid information

    auto origin = grid.getOrigin();

    auto stepsize = grid.getStepSize();

    auto numpoints = grid.getNumPoints();

    // sanity check

    std::string erridx("VisualizationDriver.compute: Invalid index of MO");

    std::string errnao("VisualizationDriver.compute: Inconsistent number of AOs");

    auto morows = nao;

    auto mocols = nmo;

    errors::assertMsgCritical(0 <= moidx && moidx < mocols, erridx);

    auto phi0 = _compPhiAtomicOrbitals(molecule, basis, origin[0], origin[1], origin[2]);

    errors::assertMsgCritical(morows == static_cast<int>(phi0.size()), errnao);

    // target MO

    #pragma omp parallel for schedule(dynamic)
    for (int ix = 0; ix < numpoints[0]; ix++)
    {
        double xp = origin[0] + stepsize[0] * ix;

        int xstride = ix * numpoints[1] * numpoints[2];

        for (int iy = 0; iy < numpoints[1]; iy++)
        {
            double yp = origin[1] + stepsize[1] * iy;

            int ystride = iy * numpoints[2];

            for (int iz = 0; iz < numpoints[2]; iz++)
            {
                double zp = origin[2] + stepsize[2] * iz;

                int zstride = iz;

                auto phi = _compPhiAtomicOrbitals(molecule, basis, xp, yp, zp);

                double psi = 0.0;

                for (int aoidx = 0; aoidx < nao; aoidx++)
                {
                    double mocoef = mocoefs[aoidx * mocols + moidx];

                    psi += mocoef * phi[aoidx];
                }

                int index = xstride + ystride + zstride;

                grid.values()[index] = psi;
            }
        }
    }
}

void
CVisualizationDriver::compute_local_grid(CCubicGrid&             grid,
                                         const CMolecule&        molecule,
                                         const CMolecularBasis&  basis,
                                         const CAODensityMatrix& density,
                                         const int               denidx,
                                         const std::string&      denspin) const
{
    // grid information

    auto origin = grid.getOrigin();

    auto stepsize = grid.getStepSize();

    auto numpoints = grid.getNumPoints();

    // sanity check

    std::string erridx("VisualizationDriver.compute: Invalid index of density matrix");

    std::string errspin("VisualizationDriver.compute: Invalid spin of density matrix");

    std::string errnao("VisualizationDriver.compute: Inconsistent number of AOs");

    auto numdens = density.getNumberOfDensityMatrices();

    errors::assertMsgCritical(0 <= denidx && denidx < numdens, erridx);

    bool alphaspin = (format::upper_case(denspin) == std::string("ALPHA"));

    bool betaspin = (format::upper_case(denspin) == std::string("BETA"));

    bool diffspin = (format::upper_case(denspin) == std::string("SPIN"));

    errors::assertMsgCritical(alphaspin || betaspin || diffspin, errspin);

    auto phi0 = _compPhiAtomicOrbitals(molecule, basis, origin[0], origin[1], origin[2]);

    const int nao = static_cast<int>(phi0.size());

    auto denrows = density.getNumberOfRows(denidx);

    auto dencols = density.getNumberOfColumns(denidx);

    errors::assertMsgCritical(denrows == nao && dencols == nao, errnao);

    // target density

    std::vector<double> rho(nao * nao);

    auto dens_alpha = density.alphaDensity(denidx);

    auto dens_beta = density.betaDensity(denidx);

    if (alphaspin)
    {
        std::memcpy(rho.data(), dens_alpha, nao * nao * sizeof(double));
    }
    else if (betaspin)
    {
        std::memcpy(rho.data(), dens_beta, nao * nao * sizeof(double));
    }
    else if (diffspin)
    {
        for (int p = 0; p < nao * nao; p++)
        {
            rho.data()[p] = dens_alpha[p] - dens_beta[p];
        }
    }

    auto rho_data = rho.data();

    // calculate densities on grid points

    #pragma omp parallel for schedule(dynamic)
    for (int ix = 0; ix < numpoints[0]; ix++)
    {
        double xp = origin[0] + stepsize[0] * ix;

        int xstride = ix * numpoints[1] * numpoints[2];

        for (int iy = 0; iy < numpoints[1]; iy++)
        {
            double yp = origin[1] + stepsize[1] * iy;

            int ystride = iy * numpoints[2];

            for (int iz = 0; iz < numpoints[2]; iz++)
            {
                double zp = origin[2] + stepsize[2] * iz;

                int zstride = iz;

                auto phi = _compPhiAtomicOrbitals(molecule, basis, xp, yp, zp);

                double psi = 0.0;

                for (int iao = 0; iao < nao; iao++)
                {
                    for (int jao = 0; jao < nao; jao++)
                    {
                        psi += phi[iao] * rho_data[iao * nao + jao] * phi[jao];
                    }
                }

                int index = xstride + ystride + zstride;

                grid.values()[index] = psi;
            }
        }
    }
}

std::vector<double>
CVisualizationDriver::getMO(const std::vector<std::vector<double>>& coords,
                            const CMolecule&                        molecule,
                            const CMolecularBasis&                  basis,
                            const int                               nao,
                            const int                               nmo,
                            const double*                           mocoefs,
                            const int                               moidx) const
{
    // sanity check

    std::string erridx("VisualizationDriver.get_mo: invalid MO index");

    errors::assertMsgCritical(0 <= moidx && moidx < nmo, erridx);

    // compute MO

    auto npoints = static_cast<int>(coords.size());

    std::vector<double> psi(npoints, 0.0);

    for (int p = 0; p < npoints; p++)
    {
        auto phi = _compPhiAtomicOrbitals(molecule, basis, coords[p][0], coords[p][1], coords[p][2]);

        for (int aoidx = 0; aoidx < nao; aoidx++)
        {
            psi[p] += mocoefs[aoidx * nmo + moidx] * phi[aoidx];
        }
    }

    return psi;
}

std::vector<double>
CVisualizationDriver::getDensity(const std::vector<std::vector<double>>& coords,
                                 const CMolecule&                        molecule,
                                 const CMolecularBasis&                  basis,
                                 const CAODensityMatrix&                 density,
                                 const std::string&                      denspin) const
{
    // sanity check

    std::string errspin("VisualizationDriver.get_density: invalid spin");

    bool alphaspin = (format::upper_case(denspin) == std::string("ALPHA"));

    bool betaspin = (format::upper_case(denspin) == std::string("BETA"));

    errors::assertMsgCritical(alphaspin || betaspin, errspin);

    std::string erridx("VisualizationDriver.get_density: multiple density matrices not supported");

    auto numdens = density.getNumberOfDensityMatrices();

    errors::assertMsgCritical(numdens == 1, erridx);

    const int denidx = 0;

    // compute density

    auto nao = density.getNumberOfRows(denidx);

    auto dens = alphaspin ? density.alphaDensity(denidx) : density.betaDensity(denidx);

    auto npoints = static_cast<int>(coords.size());

    std::vector<double> psi(npoints, 0.0);

    for (int p = 0; p < npoints; p++)
    {
        auto phi = _compPhiAtomicOrbitals(molecule, basis, coords[p][0], coords[p][1], coords[p][2]);

        for (int iao = 0; iao < nao; iao++)
        {
            for (int jao = 0; jao < nao; jao++)
            {
                psi[p] += phi[iao] * dens[iao * nao + jao] * phi[jao];
            }
        }
    }

    return psi;
}

std::vector<double>
CVisualizationDriver::getOneParticleDensity(const std::vector<std::vector<double>>& coords_1,
                                            const std::vector<std::vector<double>>& coords_2,
                                            const CMolecule&                        molecule,
                                            const CMolecularBasis&                  basis,
                                            const CAODensityMatrix&                 density,
                                            const std::string&                      spin_1,
                                            const std::string&                      spin_2) const
{
    if (format::upper_case(spin_1) != format::upper_case(spin_2))
    {
        return std::vector<double>(coords_1.size(), 0.0);
    }

    bool alphaspin = (format::upper_case(spin_1) == std::string("ALPHA"));

    // Note: getOneParticleDensity is only called by getTwoParticleDensity
    // which guarantees that density.getNumberOfDensityMatrices() == 1 and
    // we therefore use denidx == 0

    const int denidx = 0;

    // compute density

    auto nao = density.getNumberOfRows(denidx);

    auto dens = alphaspin ? density.alphaDensity(denidx) : density.betaDensity(denidx);

    auto npoints = static_cast<int>(coords_1.size());

    std::vector<double> psi(npoints, 0.0);

    for (int p = 0; p < npoints; p++)
    {
        auto phi_1 = _compPhiAtomicOrbitals(molecule, basis, coords_1[p][0], coords_1[p][1], coords_1[p][2]);

        auto phi_2 = _compPhiAtomicOrbitals(molecule, basis, coords_2[p][0], coords_2[p][1], coords_2[p][2]);

        for (int iao = 0; iao < nao; iao++)
        {
            for (int jao = 0; jao < nao; jao++)
            {
                psi[p] += phi_1[iao] * dens[iao * nao + jao] * phi_2[jao];
            }
        }
    }

    return psi;
}

std::vector<double>
CVisualizationDriver::getTwoParticleDensity(const std::vector<std::vector<double>>& coords_1,
                                            const std::vector<std::vector<double>>& coords_2,
                                            const CMolecule&                        molecule,
                                            const CMolecularBasis&                  basis,
                                            const CAODensityMatrix&                 density,
                                            const std::string&                      spin_1,
                                            const std::string&                      spin_2) const
{
    // sanity check

    std::string errspin("VisualizationDriver.get_two_particle_density: invalid spin");

    bool alphaspin_1 = (format::upper_case(spin_1) == std::string("ALPHA"));

    bool betaspin_1 = (format::upper_case(spin_1) == std::string("BETA"));

    bool alphaspin_2 = (format::upper_case(spin_2) == std::string("ALPHA"));

    bool betaspin_2 = (format::upper_case(spin_2) == std::string("BETA"));

    errors::assertMsgCritical(alphaspin_1 || betaspin_1, errspin);

    errors::assertMsgCritical(alphaspin_2 || betaspin_2, errspin);

    std::string erridx("VisualizationDriver.get_two_particle_density: multiple density matrices not supported");

    auto numdens = density.getNumberOfDensityMatrices();

    errors::assertMsgCritical(numdens == 1, erridx);

    // compute density

    auto g_11 = getOneParticleDensity(coords_1, coords_1, molecule, basis, density, spin_1, spin_1);

    auto g_22 = getOneParticleDensity(coords_2, coords_2, molecule, basis, density, spin_2, spin_2);

    auto g_12 = getOneParticleDensity(coords_1, coords_2, molecule, basis, density, spin_1, spin_2);

    auto g_21 = getOneParticleDensity(coords_2, coords_1, molecule, basis, density, spin_2, spin_1);

    auto npoints = static_cast<int>(coords_1.size());

    std::vector<double> psi(npoints, 0.0);

    for (int p = 0; p < npoints; p++)
    {
        psi[p] = g_11[p] * g_22[p] - g_12[p] * g_21[p];
    }

    return psi;
}
