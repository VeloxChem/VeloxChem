//
//                                   VELOXCHEM
//              ----------------------------------------------------
//                          An Electronic Structure Code
//
//  SPDX-License-Identifier: BSD-3-Clause
//
//  Copyright 2018-2025 VeloxChem developers
//

#include "SimdRIJGradientDriver.hpp"

#include <algorithm>
#include <numeric>
#include <string>

#include "DenseIndexFunc.hpp"
#include "ErrorHandler.hpp"
#include "SimdThreeCenterElectronRepulsionDriver.hpp"
#include "SimdThreeCenterElectronRepulsionGradientDriver.hpp"
#include "SimdTwoCenterElectronRepulsionGradientDriver.hpp"
#include "SparseTensor.hpp"
#include "TensorComponents.hpp"

auto
CSimdRIJGradientDriver::get_threshold() const -> double
{
    return _threshold;
}

auto
CSimdRIJGradientDriver::get_block_size() const -> size_t
{
    return _block_size;
}

auto
CSimdRIJGradientDriver::_omega(const std::vector<double> &fitting, const double factor) const -> CPackedMatrix
{
    const auto naux = fitting.size();

    auto omega = CPackedMatrix(naux, naux, mat_t::symmetric);

    omega.zero();

    // NOTE: two for a closed shell and a half for an open one. The energy carries
    // the metric term with a half against a Coulomb term with a whole; what the
    // factor carries on top of that is which density was fitted. The sign is taken
    // by the caller of the two-center driver rather than here.

    for (size_t p = 0; p < naux; p++)
    {
        for (size_t q = 0; q <= p; q++)
        {
            omega.data()[omega.index(p, q)] = factor * fitting[p] * fitting[q];
        }
    }

    return omega;
}

auto
CSimdRIJGradientDriver::_compute(const CMolecule           &molecule,
                                 const CMolecularBasis     &basis,
                                 const CMolecularBasis     &aux_basis,
                                 const std::vector<double> &fitting,
                                 const CPackedMatrix       &density,
                                 const double               coulomb_factor,
                                 const double               metric_factor,
                                 const std::vector<int>    &atoms,
                                 const std::vector<int>    &aux_atoms,
                                 const bool                 with_metric) const -> CPackedMatrix
{
    const auto natoms = molecule.number_of_atoms();

    const auto naux = aux_basis.dimensions_of_basis();

    errors::assertMsgCritical(fitting.size() == naux,
                              std::string("SimdRIJGradientDriver: The fitting coefficients are not of the "
                                          "auxiliary basis"));

    errors::assertMsgCritical(density.number_of_rows() == basis.dimensions_of_basis(),
                              std::string("SimdRIJGradientDriver: The density is not of the molecular basis"));

    auto gradient = CPackedMatrix(natoms, 3, mat_t::general);

    gradient.zero();

    // the atoms whose rows are filled in, as a lookup rather than a search

    auto wanted = std::vector<bool>(natoms, false);

    for (const auto iatom : atoms)
    {
        errors::assertMsgCritical((iatom >= 0) && (static_cast<size_t>(iatom) < natoms),
                                  std::string("SimdRIJGradientDriver: An atom outside the molecule was asked for"));

        wanted[static_cast<size_t>(iatom)] = true;
    }

    _three_center(gradient, molecule, basis, aux_basis, fitting, density, coulomb_factor, wanted, aux_atoms);

    // and the two-center one. The driver of it returns the sum of Omega against
    // the derivative with no sign applied, so the sign is taken here, where it is
    // known which term is being asked for.
    //
    // NOTE: the whole auxiliary basis and not a rank's share. Omega is the outer
    // product of the fitting coefficients, which every rank holds complete because
    // the fitting was solved after the ranks added their halves of its right hand
    // side. A rank which cut this term to its own auxiliary atoms would drop the
    // pairs which straddle two shares.

    if (with_metric)
    {
        const auto two_center = CSimdTwoCenterElectronRepulsionGradientDriver(_block_size);

        const auto metric_part = two_center.compute(molecule, aux_basis, _omega(fitting, metric_factor), atoms);

        for (size_t iatom = 0; iatom < natoms; iatom++)
        {
            for (size_t c = 0; c < 3; c++)
            {
                gradient.data()[gradient.index(iatom, c)] -= metric_part.at(iatom, c);
            }
        }
    }

    return gradient;
}

auto
CSimdRIJGradientDriver::compute(const CMolecule           &molecule,
                                const CMolecularBasis     &basis,
                                const CMolecularBasis     &aux_basis,
                                const std::vector<double> &fitting,
                                const CPackedMatrix       &density,
                                const std::vector<int>    &atoms,
                                const std::vector<int>    &aux_atoms,
                                const bool                 with_metric) const -> CPackedMatrix
{
    // NOTE: the factors of a closed shell. The density handed in is one spin's and
    // the fitting coefficients are of that density, so the Coulomb term carries
    // four: two from the energy being of the total density, which is twice this
    // one on each side, and two from the derivative of a product of two things
    // which are both the density.

    return _compute(molecule, basis, aux_basis, fitting, density, 4.0, 2.0, atoms, aux_atoms, with_metric);
}

auto
CSimdRIJGradientDriver::compute(const CMolecule           &molecule,
                                const CMolecularBasis     &basis,
                                const CMolecularBasis     &aux_basis,
                                const std::vector<double> &fitting,
                                const CPackedMatrix       &density) const -> CPackedMatrix
{
    auto atoms = std::vector<int>(molecule.number_of_atoms());

    std::iota(atoms.begin(), atoms.end(), 0);

    return compute(molecule, basis, aux_basis, fitting, density, atoms, {}, true);
}

auto
CSimdRIJGradientDriver::compute_open_shell(const CMolecule           &molecule,
                                           const CMolecularBasis     &basis,
                                           const CMolecularBasis     &aux_basis,
                                           const std::vector<double> &fitting,
                                           const CPackedMatrix       &density,
                                           const std::vector<int>    &atoms,
                                           const std::vector<int>    &aux_atoms,
                                           const bool                 with_metric) const -> CPackedMatrix
{
    // NOTE: the factors of an open shell. The density is the total one and the
    // coefficients are of it, so neither of the twos the closed shell picks up
    // from halving that density is here, and what is left is the plain one and the
    // half the energy carries on the metric term.

    return _compute(molecule, basis, aux_basis, fitting, density, 1.0, 0.5, atoms, aux_atoms, with_metric);
}

auto
CSimdRIJGradientDriver::compute_open_shell(const CMolecule           &molecule,
                                           const CMolecularBasis     &basis,
                                           const CMolecularBasis     &aux_basis,
                                           const std::vector<double> &fitting,
                                           const CPackedMatrix       &density) const -> CPackedMatrix
{
    auto atoms = std::vector<int>(molecule.number_of_atoms());

    std::iota(atoms.begin(), atoms.end(), 0);

    return compute_open_shell(molecule, basis, aux_basis, fitting, density, atoms, {}, true);
}

auto
CSimdRIJGradientDriver::_three_center(CPackedMatrix             &gradient,
                                      const CMolecule           &molecule,
                                      const CMolecularBasis     &basis,
                                      const CMolecularBasis     &aux_basis,
                                      const std::vector<double> &fitting,
                                      const CPackedMatrix       &density,
                                      const double               coulomb_factor,
                                      const std::vector<bool>   &wanted,
                                      const std::vector<int>    &aux_atoms) const -> void
{
    const auto nao = basis.dimensions_of_basis();

    const auto indices = denseidx::index_functions(basis);

    const auto aux_indices = denseidx::index_functions(aux_basis);

    const auto starts = denseidx::make_dense_starts(basis);

    const auto strides = denseidx::make_dense_strides(basis);

    const auto aux_starts = denseidx::make_dense_starts(aux_basis);

    const auto aux_strides = denseidx::make_dense_strides(aux_basis);

    const auto nmoms = strides.size();

    const auto aux_nmoms = aux_strides.size();

    // the density as a square, read once for every auxiliary atom below

    auto dense_d = std::vector<double>(nao * nao, 0.0);

    density.to_dense(dense_d.data());

    // every atom of the auxiliary basis, or the share the caller holds

    auto atoms = aux_atoms;

    if (atoms.empty())
    {
        atoms.resize(molecule.number_of_atoms());

        std::iota(atoms.begin(), atoms.end(), 0);
    }

    const auto eri_drv = CSimdThreeCenterElectronRepulsionDriver();

    const auto grad_drv = CSimdThreeCenterElectronRepulsionGradientDriver(_block_size);

    for (const auto iaux : atoms)
    {
        // NOTE: the pattern of this atom alone, described with the threshold the
        // calculation formed its fitting with, so the derivative indexes the same
        // atom pairs the energy was summed over.

        const auto pattern = eri_drv.make_pattern(molecule, basis, aux_basis, _threshold, {iaux});

        const auto derivative = grad_drv.compute(pattern, molecule, basis, aux_basis);

        const auto nblocks = static_cast<size_t>(pattern.number_of_blocks());

        for (size_t iblk = 0; iblk < nblocks; iblk++)
        {
            const auto &block = pattern.block(iblk);

            const auto npairs = block.number_of_pairs();

            const auto natoms = block.number_of_c_atoms();

            if ((npairs == 0) || (natoms == 0)) continue;

            const auto &a_atoms = block.a_atoms();

            const auto &b_atoms = block.b_atoms();

            const auto &c_atoms = block.c_atoms();

            const auto &a_index = indices[static_cast<size_t>(block.a_index())];

            const auto &b_index = indices[static_cast<size_t>(block.b_index())];

            const auto &c_index = aux_indices[static_cast<size_t>(block.c_index())];

            for (size_t i = 0; i < a_index.size(); i++)
            {
                for (size_t j = 0; j < b_index.size(); j++)
                {
                    for (size_t k = 0; k < c_index.size(); k++)
                    {
                        const auto [la, ia] = a_index[i];

                        const auto [lb, jb] = b_index[j];

                        const auto [lc, kc] = c_index[k];

                        // NOTE: the atom pairs this combination reaches, which is
                        // what the kernel wrote and is not the atom pairs of the
                        // block: a combination which the screening cut short holds
                        // fewer, and they are the leading ones.

                        const auto reached = block.number_of_pairs(la, ia, lb, jb, lc, kc);

                        if (reached == 0) continue;

                        const auto na = static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{la}));

                        const auto nb = static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{lb}));

                        const auto nc = static_cast<size_t>(tensor::number_of_spherical_components(std::array<int, 1>{lc}));

                        const auto ncomps = na * nb * nc;

                        const auto *values = derivative.values(iblk, la, ia, lb, jb, lc, kc);

                        const auto run = natoms * reached;

                        // NOTE: the components of a combination are the slowest
                        // index, the three of the first center then the three of
                        // the second, and the angular components sit inside them.

                        for (size_t ma = 0; ma < na; ma++)
                        {
                            for (size_t mb = 0; mb < nb; mb++)
                            {
                                for (size_t mc = 0; mc < nc; mc++)
                                {
                                    const auto angular = (ma * nb + mb) * nc + mc;

                                    for (size_t n = 0; n < natoms; n++)
                                    {
                                        const auto catom = static_cast<size_t>(c_atoms[n]);

                                        const auto q = aux_starts[catom * aux_nmoms + lc] + kc + mc * aux_strides[lc];

                                        // Gamma of this auxiliary function, which
                                        // for a Coulomb fitting is one number
                                        // times the density and is never stored.

                                        const auto cp = coulomb_factor * fitting[q];

                                        for (size_t p = 0; p < reached; p++)
                                        {
                                            const auto aatom = static_cast<size_t>(a_atoms[p]);

                                            const auto batom = static_cast<size_t>(b_atoms[p]);

                                            const auto mu = starts[aatom * nmoms + la] + ia + ma * strides[la];

                                            const auto nu = starts[batom * nmoms + lb] + jb + mb * strides[lb];

                                            auto gamma = cp * dense_d[mu * nao + nu];

                                            // NOTE: the pattern holds each unordered
                                            // pair of atoms once, so a pair of two
                                            // different atoms stands for the term of
                                            // the sum with its two orbitals the other
                                            // way round as well. A pair of one atom
                                            // holds both orders already and is not
                                            // counted twice.

                                            if (aatom != batom) gamma *= 2.0;

                                            for (size_t c = 0; c < 3; c++)
                                            {
                                                const auto at = ((c * ncomps) + angular) * run + n * reached + p;

                                                const auto bt = (((3 + c) * ncomps) + angular) * run + n * reached + p;

                                                const auto ta = gamma * values[at];

                                                const auto tb = gamma * values[bt];

                                                if (wanted[aatom]) gradient.data()[gradient.index(aatom, c)] += ta;

                                                if (wanted[batom]) gradient.data()[gradient.index(batom, c)] += tb;

                                                // NOTE: the three derivatives of
                                                // an integral sum to zero, so the
                                                // auxiliary center takes what the
                                                // other two leave.

                                                if (wanted[catom]) gradient.data()[gradient.index(catom, c)] -= ta + tb;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}
