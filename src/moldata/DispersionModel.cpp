//
//                           VELOXCHEM 1.0-RC
//         ----------------------------------------------------
//                     An Electronic Structure Code
//
//  Copyright © 2018-2021 by VeloxChem developers. All rights reserved.
//  Contact: https://veloxchem.org/contact
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
//
//  This file contains derivative work of dftd4 (v2.4.0):
//  Copyright © 2017-2019 Stefan Grimme, Sebastian Ehlert, Eike Caldeweyher

#include "DispersionModel.hpp"

#include <algorithm>
#include <cmath>

#include "CoordinationNumber.hpp"
#include "DenseLinearAlgebra.hpp"
#include "DispersionData.hpp"
#include "DispersionParameters.hpp"
#include "ErrorHandler.hpp"
#include "MathConst.hpp"
#include "MathFunc.hpp"
#include "Molecule.hpp"
#include "PartialCharges.hpp"

CDispersionModel::CDispersionModel()

    : _wf(6.0)

    , _g_a(3.0)

    , _g_c(2.0)
{
}

CDispersionModel::~CDispersionModel()
{
}

double
CDispersionModel::_zeta(const double a, const double c, const double qref, const double qmod)
{
    if (qmod < 0.0)
    {
        return std::exp(a);
    }
    else
    {
        return std::exp(a * (1.0 - std::exp(c * (1.0 - qref / qmod))));
    }
}

double
CDispersionModel::_dzeta(const double a, const double c, const double qref, const double qmod)
{
    if (qmod < 0.0)
    {
        return 0.0;
    }
    else
    {
        return -a * c * std::exp(c * (1.0 - qref / qmod)) * _zeta(a, c, qref, qmod) * qref / std::pow(qmod, 2);
    }
}

double
CDispersionModel::_cngw(const double wf, const double cn, const double cnref)
{
    return std::exp(-wf * std::pow(cn - cnref, 2));
}

CMemBlock<double>
CDispersionModel::_getWeights()
{
    CMemBlock<double> freq({0.000001, 0.050000, 0.100000, 0.200000, 0.300000, 0.400000, 0.500000, 0.600000, 0.700000, 0.800000, 0.900000, 1.000000,
                            1.200000, 1.400000, 1.600000, 1.800000, 2.000000, 2.500000, 3.000000, 4.000000, 5.000000, 7.500000, 10.00000});

    CMemBlock<double> weights(freq.size());

    for (int32_t i = 0; i < freq.size(); i++)
    {
        if (i > 0)
        {
            weights.data()[i] += 0.5 * (freq.data()[i] - freq.data()[i - 1]);
        }

        if (i < freq.size() - 1)
        {
            weights.data()[i] += 0.5 * (freq.data()[i + 1] - freq.data()[i]);
        }
    }

    return weights;
}

double
CDispersionModel::_getMaximum(const std::vector<std::vector<double>>& data)
{
    std::vector<double> values;

    for (int32_t i = 0; i < static_cast<int32_t>(data.size()); i++)
    {
        auto num = static_cast<int32_t>(data[i].size());

        if (num > 0)
        {
            values.push_back(mathfunc::max(data[i].data(), num));
        }
    }

    return mathfunc::max(values.data(), static_cast<int32_t>(values.size()));
}

int32_t
CDispersionModel::_getRound(const double value)
{
    return static_cast<int32_t>(std::round(value));
}

void
CDispersionModel::_initialize(const CMolecule& molecule)
{
    // get dispersion model parameters

    const double thopi = 3.0 / mathconst::getPiValue();

    auto refn = dispdata::getRefN();

    auto zeff = dispdata::getZeff();

    auto gam = dispdata::getChemicalHardness();

    auto refsys = dispdata::getRefSys();

    auto clsh = dispdata::getClsH();

    auto clsq = dispdata::getClsQ();

    auto refcn = dispdata::getRefCN();

    auto refcovcn = dispdata::getRefCovCN();

    auto ascale = dispdata::getAscale();

    auto hcount = dispdata::getHcount();

    auto alphaiw = dispdata::getAlphaiw();

    auto n_alpha = static_cast<int32_t>(alphaiw[1][0].size());

    auto max_elem = static_cast<int32_t>(refn.size());

    auto max_refn = mathfunc::max(refn.data(), max_elem);

    auto max_refn_2 = max_refn * max_refn;

    auto max_ref_cn = _getRound(_getMaximum(refcn)) + 1;

    // get molecular information

    auto natoms = molecule.getNumberOfAtoms();

    auto idselem = molecule.getIdsElemental();

    // initialize dispersion model data

    _atoms.clear();

    _nref.clear();

    _ncount.clear();

    _cn.clear();

    _q.clear();

    _alpha.clear();

    for (int32_t i = 0; i < max_elem; i++)
    {
        auto n_i = alphaiw[i].size();

        _atoms.push_back(0);

        _nref.push_back(0);

        _ncount.push_back(std::vector<int32_t>(n_i, 0));

        _cn.push_back(std::vector<double>(n_i, 0.0));

        _q.push_back(std::vector<double>(clsq[i]));

        _alpha.push_back(std::vector<std::vector<double>>(n_i, std::vector<double>(n_alpha, 0.0)));
    }

    _c6.clear();

    _c6.resize(max_elem * max_elem * max_refn_2, 0.0);

    // compute dispersion model data

    CMemBlock<double> alpha(n_alpha);

    std::vector<int32_t> cncount(max_ref_cn, 0);

    for (int32_t i = 0; i < natoms; i++)
    {
        cncount.assign(cncount.size(), 0);

        cncount[0] = 1;

        auto ia = idselem[i];

        if (_atoms[ia] == 0)
        {
            _nref[ia] = refn[ia];

            for (int32_t j = 0; j < refn[ia]; j++)
            {
                auto is = refsys[ia][j];

                auto iz = zeff[is];

                auto sscale = dispdata::getSscale(is);

                auto secaiw = dispdata::getSecaiw(is);

                for (int32_t k = 0; k < alpha.size(); k++)
                {
                    alpha.data()[k] = sscale * secaiw[k] * _zeta(_g_a, gam[is] * _g_c, iz, clsh[ia][j] + iz);
                }

                auto icn = _getRound(refcn[ia][j]);

                _cn[ia][j] = refcovcn[ia][j];

                cncount[icn] += 1;

                for (int32_t k = 0; k < alpha.size(); k++)
                {
                    _alpha[ia][j][k] = std::max(ascale[ia][j] * (alphaiw[ia][j][k] - hcount[ia][j] * alpha.data()[k]), 0.0);
                }
            }

            for (int32_t j = 0; j < refn[ia]; j++)
            {
                int32_t icn = cncount[_getRound(refcn[ia][j])];

                _ncount[ia][j] = icn * (icn + 1) / 2;
            }
        }

        _atoms[ia] += 1;
    }

    // integrate C6 coefficients

    auto weights = _getWeights();

    for (int32_t i = 1; i < max_elem; i++)
    {
        for (int32_t j = 1; j <= i; j++)
        {
            double* c6_ij = _c6.data() + (i * max_elem + j) * max_refn_2;

            double* c6_ji = _c6.data() + (j * max_elem + i) * max_refn_2;

            if ((_atoms[i] > 0) && (_atoms[j] > 0))
            {
                for (int32_t ii = 0; ii < _nref[i]; ii++)
                {
                    for (int32_t jj = 0; jj < _nref[j]; jj++)
                    {
                        for (int32_t k = 0; k < alpha.size(); k++)
                        {
                            alpha.data()[k] = _alpha[i][ii][k] * _alpha[j][jj][k];
                        }

                        double c6 = thopi * denblas::dot(alpha, weights);

                        c6_ij[ii * max_refn + jj] = c6;

                        c6_ji[jj * max_refn + ii] = c6;
                    }
                }
            }
        }
    }
}

void
CDispersionModel::_compWeightsAndCoefficients(const CMolecule& molecule, const std::vector<double>& covcn)
{
    // initialize weights and coefficients

    _ndim = 0;

    for (int32_t i = 0; i < static_cast<int32_t>(_atoms.size()); i++)
    {
        _ndim += _atoms[i] * _nref[i];
    }

    _gw.clear();

    _dgw.clear();

    _refc6.clear();

    _gw.resize(_ndim, 0.0);

    _dgw.resize(_ndim, 0.0);

    _refc6.resize(_ndim * _ndim, 0.0);

    auto refn = dispdata::getRefN();

    auto max_elem = static_cast<int32_t>(refn.size());

    auto max_refn = mathfunc::max(refn.data(), max_elem);

    auto max_refn_2 = max_refn * max_refn;

    // get molecular inforamtion

    auto idselem = molecule.getIdsElemental();

    auto natoms = molecule.getNumberOfAtoms();

    // compute weights

    std::vector<std::vector<int32_t>> itbl(max_refn, std::vector<int32_t>(natoms, 0));

    for (int32_t i = 0, k = 0; i < natoms; i++)
    {
        for (int32_t j = 0; j < _nref[idselem[i]]; j++, k++)
        {
            itbl[j][i] = k;
        }
    }

    for (int32_t i = 0; i < natoms; i++)
    {
        auto ia = idselem[i];

        double norm = 0.0;

        double dnorm = 0.0;

        for (int32_t j = 0; j < _nref[ia]; j++)
        {
            for (int32_t t = 0; t < _ncount[ia][j]; t++)
            {
                double twf = (t + 1) * _wf;

                double tgw = _cngw(twf, covcn[i], _cn[ia][j]);

                norm += tgw;

                dnorm += 2.0 * twf * (_cn[ia][j] - covcn[i]) * tgw;
            }
        }

        norm = 1.0 / norm;

        for (int32_t j = 0; j < _nref[ia]; j++)
        {
            auto k = itbl[j][i];

            double dexpw = 0.0;

            double expw = 0.0;

            for (int32_t t = 0; t < _ncount[ia][j]; t++)
            {
                double twf = (t + 1) * _wf;

                double tgw = _cngw(twf, covcn[i], _cn[ia][j]);

                expw += tgw;

                dexpw += 2.0 * twf * (_cn[ia][j] - covcn[i]) * tgw;
            }

            _gw[k] = expw * norm;

            if (_gw[k] != _gw[k])
            {
                auto max_cn_ia = mathfunc::max(_cn[ia].data(), _nref[ia]);

                if (max_cn_ia == _cn[ia][j])
                {
                    _gw[k] = 1.0;
                }
                else
                {
                    _gw[k] = 0.0;
                }
            }

            _dgw[k] = dexpw * norm - expw * dnorm * std::pow(norm, 2);

            if (_dgw[k] != _dgw[k])
            {
                _dgw[k] = 0.0;
            }
        }
    }

    // compute coefficients

    for (int32_t i = 0; i < natoms; i++)
    {
        auto ia = idselem[i];

        for (int32_t j = 0; j <= i; j++)
        {
            auto ja = idselem[j];

            const double* c6_jaia = _c6.data() + (ja * max_elem + ia) * max_refn_2;

            for (int32_t ji = 0; ji < _nref[ia]; ji++)
            {
                auto k = itbl[ji][i];

                for (int32_t jj = 0; jj < _nref[ja]; jj++)
                {
                    if ((j == i) && (jj > ji)) continue;

                    auto l = itbl[jj][j];

                    auto val = c6_jaia[jj * max_refn + ji];

                    _refc6[l * _ndim + k] = val;

                    _refc6[k * _ndim + l] = val;
                }
            }
        }
    }
}

double
CDispersionModel::_compTwoBodyContribution(const CMolecule&           molecule,
                                           const std::string&         xcLabel,
                                           const std::vector<double>& chg,
                                           const CDenseMatrix&        dqdr,
                                           const CDenseMatrix&        dcovcndr,
                                           CDenseMatrix&              gradient)
{
    // get functional parameters

    CDispersionParameters dispParams(xcLabel);

    auto param_s6 = dispParams.getS6();

    auto param_s8 = dispParams.getS8();

    auto param_s10 = dispParams.getS10();

    auto param_a1 = dispParams.getA1();

    auto param_a2 = dispParams.getA2();

    // get dispersion model parameters

    auto zeff = dispdata::getZeff();

    auto gam = dispdata::getChemicalHardness();

    auto r4r2 = dispdata::getR4R2();

    auto refn = dispdata::getRefN();

    auto max_elem = static_cast<int32_t>(refn.size());

    auto max_refn = mathfunc::max(refn.data(), max_elem);

    const double rthr_vdw = 4000.0;

    // get molecular inforamtion

    auto idselem = molecule.getIdsElemental();

    auto natoms = molecule.getNumberOfAtoms();

    auto xcoord = molecule.getCoordinatesX();

    auto ycoord = molecule.getCoordinatesY();

    auto zcoord = molecule.getCoordinatesZ();

    // compute two-body contribution to dispersion energy

    std::vector<std::vector<int32_t>> itbl(max_refn, std::vector<int32_t>(natoms, 0));

    for (int32_t i = 0, k = 0; i < natoms; i++)
    {
        for (int32_t j = 0; j < _nref[idselem[i]]; j++, k++)
        {
            itbl[j][i] = k;
        }
    }

    std::vector<double> zetavec(_ndim, 0.0);

    std::vector<double> dzvec(_ndim, 0.0);

    std::vector<double> dzdq(_ndim, 0.0);

    CDenseMatrix dc6dcn(natoms, 1);

    CDenseMatrix dc6dq(natoms, 1);

    for (int32_t i = 0; i < natoms; i++)
    {
        auto ia = idselem[i];

        auto iz = zeff[ia];

        for (int32_t j = 0; j < _nref[ia]; j++)
        {
            auto k = itbl[j][i];

            zetavec[k] = _zeta(_g_a, gam[ia] * _g_c, _q[ia][j] + iz, chg[i] + iz) * _gw[k];

            dzvec[k] = _zeta(_g_a, gam[ia] * _g_c, _q[ia][j] + iz, chg[i] + iz) * _dgw[k];

            dzdq[k] = _dzeta(_g_a, gam[ia] * _g_c, _q[ia][j] + iz, chg[i] + iz) * _gw[k];
        }
    }

    double ed = 0.0;

    std::string err_size("TwoBodyDispersion - Mismatch in gradient matrix size");

    errors::assertMsgCritical(gradient.getNumberOfRows() == 3, err_size);

    errors::assertMsgCritical(gradient.getNumberOfColumns() == natoms, err_size);

    for (int32_t i = 0; i < natoms; i++)
    {
        auto ia = idselem[i];

        for (int32_t j = 0; j < i; j++)
        {
            auto ja = idselem[j];

            std::vector<double> rij({xcoord[i] - xcoord[j], ycoord[i] - ycoord[j], zcoord[i] - zcoord[j]});

            double r2 = rij[0] * rij[0] + rij[1] * rij[1] + rij[2] * rij[2];

            if (r2 > rthr_vdw) continue;

            double r = std::sqrt(r2);

            double c6ij = 0.0;

            double dic6ij = 0.0;

            double djc6ij = 0.0;

            double dizij = 0.0;

            double djzij = 0.0;

            for (int32_t ji = 0; ji < _nref[ia]; ji++)
            {
                auto k = itbl[ji][i];

                for (int32_t jj = 0; jj < _nref[ja]; jj++)
                {
                    auto l = itbl[jj][j];

                    c6ij += zetavec[k] * zetavec[l] * _refc6[k * _ndim + l];

                    dic6ij += dzvec[k] * zetavec[l] * _refc6[k * _ndim + l];

                    djc6ij += zetavec[k] * dzvec[l] * _refc6[k * _ndim + l];

                    dizij += dzdq[k] * zetavec[l] * _refc6[k * _ndim + l];

                    djzij += zetavec[k] * dzdq[l] * _refc6[k * _ndim + l];
                }
            }

            auto r4r2ij = 3.0 * r4r2[ia] * r4r2[ja];

            auto r0 = param_a1 * std::sqrt(r4r2ij) + param_a2;

            auto oor6 = 1.0 / (std::pow(r2, 3) + std::pow(r0, 6));

            auto oor8 = 1.0 / (std::pow(r2, 4) + std::pow(r0, 8));

            auto oor10 = 1.0 / (std::pow(r2, 5) + std::pow(r0, 10));

            auto disp = param_s6 * oor6 + param_s8 * r4r2ij * oor8 + param_s10 * 49.0 / 40.0 * std::pow(r4r2ij, 2) * oor10;

            ed -= c6ij * disp;

            double door6 = -6.0 * std::pow(r2, 2) * r * std::pow(oor6, 2);

            double door8 = -8.0 * std::pow(r2, 3) * r * std::pow(oor8, 2);

            double door10 = -10.0 * std::pow(r2, 4) * r * std::pow(oor10, 2);

            double ddisp = param_s6 * door6 + param_s8 * r4r2ij * door8 + param_s10 * 49.0 / 40.0 * std::pow(r4r2ij, 2) * door10;

            for (int32_t d = 0; d < 3; d++)
            {
                gradient.values()[d * natoms + i] -= c6ij * ddisp * rij[d] / r;

                gradient.values()[d * natoms + j] += c6ij * ddisp * rij[d] / r;
            }

            dc6dq.values()[i] += dizij * disp;

            dc6dq.values()[j] += djzij * disp;

            dc6dcn.values()[i] += dic6ij * disp;

            dc6dcn.values()[j] += djc6ij * disp;
        }
    }

    auto prod_q = denblas::multAB(dqdr, dc6dq);

    auto prod_cn = denblas::multAB(dcovcndr, dc6dcn);

    auto prod_sum = denblas::addAB(prod_q, prod_cn, 1.0);

    for (int32_t dj = 0; dj < 3 * natoms; dj++)
    {
        gradient.values()[dj] -= prod_sum.values()[dj];
    }

    return ed;
}

double
CDispersionModel::_compThreeBodyContribution(const CMolecule&    molecule,
                                             const std::string&  xcLabel,
                                             const CDenseMatrix& dcovcndr,
                                             CDenseMatrix&       gradient)
{
    // get functional parameters

    CDispersionParameters dispParams(xcLabel);

    auto param_a1 = dispParams.getA1();

    auto param_a2 = dispParams.getA2();

    auto param_alp = dispParams.getAlp();

    // get dispersion model parameters

    auto zeff = dispdata::getZeff();

    auto gam = dispdata::getChemicalHardness();

    auto r4r2 = dispdata::getR4R2();

    auto refn = dispdata::getRefN();

    auto max_elem = static_cast<int32_t>(refn.size());

    auto max_refn = mathfunc::max(refn.data(), max_elem);

    const double rthr_mbd = 1600.0;

    const double six = 6.0;

    const double oth = 1.0 / 3.0;

    // get molecular inforamtion

    auto idselem = molecule.getIdsElemental();

    auto natoms = molecule.getNumberOfAtoms();

    auto xcoord = molecule.getCoordinatesX();

    auto ycoord = molecule.getCoordinatesY();

    auto zcoord = molecule.getCoordinatesZ();

    // prepare c6ab

    std::vector<std::vector<int32_t>> itbl(max_refn, std::vector<int32_t>(natoms, 0));

    for (int32_t i = 0, k = 0; i < natoms; i++)
    {
        for (int32_t j = 0; j < _nref[idselem[i]]; j++, k++)
        {
            itbl[j][i] = k;
        }
    }

    std::vector<double> zerovec(_ndim, 0.0);

    std::vector<double> dzvec(_ndim, 0.0);

    CDenseMatrix dc6dcn(natoms, 1);

    for (int32_t i = 0; i < natoms; i++)
    {
        auto ia = idselem[i];

        auto iz = zeff[ia];

        for (int32_t j = 0; j < _nref[ia]; j++)
        {
            auto k = itbl[j][i];

            zerovec[k] = _zeta(_g_a, gam[ia] * _g_c, _q[ia][j] + iz, iz) * _gw[k];

            dzvec[k] = _zeta(_g_a, gam[ia] * _g_c, _q[ia][j] + iz, iz) * _dgw[k];
        }
    }

    std::vector<double> c6ab(natoms * natoms, 0.0);

    std::vector<double> dc6ab(natoms * natoms, 0.0);

    for (int32_t i = 0; i < natoms; i++)
    {
        auto ia = idselem[i];

        for (int32_t j = 0; j <= i; j++)
        {
            auto ja = idselem[j];

            double c6ij = 0.0;

            double dic6ij = 0.0;

            double djc6ij = 0.0;

            for (int32_t ji = 0; ji < _nref[ia]; ji++)
            {
                auto k = itbl[ji][i];

                for (int32_t jj = 0; jj < _nref[ja]; jj++)
                {
                    auto l = itbl[jj][j];

                    c6ij += zerovec[k] * zerovec[l] * _refc6[k * _ndim + l];

                    dic6ij += dzvec[k] * zerovec[l] * _refc6[k * _ndim + l];

                    djc6ij += zerovec[k] * dzvec[l] * _refc6[k * _ndim + l];
                }
            }

            c6ab[i * natoms + j] = c6ij;

            dc6ab[i * natoms + j] = dic6ij;

            if (i != j)
            {
                c6ab[j * natoms + i] = c6ij;

                dc6ab[j * natoms + i] = djc6ij;
            }
        }
    }

    // compute three-body contribution to dispersion energy

    double eabc = 0.0;

    std::string err_size("ThreeBodyDispersion - Mismatch in gradient matrix size");

    errors::assertMsgCritical(gradient.getNumberOfRows() == 3, err_size);

    errors::assertMsgCritical(gradient.getNumberOfColumns() == natoms, err_size);

    for (int32_t i = 0; i < natoms; i++)
    {
        auto ia = idselem[i];

        for (int32_t j = 0; j < i; j++)
        {
            auto ja = idselem[j];

            std::vector<double> ijvec({xcoord[j] - xcoord[i], ycoord[j] - ycoord[i], zcoord[j] - zcoord[i]});

            double rij2 = std::pow(ijvec[0], 2) + std::pow(ijvec[1], 2) + std::pow(ijvec[2], 2);

            if (rij2 > rthr_mbd) continue;

            auto c6ij = c6ab[i * natoms + j];

            double cij = param_a1 * std::sqrt(3.0 * r4r2[ia] * r4r2[ja]) + param_a2;

            for (int32_t k = 0; k < j; k++)
            {
                auto ka = idselem[k];

                std::vector<double> ikvec({xcoord[k] - xcoord[i], ycoord[k] - ycoord[i], zcoord[k] - zcoord[i]});

                double rik2 = std::pow(ikvec[0], 2) + std::pow(ikvec[1], 2) + std::pow(ikvec[2], 2);

                if (rik2 > rthr_mbd) continue;

                std::vector<double> jkvec({xcoord[k] - xcoord[j], ycoord[k] - ycoord[j], zcoord[k] - zcoord[j]});

                double rjk2 = std::pow(jkvec[0], 2) + std::pow(jkvec[1], 2) + std::pow(jkvec[2], 2);

                if (rjk2 > rthr_mbd) continue;

                auto c6ik = c6ab[i * natoms + k];

                auto c6jk = c6ab[j * natoms + k];

                double c9 = -std::sqrt(c6ij * c6ik * c6jk);

                double cik = param_a1 * std::sqrt(3.0 * r4r2[ia] * r4r2[ka]) + param_a2;

                double cjk = param_a1 * std::sqrt(3.0 * r4r2[ja] * r4r2[ka]) + param_a2;

                double cijk = cij * cjk * cik;

                double rijk2 = rij2 * rjk2 * rik2;

                double rijk = std::sqrt(rijk2);

                double fdmp = 1.0 / (1.0 + six * std::pow(std::pow(cijk / rijk, oth), param_alp));

                double rijk3 = rijk * rijk2;

                double ang = 0.375 * (rij2 + rjk2 - rik2) * (rij2 - rjk2 + rik2) * (-rij2 + rjk2 + rik2) / (rijk3 * rijk2) + 1.0 / rijk3;

                double r9ijk = ang * fdmp;

                eabc -= r9ijk * c9;

                double r, dang, val;

                double dfdmp = -oth * six * param_alp * std::pow(std::pow(cijk / rijk, oth), param_alp) * std::pow(fdmp, 2);

                // derivatives w.r.t. r_ij

                r = std::sqrt(rij2);

                dang = -0.375 *
                       (std::pow(rij2, 3) + std::pow(rij2, 2) * (rjk2 + rik2) +
                        rij2 * (3.0 * std::pow(rjk2, 2) + 2.0 * rjk2 * rik2 + 3.0 * std::pow(rik2, 2)) -
                        5.0 * std::pow(rjk2 - rik2, 2) * (rjk2 + rik2)) /
                       (r * rijk3 * rijk2);

                val = (-dang * c9 * fdmp + dfdmp / r * c9 * ang) / r;

                for (int32_t d = 0; d < 3; d++)
                {
                    gradient.values()[d * natoms + i] -= val * ijvec[d];

                    gradient.values()[d * natoms + j] += val * ijvec[d];
                }

                // derivatives w.r.t. r_ik

                r = std::sqrt(rik2);

                dang = -0.375 *
                       (std::pow(rik2, 3) + std::pow(rik2, 2) * (rjk2 + rij2) +
                        rik2 * (3.0 * std::pow(rjk2, 2) + 2.0 * rjk2 * rij2 + 3.0 * std::pow(rij2, 2)) -
                        5.0 * std::pow(rjk2 - rij2, 2) * (rjk2 + rij2)) /
                       (r * rijk3 * rijk2);

                val = (-dang * c9 * fdmp + dfdmp / r * c9 * ang) / r;

                for (int32_t d = 0; d < 3; d++)
                {
                    gradient.values()[d * natoms + i] -= val * ikvec[d];

                    gradient.values()[d * natoms + k] += val * ikvec[d];
                }

                // derivatives w.r.t. r_jk

                r = std::sqrt(rjk2);

                dang = -0.375 *
                       (std::pow(rjk2, 3) + std::pow(rjk2, 2) * (rik2 + rij2) +
                        rjk2 * (3.0 * std::pow(rik2, 2) + 2.0 * rik2 * rij2 + 3.0 * std::pow(rij2, 2)) -
                        5.0 * std::pow(rik2 - rij2, 2) * (rik2 + rij2)) /
                       (r * rijk3 * rijk2);

                val = (-dang * c9 * fdmp + dfdmp / r * c9 * ang) / r;

                for (int32_t d = 0; d < 3; d++)
                {
                    gradient.values()[d * natoms + j] -= val * jkvec[d];

                    gradient.values()[d * natoms + k] += val * jkvec[d];
                }

                // CN derivative

                double dc9;

                dc9 = 0.5 * c9 * (dc6ab[i * natoms + j] / c6ij + dc6ab[i * natoms + k] / c6ik);

                dc6dcn.values()[i] += r9ijk * dc9;

                dc9 = 0.5 * c9 * (dc6ab[j * natoms + i] / c6ij + dc6ab[j * natoms + k] / c6jk);

                dc6dcn.values()[j] += r9ijk * dc9;

                dc9 = 0.5 * c9 * (dc6ab[k * natoms + i] / c6ik + dc6ab[k * natoms + j] / c6jk);

                dc6dcn.values()[k] += r9ijk * dc9;
            }
        }
    }

    auto prod_cn = denblas::multAB(dcovcndr, dc6dcn);

    for (int32_t dj = 0; dj < 3 * natoms; dj++)
    {
        gradient.values()[dj] -= prod_cn.values()[dj];
    }

    return eabc;
}

void
CDispersionModel::compute(const CMolecule& molecule, const std::string& xcLabel)
{
    _initialize(molecule);

    auto natoms = molecule.getNumberOfAtoms();

    CDenseMatrix dqdr(3 * natoms, natoms);

    auto chg = parchg::getPartialCharges(molecule, molecule.getCharge(), dqdr);

    CDenseMatrix dcovcndr(3 * natoms, natoms);

    auto covcn = coordnum::getCovalentCoordinationNumber(molecule, dcovcndr);

    _compWeightsAndCoefficients(molecule, covcn);

    CDenseMatrix gradient(3, natoms);

    auto ed = _compTwoBodyContribution(molecule, xcLabel, chg, dqdr, dcovcndr, gradient);

    auto eabc = _compThreeBodyContribution(molecule, xcLabel, dcovcndr, gradient);

    _energy = ed + eabc;

    _gradient = CDenseMatrix(gradient.transpose());
}

double
CDispersionModel::getEnergy() const
{
    return _energy;
}

CDenseMatrix
CDispersionModel::getGradient() const
{
    return _gradient;
}
