#
#                                   VELOXCHEM
#              ----------------------------------------------------
#                          An Electronic Structure Code
#
#  SPDX-License-Identifier: BSD-3-Clause
#
#  Copyright 2018-2025 VeloxChem developers
#
#  Redistribution and use in source and binary forms, with or without modification,
#  are permitted provided that the following conditions are met:
#
#  1. Redistributions of source code must retain the above copyright notice, this
#     list of conditions and the following disclaimer.
#  2. Redistributions in binary form must reproduce the above copyright notice,
#     this list of conditions and the following disclaimer in the documentation
#     and/or other materials provided with the distribution.
#  3. Neither the name of the copyright holder nor the names of its contributors
#     may be used to endorse or promote products derived from this software without
#     specific prior written permission.
#
#  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
#  ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
#  WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
#  DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
#  FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
#  DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
#  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
#  HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
#  LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
#  OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

from mpi4py import MPI
import numpy as np
import math
import sys
import h5py

from .veloxchemlib import (mpi_master, hartree_in_wavenumber, hartree_in_ev,
                           hartree_in_inverse_nm, fine_structure_constant,
                           extinction_coefficient_from_beta, avogadro_constant,
                           bohr_in_angstrom)
from .outputstream import OutputStream
from .distributedarray import DistributedArray
from .linearsolver import LinearSolver
from .errorhandler import assert_msg_critical

try:
    import matplotlib.pyplot as plt
    from scipy.interpolate import make_interp_spline
except ImportError:
    pass


class ComplexResponseSolverBase(LinearSolver):
    """
    Shared functionality for complex linear response solvers.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes complex linear response solver to default setup.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        super().__init__(comm, ostream)

        self.a_operator = 'electric dipole'
        self.a_components = 'xyz'
        self.b_operator = 'electric dipole'
        self.b_components = 'xyz'

        self.property = None

        self.frequencies = (0,)
        self.damping = 1000.0 / hartree_in_wavenumber()
        self.max_subspace_dim = None
        self.collapse_nvec = None

        self.collapsed_subspace = False
        self.collapsed_from_dim = None
        self.collapsed_to_dim = None

        self._input_keywords['response'].update({
            'a_operator': ('str_lower', 'A operator'),
            'a_components': ('str_lower', 'Cartesian components of A operator'),
            'b_operator': ('str_lower', 'B operator'),
            'b_components': ('str_lower', 'Cartesian components of B operator'),
            'frequencies': ('seq_range', 'frequencies'),
            'damping': ('float', 'damping parameter'),
            'max_subspace_dim':
                ('int', 'maximum reduced-space dimension before collapse'),
            'collapse_nvec':
                ('int', 'number of solution states kept after collapse'),
        })

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates response and method settings in complex liner response solver.

        :param rsp_dict:
            The dictionary of response input.
        :param method_dict:
            The dictionary of method settings.
        """

        if method_dict is None:
            method_dict = {}

        super().update_settings(rsp_dict, method_dict)

    def set_cpp_property(self, prop):
        """
        Sets CPP property (absorption or ecd or ord).

        :param prop:
            The CPP property (absorption or ecd or ord).
        """

        assert_msg_critical(prop.lower() in ['absorption', 'ecd', 'ord'],
                            f'{type(self).__name__}: invalid CPP property')

        self.property = prop.lower()

        if self.property == 'absorption':
            self.a_operator = 'electric dipole'
            self.a_components = 'xyz'
            self.b_operator = 'electric dipole'
            self.b_components = 'xyz'

        elif self.property == 'ecd':
            self.a_operator = 'magnetic dipole'
            self.a_components = 'xyz'
            self.b_operator = 'linear momentum'
            self.b_components = 'xyz'
        
        elif self.property == 'ord':
            self.a_operator = 'electric dipole'
            self.a_components = 'xyz'
            self.b_operator = 'magnetic dipole'
            self.b_components = 'xyz'        

    def _preconditioning(self, precond, v_in):
        """
        Applies preconditioner to a tuple of distributed trial vectors.

        :param precond:
            The preconditioner.
        :param v_in:
            The input trial vectors.

        :return:
            A tuple of distributed trial vectors after preconditioning.
        """

        pa = precond.data[:, 0]
        pb = precond.data[:, 1]
        pc = precond.data[:, 2]
        pd = precond.data[:, 3]

        v_in_rg = v_in.data[:, 0]
        v_in_ru = v_in.data[:, 1]
        v_in_iu = v_in.data[:, 2]
        v_in_ig = v_in.data[:, 3]

        v_out_rg = pa * v_in_rg + pb * v_in_ru + pc * v_in_iu + pd * v_in_ig
        v_out_ru = pb * v_in_rg + pa * v_in_ru + pd * v_in_iu + pc * v_in_ig
        v_out_iu = pc * v_in_rg + pd * v_in_ru - pa * v_in_iu - pb * v_in_ig
        v_out_ig = pd * v_in_rg + pc * v_in_ru - pb * v_in_iu - pa * v_in_ig

        v_mat = np.hstack((
            v_out_rg.reshape(-1, 1),
            v_out_ru.reshape(-1, 1),
            v_out_iu.reshape(-1, 1),
            v_out_ig.reshape(-1, 1),
        ))

        return DistributedArray(v_mat, self.comm, distribute=False)

    def _precond_trials(self, vectors, precond):
        """
        Applies preconditioner to distributed trial vectors.

        :param vectors:
            The set of vectors.
        :param precond:
            The preconditioner.

        :return:
            The preconditioned gerade and ungerade trial vectors.
        """

        trials_ger = []
        trials_ung = []

        for (op, w), vec in vectors.items():
            v = self._preconditioning(precond[w], vec)
            norms_2 = 2.0 * v.squared_norm(axis=0)
            vn = np.sqrt(np.sum(norms_2))

            if vn > self.norm_thresh:
                norms = np.sqrt(norms_2)
                if norms[0] > self.norm_thresh:
                    trials_ger.append(v.data[:, 0])
                if norms[1] > self.norm_thresh:
                    trials_ung.append(v.data[:, 1])
                if norms[2] > self.norm_thresh:
                    trials_ung.append(v.data[:, 2])
                if norms[3] > self.norm_thresh:
                    trials_ger.append(v.data[:, 3])

        new_ger = np.array(trials_ger).T
        new_ung = np.array(trials_ung).T

        dist_new_ger = DistributedArray(new_ger, self.comm, distribute=False)
        dist_new_ung = DistributedArray(new_ung, self.comm, distribute=False)

        return dist_new_ger, dist_new_ung

    @staticmethod
    def get_full_solution_vector(solution):
        """
        Gets a full solution vector from the distributed solution.

        :param solution:
            The distributed solution as a tuple.

        :return:
            The full solution vector.
        """

        x_realger = solution.get_full_vector(0)
        x_realung = solution.get_full_vector(1)
        x_imagung = solution.get_full_vector(2)
        x_imagger = solution.get_full_vector(3)

        if solution.rank == mpi_master():
            x_real = np.hstack((x_realger, x_realger)) + np.hstack(
                (x_realung, -x_realung))
            x_imag = np.hstack((x_imagung, -x_imagung)) + np.hstack(
                (x_imagger, x_imagger))
            return x_real + 1j * x_imag
        else:
            return None

    def _reduced_space_size(self):
        """
        Gets the total reduced-space dimension across both parity sectors.
        """

        return self._dist_bger.shape(1) + self._dist_bung.shape(1)

    def _get_max_subspace_dim(self):
        """
        Gets the maximum allowed reduced-space dimension before collapse.
        """

        if self.max_subspace_dim is not None:
            return self.max_subspace_dim

        return 20 * len(self.b_components) * max(1, len(self.frequencies))

    def _get_collapse_nvec(self):
        """
        Gets the number of solution states retained during collapse.
        """

        nreq = len(self.b_components) * max(1, len(self.frequencies))

        if self.collapse_nvec is not None:
            return max(nreq, self.collapse_nvec)

        return 2 * nreq

    def _should_collapse_subspace(self):
        """
        Checks whether the reduced-space basis should be collapsed.
        """

        if self.nonlinear:
            return False

        return self._reduced_space_size() > self._get_max_subspace_dim()

    def _orthonormalize_collapsed_space(self, vectors):
        """
        Removes linear dependence and orthonormalizes a collapsed basis.

        :param vectors:
            The distributed half-sized basis vectors.

        :return:
            The orthonormalized basis.
        """

        if vectors.data.size == 0:
            return vectors

        vectors = self._remove_linear_dependence_half_size(
            vectors, self.lindep_thresh)
        vectors = self._orthogonalize_gram_schmidt_half_size(vectors)
        vectors = self._normalize_half_size(vectors)

        return vectors

    def _print_iteration(self, relative_residual_norm, xvs):
        """
        Prints information of the iteration.

        :param relative_residual_norm:
            Relative residual norms.
        :param xvs:
            A list of tuples containing operator component, frequency, and
            property.
        """

        width = 92

        output_header = '*** Iteration:   {} '.format(self._cur_iter + 1)
        output_header += '* Residuals (Max,Min): '
        output_header += '{:.2e} and {:.2e}'.format(
            max(relative_residual_norm.values()),
            min(relative_residual_norm.values()))
        self.ostream.print_header(output_header.ljust(width))
        self.ostream.print_blank()

        if (not self.nonlinear) and (self.print_level > 1):
            output_header = 'Operator:  {} ({})'.format(self.b_operator,
                                                        self.b_components)
            self.ostream.print_header(output_header.ljust(width))
            self.ostream.print_blank()

            for op, freq, xv in xvs:
                ops_label = '<<{};{}>>_{:.4f}'.format(op, op, freq)
                rel_res = relative_residual_norm[(op, freq)]
                output_iter = '{:<15s}: {:15.8f} {:15.8f}j   '.format(
                    ops_label, -xv.real, -xv.imag)
                output_iter += 'Residual Norm: {:.8f}'.format(rel_res)
                self.ostream.print_header(output_iter.ljust(width))
            self.ostream.print_blank()

        self.ostream.flush()

    @staticmethod
    def _get_molecular_mass_amu(molecule):
        """
        Gets molecular mass in atomic mass units.

        :param molecule:
            Molecule object.

        :return:
            Molecular mass in amu.
        """

        masses = molecule.get_masses()

        assert_msg_critical(
            masses is not None,
            f'{ComplexResponseSolverBase.__name__}: unable to get atomic masses'
        )

        return float(np.sum(masses))

    def _attach_molecular_metadata(self, rsp_results, molecule):
        """
        Attaches molecular metadata to response results.

        :param rsp_results:
            Response results dictionary.
        :param molecule:
            Molecule object.

        :return:
            Updated response results dictionary.
        """

        rsp_results['molecular_mass_amu'] = self._get_molecular_mass_amu(
            molecule)

        return rsp_results

    def get_spectrum(self, rsp_results, x_unit):
        """
        Gets spectrum.

        :param rsp_results:
            The dictionary containing response results.
        :param x_unit:
            The unit of x-axis.

        :return:
            A dictionary containing the spectrum.
        """

        if self.property == 'absorption':
            return self._get_absorption_spectrum(rsp_results, x_unit)

        elif self.property == 'ecd':
            return self._get_ecd_spectrum(rsp_results, x_unit)
        
        elif self.property == 'ord':
            return self._get_ord_spectrum(rsp_results, x_unit)        

        return None

    def _get_absorption_spectrum(self, rsp_results, x_unit):
        """
        Gets absorption spectrum.
        """

        assert_msg_critical(
            x_unit.lower() in ['au', 'ev', 'nm'],
            f'{type(self).__name__}.get_spectrum: x_unit should be au, ev or nm'
        )

        au2ev = hartree_in_ev()
        auxnm = 1.0 / hartree_in_inverse_nm()

        spectrum = {'x_data': [], 'y_data': []}

        if x_unit.lower() == 'au':
            spectrum['x_label'] = 'Photon energy [a.u.]'
        elif x_unit.lower() == 'ev':
            spectrum['x_label'] = 'Photon energy [eV]'
        elif x_unit.lower() == 'nm':
            spectrum['x_label'] = 'Wavelength [nm]'

        spectrum['y_label'] = 'Absorption cross-section [a.u.]'

        freqs = rsp_results['frequencies']
        rsp_funcs = rsp_results['response_functions']

        for w in freqs:
            if w == 0.0:
                continue

            if x_unit.lower() == 'au':
                spectrum['x_data'].append(w)
            elif x_unit.lower() == 'ev':
                spectrum['x_data'].append(au2ev * w)
            elif x_unit.lower() == 'nm':
                spectrum['x_data'].append(auxnm / w)

            axx = -rsp_funcs[('x', 'x', w)].imag
            ayy = -rsp_funcs[('y', 'y', w)].imag
            azz = -rsp_funcs[('z', 'z', w)].imag

            alpha_bar = (axx + ayy + azz) / 3.0
            sigma = 4.0 * math.pi * w * alpha_bar * fine_structure_constant()

            spectrum['y_data'].append(sigma)

        return spectrum

    def _get_ecd_spectrum(self, rsp_results, x_unit):
        """
        Gets circular dichroism spectrum.
        """

        assert_msg_critical(
            x_unit.lower() in ['au', 'ev', 'nm'],
            f'{type(self).__name__}.get_spectrum: x_unit should be au, ev or nm'
        )

        au2ev = hartree_in_ev()
        auxnm = 1.0 / hartree_in_inverse_nm()

        spectrum = {'x_data': [], 'y_data': []}

        if x_unit.lower() == 'au':
            spectrum['x_label'] = 'Photon energy [a.u.]'
        elif x_unit.lower() == 'ev':
            spectrum['x_label'] = 'Photon energy [eV]'
        elif x_unit.lower() == 'nm':
            spectrum['x_label'] = 'Wavelength [nm]'

        spectrum['y_label'] = 'Molar circular dichroism '
        spectrum['y_label'] += '[L mol$^{-1}$ cm$^{-1}$]'

        freqs = rsp_results['frequencies']
        rsp_funcs = rsp_results['response_functions']

        for w in freqs:
            if w == 0.0:
                continue

            if x_unit.lower() == 'au':
                spectrum['x_data'].append(w)
            elif x_unit.lower() == 'ev':
                spectrum['x_data'].append(au2ev * w)
            elif x_unit.lower() == 'nm':
                spectrum['x_data'].append(auxnm / w)

            g_xx = -rsp_funcs[('x', 'x', w)].imag / w
            g_yy = -rsp_funcs[('y', 'y', w)].imag / w
            g_zz = -rsp_funcs[('z', 'z', w)].imag / w

            beta = -(g_xx + g_yy + g_zz) / (3.0 * w)
            delta_epsilon = beta * w**2 * extinction_coefficient_from_beta()

            spectrum['y_data'].append(delta_epsilon)

        return spectrum

    def  _get_ord_spectrum(self, rsp_results, x_unit):
        """
        Gets optical rotatory dispersion spectrum.

        :param rsp_results:
            The dictionary containing response results.
        :param x_unit:
            The unit of x-axis.

        :return:
            A dictionary containing the optical rotatory dispersion spectrum.
        """

        assert_msg_critical(
            x_unit.lower() in ['au', 'ev', 'nm'],
            f'{type(self).__name__}.get_spectrum: x_unit should be au, ev or nm'
        )

        au2ev = hartree_in_ev()
        auxnm = 1.0 / hartree_in_inverse_nm()

        spectrum = {'x_data': [], 'y_data': []}

        if x_unit.lower() == 'au':
            spectrum['x_label'] = 'Photon energy [a.u.]'
        elif x_unit.lower() == 'ev':
            spectrum['x_label'] = 'Photon energy [eV]'
        elif x_unit.lower() == 'nm':
            spectrum['x_label'] = 'Wavelength [nm]'

        molecular_mass_amu = rsp_results.get('molecular_mass_amu', None)
        convert_to_specific_rotation = (
            molecular_mass_amu is not None and molecular_mass_amu > 0.0)

        if convert_to_specific_rotation:
            spectrum['y_label'] = r'Optical rotatory dispersion [10$^3$ deg dm$^{-1}$ (g cm$^{-3}$)$^{-1}$]'
            molecular_weight = float(molecular_mass_amu)
            au2wn = hartree_in_wavenumber()
            alpha_prefactor = 1.3422940570573704e-4
        else:
            spectrum['y_label'] = 'Optical rotation parameter '
            spectrum['y_label'] += r'$\beta(\omega)$ [a.u.]'

        freqs = rsp_results['frequencies']
        rsp_funcs = rsp_results['response_functions']

        for w in freqs:
            if w == 0.0:
                continue

            if x_unit.lower() == 'au':
                spectrum['x_data'].append(w)
            elif x_unit.lower() == 'ev':
                spectrum['x_data'].append(au2ev * w)
            elif x_unit.lower() == 'nm':
                spectrum['x_data'].append(auxnm / w)

            Gxx = -rsp_funcs[('x', 'x', w)].imag
            Gyy = -rsp_funcs[('y', 'y', w)].imag
            Gzz = -rsp_funcs[('z', 'z', w)].imag

            beta = (Gxx + Gyy + Gzz) / (3.0 * w)

            if convert_to_specific_rotation:
                nu_bar = au2wn * w
                specific_rotation = alpha_prefactor * (nu_bar**2) * beta
                specific_rotation /= molecular_weight
                spectrum['y_data'].append(specific_rotation / 1000.0)
            else:
                spectrum['y_data'].append(beta)

        return spectrum

    def plot(self, cpp_results, x_unit='nm', plot_scatter=True):
        """
        Plot absorption or ECD spectrum from the CPP calculation.
        """

        assert_msg_critical(
            x_unit.lower() in ['au', 'ev', 'nm'],
            f'{type(self).__name__}.plot: x_unit should be au, ev or nm')

        assert_msg_critical('matplotlib' in sys.modules,
                            'matplotlib is required.')
        assert_msg_critical('scipy' in sys.modules, 'scipy is required.')

        cpp_spec = self.get_spectrum(cpp_results, x_unit)

        if cpp_spec is None:
            print('Nothing to plot for this complex response calculation.')
            return

        if x_unit.lower() == 'nm':
            x_data = np.array(cpp_spec['x_data'][::-1])
            y_data = np.array(cpp_spec['y_data'][::-1])
        else:
            x_data = np.array(cpp_spec['x_data'])
            y_data = np.array(cpp_spec['y_data'])

        if self.property == 'absorption':
            assert_msg_critical(
                '[a.u.]' in cpp_spec['y_label'],
                f'{type(self).__name__}.plot: In valid unit in y_label')
            na = avogadro_constant()
            a_0 = bohr_in_angstrom() * 1.0e-10
            sigma_to_epsilon = a_0**2 * 10**4 * na / (np.log(10) * 10**3)
            y_data *= sigma_to_epsilon

        spl = make_interp_spline(x_data, y_data, k=3)
        x_spl = np.linspace(x_data[0], x_data[-1], x_data.size * 10)
        y_spl = spl(x_spl)

        fig, ax = plt.subplots(figsize=(8, 5))

        if x_unit.lower() == 'nm':
            ax.set_xlabel('Wavelength [nm]')
        else:
            ax.set_xlabel(f'Excitation Energy [{x_unit.lower()}]')

        y_max = np.max(np.abs(y_data))

        if self.property == 'absorption':
            ax.set_ylabel(r'$\epsilon$ [L mol$^{-1}$ cm$^{-1}$]')
            ax.set_title('Absorption Spectrum')
            ax.set_ylim(0.0, y_max * 1.1)
        elif self.property == 'ecd':
            ax.set_ylabel(r'$\Delta \epsilon$ [L mol$^{-1}$ cm$^{-1}$]')
            ax.set_title('ECD Spectrum')
            ax.set_ylim(-y_max * 1.1, y_max * 1.1)
            ax.axhline(y=0,
                       marker=',',
                       color='k',
                       linestyle='-.',
                       markersize=0,
                       linewidth=0.2)
        elif self.property == 'ord':
            ax.set_ylabel(cpp_spec['y_label'])
            ax.set_title('Optical Rotatory Dispersion')
            ax.set_ylim(-y_max * 1.1, y_max * 1.1)
            ax.axhline(y=0,
                       marker=',',
                       color='k',
                       linestyle='-.',
                       markersize=0,
                       linewidth=0.2)

        plt.plot(x_spl, y_spl, color='black', alpha=0.8, linewidth=2.0)

        if plot_scatter:
            plt.scatter(x_data, y_data, color='darkcyan', alpha=0.9, s=15)

        plt.show()

    def _print_results(self, rsp_results, ostream=None):
        """
        Prints response results to output stream.
        """

        if self.print_level > 1:
            self._print_response_functions(rsp_results, ostream)

        if self.property == 'absorption':
            self._print_absorption_results(rsp_results, ostream)
        elif self.property == 'ecd':
            self._print_ecd_results(rsp_results, ostream)
        elif self.property == 'ord':
            self._print_ord_results(rsp_results, ostream)

    def _print_response_functions(self, rsp_results, ostream=None):
        """
        Prints response functions to output stream.
        """

        if ostream is None:
            ostream = self.ostream

        width = 92

        freqs = rsp_results['frequencies']
        rsp_funcs = rsp_results['response_functions']

        title = 'Response Functions at Given Frequencies'
        ostream.print_header(title.ljust(width))
        ostream.print_header(('=' * len(title)).ljust(width))
        ostream.print_blank()

        operator_to_name = {
            'dipole': 'Dipole',
            'electric dipole': 'Dipole',
            'electric_dipole': 'Dipole',
            'quadrupole': 'Quadru',
            'electric quadrupole': 'Quadru',
            'electric_quadrupole': 'Quadru',
            'linear_momentum': 'LinMom',
            'linear momentum': 'LinMom',
            'angular_momentum': 'AngMom',
            'angular momentum': 'AngMom',
            'magnetic dipole': 'MagDip',
            'magnetic_dipole': 'MagDip',
        }
        a_name = operator_to_name[self.a_operator]
        b_name = operator_to_name[self.b_operator]

        for w in freqs:
            title = '{:<7s} {:<7s} {:>10s} {:>15s} {:>16s}'.format(
                a_name, b_name, 'Frequency', 'Real', 'Imaginary')
            ostream.print_header(title.ljust(width))
            ostream.print_header(('-' * len(title)).ljust(width))

            for a in self.a_components:
                for b in self.b_components:
                    rsp_func_val = rsp_funcs[(a, b, w)]
                    ops_label = '<<{:>3s}  ;  {:<3s}>> {:10.4f}'.format(
                        a.lower(), b.lower(), w)
                    output = '{:<15s} {:15.8f} {:15.8f}j'.format(
                        ops_label, rsp_func_val.real, rsp_func_val.imag)
                    ostream.print_header(output.ljust(width))
            ostream.print_blank()
        ostream.flush()

    def _print_absorption_results(self, rsp_results, ostream=None):
        """
        Prints absorption results to output stream.
        """

        if ostream is None:
            ostream = self.ostream

        width = 92

        spectrum = self.get_spectrum(rsp_results, 'au')

        title = 'Linear Absorption Cross-Section'
        ostream.print_header(title.ljust(width))
        ostream.print_header(('=' * len(title)).ljust(width))
        ostream.print_blank()

        freqs = rsp_results['frequencies']

        if len(freqs) == 1 and freqs[0] == 0.0:
            text = '*** No linear absorption spectrum at zero frequency.'
            ostream.print_header(text.ljust(width))
            ostream.print_blank()
            return

        title = 'Reference: '
        title += 'J. Kauczor and P. Norman, '
        title += 'J. Chem. Theory Comput. 2014, 10, 2449-2455.'
        ostream.print_header(title.ljust(width))
        ostream.print_blank()

        assert_msg_critical(
            '[a.u.]' in spectrum['x_label'],
            f'{type(self).__name__}._print_absorption_results: In valid unit in x_label'
        )
        assert_msg_critical(
            '[a.u.]' in spectrum['y_label'],
            f'{type(self).__name__}._print_absorption_results: In valid unit in y_label'
        )

        title = '{:<20s}{:<20s}{:>15s}'.format('Frequency[a.u.]',
                                               'Frequency[eV]',
                                               'sigma(w)[a.u.]')
        ostream.print_header(title.ljust(width))
        ostream.print_header(('-' * len(title)).ljust(width))

        for w, sigma in zip(spectrum['x_data'], spectrum['y_data']):
            output = '{:<20.4f}{:<20.5f}{:>13.8f}'.format(
                w, w * hartree_in_ev(), sigma)
            ostream.print_header(output.ljust(width))

        ostream.print_blank()
        ostream.flush()

    def _print_ecd_results(self, rsp_results, ostream=None):
        """
        Prints ECD results to output stream.
        """

        if ostream is None:
            ostream = self.ostream

        width = 92

        spectrum = self.get_spectrum(rsp_results, 'au')

        title = 'Circular Dichroism Spectrum'
        ostream.print_header(title.ljust(width))
        ostream.print_header(('=' * len(title)).ljust(width))
        ostream.print_blank()

        freqs = rsp_results['frequencies']

        if len(freqs) == 1 and freqs[0] == 0.0:
            text = '*** No circular dichroism spectrum at zero frequency.'
            ostream.print_header(text.ljust(width))
            ostream.print_blank()
            return

        title = 'Reference: '
        title += 'A. Jiemchooroj and P. Norman, '
        title += 'J. Chem. Phys. 126, 134102 (2007).'
        ostream.print_header(title.ljust(width))
        ostream.print_blank()

        assert_msg_critical(
            '[a.u.]' in spectrum['x_label'],
            f'{type(self).__name__}._print_ecd_results: In valid unit in x_label'
        )
        assert_msg_critical(
            r'[L mol$^{-1}$ cm$^{-1}$]' in spectrum['y_label'],
            f'{type(self).__name__}._print_ecd_results: In valid unit in y_label'
        )

        title = '{:<20s}{:<20s}{:>28s}'.format('Frequency[a.u.]',
                                               'Frequency[eV]',
                                               'Delta_epsilon[L mol^-1 cm^-1]')
        ostream.print_header(title.ljust(width))
        ostream.print_header(('-' * len(title)).ljust(width))

        for w, delta_epsilon in zip(spectrum['x_data'], spectrum['y_data']):
            output = '{:<20.4f}{:<20.5f}{:>18.8f}'.format(
                w, w * hartree_in_ev(), delta_epsilon)
            ostream.print_header(output.ljust(width))

        ostream.print_blank()
        ostream.flush()

    def _print_ord_results(self, rsp_results, ostream=None):
        """
        Prints ORD results to output stream.
        """

        if ostream is None:
            ostream = self.ostream

        width = 92

        spectrum = self.get_spectrum(rsp_results, 'au')

        title = 'Optical Rotatory Dispersion Spectrum'
        ostream.print_header(title.ljust(width))
        ostream.print_header(('=' * len(title)).ljust(width))
        ostream.print_blank()

        freqs = rsp_results['frequencies']

        if len(freqs) == 1 and freqs[0] == 0.0:
            text = '*** No optical rotatory dispersion spectrum at zero frequency.'
            ostream.print_header(text.ljust(width))
            ostream.print_blank()
            return

        y_label = spectrum.get('y_label', 'ORD[10^3 deg dm^-1 (g cm^-3)^-1]')
        title = '{:<20s}{:<20s}{:>28s}'.format(
            'Frequency[a.u.]', 'Frequency[eV]', y_label
        )
        ostream.print_header(title.ljust(width))
        ostream.print_header(('-' * len(title)).ljust(width))

        for w, ord_value in zip(spectrum['x_data'], spectrum['y_data']):
            output = '{:<20.4f}{:<20.5f}{:>18.8f}'.format(
                w, w * hartree_in_ev(), ord_value)
            ostream.print_header(output.ljust(width))

        ostream.print_blank()
        ostream.flush()

    def write_cpp_rsp_results_to_hdf5(self, fname, rsp_results):
        """
        Writes the results of a linear response calculation to HDF5 file.
        """

        if fname and isinstance(fname, str):
            hf = h5py.File(fname, 'a')

            xlabel = self.group_label + '/frequencies'
            if xlabel in hf:
                del hf[xlabel]
            hf.create_dataset(xlabel, data=rsp_results['frequencies'])

            spectrum = self.get_spectrum(rsp_results, 'au')

            if spectrum is not None:
                y_data = np.array(spectrum['y_data'])

                if self.property == 'absorption':
                    assert_msg_critical(
                        '[a.u.]' in spectrum['y_label'],
                        f'{type(self).__name__}.write_cpp_rsp_results_to_hdf5: '
                        + 'In valid unit in y_label')
                    ylabel = self.group_label + '/sigma'
                elif self.property == 'ecd':
                    ylabel = self.group_label + '/delta-epsilon'
                elif self.property == 'ord':
                    ylabel = self.group_label + '/optical-rotation'
                if ylabel in hf:
                    del hf[ylabel]
                hf.create_dataset(ylabel, data=y_data)

            hf.close()
