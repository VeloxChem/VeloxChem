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
import sys

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .distributedarray import DistributedArray
from .linearsolver import LinearSolver
from .errorhandler import assert_msg_critical


class LinearResponseSolverBase(LinearSolver):
    """
    Thin shared scaffold for restricted and unrestricted LR solvers.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes linear response solver to default setup.
        """

        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        super().__init__(comm, ostream)

        # operators and frequencies
        self.a_operator = 'electric dipole'
        self.a_components = 'xyz'
        self.b_operator = 'electric dipole'
        self.b_components = 'xyz'
        self.property = None
        self.frequencies = (0,)

        self._input_keywords['response'].update({
            'a_operator': ('str_lower', 'A operator'),
            'a_components': ('str_lower', 'Cartesian components of A operator'),
            'b_operator': ('str_lower', 'B operator'),
            'b_components': ('str_lower', 'Cartesian components of B operator'),
            'frequencies': ('seq_range', 'frequencies'),
        })

    def update_settings(self, rsp_dict, method_dict=None):
        """
        Updates response and method settings in linear response solver.

        :param rsp_dict:
            The dictionary of response input.
        :param method_dict:
            The dictionary of method settings.
        """

        if method_dict is None:
            method_dict = {}

        super().update_settings(rsp_dict, method_dict)

    def set_lr_property(self, prop):
        """
        Sets linear-response property.

        :param prop:
            The linear-response property.
        """

        assert_msg_critical(prop.lower() in ['polarizability'],
                            f'{type(self).__name__}: invalid LR property')

        self.property = prop.lower()

        if self.property == 'polarizability':
            self.a_operator = 'electric dipole'
            self.a_components = 'xyz'
            self.b_operator = 'electric dipole'
            self.b_components = 'xyz'

    @staticmethod
    def get_full_solution_vector(solution):
        """
        Gets a full solution vector from the distributed solution.

        :param solution:
            The distributed solution as a tuple.

        :return:
            The full solution vector.
        """

        x_ger = solution.get_full_vector(0)
        x_ung = solution.get_full_vector(1)

        if solution.rank == mpi_master():
            x_ger_full = np.hstack((x_ger, x_ger))
            x_ung_full = np.hstack((x_ung, -x_ung))
            return x_ger_full + x_ung_full
        return None

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

        if self.print_level > 1:
            for op, freq, xv in xvs:
                ops_label = '<<{};{}>>_{:.4f}'.format(op, op, freq)
                rel_res = relative_residual_norm[(op, freq)]
                output_iter = '{:<15s}: {:15.8f} '.format(ops_label, -xv)
                output_iter += 'Residual Norm: {:.8f}'.format(rel_res)
                self.ostream.print_header(output_iter.ljust(width))
            self.ostream.print_blank()

        self.ostream.flush()

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

        v_in_rg = v_in.data[:, 0]
        v_in_ru = v_in.data[:, 1]

        v_out_rg = pa * v_in_rg + pb * v_in_ru
        v_out_ru = pb * v_in_rg + pa * v_in_ru

        v_mat = np.hstack((
            v_out_rg.reshape(-1, 1),
            v_out_ru.reshape(-1, 1),
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
                # gerade
                if norms[0] > self.norm_thresh:
                    trials_ger.append(v.data[:, 0])
                # ungerade
                if norms[1] > self.norm_thresh:
                    trials_ung.append(v.data[:, 1])

        new_ger = np.array(trials_ger).T
        new_ung = np.array(trials_ung).T

        dist_new_ger = DistributedArray(new_ger, self.comm, distribute=False)
        dist_new_ung = DistributedArray(new_ung, self.comm, distribute=False)

        return dist_new_ger, dist_new_ung

    def _print_polarizability_results(self, rsp_funcs, ostream):
        """
        Prints polarizability to output stream.

        :param rsp_funcs:
            The response functions.
        :param ostream:
            The output stream.
        """

        width = 92

        dipole_ops = ['dipole', 'electric dipole', 'electric_dipole']

        if self.a_operator in dipole_ops and self.b_operator in dipole_ops:

            for w in self.frequencies:
                w_str = 'Polarizability (w={:.4f})'.format(w)
                ostream.print_header(w_str.ljust(width))
                ostream.print_header(('-' * len(w_str)).ljust(width))

                valstr = '{:<5s}'.format('')
                for b in self.b_components:
                    valstr += '{:>15s}'.format(b.upper())
                ostream.print_header(valstr.ljust(width))

                for a in self.a_components:
                    valstr = '{:<5s}'.format(a.upper())
                    for b in self.b_components:
                        prop = -rsp_funcs[(a, b, w)]
                        valstr += '{:15.8f}'.format(prop)
                    ostream.print_header(valstr.ljust(width))
                ostream.print_blank()

        ostream.flush()

    def print_property(self, rsp_funcs):
        """
        Prints response property to stdout.

        :param rsp_funcs:
            The response functions.
        """

        ostream = OutputStream.create_mpi_ostream(self.comm)

        if self.property == 'polarizability':
            self._print_polarizability_results(rsp_funcs, ostream)

    def _print_results(self, rsp_funcs, ostream):
        """
        Prints response results to output stream.

        :param rsp_funcs:
            The response functions.
        :param ostream:
            The output stream.
        """

        if self.property == 'polarizability':
            self._print_polarizability_results(rsp_funcs, ostream)
