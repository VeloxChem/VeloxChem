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
from .errorhandler import assert_msg_critical

try:
    from tblite.interface import Calculator as TBLiteCalculator
    _TBLITE_AVAILABLE = True
except ImportError:
    TBLiteCalculator = None
    _TBLITE_AVAILABLE = False


class XtbDriver:
    """
    Implements xTB driver using tblite backend.
    """

    def __init__(self, comm=None, ostream=None):
        if comm is None:
            comm = MPI.COMM_WORLD

        if ostream is None:
            if comm.Get_rank() == mpi_master():
                ostream = OutputStream(sys.stdout)
            else:
                ostream = OutputStream(None)

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.ostream = ostream
        self.roots = 0

        # keep the old user-facing style for compatibility
        # accepted examples: 'gfn1', 'gfn2', 'ipea1'
        self.xtb_method = 'gfn2'

        self._xtb_calc = None
        self._xtb_res = None

    @staticmethod
    def is_available():
        """
        Returns if tblite backend is available.
        """
        return _TBLITE_AVAILABLE

    def _tblite_method_string(self):
        """
        Map legacy veloxchem method labels to tblite method names.
        """
        method = self.xtb_method.strip().lower()

        mapping = {
            'gfn1': 'GFN1-xTB',
            'gfn2': 'GFN2-xTB',
            'ipea1': 'IPEA1-xTB',
            'gfn1-xtb': 'GFN1-xTB',
            'gfn2-xtb': 'GFN2-xTB',
            'ipea1-xtb': 'IPEA1-xTB',
        }

        errmsg = f'XtbDriver: Unsupported xTB method "{self.xtb_method}".'
        assert_msg_critical(method in mapping, errmsg)

        return mapping[method]

    def set_method(self, xtb_method):
        if self.rank == mpi_master():
            self.xtb_method = xtb_method

    def get_method(self):
        if self.rank == mpi_master():
            return self.xtb_method
        return None

    def get_energy(self):
        if self.rank == mpi_master() and self._xtb_res is not None:
            return self._xtb_res.get("energy")
        return None

    def get_gradient(self):
        if self.rank == mpi_master() and self._xtb_res is not None:
            return self._xtb_res.get("gradient")
        return None

    def get_dipole(self):
        if self.rank == mpi_master() and self._xtb_res is not None:
            return self._xtb_res.get("dipole")
        return None

    def compute(self, molecule):
        """
        Performs xTB calculation via tblite backend.
        """

        errmsg = 'XtbDriver: tblite-python is not available. Please install tblite-python.'
        assert_msg_critical(self.is_available(), errmsg)

        identifiers = np.array(molecule.get_identifiers(), dtype=np.int32)
        coords_in_au = np.array(molecule.get_coordinates_in_bohr(), dtype=np.float64)

        if self.rank == mpi_master():
            self.print_title()

            charge = molecule.get_charge()
            multiplicity = molecule.get_multiplicity()
            uhf = max(multiplicity - 1, 0)

            self._xtb_calc = TBLiteCalculator(
                self._tblite_method_string(),
                identifiers,
                coords_in_au,
                charge=charge,
                uhf=uhf
            )

            # tblite verbosity is numeric in the Python API
            # 0 = quiet/minimal, 1 = default printout
            if self.ostream.is_muted:
                self._xtb_calc.set("verbosity", 0)
            else:
                self._xtb_calc.set("verbosity", 1)

            self._xtb_res = self._xtb_calc.singlepoint()

            # required quantities
            energy = self._xtb_res.get("energy")
            gradient = self._xtb_res.get("gradient")

            # optional quantities: retrieve only if available in backend/result
            result_dict = self._xtb_res.dict()

            dipole = result_dict.get("dipole")
            partial_charges = result_dict.get("charges")
            bond_orders = result_dict.get("bond-orders")
            orbital_energies = result_dict.get("orbital-energies")
            orbital_occupations = result_dict.get("orbital-occupations")

            grad2 = np.sum(gradient**2, axis=1)
            rms_grad = np.sqrt(np.mean(grad2))
            max_grad = np.max(np.sqrt(grad2))

            xtb_results = {
                'energy': energy,
                'gradient': gradient,
                'max_gradient': max_grad,
                'rms_gradient': rms_grad,
                'dipole': dipole,
                'partial_charges': partial_charges,
                'bond_orders': bond_orders,
                'orbital_energies': orbital_energies,
                'orbital_occupations': orbital_occupations,
            }

            return xtb_results

        return None

    def print_title(self):
        self.ostream.print_blank()
        self.ostream.print_header('XTB Driver')
        self.ostream.print_header(12 * '=')
        self.ostream.print_blank()

        self.ostream.print_reference('Reference:')
        self.ostream.print_reference(self.get_reference())
        self.ostream.flush()

    def get_reference(self):
        ref_str = 'C. Bannwarth, E. Caldeweyher, S. Ehlert, '
        ref_str += 'A. Hansen, P. Pracht, J. Seibert, S. Spicher, S. Grimme, '
        ref_str += 'WIREs Comput. Mol. Sci., 2020, 11, e01493'
        return ref_str