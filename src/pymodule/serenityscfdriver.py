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

from contextlib import nullcontext
from mpi4py import MPI
import hashlib
import os
import sys
import tempfile
import uuid
import numpy as np

from .veloxchemlib import mpi_master
from .outputstream import OutputStream
from .errorhandler import assert_msg_critical

from pathlib import Path
from .molecularbasis import MolecularBasis
from .resultsio import create_hdf5, write_scf_results_to_hdf5
from .inputparser import write_unparsed_input_to_hdf5

from .firstorderprop import FirstOrderProperties
from .oneeints import compute_overlap_integrals
from .resultsio import create_hdf5, write_scf_results_to_hdf5, write_results_to_hdf5

try:
    from qcserenity import serenipy as spy
    import qcserenity as qc
except ImportError:
    pass


class SerenityScfDriver:
    """
    Implements Serenity SCF driver.

    :param comm:
        The MPI communicator.
    :param ostream:
        The output stream.

    Instance variables
        - method: Electronic-structure method (`hf` or `dft`).
        - scf_mode: SCF mode (`restricted`, `unrestricted`, or `auto`).
        - basis: Basis set label for Serenity.
        - dft_functional: DFT functional label for Serenity.
        - scratch_dir: Base scratch directory for Serenity files.
        - serenity_verbose: Print Serenity output directly to stdout.
    """

    def __init__(self, comm=None, ostream=None):
        """
        Initializes Serenity SCF driver.
        """

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

        self.method = 'hf'
        self.rohf_type = None
        self.scf_mode = 'auto'
        self.basis = '6-31GS'
        self.dft_functional = 'bp86'
        
        self.dispersion = None

        #custom functional part
        self.basic_functional = None
        self.mixing_factors = None
        self.hfexchange_ratio = None
        self.lrexchange_ratio = None
        self.mu = None
        self.impl = 'xcfun'

        self.densfit_j = 'none'
        self.grid_accuracy = 7
        self.small_grid_accuracy = 7

        self.scratch_dir = None
        self.serenity_verbose = False

        self._system = None # This is the systemController object from Serenity
        self._scf_task = None # This is the Task object that is executing the given Task here: SCFTask
        self._fat_task = None
        self._gradient_task = None # Uses gradient task object to perform the scf calculation

        # Variables that store information about the SCF 
        # !!! Currently needs to be reseted within the dyanmics as BasisFunctionOnGridController is giving an Error !!!
        self._system_signature = None 
        self._active_geom_signature = None
        self._last_scf_geom_signature = None
        self._last_grad_geom_signature = None
        
        self.max_cycles = 1000

        self._energy = None
        self._gradient = None
        self._scf_results = None
        self._current_scf_mode = None
        
        # h5 part of the file
        self.filename = None
        self._skip_writing_h5 = False

    @staticmethod
    def is_available():
        """
        Returns if Serenity python driver is available.
        """

        return ('qcserenity' in sys.modules and 'spy' in globals() and
                'qc' in globals())

    def set_method(self, method_label):
        """
        Sets Serenity method or DFT functional.
        """

        if self.rank != mpi_master():
            return

        label = str(method_label).strip().lower()

        if label in ('rohf', 'roks'):
            assert_msg_critical(
                False,
                f'SerenityScfDriver: "{label}" does not select an open-shell '
                'reference constraint; use set_rohf_type("CUHF") together '
                f'with set_method("{"hf" if label == "rohf" else "dft"}").')

        hf_aliases = {'hf', 'rhf', 'uhf'}
        dft_aliases = {'dft', 'rks', 'uks'}

        if label in hf_aliases:
            self.method = 'hf'
        elif label in dft_aliases:
            self.method = 'dft'
        else:
            # Interpret unknown labels as DFT functionals.
            self.method = 'dft'
            self.dft_functional = label

        self._invalidate_cache()

    def set_basis(self, basis_label):
        """
        Sets Serenity basis label.
        """

        if self.rank == mpi_master():
            self.basis = str(basis_label)
            self._invalidate_cache()

    def set_functional(self, functional_label=None, custom_functional=None):
        """
        Sets Serenity DFT functional label.
        """

        if self.rank == mpi_master():
            self.method = 'dft'
            if functional_label:
                self.dft_functional = str(functional_label).strip().lower()
            elif custom_functional:
                if 'basicfunctional' in custom_functional:
                    self.basic_functional = custom_functional.get('basicfunctional')
                else:
                    raise ValueError('Basicfunctionals missing')
                if 'mixingfactors' in custom_functional:
                    self.mixing_factors = custom_functional.get('mixingfactors')
                else:
                    raise ValueError('Missing mixingfactors for basic functionals')
                if 'hfexchangeratio' in custom_functional:
                    self.hfexchange_ratio = custom_functional.get('hfexchangeratio')
                else:
                    raise ValueError('Missing hfexchangeratio for custom functional construction')
                if 'lrexchangeratio' in custom_functional:
                    self.lrexchange_ratio = custom_functional.get('lrexchangeratio')
                else:
                    raise ValueError('Missing lrexchangeratio for custom functional construction')
                if 'mu' in custom_functional:
                    self.mu = custom_functional.get('mu')
                else:
                    raise ValueError('Missing mu for custom functional construction')

            self._invalidate_cache()

    def set_rohf_type(self, rohf_type):
        """
        Sets the open-shell reference constraint (NONE, CUHF or SUHF).

        CUHF constrains the unrestricted high-spin reference to <S^2> = S(S+1)
        exactly, which is the ROHF-equivalent reference used by OpenQP.
        """

        if self.rank != mpi_master():
            return

        if rohf_type is None:
            self.rohf_type = None
        else:
            label = str(rohf_type).strip().upper()
            assert_msg_critical(
                label in ('NONE', 'CUHF', 'SUHF'),
                f'SerenityScfDriver: Invalid rohf_type "{rohf_type}"')
            self.rohf_type = label

        self._invalidate_cache()

    def set_scf_mode(self, scf_mode):
        """
        Sets Serenity SCF mode.
        """

        if self.rank != mpi_master():
            return

        mode = str(scf_mode).strip().lower()
        if mode in ('auto',):
            self.scf_mode = 'auto'
        elif mode in ('restricted', 'rhf', 'rks'):
            self.scf_mode = 'restricted'
        elif mode in ('unrestricted', 'uhf', 'uks'):
            self.scf_mode = 'unrestricted'
        else:
            assert_msg_critical(False,
                                f'SerenityScfDriver: Invalid scf_mode "{mode}"')

        self._invalidate_cache()

    def set_scratch_dir(self, scratch_dir):
        """
        Sets base scratch directory for Serenity files.
        """

        if self.rank == mpi_master():
            path = os.path.abspath(str(scratch_dir))
            self.scratch_dir = path + os.sep if not path.endswith(
                os.sep) else path
            self._invalidate_cache()

    def get_method(self):
        """
        Gets Serenity method.
        """

        if self.rank == mpi_master():
            if self.method == 'dft':
                return f'dft:{self.dft_functional}'
            return self.method
        return None

    def get_energy(self):
        """
        Gets Serenity energy.
        """

        if self.rank == mpi_master() and self._energy is not None:
            return self._energy
        return None

    def get_gradient(self):
        """
        Gets Serenity gradient.
        """

        if self.rank == mpi_master() and self._gradient is not None:
            return self._gradient
        return None

    def compute(self, molecule, basis=None):
        """
        Performs Serenity SCF calculation.

        :param molecule:
            The molecule.

        :return:
            The Serenity SCF results on master rank, otherwise None.
        """

        errmsg = 'SerenityScfDriver: qcserenity is not available. '
        errmsg += 'Please install/build Serenity python bindings.'
        assert_msg_critical(self.is_available(), errmsg)
        if self.rank == mpi_master():
            results = self._compute_energy_master(molecule)
            energy = float(results['scf_energy'])
        else:
            results = None
            energy = None

        self._energy = self.comm.bcast(energy, root=mpi_master())

        self.write_final_hdf5(self.get_final_hdf5_file(), molecule, basis)
        
        if self.rank == mpi_master():
            return dict(results)
        return None

    def get_final_hdf5_file(self):
        if self.filename is None:
            return None
        return f'{self.filename}.h5'

    def compute_gradient(self, molecule):
        """
        Performs Serenity gradient calculation.

        :param molecule:
            The molecule.

        :return:
            A results dictionary on master rank, otherwise None.
        """

        errmsg = 'SerenityScfDriver: qcserenity is not available. '
        errmsg += 'Please install/build Serenity python bindings.'
        assert_msg_critical(self.is_available(), errmsg)

        if self.rank == mpi_master():
            self._compute_gradient_master(molecule)
            results = {
                'energy': float(self._energy),
                'gradient': self._gradient.copy(),
            }
            gradient = self._gradient
        else:
            results = None
            gradient = None

        self._gradient = self.comm.bcast(gradient, root=mpi_master())

        if self.rank == mpi_master():
            return results
        return None

    def print_title(self):
        """
        Prints title for Serenity calculation.
        """

        self.ostream.print_blank()
        self.ostream.print_header('Serenity SCF Driver')
        self.ostream.print_header(19 * '=')
        self.ostream.print_blank()
        self.ostream.print_reference('Reference:')
        self.ostream.print_reference(self.get_reference())
        self.ostream.flush()

    def get_reference(self):
        """
        Gets reference string for Serenity.
        """

        ref_str = 'J. P. Unsleber, T. Dresselhaus, K. Klahr, D. Schnieders, '
        ref_str += 'M. Bockers, D. Barton, J. Neugebauer, '
        ref_str += 'J. Comput. Chem. 2018, 39, 788-798'
        return ref_str

    def _invalidate_cache(self):
        """
        This function is crucial for resetting the full object
        as the new initalizaiton for the dynamics is necessary
        for a correct basis function integration!
        """

        # self._system = None
        self._scf_task = None
        self._gradient_task = None
        self._system_signature = None
        self._active_geom_signature = None
        self._last_scf_geom_signature = None
        self._last_grad_geom_signature = None
        self._energy = None
        self._gradient = None
        self._scf_results = None
        self._current_scf_mode = None

    def _compute_energy_master(self, molecule):
        self._ensure_system(molecule)

        geom_signature = self._sync_geometry_if_needed(molecule)

        if self._last_scf_geom_signature != geom_signature:

            with self._serenity_output_context():
                self.print_title()
                self._scf_task.run()
                self._energy = float(self._system.getEnergy())
                ao_basis = self._veloxchem_basis(molecule)
                self._scf_results = self._collect_scf_results(molecule, ao_basis)
                self._last_scf_geom_signature = geom_signature
                self._last_grad_geom_signature = None
                self._gradient = None

        return self._scf_results

    def _compute_gradient_master(self, molecule):
        
        self._compute_energy_master(molecule)

        geom_signature = self._active_geom_signature
        if self._last_grad_geom_signature != geom_signature:
            with self._serenity_output_context():
                self._gradient_task.run()
            self._gradient = np.array(self._system.getGeometry().getGradients(),
                                      dtype=float)
            self._last_grad_geom_signature = geom_signature

    # def _collect_scf_results(self):
    #     results = {'energy': float(self._energy)}

    #     mode = self._current_scf_mode
    #     try:
    #         if mode == 'restricted':
    #             es = self._system.getElectronicStructure_R()
    #             results['orbital_energies'] = np.array(es.orbEn(), dtype=float)
    #         else:
    #             es = self._system.getElectronicStructure_U()
    #             results['orbital_energies_alpha'] = np.array(es.alphaOrbEn(),
    #                                                          dtype=float)
    #             results['orbital_energies_beta'] = np.array(es.betaOrbEn(),
    #                                                         dtype=float)
    #     except Exception:
    #         # Keep SCF results minimal if optional details are unavailable.
    #         pass

    #     return results

    def _ensure_system(self, molecule):

        """
        In Sereentiy the systemController 
        """

        mode = self._get_effective_scf_mode(molecule)
        signature = self._get_system_signature(molecule, mode)

        if (self._system is not None and self._system_signature == signature):
            # Reused system: geometry may still have changed (e.g., MD steps).
            self._sync_geometry_if_needed(molecule)
            return

        self._ensure_scratch_dir()

        settings = spy.Settings()
        settings.name = f'vlx_serenity_{uuid.uuid4().hex[:10]}'
        settings.path = self.scratch_dir
        settings.charge = int(molecule.get_charge())
        settings.spin = int(molecule.get_multiplicity() - 1)
        settings.basis.label = self.basis
        settings.basis.densFitJ = self.densfit_j
        settings.grid.accuracy = self.grid_accuracy
        settings.grid.smallGridAccuracy = self.small_grid_accuracy
        settings.scf.maxCycles = self.max_cycles

        # if mode == 'restricted':
        #     settings.scfMode = spy.SCF_MODES.RESTRICTED
        # else:
        #     settings.scfMode = spy.SCF_MODES.UNRESTRICTED
        if mode == 'restricted':
            settings.scfMode = spy.SCF_MODES.RESTRICTED
            assert_msg_critical(
                self.rohf_type in (None, 'NONE'),
                'SerenityScfDriver: rohf_type requires an unrestricted '
                'open-shell reference.')
        else:
            settings.scfMode = spy.SCF_MODES.UNRESTRICTED
            if self.rohf_type is not None and self.rohf_type != 'NONE':
                assert_msg_critical(
                    settings.spin != 0,
                    'SerenityScfDriver: rohf_type requires multiplicity > 1.')
                settings.scf.rohf = self.rohf_type

        if self.method == 'hf':
            settings.method = spy.ELECTRONIC_STRUCTURE_THEORIES.HF
        else:
            settings.method = spy.ELECTRONIC_STRUCTURE_THEORIES.DFT
            settings.dft.functional = self.dft_functional
        if self.dispersion is not None:
            settings.dft.dispersion = self.dispersion

        labels = list(molecule.get_labels())
        coords = np.array(molecule.get_coordinates_in_bohr(), dtype=float)
        geometry = spy.Geometry(labels, coords)

        with self._serenity_output_context():
            self._system = spy.System(geometry, settings)
            if mode == 'restricted':
                self._scf_task = spy.ScfTask_R(self._system)
                self._gradient_task = spy.GradientTask_R([self._system], [])
            else:
                self._scf_task = spy.ScfTask_U(self._system)
                self._gradient_task = spy.GradientTask_U([self._system], [])

        self._configure_tasks()

        self._system_signature = signature
        self._active_geom_signature = self._get_geometry_signature(molecule)
        self._last_scf_geom_signature = None
        self._last_grad_geom_signature = None
        self._energy = None
        self._gradient = None
        self._scf_results = None
        self._current_scf_mode = mode

    def _sync_geometry_if_needed(self, molecule):
        """
        Synchronizes Serenity system geometry with the current molecule.

        Returns:
            The geometry signature of the input molecule.
        """

        geom_signature = self._get_geometry_signature(molecule)
        
        if self._system is None:
            return geom_signature

        if geom_signature != self._active_geom_signature:
            coords = np.array(molecule.get_coordinates_in_bohr(), dtype=float)
            with self._serenity_output_context():
                self._system.getGeometry().setCoordinates(coords)
            self._active_geom_signature = geom_signature
            self._last_scf_geom_signature = None
            self._last_grad_geom_signature = None
            self._energy = None
            self._gradient = None
            self._scf_results = None

        return geom_signature

    def _configure_tasks(self):
        if hasattr(self._scf_task, 'generalSettings'):
            self._scf_task.generalSettings.printLevel = (
                spy.GLOBAL_PRINT_LEVELS.MINIMUM)

        if hasattr(self._gradient_task, 'generalSettings'):
            self._gradient_task.generalSettings.printLevel = (
                spy.GLOBAL_PRINT_LEVELS.MINIMUM)

        self._gradient_task.settings.gradType = spy.GRADIENT_TYPES.ANALYTICAL

    def _ensure_scratch_dir(self):
        if self.scratch_dir is None:
            self.scratch_dir = tempfile.mkdtemp(prefix='vlx_serenity_') + os.sep
        elif not self.scratch_dir.endswith(os.sep):
            self.scratch_dir += os.sep

    def _get_effective_scf_mode(self, molecule=None):
        if self.scf_mode in ('restricted', 'unrestricted'):
            return self.scf_mode
        if molecule is not None and molecule.get_multiplicity() == 1:
            return 'restricted'
        return 'unrestricted'

    def _get_system_signature(self, molecule, mode):
        labels = tuple(molecule.get_labels())
        return (
            mode,
            self.method,
            self.basis.upper(),
            self.dft_functional.lower(),
            self.rohf_type or 'NONE',
            int(molecule.get_charge()),
            int(molecule.get_multiplicity()),
            labels,
        )

    @staticmethod
    def _get_geometry_signature(molecule):
        coords = np.ascontiguousarray(molecule.get_coordinates_in_bohr(),
                                      dtype=np.float64)
        hasher = hashlib.sha1()
        hasher.update(str(coords.shape).encode('ascii'))
        hasher.update(coords.tobytes())
        return hasher.hexdigest()

    def _serenity_output_context(self):
        print(self.serenity_verbose, self.ostream.is_muted)
        if self.serenity_verbose and not self.ostream.is_muted:
            return nullcontext()
        return qc.redirectOutputToFile(os.devnull)
    
    def get_final_h5py_file(self):
        
        if self.filename is None:
            return None
        else:
            return f"{self.filename}.h5"
        
    def _veloxchem_basis(self, molecule, basis_label=None):

        if basis_label is None:
            vlx_label = self.basis
            if self.basis == "6-31GS":
                vlx_label = "6-31g*"
            return MolecularBasis.read(molecule, vlx_label)
        else:
            vlx_label = basis_label
            if self.basis == "6-31GS":
                vlx_label = "6-31g*"

            return MolecularBasis.read(molecule, vlx_label)
    
    def create_final_hdf5(self, fname, molecule, basis=None):
        if self.rank != mpi_master() or self._skip_writing_h5:
            return False
        if Path(fname).is_file():
            return True

        ao_basis = self._veloxchem_basis(molecule, basis)
        xc_label = self.dft_functional if self.method == 'dft' else 'HF'
        create_hdf5(fname, molecule, ao_basis, xc_label, '')
        return Path(fname).is_file()

    def write_final_hdf5(self, fname, molecule, basis=None):
        if fname is None:
            return False
        if not self.create_final_hdf5(fname, molecule, basis):
            return False

        if self._scf_results is not None:
            write_scf_results_to_hdf5(fname, self._scf_results)

            write_results_to_hdf5(fname, 'serenity', {
                'program': 'Serenity',
                'method': self.method,
                'scf_mode': self._current_scf_mode,
                'basis_label_serenity': self.basis,
            })

        return True

    def _collect_scf_results(self, molecule, ao_basis):
        mode = self._current_scf_mode

        if mode == 'restricted':
            es = self._system.getElectronicStructure_R()
            C_alpha = np.array(es.coeff(), dtype=float)
            C_beta = C_alpha.copy()
            E_alpha = np.array(es.orbEn(), dtype=float)
            E_beta = E_alpha.copy()
            F_alpha = np.array(es.fock(), dtype=float)
            F_beta = F_alpha.copy()
            S = np.array(es.overlap(), dtype=float)

            # Serenity totalDens is total alpha+beta density.
            D_total = np.array(es.totalDens(), dtype=float)
            D_alpha = 0.5 * D_total
            D_beta = 0.5 * D_total

            scf_type = 'restricted'

        else:
            es = self._system.getElectronicStructure_U()
            C_alpha = np.array(es.alphaCoeff(), dtype=float)
            C_beta = np.array(es.betaCoeff(), dtype=float)
            E_alpha = np.array(es.alphaOrbEn(), dtype=float)
            E_beta = np.array(es.betaOrbEn(), dtype=float)
            F_alpha = np.array(es.alphaFock(), dtype=float)
            F_beta = np.array(es.betaFock(), dtype=float)
            D_alpha = np.array(es.alphaDens(), dtype=float)
            D_beta = np.array(es.betaDens(), dtype=float)
            S = np.array(es.overlap(), dtype=float)

            scf_type = 'unrestricted'
        
        perm = self._serenity_to_veloxchem_ao_indices(molecule, ao_basis)

        S = self._reorder_serenity_ao_matrix_to_veloxchem(S, perm)
        D_alpha = self._reorder_serenity_ao_matrix_to_veloxchem(D_alpha, perm)
        D_beta = self._reorder_serenity_ao_matrix_to_veloxchem(D_beta, perm)
        F_alpha = self._reorder_serenity_ao_matrix_to_veloxchem(F_alpha, perm)
        F_beta = self._reorder_serenity_ao_matrix_to_veloxchem(F_beta, perm)

        C_alpha = self._reorder_serenity_mo_coefficients_to_veloxchem(C_alpha, perm)
        C_beta = self._reorder_serenity_mo_coefficients_to_veloxchem(C_beta, perm)

        S_vlx = compute_overlap_integrals(molecule, ao_basis)
        max_overlap_error = np.max(np.abs(S - S_vlx))

        assert_msg_critical(
            np.allclose(S, S_vlx, atol=1.0e-10, rtol=1.0e-10),
            'Serenity AO reordering failed: reordered overlap matrix does not match '
            f'VeloxChem ordering. Max error: {max_overlap_error:.3e}'
        )

        n_mo = C_alpha.shape[1]
        occ_alpha = molecule.get_aufbau_alpha_occupation(n_mo, ao_basis)
        occ_beta = molecule.get_aufbau_beta_occupation(n_mo, ao_basis)

        scf_results = {
            'eri_thresh': 1.0e-12,
            'scf_type': scf_type,
            'rohf_type': self.rohf_type or 'NONE',
            'scf_energy': float(self._energy),
            'restart': False,
            'filename': self.filename,
            'scf_history': [],
            'S': S,
            'C_alpha': C_alpha,
            'C_beta': C_beta,
            'E_alpha': E_alpha,
            'E_beta': E_beta,
            'occ_alpha': occ_alpha,
            'occ_beta': occ_beta,
            'D_alpha': D_alpha,
            'D_beta': D_beta,
            'F_alpha': F_alpha,
            'F_beta': F_beta,
            'F': (F_alpha, F_beta),
        }

        if self.method == 'dft':
            scf_results['xcfun'] = self.dft_functional.upper()

        try:
            prop = FirstOrderProperties(self.comm, self.ostream)
            prop.compute_scf_prop(molecule, ao_basis, scf_results)
            scf_results['dipole_moment'] = np.array(
                prop.get_property('dipole_moment'))
        except Exception:
            pass

        return scf_results
    
    def _serenity_to_veloxchem_ao_indices(self, molecule, ao_basis):
        serenity_order = []
    
        for atomidx in range(molecule.number_of_atoms()):
            for ang in range(ao_basis.max_angular_momentum() + 1):
                nrad = ao_basis.number_of_basis_functions([atomidx], ang)
                ncomp = 2 * ang + 1

                for radial in range(nrad):
                    for comp in range(ncomp):
                        serenity_order.append((atomidx, ang, radial, comp))

        serenity_index = {
            label: idx for idx, label in enumerate(serenity_order)
        }

        perm = []
        for ang in range(ao_basis.max_angular_momentum() + 1):
            ncomp = 2 * ang + 1

            for comp in range(ncomp):
                for atomidx in range(molecule.number_of_atoms()):
                    nrad = ao_basis.number_of_basis_functions([atomidx], ang)

                    for radial in range(nrad):
                        perm.append(serenity_index[(atomidx, ang, radial, comp)])

        nao = ao_basis.get_dimensions_of_basis()
        assert_msg_critical(
            len(perm) == nao,
            'Serenity AO reordering: inconsistent number of AO indices.'
        )

        return np.array(perm, dtype=int)


    @staticmethod
    def _reorder_serenity_ao_matrix_to_veloxchem(matrix, perm):
        matrix = np.array(matrix, dtype=float)
        return matrix[np.ix_(perm, perm)]


    @staticmethod
    def _reorder_serenity_mo_coefficients_to_veloxchem(coefficients, perm):
        coefficients = np.array(coefficients, dtype=float)
        return coefficients[perm, :]
                