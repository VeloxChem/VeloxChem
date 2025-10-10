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
from datetime import datetime, timedelta

from .veloxchemlib import mpi_master
from .mpitask import MpiTask
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .scfrestopendriver import ScfRestrictedOpenDriver
from .mmforcefieldgenerator import MMForceFieldGenerator
from .respchargesdriver import RespChargesDriver
from .excitondriver import ExcitonModelDriver
from .rixsdriver import RixsDriver
from .numerovdriver import NumerovDriver
from .mp2driver import Mp2Driver
from .peforcefieldgenerator import PEForceFieldGenerator
from .scfgradientdriver import ScfGradientDriver
from .tddftgradientdriver import TddftGradientDriver
from .tddftorbitalresponse import TddftOrbitalResponse
from .scfhessiandriver import ScfHessianDriver
from .optimizationdriver import OptimizationDriver
from .pulsedrsp import PulsedResponse
from .rsppolarizability import Polarizability
from .rspabsorption import Absorption
from .rsplinabscross import LinearAbsorptionCrossSection
from .rspcdspec import CircularDichroismSpectrum
from .rspc6 import C6
from .rspshg import SHG
from .rsptpatransition import TpaTransition
from .rspdoublerestrans import DoubleResTransition
from .rspthreepatransition import ThreePATransition
from .rsptpa import TPA
from .tdhfhessiandriver import TdhfHessianDriver
from .polarizabilitygradient import PolarizabilityGradient
from .vibrationalanalysis import VibrationalAnalysis
#from .rspcustomproperty import CustomProperty
from .visualizationdriver import VisualizationDriver
from .trajectorydriver import TrajectoryDriver
from .xtbdriver import XtbDriver
from .xtbgradientdriver import XtbGradientDriver
from .xtbhessiandriver import XtbHessianDriver
from .cli import cli
from .errorhandler import assert_msg_critical
from .localizationdriver import LocalizationDriver


def select_scf_driver(task, scf_type):
    """
    Selects SCF driver.

    :param task:
        The MPI task.
    :param scf_type:
        The type of SCF calculation (restricted, unrestricted, or
        restricted_openshell).

    :return:
        The SCF driver object.
    """

    # check number of MPI nodes
    if task.mpi_rank == mpi_master():
        n_ao = task.ao_basis.get_dimensions_of_basis()
        assert_msg_critical(task.mpi_size == 1 or task.mpi_size <= n_ao,
                            'SCF: too many MPI processes')

    nalpha = task.molecule.number_of_alpha_electrons()
    nbeta = task.molecule.number_of_beta_electrons()

    if scf_type == 'restricted' and nalpha == nbeta:
        scf_drv = ScfRestrictedDriver(task.mpi_comm, task.ostream)
    elif scf_type == 'restricted' and nalpha != nbeta:
        scf_drv = ScfUnrestrictedDriver(task.mpi_comm, task.ostream)
    elif scf_type == 'unrestricted':
        scf_drv = ScfUnrestrictedDriver(task.mpi_comm, task.ostream)
    elif scf_type == 'restricted_openshell':
        scf_drv = ScfRestrictedOpenDriver(task.mpi_comm, task.ostream)
    else:
        assert_msg_critical(False, f'SCF: invalid scf_type {scf_type}')

    return scf_drv


def select_rsp_property(task, mol_orbs, rsp_dict, method_dict):
    """
    Selects response property.

    :param task:
        The MPI task.
    :param mol_orbs:
        The molecular orbitals.
    :param rsp_dict:
        The dictionary of response input.
    :param method_dict:
        The dictionary of method settings.

    :return:
        The response property object.
    """

    # check number of MPI nodes
    if task.mpi_rank == mpi_master():
        nocc = task.molecule.number_of_alpha_electrons()
        n_ov = nocc * (mol_orbs.number_of_mos() - nocc)
        assert_msg_critical(task.mpi_size == 1 or task.mpi_size <= n_ov,
                            'Response: too many MPI processes')

    # check property type
    if 'property' in rsp_dict:
        prop_type = rsp_dict['property'].lower()
    else:
        prop_type = 'custom'

    if prop_type in [
            'polarizability',
            'dipole polarizability',
    ]:
        rsp_prop = Polarizability(rsp_dict, method_dict)

    elif prop_type in [
            'absorption',
            'uv-vis',
            'ecd',
    ]:
        rsp_prop = Absorption(rsp_dict, method_dict)

    elif prop_type in [
            'linear absorption cross-section',
            'linear absorption (cpp)',
            'linear absorption(cpp)',
            'absorption (cpp)',
            'absorption(cpp)',
    ]:
        rsp_prop = LinearAbsorptionCrossSection(rsp_dict, method_dict)

    elif prop_type in [
            'circular dichroism spectrum',
            'circular dichroism (cpp)',
            'circular dichroism(cpp)',
            'ecd (cpp)',
            'ecd(cpp)',
    ]:
        rsp_prop = CircularDichroismSpectrum(rsp_dict, method_dict)

    elif prop_type == 'c6':
        rsp_prop = C6(rsp_dict, method_dict)

    elif prop_type == 'shg':
        rsp_prop = SHG(rsp_dict, method_dict)

    elif prop_type == 'tpa transition':
        rsp_prop = TpaTransition(rsp_dict, method_dict)

    elif prop_type == 'transition dipole moment':
        rsp_prop = DoubleResTransition(rsp_dict, method_dict)

    elif prop_type == '3pa transition':
        rsp_prop = ThreePATransition(rsp_dict, method_dict)

    elif prop_type == 'tpa':
        rsp_prop = TPA(rsp_dict, method_dict)

    # elif prop_type == 'custom':
    #     rsp_prop = CustomProperty(rsp_dict, method_dict)

    else:
        assert_msg_critical(False,
                            f'Response: invalid response property {prop_type}')

    return rsp_prop


def updated_dict_with_eri_settings(settings_dict, scf_drv):
    """
    Returns an updated dictionary with ERI settings from SCF driver.

    :param settings_dict:
        The original dictionary of settings.
    :param scf_drv:
        The SCF driver.

    :return:
        An updated dictionary with updated ERI settings.
    """

    new_dict = dict(settings_dict)

    if 'eri_thresh' not in new_dict:
        new_dict['eri_thresh'] = scf_drv.eri_thresh

    if not scf_drv.restart:
        new_dict['restart'] = 'no'

    return new_dict


def main():
    """
    Runs VeloxChem with command line arguments.
    """

    program_start_time = datetime.now()

    # Parse command line

    parser = cli()
    args = parser.parse_args()

    # MPI task

    task = MpiTask([args.input_file, args.output_file], MPI.COMM_WORLD)
    task_type = task.input_dict['jobs']['task'].lower()

    # Timelimit

    if 'maximum_hours' in task.input_dict['jobs']:
        maximum_hours = float(task.input_dict['jobs']['maximum_hours'])
        program_end_time = program_start_time + timedelta(hours=maximum_hours)
    else:
        program_end_time = None

    # Method settings
    # Note: the @pe group is added to method_dict as pe_options

    method_dict = (dict(task.input_dict['method_settings'])
                   if 'method_settings' in task.input_dict else {})
    method_dict['pe_options'] = (dict(task.input_dict['pe'])
                                 if 'pe' in task.input_dict else {})

    use_xtb = ('xtb' in method_dict)

    # Exciton model

    if task_type == 'exciton':
        exciton_dict = (dict(task.input_dict['exciton'])
                        if 'exciton' in task.input_dict else {})
        exciton_dict['program_end_time'] = program_end_time
        exciton_dict['filename'] = task.input_dict['filename']

        exciton_drv = ExcitonModelDriver(task.mpi_comm, task.ostream)
        exciton_drv.update_settings(exciton_dict, method_dict)
        exciton_drv.compute(task.molecule, task.ao_basis, task.min_basis)

    # Force field generator

    if task_type == 'mm force field':
        force_field_dict = (dict(task.input_dict['force_field'])
                            if 'force_field' in task.input_dict else {})
        force_field_dict['filename'] = task.input_dict['filename']

        resp_dict = (dict(task.input_dict['resp_charges'])
                     if 'resp_charges' in task.input_dict else {})
        resp_dict['filename'] = task.input_dict['filename']

        force_field_drv = MMForceFieldGenerator(task.mpi_comm, task.ostream)
        force_field_drv.update_settings(force_field_dict, resp_dict)
        force_field_drv.compute(task.molecule, task.ao_basis)

    # Spectrum from trajectory

    if task_type == 'trajectory':
        traj_dict = (dict(task.input_dict['trajectory'])
                     if 'trajectory' in task.input_dict else {})
        spect_dict = (dict(task.input_dict['spectrum_settings'])
                      if 'spectrum_settings' in task.input_dict else {})
        rsp_dict = (dict(task.input_dict['response'])
                    if 'response' in task.input_dict else {})

        traj_dict['filename'] = task.input_dict['filename']
        traj_dict['charges'] = task.input_dict['charges']
        traj_dict['polarizabilities'] = task.input_dict['polarizabilities']

        traj_drv = TrajectoryDriver(task.mpi_comm, task.ostream)
        traj_drv.update_settings(traj_dict, spect_dict, rsp_dict, method_dict)
        traj_drv.compute(task.molecule, task.ao_basis, task.min_basis)

    # Diatomic vibronic spectrum using Numerov

    if task_type == 'numerov':
        numerov_dict = (dict(task.input_dict['numerov'])
                        if 'numerov' in task.input_dict else {})

        scf_dict = (dict(task.input_dict['scf'])
                    if 'scf' in task.input_dict else {})

        numerov_drv = NumerovDriver(task.mpi_comm, task.ostream)
        numerov_drv.update_settings(numerov_dict, scf_dict, method_dict)
        numerov_drv.compute(task.molecule, task.ao_basis, task.min_basis)

    # Self-consistent field
    run_scf = task_type in [
        'hf', 'rhf', 'uhf', 'rohf', 'scf', 'uscf', 'roscf', 'wavefunction',
        'wave function', 'mp2', 'ump2', 'romp2', 'gradient', 'uscf_gradient',
        'hessian', 'optimize', 'response', 'pulses', 'visualization', 'loprop',
        'pe force field', 'vibrational', 'polarizability_gradient', 'rixs'
    ]

    scf_type = 'restricted'
    if task_type in ['uhf', 'uscf', 'ump2', 'uscf_gradient']:
        scf_type = 'unrestricted'
    elif task_type in ['rohf', 'roscf', 'romp2']:
        scf_type = 'restricted_openshell'

    if run_scf:
        assert_msg_critical(task.molecule.number_of_atoms(),
                            'Molecule: no atoms found in molecule')

        scf_dict = (dict(task.input_dict['scf'])
                    if 'scf' in task.input_dict else {})
        scf_dict['program_end_time'] = program_end_time
        scf_dict['filename'] = task.input_dict['filename']

        if use_xtb:
            if 'potfile' in method_dict:
                errmsg = 'XtbDriver: The \'potfile\' keyword is not supported '
                errmsg += 'in XTB calculation.'
                if task.mpi_rank == mpi_master():
                    assert_msg_critical(False, errmsg)
            xtb_drv = XtbDriver(task.mpi_comm, task.ostream)
            xtb_drv.set_method(method_dict['xtb'].lower())
            xtb_drv.xtb_verbose = True
            xtb_results = xtb_drv.compute(task.molecule)
        else:
            scf_drv = select_scf_driver(task, scf_type)
            scf_drv.update_settings(scf_dict, method_dict)
            scf_results = scf_drv.compute(task.molecule, task.ao_basis,
                                          task.min_basis)

            mol_orbs = scf_drv.molecular_orbitals
            density = scf_drv.density

            if not scf_drv.is_converged:
                return

            if (scf_drv.electric_field is not None and
                    task.molecule.get_charge() != 0):
                task.finish()
                return

    # Gradient

    if task_type in ['gradient', 'uscf_gradient']:

        grad_dict = (dict(task.input_dict['gradient'])
                     if 'gradient' in task.input_dict else {})
        grad_dict['program_end_time'] = program_end_time
        grad_dict['filename'] = task.input_dict['filename']

        run_excited_state_gradient = ('response' in task.input_dict)
        run_ground_state_gradient = (not run_excited_state_gradient)

        if run_ground_state_gradient:

            if use_xtb:
                grad_drv = XtbGradientDriver(xtb_drv)
                grad_drv.update_settings(grad_dict, method_dict)
                grad_drv.compute(task.molecule)

            else:
                grad_drv = ScfGradientDriver(scf_drv)
                grad_drv.update_settings(grad_dict, method_dict)
                grad_drv.compute(task.molecule, task.ao_basis, scf_results)

        elif run_excited_state_gradient:

            rsp_dict = dict(task.input_dict['response'])
            rsp_dict['program_end_time'] = program_end_time
            rsp_dict['filename'] = task.input_dict['filename']
            rsp_dict = updated_dict_with_eri_settings(rsp_dict, scf_drv)

            orbrsp_dict = (dict(task.input_dict['orbital_response'])
                           if 'orbital_response' in task.input_dict else {})

            assert_msg_critical(
                rsp_dict['property'].lower() in ['absorption', 'uv-vis', 'ecd'],
                'Invalid response property for gradient calculation')

            rsp_prop = select_rsp_property(task, mol_orbs, rsp_dict,
                                           method_dict)
            rsp_prop.init_driver(task.mpi_comm, task.ostream)
            rsp_prop.compute(task.molecule, task.ao_basis, scf_results)

            tddftgrad_drv = TddftGradientDriver(scf_drv)
            tddftgrad_drv.update_settings(grad_dict, rsp_dict, orbrsp_dict,
                                          method_dict)
            tddftgrad_drv.compute(task.molecule, task.ao_basis, scf_drv,
                                  rsp_prop._rsp_driver, rsp_prop._rsp_property)

    # Hessian
    # TODO reconsider keeping this after introducing vibrationalanalysis class
    if task_type == 'hessian':
        hessian_dict = (dict(task.input_dict['hessian'])
                        if 'hessian' in task.input_dict else {})

        orbrsp_dict = (dict(task.input_dict['orbital_response'])
                       if 'orbital_response' in task.input_dict else {})
        orbrsp_dict['program_end_time'] = program_end_time
        orbrsp_dict['filename'] = task.input_dict['filename']

        if use_xtb:
            hessian_drv = XtbHessianDriver(xtb_drv)
            hessian_drv.update_settings(method_dict, hessian_dict)
            hessian_drv.compute(task.molecule)

        else:
            hessian_drv = ScfHessianDriver(scf_drv)
            hessian_drv.update_settings(method_dict, hessian_dict, orbrsp_dict)
            hessian_drv.compute(task.molecule, task.ao_basis)

    # Geometry optimization

    if task_type == 'optimize':

        if 'potfile' in method_dict:
            errmsg = 'OptimizationDriver: The \'potfile\' keyword is not '
            errmsg += 'supported in geometry optimization.'
            if task.mpi_rank == mpi_master():
                assert_msg_critical(False, errmsg)

        opt_dict = (dict(task.input_dict['optimize'])
                    if 'optimize' in task.input_dict else {})
        opt_dict['filename'] = task.input_dict['filename']

        run_excited_state_gradient = ('response' in task.input_dict)
        run_ground_state_gradient = (not run_excited_state_gradient)

        if run_ground_state_gradient:

            if use_xtb:
                grad_drv = XtbGradientDriver(xtb_drv)
                opt_drv = OptimizationDriver(grad_drv)
                opt_drv.keep_files = True
                opt_drv.update_settings(opt_dict)
                opt_results = opt_drv.compute(task.molecule)

            else:
                grad_drv = ScfGradientDriver(scf_drv)
                opt_drv = OptimizationDriver(grad_drv)
                opt_drv.keep_files = True
                opt_drv.update_settings(opt_dict)
                opt_results = opt_drv.compute(task.molecule, task.ao_basis,
                                              scf_results)

        elif run_excited_state_gradient:

            grad_dict = (task.input_dict['gradient']
                         if 'gradient' in task.input_dict else {})

            rsp_dict = dict(task.input_dict['response'])
            rsp_dict['program_end_time'] = program_end_time
            rsp_dict['filename'] = task.input_dict['filename']
            rsp_dict = updated_dict_with_eri_settings(rsp_dict, scf_drv)

            orbrsp_dict = (dict(task.input_dict['orbital_response'])
                           if 'orbital_response' in task.input_dict else {})

            assert_msg_critical(
                rsp_dict['property'].lower() in ['absorption', 'uv-vis', 'ecd'],
                'Invalid response property for geometry optimization')

            rsp_prop = select_rsp_property(task, mol_orbs, rsp_dict,
                                           method_dict)
            rsp_prop.init_driver(task.mpi_comm, task.ostream)
            rsp_prop.compute(task.molecule, task.ao_basis, scf_results)

            tddftgrad_drv = TddftGradientDriver(scf_drv)
            tddftgrad_drv.update_settings(grad_dict, rsp_dict, orbrsp_dict,
                                          method_dict)

            opt_drv = OptimizationDriver(tddftgrad_drv)
            opt_drv.keep_files = True
            opt_drv.update_settings(opt_dict)
            opt_results = opt_drv.compute(task.molecule, task.ao_basis, scf_drv,
                                          rsp_prop._rsp_driver,
                                          rsp_prop._rsp_property)

    # Vibrational analysis

    if task_type == 'vibrational':

        vib_dict = (task.input_dict['vibrational']
                    if 'vibrational' in task.input_dict else {})
        vib_dict['filename'] = task.input_dict['filename']

        hessian_dict = (task.input_dict['hessian']
                        if 'hessian' in task.input_dict else {})

        polgrad_dict = (task.input_dict['polarizability_gradient']
                        if 'polarizability_gradient' in task.input_dict else {})

        orbrsp_dict = (task.input_dict['orbital_response']
                       if 'orbital_response' in task.input_dict else {})
        orbrsp_dict['program_end_time'] = program_end_time
        orbrsp_dict['filename'] = task.input_dict['filename']

        rsp_dict = (task.input_dict['response']
                    if 'response' in task.input_dict else {})
        rsp_dict['filename'] = task.input_dict['filename']

        if use_xtb:
            vibrational_drv = VibrationalAnalysis(xtb_drv)
            vibrational_drv.update_settings(method_dict,
                                            vib_dict,
                                            hessian_dict=hessian_dict,
                                            cphf_dict=orbrsp_dict,
                                            rsp_dict=rsp_dict,
                                            polgrad_dict=polgrad_dict)

        else:
            vibrational_drv = VibrationalAnalysis(scf_drv)
            vibrational_drv.update_settings(method_dict,
                                            vib_dict,
                                            hessian_dict=hessian_dict,
                                            cphf_dict=orbrsp_dict,
                                            rsp_dict=rsp_dict,
                                            polgrad_dict=polgrad_dict)

        vib_results = vibrational_drv.compute(task.molecule, task.ao_basis)

    # Polarizability gradient

    if task_type == 'polarizability_gradient':

        polgrad_dict = (task.input_dict['polarizability_gradient']
                        if 'polarizability_gradient' in task.input_dict else {})
        orbrsp_dict = (task.input_dict['orbital_response']
                       if 'orbital_response' in task.input_dict else {})
        rsp_dict = (task.input_dict['response']
                    if 'response' in task.input_dict else {})

        rsp_dict['program_end_time'] = program_end_time
        rsp_dict['filename'] = task.input_dict['filename']
        rsp_dict = updated_dict_with_eri_settings(rsp_dict, scf_drv)

        rsp_prop = select_rsp_property(task, mol_orbs, rsp_dict, method_dict)
        rsp_prop.init_driver(task.mpi_comm, task.ostream)
        rsp_prop.compute(task.molecule, task.ao_basis, scf_results)

        polgrad_drv = PolarizabilityGradient(scf_drv, task.mpi_comm,
                                             task.ostream)
        polgrad_drv.update_settings(polgrad_dict, orbrsp_dict, method_dict)
        polgrad_drv.compute(task.molecule, task.ao_basis, scf_drv.scf_tensors,
                            rsp_prop._rsp_property)

    # Response

    if (task_type == 'response' and
            scf_drv.scf_type in ['restricted', 'unrestricted']):

        rsp_dict = (dict(task.input_dict['response'])
                    if 'response' in task.input_dict else {})
        rsp_dict['program_end_time'] = program_end_time
        rsp_dict['filename'] = task.input_dict['filename']
        rsp_dict = updated_dict_with_eri_settings(rsp_dict, scf_drv)
        if 'localize_mos' in rsp_dict:
            loc_drv = LocalizationDriver(task.mpi_comm, task.ostream)
            loc_drv.localize_and_write(task.molecule, task.ao_basis, scf_drv, rsp_dict['localize_mos'], write_hdf5 = True)

        rsp_prop = select_rsp_property(task, mol_orbs, rsp_dict, method_dict)
        rsp_prop.init_driver(task.mpi_comm, task.ostream,
                             method_type=scf_drv.scf_type)
        rsp_prop.compute(task.molecule, task.ao_basis, scf_results)

        if not rsp_prop.is_converged:
            return

        # Calculate the excited-state gradient if requested

        if 'property' in rsp_dict:
            prop_type = rsp_dict['property'].lower()
        else:
            prop_type = None

        if prop_type in ['absorption', 'uv-vis', 'ecd']:

            if 'orbital_response' in task.input_dict:
                orbrsp_dict = task.input_dict['orbital_response']
                if 'tamm_dancoff' in rsp_dict:
                    orbrsp_dict['tamm_dancoff'] = rsp_dict['tamm_dancoff']
                # TODO: if gradient dict is not defined, do just the
                # orbital response. This is for testing/benchmarking
                # only and should be removed or done in a better way
                # (e.g. task_type = orbital-response)
                if 'gradient' not in task.input_dict:
                    orbrsp_drv = TddftOrbitalResponse(task.mpi_comm,
                                                      task.ostream)
                    orbrsp_drv.update_settings(orbrsp_dict, method_dict)
                    orbrsp_drv.compute(task.molecule, task.ao_basis,
                                       scf_drv.scf_tensors,
                                       rsp_prop._rsp_property)
            else:
                orbrsp_dict = {}

            # Excited state gradient
            if 'gradient' in task.input_dict:
                grad_dict = task.input_dict['gradient']
                tddftgrad_drv = TddftGradientDriver(scf_drv)
                tddftgrad_drv.update_settings(grad_dict, rsp_dict, orbrsp_dict,
                                              method_dict)
                tddftgrad_drv.compute(task.molecule, task.ao_basis, scf_drv,
                                      rsp_prop._rsp_driver,
                                      rsp_prop._rsp_property)

            # Excited state Hessian and vibrational analysis
            if 'vibrational' in task.input_dict:
                freq_dict = task.input_dict['vibrational']
                tdhfhessian_drv = TdhfHessianDriver(scf_drv)
                tdhfhessian_drv.update_settings(method_dict, rsp_dict,
                                                freq_dict, orbrsp_dict)
                tdhfhessian_drv.compute(task.molecule, task.ao_basis,
                                        rsp_prop._rsp_driver)
                if task.mpi_rank == mpi_master():
                    tdhfhessian_drv.vibrational_analysis(task.molecule)

            # Excited state optimization
            if 'optimize_excited_state' in task.input_dict:
                opt_dict = task.input_dict['optimize_excited_state']
                tddftgrad_drv = TddftGradientDriver(scf_drv)
                tddftgrad_drv.update_settings(opt_dict, rsp_dict, orbrsp_dict,
                                              method_dict)
                opt_drv = OptimizationDriver(tddftgrad_drv)
                opt_drv.compute(task.molecule, task.ao_basis, scf_drv,
                                rsp_prop._rsp_driver, rsp_prop._rsp_property)

    # Pulsed Linear Response Theory

    if ((task_type == 'pulses' or 'pulses' in task.input_dict) and
            scf_drv.scf_type == 'restricted'):
        prt_dict = (task.input_dict['pulses']
                    if 'pulses' in task.input_dict else {})

        cpp_dict = updated_dict_with_eri_settings({}, scf_drv)

        pulsed_response = PulsedResponse(task.mpi_comm, task.ostream)
        pulsed_response.update_settings(prt_dict, cpp_dict, method_dict)
        pulsed_response.compute(task.molecule, task.ao_basis, scf_results)

    # Resonant inelastic X-ray scattering
    if (task_type == 'rixs') and (scf_drv.scf_type == 'restricted'):
        rixs_dict = (task.input_dict['rixs']
                     if 'rixs' in task.input_dict else {})
        rixs_dict['program_end_time'] = program_end_time
        rixs_dict['filename'] = task.input_dict['filename']
        
        rsp_dict = (dict(task.input_dict['response'])
                    if 'response' in task.input_dict else {})
        rsp_dict['program_end_time'] = program_end_time
        rsp_dict['filename'] = task.input_dict['filename']
        rsp_dict = updated_dict_with_eri_settings(rsp_dict, scf_drv)

        rsp_prop = select_rsp_property(task, mol_orbs, rsp_dict, method_dict)
        rsp_prop.init_driver(task.mpi_comm, task.ostream)
        rsp_prop.compute(task.molecule, task.ao_basis, scf_results)

        cvs_rsp_dict = (dict(task.input_dict['cvs_response'])
                    if 'cvs_response' in task.input_dict else {})
        if cvs_rsp_dict == {}:
            # Restricted subspace approach
            cvs_rsp_prop = None
        else:
            # Two-shot approach
            cvs_rsp_dict['program_end_time'] = program_end_time
            cvs_rsp_dict['filename'] = task.input_dict['filename'] + '_cvs'
            cvs_rsp_dict = updated_dict_with_eri_settings(cvs_rsp_dict, scf_drv)
            
            cvs_rsp = select_rsp_property(task, mol_orbs, cvs_rsp_dict, method_dict)
            cvs_rsp.init_driver(task.mpi_comm, task.ostream)
            cvs_rsp.compute(task.molecule, task.ao_basis, scf_results)
            cvs_rsp_prop = cvs_rsp._rsp_property

        rixs_drv = RixsDriver(task.mpi_comm, task.ostream)
        rixs_drv.update_settings(rixs_dict, method_dict)
        rixs_drv.compute(task.molecule, task.ao_basis, scf_results, rsp_prop._rsp_property, cvs_rsp_prop)
    

    # MP2 perturbation theory

    if task_type in ['mp2', 'ump2', 'romp2']:
        mp2_dict = task.input_dict['mp2'] if 'mp2' in task.input_dict else {}
        mp2_dict = updated_dict_with_eri_settings(mp2_dict, scf_drv)

        mp2_drv = Mp2Driver(task.mpi_comm, task.ostream)
        mp2_drv.update_settings(mp2_dict, method_dict)
        mp2_drv.compute(task.molecule, task.ao_basis, scf_results)

    # Cube file

    if task_type == 'visualization':
        cube_dict = (task.input_dict['visualization']
                     if 'visualization' in task.input_dict else {})

        mol_orbs = mol_orbs.broadcast(task.mpi_comm, root=mpi_master())
        density = density.broadcast(task.mpi_comm, root=mpi_master())

        vis_drv = VisualizationDriver(task.mpi_comm)
        vis_drv.gen_cubes(cube_dict, task.molecule, task.ao_basis, mol_orbs,
                          density)

    # PE force field

    if task_type in ['loprop', 'pe force field']:
        pe_ff_gen = PEForceFieldGenerator(task.mpi_comm, task.ostream)
        pe_ff_gen.compute(task.molecule, task.ao_basis, scf_results)

    # RESP and ESP charges

    if task_type in ['resp charges', 'esp charges']:
        if (task_type == 'resp charges' and 'resp_charges' in task.input_dict):
            charges_dict = task.input_dict['resp_charges']
        elif (task_type == 'esp charges' and 'esp_charges' in task.input_dict):
            charges_dict = task.input_dict['esp_charges']
        else:
            charges_dict = {}

        charges_dict['filename'] = task.input_dict['filename']

        chg_drv = RespChargesDriver(task.mpi_comm, task.ostream)
        chg_drv.update_settings(charges_dict, method_dict)

        if task_type == 'resp charges':
            # TODO: use scf_results
            chg_drv.compute(task.molecule, task.ao_basis, 'resp')
        elif task_type == 'esp charges':
            # TODO: use scf_results
            chg_drv.compute(task.molecule, task.ao_basis, 'esp')

    # All done

    task.finish()
