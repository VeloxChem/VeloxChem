#
#                           VELOXCHEM 1.0-RC3
#         ----------------------------------------------------
#                     An Electronic Structure Code
#
#  Copyright © 2018-2022 by VeloxChem developers. All rights reserved.
#  Contact: https://veloxchem.org/contact
#
#  SPDX-License-Identifier: LGPL-3.0-or-later
#
#  This file is part of VeloxChem.
#
#  VeloxChem is free software: you can redistribute it and/or modify it under
#  the terms of the GNU Lesser General Public License as published by the Free
#  Software Foundation, either version 3 of the License, or (at your option)
#  any later version.
#
#  VeloxChem is distributed in the hope that it will be useful, but WITHOUT
#  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
#  FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public
#  License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with VeloxChem. If not, see <https://www.gnu.org/licenses/>.

from mpi4py import MPI
from datetime import datetime, timedelta
import time as tm

from .veloxchemlib import mpi_initialized, mpi_master
from .mpitask import MpiTask
from .scfrestdriver import ScfRestrictedDriver
from .scfunrestdriver import ScfUnrestrictedDriver
from .scfrestopendriver import ScfRestrictedOpenDriver
from .firstorderprop import FirstOrderProperties
from .forcefieldgenerator import ForceFieldGenerator
from .respchargesdriver import RespChargesDriver
from .excitondriver import ExcitonModelDriver
from .numerovdriver import NumerovDriver
from .mp2driver import Mp2Driver
from .cnadriver import CnaAnalysisDriver
from .gopdriver import GlobalOptimizationDriver
from .loprop import LoPropDriver
from .trajectorydriver import TrajectoryDriver
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
from .rsptpa import TPA
from .tdhfhessiandriver import TdhfHessianDriver
from .polarizabilitygradient import PolarizabilityGradient
from .cphfsolver import CphfSolver
from .rspcustomproperty import CustomProperty
from .visualizationdriver import VisualizationDriver
from .xtbdriver import XtbDriver
from .xtbgradientdriver import XtbGradientDriver
from .xtbhessiandriver import XtbHessianDriver
from .veloxchemlib import DiagEriDriver
from .cli import cli
from .dftutils import print_libxc_reference
from .errorhandler import assert_msg_critical


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
        n_ao = task.ao_basis.get_dimension_of_basis(task.molecule)
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
        assert_msg_critical(False, f'SCF: invalide scf_type {scf_type}')

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

    elif prop_type == 'tpa':
        rsp_prop = TPA(rsp_dict, method_dict)

    elif prop_type == 'custom':
        rsp_prop = CustomProperty(rsp_dict, method_dict)

    else:
        assert_msg_critical(
            False, f'Response: invalide response property {prop_type}')

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
    if 'qq_type' not in new_dict:
        new_dict['qq_type'] = scf_drv.qq_type
    if not scf_drv.restart:
        new_dict['restart'] = 'no'

    return new_dict


def main():
    """
    Runs VeloxChem with command line arguments.
    """

    program_start_time = datetime.now()

    assert_msg_critical(mpi_initialized(), "MPI not initialized")

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

    if task_type == 'force field':
        force_field_dict = (dict(task.input_dict['force_field'])
                            if 'force_field' in task.input_dict else {})
        force_field_dict['filename'] = task.input_dict['filename']

        resp_dict = (dict(task.input_dict['resp_charges'])
                     if 'resp_charges' in task.input_dict else {})
        resp_dict['filename'] = task.input_dict['filename']

        force_field_drv = ForceFieldGenerator(task.mpi_comm, task.ostream)
        force_field_drv.update_settings(force_field_dict, resp_dict)
        force_field_drv.compute(task.molecule, task.ao_basis)

    # Spectrum from trajectory

    if task_type == 'trajectory':
        traj_dict = (dict(task.input_dict['trajectory'])
                     if 'trajectory' in task.input_dict else {})
        traj_dict['filename'] = task.input_dict['filename']
        traj_dict['charges'] = task.input_dict['charges']
        traj_dict['polarizabilities'] = task.input_dict['polarizabilities']

        spect_dict = (dict(task.input_dict['spectrum_settings'])
                      if 'spectrum_settings' in task.input_dict else {})

        rsp_dict = (dict(task.input_dict['response'])
                    if 'response' in task.input_dict else {})

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
        'wave function', 'mp2', 'ump2', 'romp2', 'gradient', 'hessian',
        'optimize', 'response', 'pulses', 'visualization', 'loprop',
        'vibrational', 'freq', 'cphf', 'polarizability_gradient'
    ]

    if task_type == 'visualization' and 'visualization' in task.input_dict:
        run_scf = 'read_dalton' not in task.input_dict['visualization']['cubes']

    scf_type = 'restricted'
    if task_type in ['uhf', 'uscf', 'ump2']:
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
            xtb_drv = XtbDriver(task.mpi_comm)
            xtb_drv.set_method(method_dict['xtb'].lower())
            xtb_drv.compute(task.molecule, task.ostream)
        else:
            scf_drv = select_scf_driver(task, scf_type)
            scf_drv.update_settings(scf_dict, method_dict)
            scf_results = scf_drv.compute(task.molecule, task.ao_basis,
                                          task.min_basis)

            mol_orbs = scf_drv.molecular_orbitals
            density = scf_drv.density

            if not scf_drv.is_converged:
                return

            # SCF first-order properties
            scf_prop = FirstOrderProperties(task.mpi_comm, task.ostream)
            scf_prop.compute_scf_prop(task.molecule, task.ao_basis, scf_results)
            if task.mpi_rank == mpi_master():
                scf_prop.print_properties(task.molecule)

            if (scf_drv.electric_field is not None and
                    task.molecule.get_charge() != 0):
                task.finish()
                return

    # Gradient

    if task_type == 'gradient':

        run_excited_state_gradient = ('response' in task.input_dict)
        run_ground_state_gradient = (not run_excited_state_gradient)

        if run_ground_state_gradient:

            if use_xtb:
                grad_drv = XtbGradientDriver(task.mpi_comm, task.ostream)
                grad_drv.compute(task.molecule, xtb_drv)
            elif scf_drv.scf_type == 'restricted':
                grad_drv = ScfGradientDriver(task.mpi_comm, task.ostream)
                grad_drv.compute(task.molecule, task.ao_basis, scf_drv)

        elif run_excited_state_gradient:

            grad_dict = (task.input_dict['gradient']
                         if 'gradient' in task.input_dict else {})

            rsp_dict = dict(task.input_dict['response'])
            rsp_dict['program_end_time'] = program_end_time
            rsp_dict['filename'] = task.input_dict['filename']
            rsp_dict = updated_dict_with_eri_settings(rsp_dict, scf_drv)

            assert_msg_critical(
                rsp_dict['property'].lower() in ['absorption', 'uv-vis', 'ecd'],
                'Invalid response property for gradient calculation')

            rsp_prop = select_rsp_property(task, mol_orbs, rsp_dict,
                                           method_dict)
            rsp_prop.init_driver(task.mpi_comm, task.ostream)

            tddftgrad_drv = TddftGradientDriver(task.mpi_comm, task.ostream)
            tddftgrad_drv.update_settings(grad_dict, rsp_dict, method_dict)
            tddftgrad_drv.compute(task.molecule, task.ao_basis, scf_drv,
                                  rsp_prop.rsp_driver)
    # Hessian

    if task_type == 'hessian':
        hessian_dict = (task.input_dict['hessian']
                        if 'hessian' in task.input_dict else {})

        if use_xtb:
            hessian_drv = XtbHessianDriver(task.mpi_comm, task.ostream)
            hessian_drv.update_settings(method_dict, hessian_dict)
            hessian_drv.compute(task.molecule, xtb_drv)

        elif scf_drv.scf_type == 'restricted':
            hessian_drv = ScfHessianDriver(task.mpi_comm, task.ostream)
            hessian_drv.update_settings(method_dict, hessian_dict)
            hessian_drv.compute(task.molecule, task.ao_basis, scf_drv)

        if task.mpi_rank == mpi_master():
            hessian_drv.vibrational_analysis(task.molecule)

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
                grad_drv = XtbGradientDriver(task.mpi_comm, task.ostream)
                opt_drv = OptimizationDriver(grad_drv)
                opt_drv.update_settings(opt_dict)
                opt_drv.compute(task.molecule, xtb_drv)

            elif scf_drv.scf_type == 'restricted':
                grad_drv = ScfGradientDriver(task.mpi_comm, task.ostream)
                opt_drv = OptimizationDriver(grad_drv)
                opt_drv.update_settings(opt_dict)
                opt_drv.compute(task.molecule, task.ao_basis, scf_drv)

        elif run_excited_state_gradient:

            grad_dict = (task.input_dict['gradient']
                         if 'gradient' in task.input_dict else {})

            rsp_dict = dict(task.input_dict['response'])
            rsp_dict['program_end_time'] = program_end_time
            rsp_dict['filename'] = task.input_dict['filename']
            rsp_dict = updated_dict_with_eri_settings(rsp_dict, scf_drv)

            assert_msg_critical(
                rsp_dict['property'].lower() in ['absorption', 'uv-vis', 'ecd'],
                'Invalid response property for geometry optimization')

            rsp_prop = select_rsp_property(task, mol_orbs, rsp_dict,
                                           method_dict)
            rsp_prop.init_driver(task.mpi_comm, task.ostream)

            tddftgrad_drv = TddftGradientDriver(task.mpi_comm, task.ostream)
            tddftgrad_drv.update_settings(grad_dict, rsp_dict, method_dict)

            opt_drv = OptimizationDriver(tddftgrad_drv)
            opt_drv.update_settings(opt_dict)
            opt_drv.compute(task.molecule, task.ao_basis, scf_drv,
                            rsp_prop.rsp_driver, rsp_prop._rsp_property)

    # Ground state Hessian / Vibrational analysis

    #if task_type in ['freq', 'frequencies']:
    if task_type in ['vib', 'vibrational']:

        #if 'frequencies' in task.input_dict:
        #    freq_dict = task.input_dict['frequencies']
        #else:
        #    freq_dict = {}

        vib_dict = (task.input_dict['vibrational']
                    if 'vibrational' in task.input_dict else {})
        polgrad_dict = (task.input_dict['polarizability_gradient'] 
                        if 'polarizability_gradient' in task.input_dict else {})
        orbrsp_dict = (task.input_dict['orbital_response']
                       if 'orbital_response' in task.input_dict else {})
        rsp_dict = (dict(task.input_dict['response'])
                    if 'response' in task.input_dict else {})

        hessian_drv = ScfHessianDriver(task.mpi_comm, task.ostream)
        #hessian_drv.update_settings(method_dict, freq_dict)
        hessian_drv.update_settings(method_dict, hess_dict = vib_dict, 
                                    cphf_dict = orbrsp_dict, rsp_dict = rsp_dict,
                                    polgrad_dict = polgrad_dict)
        hessian_drv.compute(task.molecule, task.ao_basis, scf_drv)
        # TODO: add output file name for geomeTRIC vibrational analysis
        if task.mpi_rank == mpi_master():
            hessian_drv.vibrational_analysis(task.molecule)

    # TODO: maybe remove (this was included for testing the MPI-parallelization)
    if task_type in ['cphf']:
        if 'cphf_settings' in task.input_dict:
            cphf_dict = task.input_dict['cphf_settings']
        else:
            cphf_dict = {}

        cphf_drv = CphfSolver(task.mpi_comm, task.ostream)
        cphf_drv.update_settings(cphf_dict, method_dict)
        cphf_drv.compute(task.molecule, task.ao_basis, scf_drv.scf_tensors) 


    # Polarizability gradient

    if task_type in ['polarizability_gradient']:

        polgrad_dict = (task.input_dict['polarizability_gradient']
                     if 'polarizability_gradient' in task.input_dict else {})
        orbrsp_dict = (task.input_dict['orbital_response']
                       if 'orbital_response' in task.input_dict else {})

        rsp_dict = (dict(task.input_dict['response'])
                    if 'response' in task.input_dict else {})
        rsp_dict['program_end_time'] = program_end_time
        rsp_dict['filename'] = task.input_dict['filename']
        rsp_dict = updated_dict_with_eri_settings(rsp_dict, scf_drv)

        rsp_prop = select_rsp_property(task, mol_orbs, rsp_dict, method_dict)
        rsp_prop.init_driver(task.mpi_comm, task.ostream)
        rsp_prop.compute(task.molecule, task.ao_basis, scf_results)
        
        polgrad_drv = PolarizabilityGradient(task.mpi_comm, task.ostream)
        polgrad_drv.update_settings(polgrad_dict, orbrsp_dict, method_dict)
        polgrad_drv.compute(task.molecule, task.ao_basis, 
                            scf_drv.scf_tensors, rsp_prop._rsp_property)
    

    # Response

    if task_type == 'response' and scf_drv.scf_type == 'restricted':
        rsp_dict = (dict(task.input_dict['response'])
                    if 'response' in task.input_dict else {})
        rsp_dict['program_end_time'] = program_end_time
        rsp_dict['filename'] = task.input_dict['filename']
        rsp_dict = updated_dict_with_eri_settings(rsp_dict, scf_drv)

        rsp_prop = select_rsp_property(task, mol_orbs, rsp_dict, method_dict)
        rsp_prop.init_driver(task.mpi_comm, task.ostream)
        rsp_prop.compute(task.molecule, task.ao_basis, scf_results)

        if not rsp_prop.is_converged:
            return

        # Calculate the excited-state gradient

        if 'property' in rsp_dict:
            prop_type = rsp_dict['property'].lower()
        else:
            prop_type = None

        if prop_type in ['absorption', 'uv-vis', 'ecd']: 

            if 'orbital_response' in task.input_dict:
                orbrsp_dict = task.input_dict['orbital_response']
                orbrsp_drv = TddftOrbitalResponse(task.mpi_comm, task.ostream)
                if 'tamm_dancoff' in rsp_dict:
                    orbrsp_dict['tamm_dancoff'] = rsp_dict['tamm_dancoff']
                orbrsp_drv.update_settings(orbrsp_dict, method_dict)
                orbrsp_drv.compute(task.molecule, task.ao_basis,
                                   scf_drv.scf_tensors, rsp_prop._rsp_property)

            else:
                orbrsp_dict = {}

            # Excited state gradient
            if 'gradient' in task.input_dict:
                grad_dict = task.input_dict['gradient']
                tddftgrad_drv = TddftGradientDriver(task.mpi_comm, task.ostream)
                tddftgrad_drv.update_settings(grad_dict, rsp_dict,
                                             orbrsp_dict, method_dict)
                tddftgrad_drv.compute(task.molecule, task.ao_basis,
                                     scf_drv,
                                     rsp_prop._rsp_driver,
                                     rsp_prop._rsp_property)

            # Excited state Hessian and vibrational analysis
            if 'vibrational' in task.input_dict:
                freq_dict = task.input_dict['vibrational']
                tdhfhessian_drv = TdhfHessianDriver(scf_drv, 
                                                    task.mpi_comm,
                                                    task.ostream)
                tdhfhessian_drv.update_settings(method_dict, rsp_dict,
                                                freq_dict, orbrsp_dict)
                tdhfhessian_drv.compute(task.molecule, task.ao_basis,
                                        rsp_prop._rsp_driver)
                if task.mpi_rank == mpi_master():
                    tdhfhessian_drv.vibrational_analysis(task.molecule)

            # Excited state optimization
            if 'optimize_excited_state' in task.input_dict:
                opt_dict = task.input_dict['optimize_excited_state']
                tddftgrad_drv = TddftGradientDriver(task.mpi_comm, task.ostream)
                tddftgrad_drv.update_settings(opt_dict, rsp_dict,
                                             orbrsp_dict, method_dict)
                opt_drv = OptimizationDriver(tddftgrad_drv)
                opt_drv.compute(task.molecule, task.ao_basis,
                                scf_drv,
                                rsp_prop._rsp_driver,
                                rsp_prop._rsp_property)

        else:
            task.ostream.print_blank()
            info_msg = 'The excited state derivatives '
            info_msg += 'can only be computed if the response '
            info_msg += 'property is "absorption", "uv-vis", or "ecd"'
            task.ostream.print_info(info_msg)
            info_msg = 'Computation of gradient/Hessian will be skipped.'
            task.ostream.print_info(info_msg)
            task.ostream.print_blank()
            task.ostream.flush()

    # Pulsed Linear Response Theory

    if ((task_type == 'pulses' or 'pulses' in task.input_dict) and
            scf_drv.scf_type == 'restricted'):
        prt_dict = (task.input_dict['pulses']
                    if 'pulses' in task.input_dict else {})

        cpp_dict = updated_dict_with_eri_settings({}, scf_drv)

        pulsed_response = PulsedResponse(task.mpi_comm, task.ostream)
        pulsed_response.update_settings(prt_dict, cpp_dict, method_dict)
        pulsed_response.compute(task.molecule, task.ao_basis, scf_results)

    # MP2 perturbation theory

    if task_type in ['mp2', 'ump2', 'romp2']:
        mp2_dict = task.input_dict['mp2'] if 'mp2' in task.input_dict else {}
        mp2_dict = updated_dict_with_eri_settings(mp2_dict, scf_drv)

        mp2_drv = Mp2Driver(task.mpi_comm, task.ostream)
        mp2_drv.update_settings(mp2_dict, method_dict)
        mp2_drv.compute(task.molecule, task.ao_basis, mol_orbs,
                        scf_drv.scf_type)

    # Cube file

    if task_type == 'visualization':
        cube_dict = (task.input_dict['visualization']
                     if 'visualization' in task.input_dict else {})

        if 'read_dalton' not in task.input_dict['visualization']['cubes']:
            mol_orbs.broadcast(task.mpi_rank, task.mpi_comm)
            density.broadcast(task.mpi_rank, task.mpi_comm)
        else:
            mol_orbs = None
            density = None

        vis_drv = VisualizationDriver(task.mpi_comm)
        vis_drv.gen_cubes(cube_dict, task.molecule, task.ao_basis, mol_orbs,
                          density)

    # LoProp

    if task_type == 'loprop':
        loprop_driver = LoPropDriver(task.mpi_comm, task.ostream)
        loprop_driver.compute(task.molecule, task.ao_basis, scf_results)

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
            chg_drv.compute(task.molecule, task.ao_basis, 'resp')
        elif task_type == 'esp charges':
            chg_drv.compute(task.molecule, task.ao_basis, 'esp')

    # Test of electron repulstion integrals

    if task_type == 'eritest':
        print('*** Testing Two Electron Implementation ***')
        tm0 = tm.time()
        eri_driver = DiagEriDriver()
        pgblock = eri_driver.compute(task.molecule, task.ao_basis)
        print('Diagonal Eri Driver: ', tm.time() - tm0, ' sec.')

    # CNA correlation analysis

    if task_type == 'cna':
        if 'cna' in task.input_dict:
            cna_dict = task.input_dict['cna']
        else:
            cna_dict = {}

        cna_drv = CnaAnalysisDriver(task.mpi_comm, task.ostream)
        cna_drv.update_settings(cna_dict)
        cna_drv.compute()

    # Global optimization with tree-growth scheme

    if task_type == 'gop':
        if 'gop' in task.input_dict:
            gop_dict = task.input_dict['gop']
        else:
            gop_dict = {}

        gop_drv = GlobalOptimizationDriver(task.mpi_comm, task.ostream)
        gop_drv.update_settings(gop_dict)
        gop_drv.compute(task.input_dict['filename'], task.molecule)

    # All done

    print_libxc_reference(method_dict.get('xcfun', None), task.ostream)
    task.finish()
