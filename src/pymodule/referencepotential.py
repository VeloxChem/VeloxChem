# SPDX-License-Identifier: BSD-3-Clause
"""Serial, isolated-molecule reference energies and Cartesian gradients."""

import copy
import hashlib
import importlib.metadata
import inspect
import json
from pathlib import Path

import numpy as np

from .potential import PotentialResult


def _json_value(value):
    """Convert settings, including VeloxChem functional objects, to JSON."""
    if hasattr(value, 'get_func_label'):
        return value.get_func_label()
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, dict):
        return {str(k): _json_value(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_json_value(v) for v in value]
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    raise TypeError(f'Cannot record setting of type {type(value).__name__}.')


def _identity(value):
    return hashlib.sha256(json.dumps(_json_value(value), sort_keys=True,
                                    allow_nan=False).encode()).hexdigest()


def _versions():
    versions = {}
    for package in ('veloxchem', 'numpy', 'tblite', 'xtb', 'dftd4', 'openmm'):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    return versions


class ReferencePotential:
    """Adapt xTB, SCF or an explicitly described MM driver to PotentialResult.

    The target is frozen at construction. Each evaluation checks settings,
    atom order, charge and multiplicity. SCF is recomputed before *every*
    gradient; no cached gradient is paired with a later geometry. An optional
    ``evaluation_callback(molecule, result)`` can collect successful reference
    evaluations during native potential dynamics, without evaluating twice.

    :param driver: XtbDriver, SCF driver, MMDriver or OpenMMDriver.
    :param molecule: The isolated molecular system to bind.
    :param basis_set_label: Required for SCF.
    :param target: Required explicit parameter provenance for MM.
    """

    capabilities = frozenset(('energy', 'gradient'))
    model_role = 'full'

    def __init__(self, driver, molecule, basis_set_label=None, target=None,
                 evaluation_callback=None):
        from .xtbdriver import XtbDriver
        from .scfdriver import ScfDriver
        from .mmdriver import MMDriver
        from .openmmdriver import OpenMMDriver

        if driver.comm.Get_size() != 1:
            raise ValueError('ReferencePotential supports one MPI rank only.')
        self.driver = driver
        self.comm = driver.comm
        self.basis_set_label = basis_set_label
        self.evaluation_callback = evaluation_callback
        self._numbers = list(molecule.get_element_ids())
        self._charge = molecule.get_charge()
        self._multiplicity = molecule.get_multiplicity()
        self._gradient_driver = None
        if isinstance(driver, XtbDriver):
            if not driver.is_available():
                raise ImportError('Reference labeling requires xtb-python or tblite-python.')
            self.kind = 'xtb'
        elif isinstance(driver, ScfDriver):
            from .scfgradientdriver import ScfGradientDriver
            if not basis_set_label:
                raise ValueError('SCF labeling requires basis_set_label.')
            if driver.ri_jk or driver.electric_field is not None or driver.solvation_model:
                raise ValueError('This adapter requires isolated SCF without RI-JK or applied field.')
            if (getattr(driver, 'potfile', None) is not None or
                    getattr(driver, 'point_charges', None) is not None or
                    getattr(driver, 'pressure', 0.0) not in (None, 0.0)):
                raise ValueError('Embedding is outside the isolated dataset target.')
            self.kind = 'scf'
            self._gradient_driver = ScfGradientDriver(driver)
        elif isinstance(driver, (MMDriver, OpenMMDriver)):
            if not target or not target.get('parameters_sha256'):
                raise ValueError('MM labels require target with parameters_sha256 and method.')
            if not target.get('method'):
                raise ValueError('MM labels require an explicit method description.')
            if isinstance(driver, OpenMMDriver):
                if not hasattr(driver, 'simulation'):
                    raise ValueError('Load the OpenMM topology before creating a reference target.')
                system = driver.simulation.system
                numbers = [atom.element.atomic_number for atom in driver.simulation.topology.atoms()]
                if (system.usesPeriodicBoundaryConditions() or numbers != self._numbers or
                        any(system.isVirtualSite(i) for i in range(system.getNumParticles()))):
                    raise ValueError('OpenMM labels require the same isolated atoms without virtual sites.')
            self.kind = 'mm'
        else:
            raise TypeError('Use XtbDriver, an SCF driver, MMDriver or OpenMMDriver.')
        self._settings = self._driver_settings()
        self._target = dict(method=self._settings, charge=self._charge,
                            multiplicity=self._multiplicity, state='ground state',
                            environment='isolated molecule', model_role='full',
                            backend=type(driver).__name__, versions=_versions(),
                            energy_unit='hartree', gradient_unit='hartree/bohr')
        if target:
            for key, value in target.items():
                if key in self._target and key != 'method' and value != self._target[key]:
                    raise ValueError(f'Conflicting target field: {key}.')
                self._target[key] = _json_value(value)
        # Always retain the actual settings, even with a human-readable method.
        self._target['driver_settings'] = self._settings
        implementations = [type(driver)]
        if self._gradient_driver is not None:
            implementations.append(type(self._gradient_driver))
            # SCF algorithm lives in the base class, not only the selected spin driver.
            implementations.append(ScfDriver)
        self._target['implementation_sha256'] = {
            cls.__name__: hashlib.sha256(Path(inspect.getfile(cls)).read_bytes()).hexdigest()
            for cls in implementations}
        self.model_id = _identity(self._target)

    def _driver_settings(self):
        if self.kind == 'xtb':
            return dict(method=self.driver.get_method(),
                        backend='tblite' if self.driver.use_tblite() else 'xtb')
        if self.kind == 'mm':
            # The caller must freeze the parameterized MM system.
            if hasattr(self.driver, 'simulation'):
                import openmm
                return dict(system_sha256=hashlib.sha256(openmm.XmlSerializer.serialize(
                    self.driver.simulation.system).encode()).hexdigest(),
                            context_parameters=dict(self.driver.simulation.context.getParameters()))
            if not getattr(self.driver, 'force_field', None):
                raise ValueError('Load the MM force field before creating a reference target.')
            return _json_value(self.driver.force_field)
        settings = {}
        # Record all exposed method/SCF settings except output/restart controls.
        omitted = {'filename', 'checkpoint_file', 'restart', 'timing', 'profiling'}
        for group in self.driver._input_keywords.values():
            for key in group:
                if key not in omitted and hasattr(self.driver, key):
                    settings[key] = _json_value(getattr(self.driver, key))
        settings['basis_set_label'] = self.basis_set_label
        # Resolve defaults which SCF converts to objects during compute.
        settings['xcfun'] = (settings.get('xcfun') or 'HF').upper()
        return settings

    @property
    def manifest(self):
        return dict(model_id=self.model_id, target=copy.deepcopy(self._target),
                    model_role='full')

    def compute(self, molecule, properties=('energy', 'gradient')):
        """Return full energy in Eh and unweighted Cartesian gradient in Eh/bohr."""
        if not set(properties).issubset(self.capabilities):
            raise ValueError('ReferencePotential provides energy and gradient only.')
        if (list(molecule.get_element_ids()) != self._numbers or
                molecule.get_charge() != self._charge or
                molecule.get_multiplicity() != self._multiplicity):
            raise ValueError('Reference molecule differs in atom order or electronic state.')
        if not np.all(np.isfinite(molecule.get_coordinates_in_angstrom())):
            raise ValueError('Reference coordinates must be finite.')
        if self._driver_settings() != self._settings:
            raise ValueError('Reference settings changed; create a new target.')
        if self.kind == 'xtb':
            values = self.driver.compute(molecule)
            result = PotentialResult(values['energy'], values['gradient'], self.model_id)
        elif self.kind == 'scf':
            from .molecularbasis import MolecularBasis
            basis = MolecularBasis.read(molecule, self.basis_set_label)
            self.driver.compute(molecule, basis)
            if not self.driver.is_converged:
                raise RuntimeError('SCF did not converge; no label was accepted.')
            self._gradient_driver.compute(molecule, basis)
            result = PotentialResult(self.driver.get_scf_energy(),
                                     self._gradient_driver.get_gradient(), self.model_id)
        else:
            self.driver.compute(molecule)
            result = PotentialResult(self.driver.get_energy(), self.driver.get_gradient(),
                                     self.model_id)
        if result.gradient.shape != (len(self._numbers), 3):
            raise ValueError('Reference gradient shape differs from the molecule.')
        if self.evaluation_callback is not None:
            self.evaluation_callback(molecule, result)
        return result
