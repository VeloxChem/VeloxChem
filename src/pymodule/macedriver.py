# SPDX-License-Identifier: BSD-3-Clause
"""Optional, serial MACE energy/gradient evaluation for isolated molecules."""

import copy
import hashlib
import importlib.metadata
import json
from pathlib import Path
from time import perf_counter

import numpy as np
from mpi4py import MPI

from .potential import PotentialResult
from .veloxchemlib import bohr_in_angstrom, hartree_in_ev

class MACEDriver:
    """Evaluate one explicitly selected local MACE checkpoint.

    Only CPU float64, one MPI rank, and a fixed isolated electronic target
    are supported. No weights are downloaded. Charge/multiplicity describe
    the fixed target; they do not add conditioning to an unconditioned model.
    ``metadata_file`` is an optional JSON object containing ``target``,
    ``model_role`` (full or residual), and/or ``checkpoint_sha256``. Unknown
    training provenance remains unknown. No external corrections are added.
    """

    capabilities = frozenset(('energy', 'gradient'))

    def __init__(self, checkpoint, device='cpu', dtype='float64', head=None,
                 metadata_file=None, target=None, charge=None,
                 multiplicity=None, model_role=None, comm=None, ostream=None):
        
        self.comm = MPI.COMM_WORLD if comm is None else comm
        if self.comm.Get_size() != 1:
            raise ValueError('MACEDriver supports one MPI rank only.')

        if device != 'cpu' or dtype != 'float64':
            raise ValueError('MACEDriver currently requires device="cpu", dtype="float64".')
        path = Path(checkpoint).expanduser().resolve(strict=True)
        if not path.is_file():
            raise ValueError('MACE checkpoint must be one local file.')
        
        metadata = {} if metadata_file is None else json.loads(
            Path(metadata_file).read_text(encoding='utf-8'))
        
        if not isinstance(metadata, dict):
            raise ValueError('MACE metadata must be a JSON object.')
        
        declared = copy.deepcopy(metadata.get('target', metadata.get('electronic_target', {})))
        if not isinstance(declared, dict):
            raise ValueError('MACE target must be a JSON object.')
        if target is not None:
            if not isinstance(target, dict):
                raise ValueError('MACE target must be a JSON object.')
            for key, value in target.items():
                if key in declared and declared[key] is not None and declared[key] != value:
                    raise ValueError(f'Conflicting target declaration: {key}.')
                declared[key] = copy.deepcopy(value)
        
        for key, value, default in [('charge', charge, 0), ('multiplicity', multiplicity, 1)]:
            if value is not None and key in declared and declared[key] != value: # making sure that the provided charge/multiplicity does not conflict with the metadata
                raise ValueError(f'Conflicting target declaration: {key}.')
            declared[key] = declared.get(key, default) if value is None else value
        
        q, mult = declared['charge'], declared['multiplicity']
        if not np.isfinite(q) or int(q) != q:
            raise ValueError('An integer molecular charge is required.')
        if not np.isfinite(mult) or int(mult) != mult or mult < 1:
            raise ValueError('Multiplicity must be a positive integer (2S+1).')
        declared.setdefault('state', 'unknown')
        declared.setdefault('method', None)
        declared.setdefault('environment', 'isolated molecule')
        if declared['environment'] != 'isolated molecule':
            raise ValueError('Only the isolated molecule environment is supported.')
        role = metadata.get('model_role', 'full') if model_role is None else model_role
        if role not in ('full', 'residual') or (
                'model_role' in metadata and role != metadata['model_role']):
            raise ValueError('Invalid or conflicting model_role (full or residual).')
        for key, expected in [('native_energy_unit', 'eV'), ('native_length_unit', 'angstrom')]:
            if metadata.get(key, expected) != expected:
                raise ValueError(f'Unsupported {key}: {metadata[key]}.')
        checksum = hashlib.sha256()
        with path.open('rb') as checkpoint_stream:
            for block in iter(lambda: checkpoint_stream.read(1024 * 1024), b''):
                checksum.update(block)
        digest = checksum.hexdigest()
        if metadata.get('checkpoint_sha256', digest) != digest:
            raise ValueError('MACE checkpoint SHA-256 does not match metadata.')
        start = perf_counter()
        # Optional ML dependencies are needed only when constructing this driver.
        try:
            from ase import Atoms
            from mace.calculators import MACECalculator
            import torch
        except ImportError as exc:
            raise ImportError('MACEDriver requires ASE, PyTorch, and mace-torch; '
                              'install them in the active VeloxChem environment.') from exc
        try:
            calculator = MACECalculator(model_paths=[str(path)], device=device,
                                        default_dtype=dtype, head=head)
        except Exception as exc:
            raise RuntimeError(f'Cannot load MACE checkpoint {path}: {exc}') from exc
        heads = list(calculator.available_heads)
        if head is None and len(heads) != 1:
            raise ValueError(f'Multiple model heads {heads}; select one explicitly.')
        effective = calculator.head
        if head is not None and (head not in heads or effective != head):
            raise ValueError(f'Requested MACE head {head!r} was not selected; available: {heads}.')
        if calculator.num_models != 1:
            raise ValueError('Exactly one MACE model is required.')
        model = calculator.models[0]
        family = type(model).__name__
        if family not in ('MACE', 'ScaleShiftMACE'):
            raise ValueError(f'Model family {family} is outside the initial MACE energy/force scope.')
        parameter = next(model.parameters())
        if parameter.dtype != torch.float64 or parameter.device.type != 'cpu':
            raise ValueError('MACE did not initialize on CPU in float64.')
        elements = [int(z) for z in calculator.z_table.zs]
        versions = {name: importlib.metadata.version(name)
                    for name in ('mace-torch', 'ase', 'torch', 'numpy', 'veloxchem')}
        
        identity = dict(schema_version=1, checkpoint_sha256=digest,
                        model_family=family, head=effective, device=device,
                        dtype=dtype, elements=elements, target=declared,
                        model_role=role, native_energy_unit='eV',
                        native_length_unit='angstrom', energy_convention='full model output',
                        external_corrections=[], packages=versions)
        
        self._model_id = hashlib.sha256(json.dumps(identity, sort_keys=True,
                                                  allow_nan=False).encode()).hexdigest()
        
        self._manifest = dict(identity, model_id=self._model_id, checkpoint=str(path),
                              requested_head=head, available_heads=heads,
                              supplied_metadata=metadata,
                              target_source='declaration; not inferred from checkpoint',
                              load_seconds=perf_counter() - start)
        
        self._calculator = calculator
        self._Atoms = Atoms
        self._target = declared
        self._elements = frozenset(elements)
        self._role = role
        self._atoms = None
        self._result = None
        self.evaluation_count = 0
        self.ostream = ostream

    @property
    def model_id(self):
        return self._model_id

    @property
    def model_role(self):
        return self._role

    @property
    def manifest(self):
        """A detached JSON-serializable model configuration/provenance record."""
        return copy.deepcopy(self._manifest)

    
    def compute(self, molecule, properties=('energy', 'gradient'), *, pbc=False, cell=None):
        """Return energy and physical Cartesian gradient in atomic units.

        Every call refreshes ASE results, including model information inputs.
        The model itself is retained. Returned arrays own their storage.
        """
        self._result = None
        
        if not properties or set(properties) - self.capabilities:
            raise NotImplementedError('MACEDriver supports energy and gradient only.')
        
        if np.any(pbc) or cell is not None:
            raise NotImplementedError('MACEDriver does not support cells or periodicity.')
        numbers = np.asarray(molecule.get_element_ids())
        positions = np.asarray(molecule.get_coordinates_in_angstrom(), dtype=float)
        if (numbers.ndim != 1 or not numbers.size or not np.all(np.isfinite(numbers))
                or not np.all(numbers == np.floor(numbers))):
            raise ValueError('Atomic numbers must be a nonempty integer vector.')
        numbers = numbers.astype(int)
        if positions.shape != (len(numbers), 3) or not np.all(np.isfinite(positions)):
            raise ValueError('Coordinates must be finite with shape (N, 3).')
        missing = set(numbers) - self._elements
        if missing:
            raise ValueError(f'MACE checkpoint does not support atomic numbers {sorted(missing)}.')
        if (molecule.get_charge() != self._target['charge'] or
                molecule.get_multiplicity() != self._target['multiplicity']):
            raise ValueError('Molecule charge/multiplicity differs from the fixed MACE target.')
        if self._atoms is None or not np.array_equal(self._atoms.numbers, numbers):
            self._atoms = self._Atoms(numbers=numbers, positions=positions, pbc=False)
        else:
            self._atoms.set_positions(positions)
        self._atoms.info.update(charge=self._target['charge'], spin=self._target['multiplicity'])
        try:
            self._calculator.reset()
            self._calculator.calculate(self._atoms, properties=['energy', 'forces'])
            energy = self._calculator.results['energy']
            forces = np.asarray(self._calculator.results['forces'])
            if np.asarray(energy).shape != () or forces.shape != positions.shape:
                raise ValueError('MACE returned an invalid energy/force shape.')
            result = PotentialResult(float(energy) / hartree_in_ev(),
                                     -forces * (bohr_in_angstrom() / hartree_in_ev()),
                                     self.model_id)
        except Exception as exc:
            raise RuntimeError(f'MACE evaluation {self.evaluation_count + 1} failed '
                               f'for {self.model_id} ({len(numbers)} atoms): {exc}') from exc
        self.evaluation_count += 1
        self._result = result
        return result

    
    def get_energy(self):
        if self._result is None:
            raise RuntimeError('No successful MACE evaluation is available.')
        return self._result.energy

    
    def get_gradient(self):
        if self._result is None:
            raise RuntimeError('No successful MACE evaluation is available.')
        return self._result.gradient.copy()
