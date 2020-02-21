import os

import h5py
import numpy as np

from loprop.veloxchem import MolFragVeloxChem as MolFrag

from .veloxchemlib import (
    OverlapIntegralsDriver,
    ao_matrix_to_dalton,
    DenseMatrix,
)


class LoPropDriver:

    def __init__(self, task):
        self.task = task
        self.settings = task.input_dict.get('loprop', {})

        scf_checkpoint_file = task.input_dict['scf']['checkpoint_file']
        self.settings['scf_checkpoint_file'] = scf_checkpoint_file

        self.checkpoint = self.settings['checkpoint_file']
        self.tmpdir = os.path.dirname(self.checkpoint)

    def compute(self):
        self.save_overlap()
        self.save_orbital_info()
        self.save_density()
        self.save_coordinates()

        molfrag = MolFrag(self.tmpdir, **self.settings)

        pot_output = molfrag.output_potential_file(0, 0, 0)
        self.task.ostream.print_line(pot_output)

    def save_overlap(self):

        task = self.task
        molecule = task.molecule
        basis = task.ao_basis

        ovldrv = OverlapIntegralsDriver(task.mpi_comm)
        S = ovldrv.compute(molecule, basis).to_numpy()

        S = ao_matrix_to_dalton(
            DenseMatrix(S),
            task.ao_basis,
            task.molecule
        )

        with h5py.File(self.checkpoint, 'w') as f:
            f['ao_overlap_matrix'] = S.to_numpy()

    def save_density(self):

        task = self.task

        da = self.settings['density'].alpha_to_numpy(0)
        db = self.settings['density'].beta_to_numpy(0)

        D = ao_matrix_to_dalton(
            DenseMatrix(da + db),
            task.ao_basis,
            task.molecule,
        )

        with h5py.File(self.checkpoint, 'a') as f:
            f['ao_density_matrix'] = D.to_numpy()

    def save_coordinates(self):
        x = self.task.molecule.x_to_numpy()
        y = self.task.molecule.y_to_numpy()
        z = self.task.molecule.z_to_numpy()
        r = np.array([x, y, z]).T
        with h5py.File(self.checkpoint, 'a') as f:
            f['nuclear_coordinates'] = r

    def save_orbital_info(self):
        with h5py.File(self.checkpoint, 'a') as f:
            f['contracted_per_atom'] = self.get_cpa()
            for i, occ in enumerate(self.get_opa()):
                f[f'occupied_per_atom/{i}'] = occ

    def get_cpa(self):
        ...

    def get_opa(self):
        ...

    def get_coordinates(self):
        with h5py.File(self.checkpoint, 'r') as f:
            r = f['nuclear_coordinates'][...]
        return r
