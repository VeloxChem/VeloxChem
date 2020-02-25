from collections import Counter
import os

import h5py
import numpy as np

from loprop.veloxchem import MolFragVeloxChem as MolFrag

from .veloxchemlib import (
    OverlapIntegralsDriver,
    ao_matrix_to_dalton,
    DenseMatrix,
)

from .inputparser import InputParser


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
        elements = self.task.molecule.elem_ids_to_numpy()
        """
        bp = InputParser(current basis)
        for e in elements
            get e:th entry in bp.input_dict as e
                 given 3s 2p 1d dict sum up 3 + 2*3 + 1 * 5
        """

        basis = self.task.ao_basis.get_label()
        basis_file = f'basis/{basis}'
        bp = InputParser(basis_file)
        basis = bp.get_dict()
        keys = list(basis.keys())

        cpa = []
        for e in elements:
            k = keys[e - 1]
            atoms_data = basis[k]
            count_per_angmom = count_contracted(atoms_data)
            cpa.append(count_contracted_on_atom(count_per_angmom))

        return cpa

    def count_contracted(self):
        return 1

    def get_opa(self):
        """
        Returns list of occupied for case of single atom
        i.e. first row 1s, second row 1s, 2s, 2p, whole shell
        """
        elements = self.task.molecule.elem_ids_to_numpy()

        basis = self.task.ao_basis.get_label()
        basis_file = f'basis/{basis}'
        bp = InputParser(basis_file)
        basis = bp.get_dict()
        keys = list(basis.keys())

        opa = []
        for e in elements:
            opa.append([])
            k = keys[e - 1]
            atoms_data = basis[k]
            count_per_angmom = count_contracted(atoms_data)

            # For H and He: 1s
            opa[-1].append(0)

            # For Li-Ne: + 2s 2p
            if e>= 3:
                opa[-1].append(1)
                offset_p = count_per_angmom['S']
                opa[-1].append(offset_p + 0)
                opa[-1].append(offset_p + 1)
                opa[-1].append(offset_p + 2)

        return opa

    def get_coordinates(self):
        with h5py.File(self.checkpoint, 'r') as f:
            r = f['nuclear_coordinates'][...]
        return r

    @classmethod
    def verify_input(settings):
        """
        Verify input consistent with supported definitions
        """
        if 'localize' in settings:
            valid_options = ['charges']
            undefined = set(settings['localize']) ^ set(valid_options)
            if undefined:
                raise NotImplementedError(str(undefined))


def count_contracted(atombasis: list) -> dict:
    """
    Given atomic block in basis filei format return dict
    which maps angular momentum to number of contracted
    >>> count_contracted('S: 1 1\n1.0 1.9\n')
    {'S': 1}
    """
    c = Counter(
        line[0]
        for line in atombasis
        if line and line[0] in "SPDFGHI"
    )
    return dict(c)


def count_contracted_on_atom(atombasis: dict) -> int:
    """
    Returns total contracted given angular momentum count
    >>> count_contracted_on_atom({'S': 1})
    1
    >>> count_contracted_no_atom({'S': 2, 'P': 1})
    5
    """
    multiplicity = dict(S=1, P=3, D=5, F=7, G=9, H=11, I=13)
    return sum(multiplicity[k]*v for k, v in atombasis.items())
