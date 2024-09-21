from .veloxchemlib import NuclearPotentialDriver
from .veloxchemlib import ElectricDipoleMomentDriver


def compute_nuclear_potential_integrals(molecule, basis):
    """
    Computes nuclear potential integrals.

    :param molecule:
        The molecule.
    :param basis:
        The molecular basis set.

    :return:
        The nuclear potential integral matrix.
    """

    npot_drv = NuclearPotentialDriver()
    npot_mat = npot_drv.compute(molecule, basis)

    # Note: factor -1.0 for electron charge
    return -1.0 * npot_mat.full_matrix().to_numpy()


def compute_electric_dipole_integrals(molecule, basis, origin=(0.0, 0.0, 0.0)):
    """
    Computes electric dipole integrals.

    :param molecule:
        The molecule.
    :param basis:
        The molecular basis set.

    :return:
        A tuple containint the electric dipole integral matrices.
    """

    dip_drv = ElectricDipoleMomentDriver()
    dip_mats = dip_drv.compute(molecule, basis, list(origin))

    return tuple([
        # Note: factor -1.0 for electron charge
        -1.0 * dip_mats.matrix('X').full_matrix().to_numpy(),
        -1.0 * dip_mats.matrix('Y').full_matrix().to_numpy(),
        -1.0 * dip_mats.matrix('Z').full_matrix().to_numpy(),
    ])
