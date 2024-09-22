from os import environ, cpu_count
from pathlib import Path

from .mklconf import configure_mkl_rt


def get_basis_path():
    """
    Returns location of basis files within module.

    :return:
        The location of basis files within module.
    """

    return Path(__file__).parent / "basis"


def set_vlxbasispath():
    """
    Sets location of basis files.
    """

    if 'VLXBASISPATH' not in environ:
        environ['VLXBASISPATH'] = str(get_basis_path())


def set_omp_num_threads(ncores=None):
    """
    Sets number of OpenMP threads.

    :param ncores:
        The number of cores available for OpenMP threads.
    """

    if ncores is not None:
        environ['OMP_NUM_THREADS'] = str(ncores)
    else:
        if 'OMP_NUM_THREADS' not in environ:
            ncores = cpu_count()
            if ncores is None:
                ncores = 1
            environ['OMP_NUM_THREADS'] = str(ncores)
