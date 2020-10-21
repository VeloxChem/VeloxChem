from multiprocessing import cpu_count
from os import environ
from pathlib import Path
from sys import stdout


def set_vlxbasispath():
    if 'VLXBASISPATH' not in environ:
        module_path = Path(__file__).parent
        environ['VLXBASISPATH'] = str(Path(module_path) / 'basis')


def set_omp_num_threads(ncores=None):
    if ncores is not None:
        environ['OMP_NUM_THREADS'] = ncores
    else:
        if 'OMP_NUM_THREADS' not in environ:
            ncores = cpu_count()
            environ['OMP_NUM_THREADS'] = str(ncores)
            print('* Warning * Environment variable OMP_NUM_THREADS not set.',
                  file=stdout)
            print('* Warning * Setting OMP_NUM_THREADS to {:d}.'.format(ncores),
                  file=stdout)
