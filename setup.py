from setuptools import setup
from setuptools.command.build_py import build_py as SetuptoolsBuildPy
from setuptools.command.install import install as SetuptoolsInstall
from config.generate_setup import generate_setup
import multiprocessing
import subprocess
import sys
from pathlib import PurePath
from pathlib import Path
from os import environ


class MyBuildPy(SetuptoolsBuildPy):

    def run(self):
        self.update_setup_file()
        python_exe = 'python{}.{}'.format(sys.version_info[0],
                                          sys.version_info[1])
        self.update_vlx_script(python_exe)
        self.make_veloxchem()
        SetuptoolsBuildPy.run(self)
        self.update_vlx_script('python3')

    def make_veloxchem(self):
        cmd = 'make -C src -s -j '
        if 'OMP_NUM_THREADS' in environ:
            cmd += environ['OMP_NUM_THREADS']
        else:
            cmd += str(multiprocessing.cpu_count())
        cmd += ' release'
        process = subprocess.Popen(cmd.split(),
                                   stdout=sys.stdout,
                                   stderr=subprocess.STDOUT)
        if process.wait() != 0:
            print('Error: failed to build veloxchem')
            sys.exit(1)

    def update_setup_file(self):
        setup_file = Path('src', 'Makefile.setup')
        if not setup_file.is_file():
            template_file = Path('config', 'Setup.template')
            generate_setup(str(template_file), str(setup_file))

    def update_vlx_script(self, python_exe):
        vlx_file = Path('src', 'vlx')
        with open(str(vlx_file), 'w', encoding='utf-8') as f_vlx:
            print('#!/usr/bin/env {}'.format(python_exe), file=f_vlx)
            print('from veloxchem.main import main', file=f_vlx)
            print('main()', file=f_vlx)


class MyInstall(SetuptoolsInstall):

    def run(self):
        SetuptoolsInstall.run(self)


setup(
    name='veloxchem',
    version='1.0rc1.post1',
    packages=[
        'veloxchem',
    ],
    package_dir={
        'veloxchem': str(PurePath('src', 'pymodule')),
    },
    package_data={
        'veloxchem': [
            'veloxchemlib.so',
            str(PurePath('basis', '*')),
        ],
    },
    scripts=[
        str(PurePath('src', 'vlx')),
    ],
    python_requires='>=3.5',
    install_requires=[
        'pybind11>=2.5.0',
        'mpi4py>=3.0',
        'numpy>=1.13',
        'h5py>=2.8',
        'loprop>=0.2.3',
        'psutil>=5.6.7',
    ],
    tests_require=[
        'pytest>=5.1.2',
        'pytest-cov>=2.7.1',
    ],
    test_suite='python_tests',
    cmdclass={
        'build_py': MyBuildPy,
        'install': MyInstall,
    },
    description='VeloxChem: an electronic structure code',
    author='VeloxChem developers',
    author_email='info@veloxchem.org',
    url='http://veloxchem.org/',
)
