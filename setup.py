from setuptools import setup
from setuptools.command.build_py import build_py as SetuptoolsBuildPy
from setuptools.command.install import install as SetuptoolsInstall
from config.generate_setup import generate_setup
import multiprocessing
import subprocess
import sys
import os


class MyBuildPy(SetuptoolsBuildPy):

    def run(self):
        self.check_setup_file()
        python_exe = 'python{}.{}'.format(sys.version_info[0],
                                          sys.version_info[1])
        self.update_vlx_script(python_exe)
        self.make_veloxchem()
        SetuptoolsBuildPy.run(self)
        self.update_vlx_script('python3')

    def make_veloxchem(self):
        cmd = 'make -C src -s -j '
        if 'OMP_NUM_THREADS' in os.environ:
            cmd += os.environ['OMP_NUM_THREADS']
        else:
            cmd += str(multiprocessing.cpu_count())
        process = subprocess.Popen(cmd.split(),
                                   stdout=sys.stdout,
                                   stderr=subprocess.STDOUT)
        if process.wait() != 0:
            print('Error: failed to build veloxchem')
            sys.exit(1)

    def check_setup_file(self):
        setup_file = os.path.join('src', 'Makefile.setup')
        if not os.path.isfile(setup_file):
            template_file = os.path.join('config', 'Setup.template')
            print('*** Generating Makefile.setup...')
            generate_setup(template_file, setup_file)

    def update_vlx_script(self, python_exe):
        vlx_file = os.path.join('src', 'vlx')
        with open(vlx_file, 'w', encoding='utf-8') as f_vlx:
            print('#!/usr/bin/env {}'.format(python_exe), file=f_vlx)
            print('from veloxchem.main import main', file=f_vlx)
            print('main()', file=f_vlx)


class MyInstall(SetuptoolsInstall):

    def run(self):
        SetuptoolsInstall.run(self)


setup(
    name='veloxchem',
    version='1.0rc1post1',
    packages=[
        'veloxchem',
    ],
    package_dir={
        'veloxchem': os.path.join('src', 'pymodule'),
    },
    package_data={
        'veloxchem': [
            'veloxchemlib.so',
            os.path.join('basis', '*'),
        ],
    },
    scripts=[
        os.path.join('src', 'vlx'),
    ],
    python_requires='>=3.5',
    install_requires=[
        'pybind11>=2.3.0',
        'mpi4py>=3.0',
        'numpy>=1.13',
        'h5py>=2.8',
        'loprop>=0.2.3',
        'psutil>=5.7.0',
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
