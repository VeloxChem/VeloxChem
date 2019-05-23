from setuptools import setup
from setuptools.command.build_py import build_py as SetuptoolsBuildPy
from setuptools.command.install import install as SetuptoolsInstall
from config.generate_setup import generate_setup
import subprocess
import sys
import os


class MyBuildPy(SetuptoolsBuildPy):

    def run(self):
        self.check_setup_file()
        self.make_veloxchem()
        SetuptoolsBuildPy.run(self)

    def make_veloxchem(self):
        process = subprocess.Popen('make -C src'.split(),
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


class MyInstall(SetuptoolsInstall):

    def run(self):
        package_dir = os.path.join('build', 'python', 'veloxchem')
        if not os.path.isdir(package_dir):
            try:
                os.makedirs(package_dir)
            except IOError:
                raise IOError(
                    'Unable to create package dir {}'.format(package_dir))
        SetuptoolsInstall.run(self)


setup(
    name='VeloxChem',
    version='0.0',
    packages=[
        'veloxchem',
    ],
    package_dir={
        'veloxchem': os.path.join('build', 'python', 'veloxchem'),
    },
    package_data={
        'veloxchem': [
            'veloxchemlib.so',
            os.path.join('basis', '*'),
        ],
    },
    scripts=[
        os.path.join('build', 'bin', 'VeloxChemMain.py'),
    ],
    python_requires='>=3.5',
    install_requires=[
        'mpi4py>=3.0',
        'numpy>=1.13',
        'h5py>=2.8',
    ],
    cmdclass={
        'build_py': MyBuildPy,
        'install': MyInstall,
    },
    description='VeloxChem: an electronic structure code',
)
