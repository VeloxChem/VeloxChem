from setuptools import setup
from setuptools.command.build_py import build_py as SetuptoolsBuildPy
from setuptools.command.install import install as SetuptoolsInstall
import subprocess
import sys


class MyBuildPy(SetuptoolsBuildPy):

    def run(self):
        if run_process('make -C src'):
            SetuptoolsBuildPy.run(self)
        else:
            print('Error: failed to build veloxchem')
            sys.exit(1)


class MyInstall(SetuptoolsInstall):

    def run(self):
        if run_process('mkdir -p build/python/veloxchem'):
            SetuptoolsInstall.run(self)
        else:
            print('Error: failed to create folder build/python/veloxchem')
            sys.exit(1)


def run_process(command):
    process = subprocess.Popen(command.split(), stderr=subprocess.STDOUT)
    return process.wait() == 0


setup(
    name='VeloxChem',
    version='0.0',
    packages=['veloxchem'],
    package_dir={'veloxchem': 'build/python/veloxchem'},
    package_data={'veloxchem': ['veloxchemlib.so', 'basis/*']},
    scripts=['build/bin/VeloxChemMain.py'],
    cmdclass={
        'build_py': MyBuildPy,
        'install': MyInstall
    },
    description='VeloxChem: an electronic structure code',
)
