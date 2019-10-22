Install
=======


Conda
^^^^^

It is recommended to install veloxchem in a separate conda environment, e.g. 
::

    $ conda create -n vlx
    $ conda activate vlx
    (vlx) $ conda -c conda-forge install veloxchem

Conda binary distributions have been generated with 
 
* Ubuntu 18.04 LTS.

From source
^^^^^^^^^^^

Ubuntu (Intel CPU)
++++++++++++++++++

.. code-block:: bash

    # Install MKL
    # https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-apt-repo

    wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
    apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB

    sudo sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
    sudo apt-get update
    sudo apt-get install intel-mkl-2018.2-046

    # Install MPICH and Python

    sudo apt-get install mpich python3 python3-dev python3-pip

    # Install numpy, h5py, pybind11, pytest and mpi4py

    pip3 install numpy h5py pybind11 pytest
    pip3 install --no-binary=mpi4py mpi4py

    # Install VeloxChem

    git clone https://gitlab.com/rinkevic/VeloxChemMP.git
    cd VeloxChemMP
    python3 config/generate_setup.py
    python3 setup.py install

MacOS (Intel CPU)
+++++++++++++++++

.. code-block:: bash

    # Install MKL
    # https://software.intel.com/en-us/mkl

    # Install libomp, MPICH and Python

    brew install libomp
    brew install mpich
    brew install python

    # Install numpy, h5py, pybind11, pytest and mpi4py

    pip3 install numpy h5py pybind11 pytest
    pip3 install --no-binary=mpi4py mpi4py

    # Install VeloxChem

    git clone https://gitlab.com/rinkevic/VeloxChemMP.git
    cd VeloxChemMP
    python3 config/generate_setup.py
    python3 setup.py install
