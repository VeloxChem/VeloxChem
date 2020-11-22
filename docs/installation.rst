Installation
============

Obtaining the source code
^^^^^^^^^^^^^^^^^^^^^^^^^

The source code can be downloaded from this repository

.. code-block:: bash

    https://gitlab.com/veloxchem/veloxchem

Installing from source
^^^^^^^^^^^^^^^^^^^^^^

With Anaconda (x86-64 or POWER processor)
+++++++++++++++++++++++++++++++++++++++++

- Create virtual environment ::

    conda create --name vlx python=3.7
    conda activate vlx

- Install compiler and math library

  - on Linux (x86-64 processor) ::

        conda install gxx_linux-64 mkl mkl-include

  - on Linux (POWER processor) ::

        conda install gxx_linux-ppc64le openblas=0.3.3
        export OPENBLASROOT=/path/to/your/anaconda3/envs/vlx

  - on MacOS (x86-64 processor) ::

        conda install clangxx_osx-64 mkl mkl-include

- Install MPI, mpi4py and other Python modules ::

    conda install openmpi
    python3 -m pip install mpi4py --no-binary=mpi4py

    conda install pybind11 psutil pytest numpy h5py
    python3 -m pip install loprop

- Compile VeloxChem ::

    cd VeloxChem
    python3 config/generate_setup.py
    python3 setup.py install

- Use VeloxChem in Jupyter notebook ::

    conda install jupyter
    ipython kernel install --name vlx
    jupyter notebook
    # choose "vlx" in the drop-down menu under "New"

Cray platform (x86-64 or ARM processor)
+++++++++++++++++++++++++++++++++++++++

- Load cray modules ::

    module swap PrgEnv-cray PrgEnv-gnu
    module load cray-hdf5
    module load cray-python/3.6.5.7

- Install mpi4py

    1. Download ``mpi4py-3.0.3.tar.gz`` from https://pypi.org/project/mpi4py/#files
    2. ``tar xf mpi4py-3.0.3.tar.gz && cd mpi4py-3.0.3``
    3. Append the following lines to ``mpi.cfg`` ::

        [cray]
        mpicc         = cc
        mpicxx        = CC
        extra_compile_args   = -shared
        extra_link_args      = -Wl,-rpath,/opt/cray/pe/mpt/7.7.9/gni/mpich-gnu/8.2/lib

    4. Build and install mpi4py ::

        python3 setup.py build --mpi=cray
        python3 setup.py install --prefix=/path/to/your/mpi4py
        export PYTHONPATH=/path/to/your/mpi4py/lib/python3.6/site-packages:$PYTHONPATH

- Install other Python modules ::

    python3 -m pip install pybind11 h5py psutil loprop --user

- Use compiler wrapper to compile VeloxChem ::

    export CXX=CC

    cd VeloxChem
    python3 config/generate_setup.py
    python3 setup.py install

Debian based Linux
++++++++++++++++++

- Install MKL from https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-apt-repo ::

    wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
    apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
    sudo sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
    sudo apt-get update
    sudo apt-get install intel-mkl-2018.2-046

- Install MPI, Python, mpi4py and other Python modules ::

    sudo apt-get install mpich python3 python3-dev python3-pip
    sudo pip3 install --no-binary=mpi4py mpi4py
    sudo pip3 install numpy h5py pybind11 pytest loprop psutil

- Install VeloxChem ::

    cd VeloxChem
    python3 config/generate_setup.py
    python3 setup.py install

RPM based Linux
+++++++++++++++

- Install MKL from https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-yum-repo ::

    sudo yum-config-manager --add-repo https://yum.repos.intel.com/mkl/setup/intel-mkl.repo
    sudo rpm --import https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
    sudo yum install intel-mkl-2018.2-046

- Install MPI, Python, mpi4py and other Python modules ::

    sudo yum install mpich-3.2-devel python3-devel
    sudo ln -s /usr/lib64/mpich-3.2/bin/mpirun /usr/bin/mpirun
    sudo ln -s /usr/lib64/mpich-3.2/bin/mpicxx /usr/bin/mpicxx
    sudo ln -s /usr/lib64/mpich-3.2/bin/mpicc /usr/bin/mpicc
    sudo pip3 install --no-binary=mpi4py mpi4py
    sudo pip3 install numpy h5py pybind11 pytest loprop psutil

- Install VeloxChem ::

    cd VeloxChem
    python3 config/generate_setup.py
    python3 setup.py install

PowerLinux
++++++++++

- See installation instructions `With Anaconda (x86-64 or POWER processor)`_

MacOS
+++++

- See installation instructions `With Anaconda (x86-64 or POWER processor)`_

Windows
+++++++

- Soon to come!


Installing binaries
^^^^^^^^^^^^^^^^^^^

Docker
++++++

A docker image with pre-compiled veloxchem based on ubuntu:18.04 is available
on `Docker Hub <https://hub.docker.com/r/veloxchem/veloxchem>`_.

.. code-block:: bash

    $ docker run -it veloxchem/veloxchem:1.0rc1
    # root@fcc794d899c7:/veloxchem# which vlx
    /usr/local/bin/vlx

Dependencies
^^^^^^^^^^^^

The CPPE library for polarizable embedding
++++++++++++++++++++++++++++++++++++++++++

If you are interested in using the CPPE library for polarizable embedding,
please install it according to the instructions below. Note that cmake is
needed to build the CPPE library.

.. code-block:: bash

    # Build CPPE
    git clone https://github.com/maxscheurer/cppe
    cd cppe; mkdir build; cd build
    cmake -DENABLE_PYTHON_INTERFACE=ON ..
    make

    # Set up python path
    export PYTHONPATH=/path/to/your/cppe/build/stage/lib:$PYTHONPATH

    # Make sure that cppe can be imported
    python3 -c 'import cppe'

