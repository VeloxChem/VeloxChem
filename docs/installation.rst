Install
=======

From source
^^^^^^^^^^^

Get the source code
+++++++++++++++++++

The source code can be downloaded from this repository

.. code-block:: bash

    https://gitlab.com/veloxchem/veloxchem

Debian based Linux
++++++++++++++++++

.. code-block:: bash

    # Install MKL from
    # https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-apt-repo

    # For example:
    # wget https://apt.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
    # apt-key add GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
    # sudo sh -c 'echo deb https://apt.repos.intel.com/mkl all main > /etc/apt/sources.list.d/intel-mkl.list'
    # sudo apt-get update
    # sudo apt-get install intel-mkl-2018.2-046

    # Install MPICH (or OpenMPI) and Python

    sudo apt-get install mpich python3 python3-dev python3-pip

    # Install python modules

    sudo pip3 install numpy h5py pybind11 pytest loprop psutil
    sudo pip3 install --no-binary=mpi4py mpi4py

    # Install VeloxChem

    cd VeloxChem
    python3 config/generate_setup.py
    python3 setup.py install

RPM based Linux
+++++++++++++++

.. code-block:: bash

    # Install MKL from
    # https://software.intel.com/en-us/articles/installing-intel-free-libs-and-python-yum-repo

    # For example:
    # sudo yum-config-manager --add-repo https://yum.repos.intel.com/mkl/setup/intel-mkl.repo
    # sudo rpm --import https://yum.repos.intel.com/intel-gpg-keys/GPG-PUB-KEY-INTEL-SW-PRODUCTS-2019.PUB
    # sudo yum install intel-mkl-2018.2-046

    # Install MPICH (or OpenMPI) and Python
    # Note: need to manually link MPI executables

    sudo yum install mpich-3.2-devel python3-devel
    sudo ln -s /usr/lib64/mpich-3.2/bin/mpirun /usr/bin/mpirun
    sudo ln -s /usr/lib64/mpich-3.2/bin/mpicxx /usr/bin/mpicxx
    sudo ln -s /usr/lib64/mpich-3.2/bin/mpicc /usr/bin/mpicc

    # Install python modules

    sudo pip3 install numpy h5py pybind11 pytest loprop psutil
    sudo pip3 install --no-binary=mpi4py mpi4py

    # Install VeloxChem

    cd VeloxChem
    python3 config/generate_setup.py
    python3 setup.py install

PowerLinux
++++++++++

.. code-block:: bash

    # Install OpenBLAS

    wget https://github.com/xianyi/OpenBLAS/archive/v0.3.4.tar.gz
    tar xf v0.3.4.tar.gz
    cd OpenBLAS-0.3.4
    make TARGET=POWER8 CC=gcc FC=gfortran USE_OPENMP=1
    make PREFIX=<path-to-your-openblas> install
    export OPENBLASROOT=<path-to-your-openblas>

    # Install MPICH (or OpenMPI)

    wget http://www.mpich.org/static/downloads/3.2.1/mpich-3.2.1.tar.gz
    tar xf mpich-3.2.1.tar.gz
    cd mpich-3.2.1
    ./configure --disable-fortran --prefix=<path-to-your-mpich> CC=gcc CXX=g++
    make && make install
    export PATH=<path-to-your-mpich>/bin:$PATH
    export LD_LIBRARY_PATH=<path-to-your-mpich>/lib:$LD_LIBRARY_PATH

    # Install Anaconda (Python 3.7 version) for Power8 and Power9 from
    # https://www.anaconda.com/distribution/

    # For example:
    # wget https://repo.anaconda.com/archive/Anaconda3-2019.10-Linux-ppc64le.sh
    # bash Anaconda3-2019.10-Linux-ppc64le.sh

    # Install python modules

    pip install numpy h5py pybind11 pytest loprop psutil
    pip install --no-binary=mpi4py mpi4py

    # Install VeloxChem

    cd VeloxChem
    python3 config/generate_setup.py
    python3 setup.py install

Arm (Isambard Cray XC50 system)
+++++++++++++++++++++++++++++++

.. code-block:: bash

    # Load cray modules

    module swap PrgEnv-cray PrgEnv-gnu
    module load cray-hdf5
    module load cray-python/3.6.5.7

    # Install python modules

    python3 -m pip install pybind11 --user
    python3 -m pip install h5py --user
    python3 -m pip install loprop psutil --user

    # Manually install mpi4py

    # 1. Download mpi4py-3.0.3.tar.gz from https://pypi.org/project/mpi4py/#files
    # 2. tar xf mpi4py-3.0.3.tar.gz && cd mpi4py-3.0.3
    # 3. Add the following lines to mpi.cfg
    # [cray]
    # mpicc         = cc
    # mpicxx        = CC
    # extra_compile_args   = -shared
    # extra_link_args      = -Wl,-rpath,/opt/cray/pe/mpt/7.7.9/gni/mpich-gnu/8.2/lib
    python3 setup.py build --mpi=cray
    python3 setup.py install --prefix=<path-to-your-mpi4py>
    export PYTHONPATH=<path-to-your-mpi4py>/lib/python3.6/site-packages:$PYTHONPATH

    # Setup compiler wrapper

    export CXX=CC

    # Install VeloxChem

    cd VeloxChem
    python3 config/generate_setup.py
    python3 setup.py install

MacOS
+++++

.. code-block:: bash

    # Download and install MKL from
    # https://software.intel.com/en-us/mkl

    # Install libomp, MPICH and Python

    brew install libomp
    brew install mpich
    brew install python

    # Install python modules

    pip3 install numpy h5py pybind11 pytest loprop psutil
    pip3 install --no-binary=mpi4py mpi4py

    # Install VeloxChem

    cd VeloxChem
    python3 config/generate_setup.py
    python3 setup.py install


Windows
+++++++

Soon to come!


Binaries
^^^^^^^^

Docker
++++++

A docker image with pre-compiled veloxchem based on ubuntu:18.04 is available
on `Docker Hub <https://hub.docker.com/r/veloxchem/veloxchem>`_.

.. code-block:: bash

    $ docker run -it veloxchem/veloxchem:1.0rc1
    # root@fcc794d899c7:/veloxchem# which vlx
    /usr/local/bin/vlx

The CPPE library for polarizable embedding
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
    export PYTHONPATH=<path-to-your-cppe>/build/stage/lib:$PYTHONPATH

    # Make sure that cppe can be imported
    python3 -c 'import cppe'

