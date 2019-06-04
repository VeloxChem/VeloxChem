Getting started
===============

Sample input
------------

`water.inp`
::

    @jobs
    task: hf
    @end

    @method settings
    basis: def2-svp
    basis path: ../basis
    @end

    @molecule
    charge: 0
    multiplicity: 1
    units: au
    xyz:  
    O   0.0   0.0   0.0
    H   0.0   1.4   1.1
    H   0.0  -1.4   1.1
    @end 

Running
-------

::

    $ VeloxChem water.inp water.out



