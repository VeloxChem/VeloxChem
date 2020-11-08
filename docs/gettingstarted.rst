Getting started
===============

Structure of input file
-----------------------

The VeloxChem input file consists of multiple groups marked with
``@group name`` and ``@end``. For example, the following input file
has three groups: ``jobs``, ``method settings``, and ``molecule``.

.. code-block:: bash

  @jobs
  task: scf
  @end

  @method settings
  xcfun: b3lyp
  basis: def2-svp
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

.. toctree::
   :maxdepth: 2

   inputs/keywords.rst

Example input files
-------------------

.. toctree::
   :maxdepth: 2

   inputs/scf.rst
   inputs/pe.rst
   inputs/rsp.rst
   inputs/mp2.rst
   inputs/cube.rst

Running with input file
-----------------------

::

    $ export OMP_NUM_THREADS=6
    $ mpirun -n 2 python3 -m veloxchem water.inp water.out

Running in Jupyter notebook
---------------------------

::

    from mpi4py import MPI
    import veloxchem as vlx
    import sys

    molecule_string = """
        O 0 0 0
        H 0 0 1.795239827225189
        H 1.693194615993441 0 -0.599043184453037"""

    basis_set_label = 'def2-svp'

    scf_settings = {'conv_thresh': 1.0e-6}
    method_settings = {'xcfun': 'b3lyp', 'grid_level': 4}

    comm = MPI.COMM_WORLD
    ostream = vlx.OutputStream(sys.stdout)

    molecule = vlx.Molecule.read_str(molecule_string, units='angs')
    basis = vlx.MolecularBasis.read(molecule, basis_set_label)

    ostream.print_block(molecule.get_string())
    ostream.print_block(basis.get_string('Atomic Basis', molecule))
    ostream.flush()

    scfdrv = vlx.ScfRestrictedDriver(comm, ostream)
    scfdrv.update_settings(scf_settings, method_settings)
    scfdrv.compute(molecule, basis)
