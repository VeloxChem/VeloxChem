Linear Response
---------------

+ Polarizability

  .. literalinclude:: polar.inp
     :language: bash


+ Random phase approximation (RPA)

  .. literalinclude:: rpa.inp
     :language: bash


+ Complex polarization propagator (CPP)

  .. literalinclude:: cpp.inp
     :language: bash


+ Linear response with pulsed external field: Correction to time-domain dipole moment

The current version of VeloxChem allows for calculation of the linear electric dipole response for the frequency region of a single Gaussian-envelope pulse deemed sufficiently large (using the `@pulses` module). The time-domain shape parameters for this pulse may be specified by the user, optionally storing a collection of pertinent results in an HDF5-formatted file or a plaintext ASCII file whose name may be specified by the user. An example of an input file that when run will carry out such a calculation is given below. For more documentation about the available keywords, please consult the source file whose path from the VeloxChem root folder is `/src/pymodule/pulsedrsp.py`. Note in particular that the default of carrier envelope phase may need adjustment to match your desired setup.

Sample input:

  .. literalinclude:: linrsp_pulsed.inp
     :language: bash

If HDF5-formatted data was produced during this calculation, that data may used for plot generation using the script located at `utils/pulsed_response_plot.py` from the VeloxChem root folder. Please note that Python version 3 is required to run this script. Please also note that other standard python modules such as `matplotlib` must be installed on the system from which this script is run. The script will take the HDF5-formatted data produced during the VeloxChem calculation and generate a plot of the real and imaginary frequency-domain electric dipole polarizability, a representation of the perturbing field in the frequency domain, the resulting (real-valued) first-order dipole moment correction in the time domain and a representation of the perturbing field in the time domain. For more information and further description of how to run this script, please consult the documentation written inside it.
