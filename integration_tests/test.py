from mpi4py import MPI
from VeloxChemMP import *

app_manager = AppManager.create("test.inp", "test.out")

app_manager.execute()

assert app_manager.get_state() == True
