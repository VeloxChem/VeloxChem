from VeloxChemMP import *

inpdata   = CInputData()
xyzreader = CMolXYZReader()
mol       = CMolecule()

cout = COutputStream("test.out")
cin  = CInputStream("test.inp", cout)

cin.read(inpdata, cout)
xyzreader.parse(mol, inpdata, cout)

mol.print_geometry(cout)
cout.flush()

assert xyzreader.get_state() == True
