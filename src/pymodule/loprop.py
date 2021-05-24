from mpi4py import MPI
import sys
import numpy as np
import os
from collections import Counter

from .scfrestdriver import ScfRestrictedDriver
from .lrsolver import LinearResponseSolver
from .inputparser import InputParser
from .outputstream import OutputStream
from .veloxchemlib import get_basis_function_indices_for_atom
from .veloxchemlib import ElectricDipoleIntegralsDriver



class LoPropDriver:    
     def __init__(self, task, comm=None, ostream=None):
        """
        Initializes trajectory driver.
        """
        self.task = task
        if comm is None:
            comm = MPI.COMM_WORLD
        if ostream is None:
            ostream = OutputStream(sys.stdout)

        self.comm = comm
        self.rank = self.comm.Get_rank()
        self.nodes = self.comm.Get_size()

        self.ostream = ostream
        
        
     def compute(self):
         #for checking
         np.set_printoptions(precision=4, suppress=True, linewidth=132)    
         """
         param: S:
             overlap matrix
         param: C:
             MO coefficient matrix
         param: D:
             density matrix
         param: natoms:
             the number of atoms in molecule/fragment
         param: norb:
             number of orbitals 
         param: ao_occ:
             occupied orbitals, ordered according to atoms then principal quantum number
         param: ao_vir:
             virtual orbitals
             
         """
         task = self.task
         molecule = task.molecule
         basis = task.ao_basis        
         natoms = molecule.number_of_atoms()
         
         scf_drv = ScfRestrictedDriver(self.comm, self.ostream)
         scf_drv.compute(molecule,basis)
         S = scf_drv.scf_tensors['S']
         C = scf_drv.scf_tensors['C']
         D = scf_drv.scf_tensors['D'][0]+ scf_drv.scf_tensors['D'][1]

         norb = np.shape(S)[0]
         
         #re-arrange AOs according to principal 
         re_arranged_indices=[]
         re_angmom=[]
         for i in range(molecule.number_of_atoms()):
             indices, angmoms = get_basis_function_indices_for_atom(
                molecule, basis, i)
             re_arranged_indices.append(indices)
             re_angmom.append(angmoms)
         re_arranged_indices=np.concatenate((re_arranged_indices))

        #obtain occupied&virtual orbital lists
         orb_per_atom,ao_occ, ao_vir=self.get_ao_indices()
         
         #TO transforation:re-arrange S     
         S0=S[re_arranged_indices,:][:,re_arranged_indices]

         #T1 transformation:Gram-Schmidt
         T1=np.zeros((norb,norb))
         ri=0
         for a in range(natoms):
             rf=ri+orb_per_atom[a]
             s = S0[ri:rf,ri:rf]
             L = np.linalg.cholesky(s)
             t = np.linalg.inv(L.T)  
             T1[ri:rf,ri:rf] = t
             ri += orb_per_atom[a]
         S1=T1.T@S0@T1
         #print('\n',T1,'\n')
         #print(S1)
         
         
        #T2 transformation:Lodwin occupied and virtual
         T2 = np.zeros((norb,norb))
         n_occ = len(ao_occ)
         n_vir = len(ao_vir)

         s = S1[ao_occ,:][:,ao_occ]
         t = self.lowdin_ortho(s)

         for r in range(n_occ):
            for c in range(n_occ):
                T2[ao_occ[r],ao_occ[c]] = t[r,c]

         s = S1[ao_vir,:][:,ao_vir]
         t = self.lowdin_ortho(s)

         for r in range(n_vir):
            for c in range(n_vir):
                T2[ao_vir[r],ao_vir[c]] = t[r,c]
         S2 = T2.T@S1@T2     
         
         #T3: projection to virtual
         T3 = np.identity((norb))
        #selected the overlap between occupied and virtual
         T3_ov = S2[ao_occ,:][:,ao_vir]
        #projektion
         ao_occupied=0
         for ao_o in ao_occ:
             ao_virtual=0
             for ao_v in ao_vir:
                 T3[ao_o, ao_v]= -T3_ov[ao_occupied, ao_virtual]
                 ao_virtual+=1
             ao_occupied+=1
         S3 = T3.T@S2@T3
         
         #T4: lodwin virtural
         T4_virtual=S3[ao_vir,:][:,ao_vir]
         T4_virtual=self.lowdin_ortho(T4_virtual)
         T4 = np.identity((norb))
         
         ao_virtual=0
         for ao_v in ao_vir:
            ao_virtual_1=0
            for ao_v_1 in ao_vir:
                T4[ao_v, ao_v_1]=T4_virtual[ao_virtual,ao_virtual_1]
                ao_virtual_1+=1
            ao_virtual+=1    
         #S4  = T4.T@S3@T4
    
        #total transformation T becomes:
         T=T1@T2@T3@T4
                       
         #obtain density matrix D in loprop basis set
         T_inv=np.linalg.inv(T)
         D_loprop= D[re_arranged_indices,:][:,re_arranged_indices]
         D_loprop=T_inv@D_loprop@T_inv.T
  
         #calculated localised charges
         Qab = np.zeros((natoms,natoms))
         start_point=0
         for atom in range(natoms):
             select=np.arange(0+start_point,orb_per_atom[atom]+start_point) 
             Qab[atom,atom]= -np.trace(D_loprop[select,:][:,select])
             start_point += orb_per_atom[atom]

         #nuclear charge
         nuclear_charge = self.task.molecule.elem_ids_to_numpy()
         for a in range(natoms):
               Qab[a,a]= Qab[a,a]+nuclear_charge[a]
               
        #solve linear response
         lrs_drv = LinearResponseSolver(self.comm, self.ostream)
         lrs_out=lrs_drv.compute(molecule,basis,scf_drv.scf_tensors)         
         
         #obtain response vectors
         Nx = lrs_out['solutions'][('x', 0)]
         Ny = lrs_out['solutions'][('y', 0)]
         Nz = lrs_out['solutions'][('z', 0)]
         
         nocc=task.molecule.number_of_alpha_electrons()
         #unpact response vectors to matrix form
         kappa_x = self.lrvec2mat(Nx, nocc, norb)
         kappa_y = self.lrvec2mat(Ny, nocc, norb)
         kappa_z = self.lrvec2mat(Nz, nocc, norb)
         
         #perturbed densities
         # factor of 2 from spin-adapted excitation vectors
         Dx = 2 * np.einsum('PQ,aP,bQ->ab', kappa_x, C, C)
         Dy = 2 * np.einsum('PQ,aP,bQ->ab', kappa_y, C, C)
         Dz = 2 * np.einsum('PQ,aP,bQ->ab', kappa_z, C, C)
         #convert to loprop basis
         Dx=Dx[re_arranged_indices,:][:,re_arranged_indices]
         Dy=Dy[re_arranged_indices,:][:,re_arranged_indices]
         Dz=Dz[re_arranged_indices,:][:,re_arranged_indices]
         Dx_loprop =T_inv@Dx@T_inv.T
         Dy_loprop =T_inv@Dy@T_inv.T
         Dz_loprop =T_inv@Dz@T_inv.T
         Dk_loprop=np.concatenate(([Dx_loprop],[Dy_loprop],[Dz_loprop]))
         
         #dipole
         dipole_drv = ElectricDipoleIntegralsDriver()
         dipole_mats = dipole_drv.compute(molecule, basis)
         x_ao = dipole_mats.x_to_numpy()
         y_ao = dipole_mats.y_to_numpy()
         z_ao = dipole_mats.z_to_numpy()
         #convert to loprop basis
         x_ao=x_ao[re_arranged_indices,:][:,re_arranged_indices]
         y_ao=y_ao[re_arranged_indices,:][:,re_arranged_indices]
         z_ao=z_ao[re_arranged_indices,:][:,re_arranged_indices]
         x_ao_loprop= T.T@x_ao@T
         y_ao_loprop= T.T@y_ao@T
         z_ao_loprop= T.T@z_ao@T
         dipole=np.concatenate(([x_ao_loprop],[y_ao_loprop],[z_ao_loprop]))
         
         '''
         localised polarisabilities:
             aAB = -rAB delta DAB + dQab (Ra-Rb)
             
             where dQab is obtained by solving the following Lagragian with 
             minimal charge transfer, here a penaty function is also introduced
             
             contribution from localised dipoles: rAB*delta DAB = Aab
                                 
         '''
         
         #dQa:charge shift per atom
         dQa = np.zeros((natoms,3))
         #dQa array as [[Ax,Ay,Az],[Bx,By,Bz]]
         it=0
         for a in range(natoms):
             select=np.arange(0+it,orb_per_atom[a]+it)
             dQa[a][0]= np.trace(Dx_loprop[select,:][:,select])
             dQa[a][1]= np.trace(Dy_loprop[select,:][:,select])
             dQa[a][2]= np.trace(Dz_loprop[select,:][:,select])
             it= it + orb_per_atom[a]
             
         #r
         molecule_coord=molecule.get_coordinates()
         r=np.zeros((natoms,natoms,3))
         for i in range(natoms):
             for j in range(natoms):
                 if i == j:
                     #if a=b rab=ra
                     r[i][j]= molecule_coord[i]
                 else:
                     #if a =!b, rab = (ra-rb)/2 
                     a =np.absolute(molecule_coord[i]-molecule_coord[j])
                     r[i][j]=a/2


         #contribution from localised dipole
         Aab=np.zeros((3,3, natoms,natoms))

        # 3 stands for x y z directions
         for i in range(3):
            for j in range(3):
                it_a=0
                for a in range(natoms):
                    select_a=np.arange(0+it_a,orb_per_atom[a]+it_a)             
                    it_b=0
                    for b in range(natoms):              
                        #select the subblock[a][b] region in dks_lp
                        select_b=np.arange(0+it_b ,orb_per_atom[b]+it_b)     
                        #selected the lp basis for subblock[A][B] in purterbed density matrix             
                        D_AB = Dk_loprop[i][select_a,:][:,select_b]              
                        #selected the dipole matrice for subblock[A][B] in purterbed density matrix    
                        dipole_select= dipole[j][select_a,:][:,select_b]
                        dipole_select=dipole_select.transpose()           
                        #-delta(rc) = -rAB * delta D_AB          
                        Aab[i,j,a,b]=  np.trace(dipole_select @ D_AB)                
                        it_b = it_b + orb_per_atom[b]
                
                    Aab[i,j,a,a] -= dQa[a,j]* r[a,a,i]               
                    it_a=it_a + orb_per_atom[a]

        #Lagragian
         Fab=np.zeros((natoms,natoms))
         for a in range(natoms):
             za= nuclear_charge[a]
             Ra=molecule_coord[a]
             for b in range(a):
                zb=nuclear_charge[b]
                ##read the coordinate
                Rb=molecule_coord[b]
                Fab[a,b]=self.pf(za,Ra,zb,Rb)
                Fab[b,a]=Fab[a][b]
             for a in range(natoms):
                Fab[a,a] += -sum(Fab[a,:])
        
         Lab=Fab+ self.sf(Fab)
         
         dQa=np.swapaxes(dQa,0,1)      
         la = [np.linalg.solve(Lab,rhs) for rhs in dQa]
     
        #dQab
         dQab=np.zeros((natoms,natoms,3))
         for i in range(3):
             for a in range(natoms):
                 za =nuclear_charge[a]
                 Ra = molecule_coord[a]
                 for b in range(a):
                     zb=nuclear_charge[b]
                     Rb=molecule_coord[b]
                     dQab[a,b,i]= -(la[i][a]-la[i][b])*self.pf(za,Ra,zb,Rb)
                     dQab[b,a,i]-= dQab[a,b,i]

         # dRab matrix: mid point of Ra-Rb
         dRab=np.zeros((natoms,natoms,3))

         for a in range(natoms):
            for b in range(natoms):  
                for i in range(3):
                    dRab[a][b][i] =(molecule_coord[a][i]- molecule_coord[b][i])/2
                    dRab[b][a][i] =- dRab[a][b][i]

        #dAab charge transfer from bond polarisability
         dAab = np.zeros((3,3,natoms,natoms))
         for a in range(natoms):
             for b in range(natoms):
                 for i in range(3):
                     for j in range(3): 
                         dAab[i,j,a,b]=dRab[a,b,i]*dQab[a,b,j]+dRab[a,b,j]*dQab[a,b,i]

         local_polarisabilities= Aab + 0.5 * dAab
         #molecular polarisabilities
         Am = (Aab + 0.5 * dAab)  .sum(axis=3).sum(axis=2)
         
         self.print_results(natoms,Qab, local_polarisabilities, Am)

    #shift function
     def sf(self,F:np.ndarray)->float:
        return 2*np.max(np.abs(F))                
     #penaty function 
     def pf(self,za,Ra,zb,Rb):
         #Ra/Rb are the position vectors     
        AUG2AU=1.8897259886
        RBS=np.array([0,0.25,0.25,1.45,1.05,0.85,0.7,0.65,0.5,0.43,1.8,1.5,1.25,1.1,1.0,1.0,1.0,1.0])*AUG2AU
        ra=RBS[za]
        rb=RBS[zb]
        #Ra and Rb taking from r_mat: 0,1,2 represents x y z directions
        rab2=(Ra[0]-Rb[0])**2 + (Ra[1]-Rb[1])**2+(Ra[2]-Rb[2])**2
        #alpha set to 2
        #a minus sign is already in the expression
        f=0.5*np.exp(-2*(rab2/(ra+rb)**2))
        return f   
         
    #unpact response vector to matrix form
     def lrvec2mat(self,vec, nocc, norb):
         kappa = np.zeros((norb,norb))    
         i = 0
         for r in range(nocc):
            for c in range(nocc, norb):
                kappa[r,c] = vec[i]
                i+=1            
         kappa = kappa + kappa.T    
         return kappa     
  
    #lowddin orthonomalisation
     def lowdin_ortho(self,A):
        eigs, U = np.linalg.eigh(A)
        X = np.einsum("Pk,k,Qk->PQ", U, 1.0/np.sqrt(eigs), U)
        return X     
 
     def get_ao_indices(self):

         elements = self.task.molecule.elem_ids_to_numpy()
         basis = self.task.ao_basis.get_label()
       
         #get basis file         
         local_name = basis
         basis_file=[]
         lib_member = os.path.join(os.environ['VLXBASISPATH'], basis)
         if os.path.exists(local_name):
              basis_file= os.path.abspath(local_name)
         elif os.path.exists(lib_member):
              basis_file = lib_member
        
              
         bp = InputParser(basis_file)
         basis = bp.get_dict()
         keys = list(basis.keys())
          
         ao_occ = []
         ao_vir = []
         iterr=0
         multiplicity = dict(S=1, P=3, D=5, F=7, G=9, H=11, I=13)
         orb_per_atom=[]
         for e in elements:
             k = keys[e - 1]
             atoms_data = basis[k]
             c = Counter(line[0] for line in atoms_data if line and line[0] in "SPDFGHI")
          
            # For H and He: 1s            
             ao_occ.append(iterr)

            # For Li-Ne: + 2s 2p
             if e >= 3:
                ao_occ.append(iterr+1)
                offset_p = c['S']+iterr
                ao_occ.append(offset_p + 0)
                ao_occ.append(offset_p + 1)
                ao_occ.append(offset_p + 2)    
            
             #sum no. orbitals. used to calculate virtual orbitals later
             orb=sum(multiplicity[k] * v for k, v in c.items())
             iterr=iterr+orb
             orb_per_atom.append(orb)
             
             
         total_orb_array=np.arange(iterr)
            
         for i in total_orb_array:
             if i not in ao_occ:
                 ao_vir.append(i)    
         
         return orb_per_atom,ao_occ, ao_vir


         
     def print_results(self,natoms,Qab, local_polarisabilities, Am):
                 
        """
        Prints results for loprop calculation.
        """
        element_names=self.task.molecule.get_labels()    
        width=92
        
        self.ostream.print_blank()
        self.ostream.print_header('Local Properties(loprop) Calculations')
        self.ostream.print_header(19 * '=')
        self.ostream.print_blank()
        
        output_header ='*** Molecular Polarisabilities *** '
        self.ostream.print_header(output_header.ljust(width))  
        
        direction=["x","y","z"]
        for a in range(3):
            output_iter = '<<{};{}>>: {:15.8f} '.format(direction[a],direction[a],Am[a,a])
            self.ostream.print_header(output_iter.ljust(width))
        self.ostream.print_blank()
      
        #print localised chagres
        output_header ='*** Loprop localised charges *** '
        self.ostream.print_header(output_header.ljust(width))       
        self.ostream.print_blank()
        for a in range(natoms):
            output_iter = '{:<5s}: {:15.8f} '.format(element_names[a], Qab[a,a])
            self.ostream.print_header(output_iter.ljust(width))
            
        output_header ='*** Localised polarisabilities *** '
        self.ostream.print_header(output_header.ljust(width))       
        self.ostream.print_blank()

        self.ostream.print_blank()
        
        
        self.ostream.flush()
