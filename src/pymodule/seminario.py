import numpy as np
class Seminario:
    def __init__(self,H,coords):
        self.H=H
        self.coords = coords

    def eig_sum(self,a,b):
        """
        Calculate the eigenvalues and eigenvectors of the 3x3 submatrix of H starting with the left corner at (3*a,3*b), and return the sum of the vectors multiplied by their respective eigenvalue
        """
        k=-self.H[3*a:3*(a+1),3*b:3*(b+1)]
        eig = np.linalg.eig(k)
        l = eig.eigenvalues
        v = eig.eigenvectors
        v = np.transpose(v)
        return np.sum(v * l[:, None], axis=0)
    
    def u(self,a,b):
        disp = self.coords[a]-self.coords[b]
        return disp/np.linalg.norm(disp)
    
    #Equation (11) in the seminario paper
    def u_N(self,a,b,c):
        uab=self.u(a,b)
        ubc=self.u(b,c)
        cross=np.cross(uab,ubc) #Signs don't matter because all normal-vectors only appear inside |x| terms
        return cross/np.linalg.norm(cross)

    def Rsq(self,a,b):
        return np.sum(np.square(self.coords[a]-self.coords[b]))

    #equation (17)
    def dihed_fc(self,a,b,c,d):
        """
        Calculate the dihedral force constant

        a,b(int): atom indices indicating the dihedral, starting at 0 and corresponding to the coords
        """
        u_ab = self.u(a,b)
        u_bc = self.u(b,c)
        u_cd = self.u(c,d)
        prefac1 = self.Rsq(a,b)*np.sum(np.square(np.cross(u_ab,u_bc)))
        prefac2 = self.Rsq(c,d)*np.sum(np.square(np.cross(u_bc,u_cd)))
        denom1=prefac1*np.linalg.norm(self.u_N(a,b,c)*self.eig_sum(a,b))
        denom2=prefac2*np.linalg.norm(self.u_N(b,c,d)*self.eig_sum(c,d))
        return 1/(1/denom1+1/denom2)

    #equation (14) 
    def angle_fc(self,a,b,c):
        """
        Calculate the angle force constant
        a,b,c(int): atom indices indicating the angle, starting at 0 and corresponding to the coords"""
        u_N=self.u_N(a,b,c)
        u_P_a=np.cross(u_N,self.u(a,b))
        u_P_c=np.cross(self.u(c,b),u_N)

        denom1=self.Rsq(a,b)*np.linalg.norm(u_P_a*self.eig_sum(a,b))
        denom2=self.Rsq(c,b)*np.linalg.norm(u_P_c*self.eig_sum(c,b))
        return 1/(1/(denom1)+1/(denom2))
        
    #equation (10)
    def bond_fc(self,a,b):
        """
        Calculate the bond force constant

        a,b(int): atom indices indicating the bond, starting at 0 and corresponding to the coords
        """
        u_ab = self.u(a,b)
        return np.linalg.norm(u_ab*self.eig_sum(a,b))