import numpy as np
class Seminario:
    
    def __init__(self,H,coords):
        self.H=H
        self.coords = coords


    def eig_sum(self,a,b,vec):
        """
        Calculate the eigenvalues and eigenvectors of the 3x3 submatrix of H starting with the left corner at (3*a,3*b), and return the sum of the vectors multiplied by their respective eigenvalue
        """
        k=-self.H[3*a:3*(a+1),3*b:3*(b+1)]
        eig = np.linalg.eig(k)
        l = eig.eigenvalues
        v = eig.eigenvectors
        v = np.transpose(v)

        sum = 0
        for i in range(3):
            sum+=l[i]*abs(np.dot(vec,v[i]))
        return sum
    
    def u(self,a,b):
        disp = self.coords[a]-self.coords[b]
        return disp/np.linalg.norm(disp)
    
    #Equation (11) in the seminario paper
    def u_N(self,a,b,c):
        uab=self.u(a,b)
        ucb=self.u(c,b)
        cross=np.cross(ucb,uab) #Signs don't matter because all normal-vectors only appear inside |x| terms
        return cross/np.linalg.norm(cross)

    def Rsq(self,a,b):
        return np.sum(np.square(self.coords[a]-self.coords[b]))

    #equation (17)
    def dihed_fc(self,a,b,c,d):
        """
        Calculate the dihedral force constant

        a,b(int): atom indices indicating the dihedral, starting at 0 and corresponding to the coords
        """
        def calc_par(i,j,k,l):
            u_ij = self.u(i,j)
            u_jk = self.u(j,k)
            u_kl = self.u(k,l)
            prefac1 = self.Rsq(i,j)*np.sum(np.square(np.cross(u_ij,u_jk)))
            prefac2 = self.Rsq(k,l)*np.sum(np.square(np.cross(u_jk,u_kl)))

            denom1=prefac1*self.eig_sum(i,j,self.u_N(i,j,k))
            denom2=prefac2*self.eig_sum(k,l,self.u_N(j,k,l))
            
            return (1/(1/denom1+1/denom2)).real
        
        fc = max(0,0.5*(calc_par(a,b,c,d)+calc_par(a,c,b,d)))

        return fc

    #equation (14) 
    def angle_fc(self,a,b,c):
        """
        Calculate the angle force constant
        a,b,c(int): atom indices indicating the angle, starting at 0 and corresponding to the coords"""
        def calc_par(i,j,k):
            u_N=self.u_N(i,j,k)
            u_P_i=np.cross(u_N,self.u(i,j))
            u_P_k=np.cross(self.u(k,j),u_N)

            denom1=self.Rsq(i,j)*self.eig_sum(i,j,u_P_i)
            denom2=self.Rsq(k,j)*self.eig_sum(k,j,u_P_k)

            return 1/(1/(denom1)+1/(denom2)).real
        return max(0,0.5*(calc_par(a,b,c)+calc_par(c,b,a)))
    #equation (10)
    def bond_fc(self,a,b):
        """
        Calculate the bond force constant

        a,b(int): atom indices indicating the bond, starting at 0 and corresponding to the coords
        """
        def calc_par(i,j):
            u_ab = self.u(i,j)
            return self.eig_sum(i,j,u_ab).real
        
        return max(0,0.5*(calc_par(a,b)+calc_par(b,a)))