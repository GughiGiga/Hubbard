import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.special
import matplotlib.pylab as pl
from scipy.optimize import curve_fit
import glob




def vecpathnew(x, y, BC, tp, ts = False):
    if ts == True:
        return glob.glob(f'/export/ggigante/Hubbard_Sim/{x}x{y}/' + BC + f'_tp{tp}_ts/*/eigvecs.txt')
    else:
        return glob.glob(f'/export/ggigante/Hubbard_Sim/{x}x{y}/' + BC + f'_tp{tp}/*/eigvecs.txt')


def valpathnew(x, y, BC, tp, ts = False):
    if ts == True:
        return glob.glob(f'/export/ggigante/Hubbard_Sim/{x}x{y}/' + BC + f'_tp{tp}_ts/*/eigvals.h5')
    else:
        return glob.glob(f'/export/ggigante/Hubbard_Sim/{x}x{y}/' + BC + f'_tp{tp}/*/eigvals.h5')





class Data:

    def __init__(self, x, y, BC, tp, ts = False):
        
        self.x = x
        self.y = y
        self.BC = BC
        self.tp = tp
        self.ts = ts
        self.valfile = valpathnew(x, y, BC, tp, ts)
        self.vecfile = vecpathnew(x, y, BC, tp, ts)
        self.U = []
        self.ts = []
        for i in self.vecfile:
            n = i.index('U')
            p = i[n:].index('e')
            self.U.append(float(i[n:][1:p - 1]))
            if ts == True:
                n = i.index('U')
                q = i[n:].index('t')
                p = i[n:].index('e')
                self.U.append(float(i[n:][1:q]))
                self.ts.append(float(i[n:][q + 2:p]))

        self.U.sort()
        self.ts.sort()
                
                
    def __repr__(self):
        return f'{self.x}x{self.y} with ' + self.BC + ' boundary conditions'
        
    def get_static_obs(self):
        '''Returns the static_observables file of the object'''
        return self.vecfile
    
    def get_Uvals(self):
        '''Returns the values of U of the object'''
        return self.U

    def get_ts(self):
        '''Returns the values of ts of the object'''
        return self.ts
    
    def get_tsvals(self):
        '''Returns the values of U of the object'''
        return self.ts

    def get_nupsndowns(self):
        '''Returns the values of U of the object'''
        with h5py.File(self.valfile[o], 'r') as fp:
            N = int(np.array(fp['results/eigenvalues/N']))
            nup = [int(np.array(fp[f'results/eigenvalues/sectors/{i}/nup'])) for i in range(N)]
            ndown = [int(np.array(fp[f'results/eigenvalues/sectors/{i}/ndown'])) for i in range(N)]
        return nup, ndown
    
    def get_eigvals_all(self, U):
        '''Returns all the eigenvalues of the object in ascending order'''
        o = self.U.index(float(U))
        val = 0
        with h5py.File(self.valfile[o], 'r') as fp:
            val = np.array(fp['results/eigenvalues/data'])
        return val
    
    def get_sectors(self, U):
        '''Returns all the sectors in the same order as the eigenvalues'''
        o = self.U.index(float(U))
        sectors = 0
        with h5py.File(self.valfile[o], 'r') as fp:
            N = int(np.array(fp['results/eigenvalues/N']))
            nup = [int(np.array(fp[f'results/eigenvalues/sectors/{i}/nup'])) for i in range(N)]
            ndown = [int(np.array(fp[f'results/eigenvalues/sectors/{i}/ndown'])) for i in range(N)]
            sectors = list(zip(nup, ndown))
        return sectors

    def get_sectors_with_Ntot(self, U, Ntot):
        '''Returns all the sectors for a specific value of the total particle number Ntot'''
        secs = self.get_sectors(U)
        secsnew = []
        for i in secs:
            if i[0] + i[1] == Ntot:
                secsnew.append(i)
        return secsnew

    def get_eigvecs(self, U, nup, ndown, eigvecnumber):

        '''Takes as input the spin occupation numbers and the number of the state (eigvecnumber) and returns the corresponding wavefunction stored in a dictionary'''
        o = self.U.index(float(U))
        
        coeff = []
        eigkets = []
            
        with open(self.vecfile[o], 'r') as fp:

            lines = []
            
            if eigvecnumber > 3:
                raise ValueError("Invalid eigenvector number")
            counter = 0
            for n, line1 in enumerate(fp):
                a = line1.split()
                if a[0] == 'Eigenvector' and int(a[5]) == nup and int(a[6]) == ndown:
                    counter += 1
                if counter - 1 == eigvecnumber:
                    lines.append(a)
                if a[0] == 'Contributions':
                    break       
            for i in lines[1:-1]:
                coeff.append(float(i[0]))
                eigkets.append(i[2])

        Z = dict(zip(eigkets, coeff))

        return Z
    
    def get_eigvecs_all(self, U, tol):

        o = self.U.index(float(U))
        
        '''Returns all the wavefunctions in the same order as the eigenvalues'''

        Total = []

        eigvals = self.get_eigvals_all(U)
        nup, ndown = self.get_nupsndowns(U)
        
        for i, j, k in zip(nup, ndown, eigvals):
        
            coeff = []
            eigkets = []

            with open(self.vecfile, 'r') as fp:

                if eigvecnumber > 3:
                    raise ValueError("Invalid eigenvector number")
                counter = 0
                lines = []
                for n, line1 in enumerate(fp):
                    a = line1.split()
                    if a[0] == 'Eigenvector' and int(a[5]) == i and int(a[6]) == j and abs(float(a[4]) - k) < 10**(-4):
                        counter += 1
                    if counter - 1 == eigvecnumber:
                        lines.append(a)
                    if a[0] == 'Contributions':
                        break       
                for i in lines[1:-1]:
                    coeff.append(float(i[0]))
                    eigkets.append(i[2])

            Z = dict(zip(eigkets, coeff))
            Total.append(Z)
        
        return Total
    
    
    
    def get_eigvecs_with_eigvals(self, U, ts, eigval):

        '''Returns the wavefunction associated to a specific eigenvalue'''
        o = self.U.index(float(U))
        
        coeff = []
        eigkets = []
            
        with open(self.vecfile[o], 'r') as fp:

            lines = []
            counter = 0
            for n, line1 in enumerate(fp):
                a = line1.split()
                if a[0] == 'Eigenvector' and abs(float(a[4]) - eigval) < 10**(-4):
                    counter += 1
                if counter == 1:
                    lines.append(a)
                if a[0] == 'Contributions':
                    break       
            for i in lines[1:-1]:
                coeff.append(float(i[0]))
                eigkets.append(i[2])

        Z = dict(zip(eigkets, coeff))

        return Z
        
    
    
    def get_eigvals(self, U, nup, ndown):

        '''Returns the eigenvalue for a specific sector stored in a list'''
        
        f = self.valfile
        eigvalss = []

        for i in range(self.N):

            if self.nups[i] == nup and self.ndowns[i] == ndown:
                eigval = self.eigvals[i]
                eigvals.append(float(eigval))

        return self.eigvalss

    
    def SpinEig(self, U, eigval, Ntot):

        '''Takes as input a specific eigenvalue, the total particle number and returns the corresponding value of the total spin S'''

        if Ntot%2 == 0:

            NUP = 0
            NDOWN = 0

            for i in range(int(Ntot/2) + 1):

                nup = Ntot/2 - i
                ndown = Ntot/2 + i

                for l, m in zip(self.get_eigvals(nup, ndown), self.get_eigvals(ndown, nup)):
                    if abs(eigval - l) < 10**(-4) or abs(eigval - m) < 10**(-4):
                        NUP = nup
                        NDOWN = ndown

        if Ntot%2 == 1:

            for i in range(int((Ntot + 1)/2)):

                nup = (Ntot + 1)/2 - i
                ndown = (Ntot - 1)/2 + i
                NUP = 0
                NDOWN = 0

                for l, m in zip(self.get_eigvals(nup, ndown), self.get_eigvals(ndown, nup)):
                    if abs(eigval - l) < 10**(-4) or abs(eigval - m) < 10**(-4):
                        NUP = nup
                        NDOWN = ndown

        return abs(NUP - NDOWN)/2
    

    def get_contributions(self, U, nup, ndown, eigvecnumber):

        '''Returns the number of printed terms in the wavefunction expansion and the fraction of the total state covered by the printed one (squared norm of the printed state)'''
        
        N = []
        contributions = []
        
        with open(self.vecfile, 'r') as fp:
        
            if eigvecnumber > 3:
                raise ValueError("Invalid eigenvector number")
            L = fp.readlines()
            counter = 0
            line2 = []
            for line1 in L:
                a = line1.split()
                if a[0] == 'Contributions' and int(a[-2]) == nup and int(a[-1]) == ndown:
                    counter += 1
                else:
                    continue
                if counter - 1 == eigvecnumber:
                    break
            n = L.index(line1)
            if n == 0:
                raise ValueError('No eigenvalues available')
            for line2 in L[n:]:
                b = line2.split()
                if b[0] == 'Eigenvector':
                    break

            m = L.index(line2)
            K = L[n + 1:m]

            for i in K:
                S = i.split()
                N.append(float(S[0]))
                contributions.append(float(S[-1]))


        return N, contributions
