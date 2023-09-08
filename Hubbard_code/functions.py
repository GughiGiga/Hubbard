import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import h5py
import scipy.special
import matplotlib.pylab as pl
from scipy.optimize import curve_fit


def Fidelity(A, B):

    '''Takes two pure states stored in dictionaries and returns the fidelity between them'''
    
    def Inner(A, B):

        F = 0

        for key in A:
        
            if key not in B:
                continue
        
            F += A[key]*B[key]

        return F
    
    return Inner(A, B)**2


def Norm(state):

    '''Takes as input a state and normalizes it'''
    
    vec = state
    norm = Fidelity(state, state)**(1/4)
    for key in state:
        a = vec[key]/norm
        vec[key] = a
    return vec



def Tensor(vecs):
    '''Takes a list of states stored in dictionaries and returns their normalized tensor product as a dictionary'''

    def Tensor(A, B):

        '''Takes two pure states stored in dictionaries and returns their tensor product as a dictionary'''
        
        Tens = []
        Coeff = []
        Coeff2 = []
        Tens2 = []
        
        for key1 in A:
            for key2 in B:
                
                n1 = key1.replace('|', ' ', 1)
                nn1 = n1.strip()
                m1 = nn1.index('|')
                L11 = nn1[:m1]
                L12 = nn1[m1:-1]
                n2 = key2.replace('|', ' ', 1)
                nn2 = n2.strip()
                m2 = nn2.index('|')
                L21 = nn2[:m2]
                L22 = nn2[m2 + 1:]
                final = '|' + L11 + L21 + L12 + L22
                Tens.append(final)
                Coeff.append(A[key1]*B[key2])
                
                        
        return dict(zip(Tens, Coeff))
    
    Tens = Tensor(vecs[0], vecs[1])
    i = 2
    while i <= len(vecs) - 1:
        Tens = Tensor(Tens, vecs[i])
        i = i + 1

    return Norm(Tens)




def Eexpt(x, y, U, diag = False):

    '''Takes as input the dimesions of the lattice (x and y), a specific value of U and the charateristics of the lattice and returns the expectation value of the energy'''
    
    A = Data(x, y, U, diag)
    fp = A.get_static_obs()
    line = fp.readlines()[1]
    linee = line.split()[2].replace("{", " ")
    lineee = linee.replace("}", " ")
    
    return float(lineee)










def D_occ(x, y, site, U):

    '''Takes as input the dimensions of the lattice (x and y), a specific value of U and the charateristics of the lattice and returns the expectation value of the energy'''
    
    A = Data(x, y, U, diag)
    fp = A.get_static_obs()
    
    line = fp.readlines()[0].split()[1 + site]
    
    if site == 1:
        a = line.replace("{", " ")
        b = a.replace(",", " ")
        return float(b)
    
    elif site == x*y:
        a = line.replace("{", " ")
        b = a.replace(",", " ")
        return float(b)
    
    else:
        a = line.replace(",", " ")
        return float(a)









def RelAlign(x, y, site1, site2, U):

    '''Takes as input the dimensions of the lattice (x and y), the coordinate number of two sites, a specific value of U and the charateristics of the lattice and returns the expectation value of the relative spin alignment between the two chosen sites'''
    
    A = Data(x, y, U, diag)
    fp = A.get_static_obs()
    line = fp.readlines()[3].split()[site1 + 6*(site1 - 1) + site2 + 1]
    
    if site1 + 6*(site1 - 1) + site2 == 2:
        a = line.replace("{", " ")
        b = a.replace(",", " ")
        return float(b)
    
    elif site1 + 6*(site1 - 1) + site2 == (x*y)**2 - 1:
        a = line.replace("}", " ")
        b = a.replace(",", " ")
        return float(b)
    
    else:
        a = line.replace(",", " ")
        return float(a)







def Density(x, y, site, U, Diag = False):

    '''Takes as input the dimensions of the lattice (x and y), the coordinate number of a site, a specific value of U and the charateristics of the lattice and returns the expectation value of the fermion density on the chosen site.'''
    
    A = Data(x, y, U, diag)
    fp = A.get_static_obs()
    line = fp.readlines()[4].split()[site + 1]
    
    if site + 1 == 2:
        a = line.replace("{", " ")
        b = a.replace(",", " ")
        return float(b)
    
    elif site == (x*y) - 1:
        a = line.replace("}", " ")
        b = a.replace(",", " ")
        return float(b)
    
    else:
        a = line.replace(",", " ")
        return float(a)




def Symmetry(ket, symm):

    '''Takes as input an eigenket stored in a string and a tuple (symm) that specifies the new ordering of the sites and returns the eigenket corresponding to the transformation.'''
    
    s = ket.replace("|", " ", 1)
    sl = s.replace(">", " ")
    sll = sl.strip()
    n = sll.index('|')
    L = sll[0:n]
    L2 = sll[n + 1:]
    newL = np.zeros(len(L))
    newL2 = np.zeros(len(L2))
    
    for i, j in zip(range(len(L)), symm):
        newL[j] = L[i]
        newL2[j] = L2[i]
        
    newLL = [str(int(i)) for i in list(newL)]
    newLL2 = [str(int(i)) for i in list(newL2)]
    newstr = ''.join(newLL)
    newstr2 = ''.join(newLL2)
    final = '|' + newstr + '|' + newstr2 + ">"
    
    return final




def TotalSymmetry(state, symm):

    '''Takes as input a state stored in a dictionary and a tuple (symm) that specifies the new ordering of the sites and returns the state corresponding to the transformation.'''
    vec = state
    
    for key in vec:
        
        a = Symmetry(key, symm)
        
        if a not in vec:
            continue
        
        vec[key] = vec[a]
        
    return vec




def SymmetryCheck(state, symm, tol):

    '''Takes as input a state stored in a dictionary and a tuple (symm) that specifies the new ordering of the sites and returns whether the state is symmetric, antisymmetric or neither with respect to the transformation.'''
    
    L = []
    
    for key in state:
        a = Symmetry(key, symm)
        
        if a not in state:
            continue
                
        if abs(state[key] - state[a]) <= tol:
            L.append(1)
            
        elif abs(state[key] + state[a]) <= tol:
            L.append(-1)
            
        else:
            L.append(0)
            
    LL = np.array(L)
    
    if np.all(LL == np.full(len(LL), 1)):
        return 1
    
    elif np.all(LL == np.full(len(LL), -1)):
        return -1
    
    else:
        return 0










def Phase(file, mulower, muupper, nmu, symm):

    '''Takes as input a list of Data objects (file), a lower and an upper bound for the chemical potential mu and a different ordering 
    for the sites and returns the values of U and mu stored in a meshgrid, and the values of the particle number N of the ground-state 
    for each pair (U, mu), the truth-values of the symmetry of the ground-state, the truth-values of the type of the first excitation (spin-flip or N-change), the value of the first excitation energy and the boundaries between different regions stored in arrays  '''
    
    u = np.array(file.get_Uvals())
    nu = len(list(u))
    mu = np.linspace(mulower, muupper, nmu)
    MU, U = np.meshgrid(mu, u)
    Z = np.zeros((nu, nmu))
    sy = np.zeros((nu, nmu))
    ecc = np.zeros((nu, nmu))
    diff = np.zeros((nu, nmu))
    boundaries = np.zeros((nu, nmu))

    for i, U in enumerate(u):
        print('U =', U)
        eigvals = file.get_eigvals_all(U)
        sectors = file.get_sectors(U)
        Ntot = np.array([k + m for k, m in sectors])
        eigvecs = [file.get_eigvecs_with_eigvals(U, k) for k in eigvals]
        
        for j, Mu in enumerate(mu):
            print('mu =', Mu)
            En = eigvals - Ntot*Mu
            minn = min(En)
            n = list(En).index(minn)
            minn2 = min([i for i in En if i != minn])
            m = np.where(En == minn2)[0][0]
            if list(Ntot)[m] == list(Ntot)[n]:
                ecc[i, j] = 0
            else:
                ecc[i, j] = 1
            vec = eigvecs[n]
            sy[i, j] = SymmetryCheck(vec, symm, 0.0001)
            Z[i, j] = Ntot[n]
            if Z[i, j - 1] != Z[i, j]:
                boundaries[i, j] = mu[j]
            diff[i, j] = minn2 - minn

    return MU, U, Z, sy, ecc, diff, boundaries



def TensorPhase44(file, mulower, muupper, nmu, plaquette = '22'):

    '''Takes as input a list of Data objects (file) of the separate plaquettes and their dimensions (the plaquette argument), a lower and an upper bound for the chemical potential mu and a different ordering 
    for the sites and returns the values of U and mu stored in a meshgrid and the values of the particle number N of the ground-state for each pair (mu, U) stored in arrays.'''

    u = np.array(file.get_Uvals())
    nu = len(list(u))
    mu = np.linspace(mulower, muupper, nmu)
    MU, U = np.meshgrid(mu, u)
    Z = np.zeros((nu, nmu))
    sy = np.zeros((nu, nmu))
    ecc = np.zeros((nu, nmu))
    diff = np.zeros((nu, nmu))
    boundaries = np.zeros((nu, nmu))

    for i in range(len(file)):
        eigvals = file.get_eigvals_all(U)
        sectors = file.get_sectors(U)
        Ntot = np.array([k[0] + k[1] for k in sectors])
        print(sectors)
        Ntotnew = []
        Ennew = []
        Sectnew = []
        for s in range(len(Ntot)):
            a = sectors[s]
            if a in Sectnew:
                continue
            Sectnew.append(a)
            Ntotnew.append(Ntot[s])
            Ennew.append(eigvals[s])

        for j in range(len(mu)):

            En = []
            Sects = []
            Ntott = []

            if plaquette == '22':

                for k in range(len(Ennew)):
                    for l in range(len(Ennew)):
                        for m in range(len(Ennew)):
                            for p in range(len(Ennew)):
                                eig = Ennew[k] + Ennew[l] + Ennew[m] + Ennew[p] - (Ntotnew[k] + Ntotnew[l] + Ntotnew[m] + Ntotnew[p])*mu[j]
                                
                                if eig in En:
                                    continue

                                sec = (Sectnew[k][0] + Sectnew[l][0] + Sectnew[m][0] + Sectnew[p][0], Sectnew[k][1] + Sectnew[l][1] + Sectnew[m][1] + Sectnew[p][1])
                                N = Ntotnew[k] + Ntotnew[l] + Ntotnew[m] + Ntotnew[p]
                                En.append(eig)
                                
                                Sects.append(sec)
                                Ntott.append(N)
                minn = min(En)
                n = list(En).index(minn)
                print(j)
    #            minn2 = min([i for i in En if i != minn])
                Z[i, j] = Ntott[n]
    #            diff[i, j] = minn2 - minn

            if plaquette == '24' or plaquette == '42':

                for k in range(len(Ennew)):
                    for l in range(len(Ennew)):
                        eig = Ennew[k] + Ennew[l] - (Ntotnew[k] + Ntotnew[l])*mu[j]     
                        if eig in En:
                            continue
                        sec = (Sectnew[k][0] + Sectnew[l][0], Sectnew[k][1] + Sectnew[l][1])
                        N = Ntotnew[k] + Ntotnew[l]
                        En.append(eig)      
                        Sects.append(sec)
                        Ntott.append(N)
                minn = min(En)
                n = list(En).index(minn)
                print(j)
    #            minn2 = min([i for i in En if i != minn])
                Z[i, j] = Ntott[n]
    #            diff[i, j] = minn2 - minn
        
    return MU, U, Z, sy, ecc, diff, boundaries


def EnDiff(x, y, file, mulower, muupper, nmu, tol):

    '''Takes as input a list of Data objects (file) and the dimenisons of the lattice, a lower and an upper bound for the chemical potential mu and a different ordering 
    for the sites and returns the values of U and mu stored in a meshgrid and the values of the energy difference between any given state and the ground-state
    for each pair (mu, U) stored in arrays.'''
    
    nu = len(file.get_Uvals())
    u = file.get_Uvals()
    mu = np.linspace(mulower, muupper, nmu)
    MU, U = np.meshgrid(mu, u)
    track = np.zeros((nu, nmu, (x*y + 1)**2))
    MU = MU
    secs = []
    

    for p in range(x*y + 1):
        for q in range(x*y + 1):
            secs.append((p, q))
    
    for i, U in enumerate(list(u)):
        print('U = ', U)
        eigvals = file.get_eigvals_all(U)
        sectors = list(file.get_sectors(U))
        Ntot = np.array([i[0] + i[1] for i in file.get_sectors(U)])
        secss = sectors

        for j in range(len(mu)):
            print('mu = ', mu[j])          
            En = eigvals - Ntot*(MU[i,j])
            minn = min(En)
            n = list(En).index(minn)

            for l in range(len(secss)):
                diff = En[l] - minn
                print(diff)
                if secs[l] not in secss or diff < tol:
                    continue
                q = secss.index(secs[l])
                track[i, j, l] = diff

    return MU, U, track



def RDM(state, j):

	empty = 0
	up = 0
	down = 0
	double = 0

	for key in state:
		n = key[1:].index('|')
		if key[1:][j] == '0' and key[1:][n + 1 + j] == '0':
			empty += state[key]**2
		if key[1:][j] == '1' and key[1:][n + 1 + j] == '0':
			up += state[key]**2
		if key[1:][j] == '0' and key[1:][n + 1 + j] == '1':
			down += state[key]**2
		if key[1:][j] == '1' and key[1:][n + 1 + j] == '1':
			double += state[key]**2

	A = np.zeros((4, 4))
	A[0, 0] = empty
	A[1, 1] = up
	A[2, 2] = down
	A[3, 3] = double

	return A

from scipy.linalg import logm

def Entropy(state, j):

	
	from scipy.linalg import logm

	A = RDM(state, j)

	return -np.trace(np.matmul(A, logm(A)))/np.log(2)



