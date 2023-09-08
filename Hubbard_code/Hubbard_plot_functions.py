import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import h5py
import scipy.special
import matplotlib.pylab as pl
from scipy.optimize import curve_fit
from functions import *


'''Plotting functions associated to each general purpose function'''


def FidPlot(x, y, nup, ndown, Ufinal, nsteps, symm, eigvecnumber, eigvecnumber0, U0 = 0.001, diag = False):
    
    A = Data(x, y, U, diag)
    G0 = A.get_eigvecs(nup, ndown, eigvecnumber0)
    
    u = np.linspace(U0, Ufinal, nsteps)
    Fid = []
    Symmetry = []
    
    for i in u:
        H = Data(x, y, U, diag)
        h = H.get_eigvecs(nup, ndown, eigvecnumber)
        f = Fidelity(h, G0)
        Fid.append(f)
        s = SymmetryCheck(h, symm, 0.001)
        Symmetry.append(s)
    
    return Fid, Symmetry



def D_occPlot(x, y, site, Diag = False):
    
    Docc = [D_occ(x, y, site, k, Diag = False) for k in u]
    
    return Docc



def RelAlignPlot(x, y, site1, site2, Diag = False):
    
    Rel = [RelAlign(x, y, site1, site2, k, Diag = False) for k in u]
    
    return Rel



def N_effPlot(x, y, Diag = False):
    
    Neff = [N_eff(x, y, k, Diag = False) for k in u]
    
    return Neff


def EigvalsPlot(x, y, nup, ndown, eigvalnumber, diag = False):
    
    Eigg = [Eigvals(x, y, nup, ndown, k, diag = False)[eigvalnumber] for k in u]
    
    return Eigg


def colors(a):
    
    if a == 1:
        return 'blue'
    
    elif a == -1:
        return 'red'
    
    else:
        return 'magenta'


def SpinPlot(x, y, U, diag = False):
    
    A = Data(x, y, U, diag)
    eigvals = A.get_eigvals_all()
    sectors = A.get_sectors()
    Ntot = np.array([i + j for i, j in sectors])
    spincolor = [A.SpinEig(i, x, y, j, U) for i, j in zip(eigvals, Ntot)]
    En = eigvals - Ntot*10
    
    return Ntot, En, spincolor


def PhasePlot(file, mulower, muupper, nmu, symm):
    
    MU, U, Z, symm, ecc, diff, boundaries = Phase(file, mulower, muupper, nmu, symm)
    maxx = int(np.amax(Z))
    minn = int(np.amin(Z))
    L = [f'C{n}' for n in range(1, maxx + 2)]
    cmap = mpl.colors.ListedColormap(L)
    bounds = np.arange(minn-0.5, maxx + 1, 1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    
    return MU, U, Z, cmap, norm, symm, ecc, diff, boundaries


def TensorPhasePlot44(file, mulower, muupper, nmu, plaquette = '22'):
    
    MU, U, Z = TensorPhase44(file, mulower, muupper, nmu, plaquette)
    maxx = int(np.amax(Z))
    minn = int(np.amin(Z))
    cmap = plt.cm.tab20
    cmaplist = [cmap(i) for i in range(cmap.N)][:maxx + 1]
    cmap = mpl.colors.LinearSegmentedColormap.from_list('Custom cmap', cmaplist, cmap.N)
    bounds = np.arange(-0.5, maxx + 1.5, 1)
    norm = mpl.colors.BoundaryNorm(bounds, cmap.N)
    
    return MU, U, Z, cmap, norm