import h5py
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.special
import matplotlib.pylab as pl
from scipy.optimize import curve_fit
import glob
from functions import *
from Hubbard_Data import *
from Hubbard_plot_functions import *

'''A44 = Data(4, 4, 'periodic', 0)

MU, U, Z, cmap, norm, symm, ecc, diff, boundaries = PhasePlot(A44, -10, 0, 1000, (15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3 ,2 ,1 ,0))

for i, obj in enumerate([MU, U, Z, cmap, norm, symm, ecc, diff, boundaries]):
    np.save(f'/export/ggigante/arr{i}', obj)


A44 = Data(4, 4, 'periodic', 0)

MU, U, Z = EnDiff(4, 4, A44, -10, 0, 1000, 0.25)

for i, obj in enumerate([MU, U, Z]):
    np.save('/export/ggigante/Endiff_arr{i}.npy', obj)'''

'''A = []

for i in range(9):
    A.append(np.load(f'/export/ggigante/arr{i}.npy', allow_pickle=True))

print(A[4])

MU, U = np.meshgrid(np.array([0.001, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]), np.linspace(-10, 0, 1000))

plt.pcolormesh(MU, U, A[2].T, cmap = 'tab20' , shading = 'nearest')

plt.colorbar(ticks = np.linspace(0, 16, 17))

plt.show()

plt.pcolormesh(MU, U, A[5].T, cmap = 'jet' , shading = 'nearest')

plt.colorbar(ticks = np.linspace(0, 16, 17))

plt.show()

for i in range(13):

    fig, ax1 = plt.subplots()

    ax1.scatter(np.linspace(-10, 0, 1000), A[-2][i,:], c = A[-3][i,:])

    ax2 = ax1.twiny()

    ax2.set_xticks(A[2][i,:])

    for j in range(1000):
        ax1.axvline(A[-1][i, j])
    plt.show()'''
    
A44 = Data(4, 4, 'periodic', 0)

A = EnDiff(4, 4, A44, -10, 0, 1000, 0.25)

for i in range(13):
    for j in range(len(A[2][0,0,:])):
        plt.scatter(np.linspace(-10, 0, 1000), A[2][i,:,j])
    plt.show()