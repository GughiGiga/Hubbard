from Hubbard_Data import *
from Hubbard_plot_functions import *
from functions import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import h5py
import scipy.special

s = (0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15)



A44 = [Data(4, 4, 'periodic', 0).get_eigvecs(u, 8, 8, 0) for u in [0.001, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]]

S44 = [TotalSymmetry(i, s) for i in A44]

A22 = [Data(2, 2, 'periodic', 0).get_eigvecs(u, 2, 2, 0) for u in [0.001, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]]

T44 = [Tensor([i, i, i, i]) for i in A22]

Fid44 = [Fidelity(i, j) for i, j in zip(S44, T44)]

plt.scatter([0.001, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12], Fid44, markers = 'x')

plt.xlabel('U', fontsize = 15)
plt.ylabel('Fidelity', fontsize = 15)

plt.title('Fidelity between tensor product of four 2x2 plaquettes and 4x4 lattice', fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)

plt.savefig('/export/ggigante/Hubbard_figures/4x4/2x2tensx4.png')

plt.show()

