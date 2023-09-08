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

A44 = Data(4, 4, 'periodic', 0)

MU, U, Z = EnDiff(4, 4, A44, -10, 0, 1000, 0.25)

for i, obj in enumerate([MU, U, Z]):
    np.save('/export/ggigante/Endiff_arr{i}.npy', obj)