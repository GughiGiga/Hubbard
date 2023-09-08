from Hubbard_Data import *
from Hubbard_plot_functions import *
from functions import *
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import h5py
import scipy.special



U22 = Data(2, 2, 'open', 0).get_Uvals()

Vecs22 = [Data(2, 2, 'open', 0).get_eigvecs(u, 2, 2, 0) for u in U22]

plt.scatter(U22, [Entropy(i, 0) for i in Vecs22])

plt.show()
'''A24 = [Data(2, 4, i, 0.0, 0, periodic = False, u100 = True) for i in np.linspace(0.001, 100, 100)[:11]]
A24periodic = [Data(2, 4, i, 0.0, 0, periodic = True, u100 = True) for i in np.linspace(0.001, 100, 100)[:11]]

A22 = [Data(2, 2, i, 0.0, 0, periodic = False, u100 = True) for i in np.linspace(0.001, 100, 100)[:11]]
A22periodic = [Data(2, 2, i, 0.0, 0, periodic = True, u100 = True) for i in np.linspace(0.001, 100, 100)[:11]]'''

'''def Histogram(File):
    A = File.get_eigvecfile()
    S = File.get_sectors()
    Contr = []
    with open(A, 'r') as f:
        read = f.readlines()
        for line in read:
            a = line.split()
            if a[0] == 'Contributions':
                n = read.index(line)
                for line2 in read[n:]:
                    if line2[0] == 'E':
                        break
                m = read[n:].index(line2)
                Contr.append(float(read[n:][m - 1].split()[-1]))
            print(Contr)
    Contr[-1] = 1.0
    return dict(zip([f'{i[0]} {i[1]}' for i in S], Contr))

from matplotlib.ticker import FuncFormatter




def Contr(x, pos):
    'The two args are the value and tick position'
    #return '$%1.1fM' % (x * 1e-6)
    return f'{x}'
formatter = FuncFormatter(Contr)




def HistoPlot(file, U):

    A = Histogram(file)
    Contributions = np.array([A[key] for key in A])
    Secs = np.array([key for key in A])

    np.save(f'/export/Contributions_U={U}.npy', Contributions)
    np.save(f'/export/Contributions_U={U}.npy', Secs)

    return Secs, Contributions



for i in [0.001, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12]:

    filee = Data(4, 4, i, 0.0, 0, periodic = True)
    Secs, Contributions = HistoPlot(filee, i)
    
    fig, ax = plt.subplots()
    ax.yaxis.set_major_formatter(formatter)
    ax.bar(x=np.arange(len(Contributions)), # The x coordinates of the bars. 
        height=Contributions, # the height(s) of the vars 
        color="green", 
        align="center",
        tick_label=Secs)
    ax.set_ylabel('Contribution', fontsize = 15)
    ax.set_title('Coverage of the Hilbert space for each sector for 300000 prints', fontsize = 15)
    ax.set_xlabel('Sectors', fontsize = 15)
    plt.show()
'''


'''print(list(A))





PL = dict(zip([f'{i[0]} {i[1]}' for i in A44.sectors], [0.918681, 0.905782, 0.905782, 0.923996, 0.923995, 0.96157, 0.96157, 0.99298, 0.99298, 1.0, 1.0, 0.836766, 0.833463, 0.836778, 0.834973, 0.840172, 0.838818, 0.837577, 0.839378, 0.904894, 0.906168, 0.904426, 0.903994, 1.0, 1.0, 0.957867, 0.957997, 0.957751, 0.958152, 0.99726, 0.997282, 0.997468, 0.997382, 0.999934, 0.999868, 1.0, 1.0, 1.0, 1.0, 0.871882, 0.87282, 0.873024, 0.872379, 0.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.974994, 0.976616, 0.976321, 0.974954, 0.9375, 0.9375, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.775919, 0.77068, 0.7702, 0.76824, 0.914475, 0.911322, 0.916576, 0.910667, 0.999194, 0.999179, 0.99914, 0.99917, 0.999999, 1.0, 0.999987, 0.999987, 1.0, 1.0, 1.0, 1.0, 0.632457, 0.70859, 0.700901, 0.706044, 0.70166, 0.832892, 0.832892, 0.832891, 0.832892, 1.0, 1.0, 1.0, 1.0, 0.985676, 0.984776, 0.985615, 0.986777, 1.0, 1.0, 1.0, 1.0, 0.582835, 0.570763, 0.561845, 0.573524, 1.0, 1.0, 1.0, 1.0, 0.655778, 0.656545, 0.639533, 0.661673, 0.965445, 0.964408, 0.962947, 0.964373, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.515972, 0.515972, 0.620637, 0.610286, 0.636625, 0.625623, 0.729587, 0.729587, 0.729586, 0.729587, 0.950872, 0.953476, 0.953096, 0.949892, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.678585, 0.678161, 0.672138, 0.67645, 1.0, 1.0, 1.0, 1.0, 0.786351, 0.761655, 0.867625, 0.847432, 0.855761, 0.867574, 0.959, 0.966344, 0.95637, 0.961426, 1.0, 1.0, 1.0, 1.0, 0.953528, 0.94923, 0.949064, 0.997901, 0.997777, 0.947174, 0.99771, 0.997657, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 0.989497, 0.989497, 0.999987, 0.999999, 1.0, 0.999987, 0.875, 0.9375, 0.0, 1.0]))



Sectors = list(PL.keys())


Contributions = [PL[i] for i in Sectors if PL[i] < 0.95]

Secs = [key for key in PL if PL[key] in Contributions]

fig, ax = plt.subplots()
ax.yaxis.set_major_formatter(formatter)
ax.bar(x=np.arange(len(Contributions)), # The x coordinates of the bars. 
       height=Contributions, # the height(s) of the vars 
       color="green", 
       align="center",
       tick_label=Secs)
ax.set_ylabel('Contribution', fontsize = 15)
ax.set_title('Coverage of the Hilbert space for each sector for 300000 prints', fontsize = 15)
ax.set_xlabel('Sectors', fontsize = 15)
plt.show()'''

'''MU, U, Z, cmap, norm, symm, ecc, diff, boundaries  = PhasePlot(A44, -10, 0, 10, (15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3, 2, 1, 0))''''''

for i in len([MU, U, Z, symm, ecc, diff, boundaries]):
    np.save(f'/export/Phase_{i}.npy', [MU, U, Z, symm, ecc, diff, boundaries][i])



plt.pcolormesh(A[0], A[1], A[2], cmap = A[3], norm = A[4])
plt.colorbar(ticks = np.linspace(0, 16, 17))
plt.title('Phase diagram of the tensor product of two 2x4 plaquettes (open BC)', fontsize = 10)
plt.xlabel(r'$/mu$', fontsize = 15)
plt.ylabel('U', fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.xlim(-10, 0)


A = TensorPhasePlot44(A24periodic, -10, 0, 1000, plaquette = '24')

plt.pcolormesh(A[0], A[1], A[2], cmap = A[3], norm = A[4])

plt.title('Phase diagram of the tensor product of two 2x4 plaquettes (periodic BC)', fontsize = 10)
plt.xlabel(r'$/mu$', fontsize = 15)
plt.ylabel('U', fontsize = 15)
plt.xticks(fontsize = 15)
plt.yticks(fontsize = 15)
plt.xlim(-10, 0)
plt.savefig('/export/TensorPhase2424_periodic')'''



'''


S = [Data(4, 4, i, 0.0, 0, periodic = True).get_sectors() for i in [0.001, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]

arr = [Data(4, 4, i, 0.0, 0, periodic = True).get_eigvals_all() for i in [0.001, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]

R = [list(zip(S[i], arr[i])) for i in range(len(S))]








for i in range(len(R)):
    file = open(f'/export/eig-sec_small_U={i}.txt','w')
    for j in range(len(R[i])):
	    file.write(f'{R[i][j]}'+"\n")
    file.close()
'''


'''MU, U, Z, cmap, norm, symm, ecc, diff, boundaries = PhasePlot(A44, -10, 0, 1000, (15, 14, 13, 12, 11, 10, 9, 8, 7, 6, 5, 4, 3 ,2 ,1 ,0))

for i, obj in enumerate([MU, U, Z, cmap, norm, symm, ecc, diff, boundaries]):
    np.save(f'/export/ggigante/arr{i}', obj)


A = [np.load(f'/export/ggigante/arr{i}.npy', allow_pickle=True) for i in range(9)]

print(A[5])

'''






