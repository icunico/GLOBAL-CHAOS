import ordpy as od
import os, sys
import numpy as np
import scipy.integrate
import pyfabm
import netCDF4 as nc
from scipy.fft import fft, fftfreq
from scipy.signal import find_peaks
import math
import subprocess
from scipy.stats import variation
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter

# Inizializza figura e assi
fig, ax = plt.subplots()

# Caricamento modello
model = pyfabm.Model('fabm_netcom_tot.yaml')
model.cell_thickness = 1.
assert model.checkReady(), 'Model not ready'

# Time derivative
def dy(t0, y):
        model.state[:] = y
        return model.getRates()

# Integrazione
t_ex =  3600 * 100000
t_eval = np.linspace(0, t_ex, 50000)
#t_eval = np.linspace(0, t_ex, 50000)
sol = scipy.integrate.solve_ivp(dy, [0., t_ex], model.state, t_eval=t_eval, max_step=5000)

t = sol.t
y = sol.y.T
y = y[40000:]
t_chaos =t[40000:]
#y = y[40000:]
#t_chaos =t[40000:]

# Estrazione  specie
specie_da_controllare_nomi = ['P1/DW1','P1/DW2','P1/DW3','C1/DWC1','C1/DWC2'] #'X/DWX1'
variable_names = [var.path for var in model.state_variables]
print(variable_names)
specie_da_controllare = [i for i, name in enumerate(variable_names) if name in specie_da_controllare_nomi]
print(specie_da_controllare_nomi)
# Flag negatività
flag_neg = False
for iv in specie_da_controllare:
    if np.any(y[:, iv] < -0.01):  
        flag_neg = True

# Statistiche
y_controllo = y[:, specie_da_controllare]
threshold = 0.01
mean_y = np.mean(y_controllo, axis=0)
std_y = np.std(y_controllo, axis=0)

# Chaos detection
taux = np.logspace(np.log(1), np.log10(len(t_chaos)//2), 50, dtype=int)
flag_chaos = False
cmax = 0.5

for iv in specie_da_controllare:
    T = y[:, iv]
    H = np.zeros(len(taux))
    C = np.zeros(len(taux))

    for it, tx in enumerate(taux):
        try:
            H[it], C[it] = od.complexity_entropy(T, taux=tx, dx=6)
        except:
            H[it], C[it] = np.nan, np.nan
             # Stampa i valori appena calcolati
    # Stampa H e C per la specie corrente
    # Se TUTTI sono NaN, allora non c'è caos
    for it in range(len(taux)):
         if not np.isnan(H[it]) and not np.isnan(C[it]):
             if 0.45 <= H[it] <= 0.7 and abs(C[it] - cmax) < 0.1 * cmax:
                 flag_chaos = True
                 break
    if flag_chaos:
        break

# Steady-state detection
ratio = np.where(mean_y >= threshold, std_y / mean_y, 0)
flag_ss = np.max(ratio) <= 0.01

# Extinction check
#flag_extin = False
#for iv in specie_da_controllare: 
#    if np.max(y[:, iv]) < 0.001:
#        print(f"Specie {iv} non ha mai superato la soglia (max = {np.max(y[:, iv]):.4f})")
#        flag_extin = True
#print(flag_extin)
# Nomi delle variabili nei gruppi
group_XYZ = ['P1/DW1']
group_C = ['C1/DWC1', 'C1/DWC2']
group_P = ['P1/DW1', 'P1/DW2', 'P1/DW3']
flag_extin = False
# Mappa: nome → indice in y (usando il path delle variabili del modello)
variable_names = [var.path for var in model.state_variables]
name_to_index = {name: i for i, name in enumerate(variable_names)}
threshold_ext=0.001
# Funzione per controllare estinzione: restituisce True se si estingue
def is_extinct(name):
    idx = name_to_index[name]
    return np.max(y[:, idx]) < threshold_ext

# Condizione 1: almeno uno tra X, Y, Z estinto
if any(is_extinct(name) for name in group_XYZ):
    flag_extin = True

# Condizione 2: entrambi C1 e C2 estinti
elif all(is_extinct(name) for name in group_C):
    flag_extin = True

# Condizione 3: tutti P1, P2, P3 estinti
elif all(is_extinct(name) for name in group_P):
    flag_extin = True

#print(f"flag_extin = {flag_extin}")
#	print(flag_chaos)
# Plot
#if flag_extin or flag_neg:
#    ax.scatter(dc, dx, c='k')
#elif flag_ss:
#    ax.scatter(dc, dx, c='g')
#elif flag_chaos:
#    ax.scatter(dc, dx, c='r')
#else:
#    ax.scatter(dc, dx, c='gold')

# Axis and labels
#ax.set_xlabel("dC1", fontsize=14)  
#ax.set_ylabel("tau", fontsize=14)  
#formatter = FuncFormatter(lambda x, _: f"{x:.2f}")
#ax.xaxis.set_major_formatter(formatter)
#ax.yaxis.set_major_formatter(formatter)

# Save figure
#fig.savefig("bubblef3.png")

# Frequenze per analisi (non usato ma lasciato)
Nt = t.shape[0]
deltaT = t[1] - t[0]
laststeps = 1000
freq = (1 / deltaT) * np.linspace(0, laststeps / 2, int(laststeps / 2)) / laststeps

# Save results
fileoutput = 'result.nc'
f = nc.Dataset(fileoutput, mode='w')
flag_dim = f.createDimension('flag', 1)

flag = f.createVariable('flag', np.int8, ('flag',))
flag.units = '1=True, 0=False'
flag.long_name = 'flag'
f.variables['flag'][:] = 1  

chaos = f.createVariable('flag_chaos', np.int8, ('flag',))
chaos.units = '1=True, 0=False'
chaos.long_name = 'flag_chaos'
f.variables['flag_chaos'][:] = int(flag_chaos)

ss = f.createVariable('flag_ss', np.int8, ('flag',))
ss.units = '1=True, 0=False'
ss.long_name = 'flag_ss'
f.variables['flag_ss'][:] = int(flag_ss)

extin = f.createVariable('flag_extin', np.int8, ('flag',))
extin.units = '1=True, 0=False'
extin.long_name = 'flag_extin'
f.variables['flag_extin'][:] = int(flag_extin)

neg = f.createVariable('flag_neg', np.int8, ('flag',))
neg.units = '1=True, 0=False'
neg.long_name = 'flag_neg'
f.variables['flag_neg'][:] = int(flag_neg)

for iv in specie_da_controllare:
    plt.plot(t_chaos, y[:, iv], label=variable_names[iv])  # `t` è il vettore dei tempi

plt.xlabel('Tempo')
plt.ylabel('Valore della specie')
plt.title('Dinamiche temporali delle specie selezionate')
plt.legend()
plt.grid(True)
plt.tight_layout()
plt.show()

f.close()

