import ordpy as od
import nolds as nd
import numpy as np
import scipy.integrate
import pyfabm
import netCDF4 as nc
from scipy.stats import variation
import matplotlib.pyplot as plt

# Caricamento modello
model = pyfabm.Model('fabm_netcom_tot.yaml')
model.cell_thickness = 1.
assert model.checkReady(), 'Model not ready'

# Inizializza flag
flag_extin = False
flag_neg = False
flag_chaos = False
flag_ss = False

# Estrai le specie da controllare
specie_da_controllare_nomi = ['P1/DW1', 'C1/DWC1', 'X/DWX1']
variable_names = [var.path for var in model.state_variables]
specie_da_controllare = [i for i, name in enumerate(variable_names) if name in specie_da_controllare_nomi]

# Contatori eventi
event_counters = {iv: 0 for iv in specie_da_controllare}

# Funzione derivata temporale
def dy(t0, y):
    model.state[:] = y
    return model.getRates()

# Funzione evento per estinzione
def event(t, y):
    for iv in specie_da_controllare:
        if y[iv] < 0.001:
            event_counters[iv] += 1
            if event_counters[iv] >= 20000:
                return 0  # Ferma l'integrazione
        else:
            event_counters[iv] = 0
    return 1

event.terminal = True
event.direction = -1

# Integrazione temporale
t_ex = 3600 * 100000
t_eval = np.linspace(0, t_ex, 50000)
sol = scipy.integrate.solve_ivp(dy, [0., t_ex], model.state, t_eval=t_eval, max_step=5000, events=event)

# Estinzione?
if any(ev.size > 0 for ev in sol.t_events):
    flag_extin = True
else:
    # Estrai dati
    t = sol.t
    y = sol.y.T
    y = y[40000:]
    t_chaos = t[40000:]

    # Controllo negatività
    for iv in specie_da_controllare:
        if np.any(y[:, iv] < -0.01):
            flag_neg = True

    # Statistiche
    y_controllo = y[:, specie_da_controllare]
    threshold = 0.01
    mean_y = np.mean(y_controllo, axis=0)
    std_y = np.std(y_controllo, axis=0)

    # Chaos detection
#    taux = np.logspace(np.log(1), np.log10(len(t_chaos)//2), 20, dtype=int)
#    cmax = 0.5

#    for iv in specie_da_controllare:
#        T = y[:, iv]
#        H = np.zeros(len(taux))
#        C = np.zeros(len(taux))

#        for it, tx in enumerate(taux):
#            try:
#                H[it], C[it] = od.complexity_entropy(T, taux=tx, dx=6)
#            except:
#                H[it], C[it] = np.nan, np.nan

#        for it in range(len(taux)):
#            if not np.isnan(H[it]) and not np.isnan(C[it]):
#                if 0.45 <= H[it] <= 0.7 and abs(C[it] - cmax) < 0.1 * cmax:
#                    flag_chaos = True
#                    break
 #       if flag_chaos:
 #           break
    # Lyap exponent
    flag_chaos = False
    lyap_threshold = 0.0001  # Soglia minima per considerare la dinamica caotica

    for iv in specie_da_controllare:
        T = y[:, iv]
        try:
            max_lyap = nd.lyap_r(T, emb_dim=6)
            print(f"Lyapunov max per variabile {iv}: {max_lyap:.4f}")
            if max_lyap > lyap_threshold:
                flag_chaos = True
                break
        except Exception as e:
            print(f"Errore nel calcolo Lyapunov per variabile {iv}: {e}")
            print(flag_chaos)    
    # Steady-state detection
    ratio = np.where(mean_y >= threshold, std_y / mean_y, 0)
    flag_ss = np.max(ratio) <= 0.01

# Salvataggio su NetCDF
fileoutput = 'result.nc'
f = nc.Dataset(fileoutput, mode='w')
f.createDimension('flag', 1)

f.createVariable('flag', np.int8, ('flag',))[:] = 1
f.createVariable('flag_chaos', np.int8, ('flag',))[:] = int(flag_chaos)
f.createVariable('flag_ss', np.int8, ('flag',))[:] = int(flag_ss)
f.createVariable('flag_extin', np.int8, ('flag',))[:] = int(flag_extin)
f.createVariable('flag_neg', np.int8, ('flag',))[:] = int(flag_neg)
f.close()

# Output a schermo
if flag_ss:
    print("flag_ss = True")
elif flag_chaos:
    print("flag_chaos = True")
elif not flag_extin and not flag_neg:
    print("gold")
else:
    print("Condizione non classificata (estinzione o negatività)")

