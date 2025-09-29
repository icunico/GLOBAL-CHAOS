import pickle
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import itertools
from mpl_toolkits.mplot3d import Axes3D


pkname = 'fussmann.pickle'
infile = open(pkname,'rb')
new_dict = pickle.load(infile)
infile.close()

for key, value in new_dict.items():
    print(f"Chiave: {key}")
    print(f"Tipo: {type(value)}")

    # Se è un array numpy o una lista, stampiamo i primi elementi
    if isinstance(value, (list, np.ndarray)):
        print(f"Valori (primi 10): {value[:10]}")
    else:
        print(f"Valore: {value}")

    print("-" * 40)

y = new_dict['Y']  #matrice delle flag  
x = new_dict['X']  #parametro di sensibilta'

flag_chaos = y[:, 0]
flag_neg = y[:, 1]
flag_ss = y[:, 2]
flag_extin = y[:, 3]

#print("flag_chaos =", flag_chaos)
#print("flag_neg =", flag_neg)
#print("flag_ss =", flag_ss)
#print("flag_extin =", flag_extin)

# Creazione delle combinazioni con itertools.product
#combinations = x[:, 0] 
combinations = x[:, [0, 1, 2]]
print(combinations)

fig, ax = plt.subplots(figsize=(8, 6))  # Puoi anche modificare la dimensione della figura
ax = fig.add_subplot(111, projection='3d')

for i, (dc, dx, dy) in enumerate(combinations):
    # Prendi il valore per ogni flag associato alla combinazione corrente
    current_flag_chaos = flag_chaos[i]
    current_flag_neg = flag_neg[i]
    current_flag_extin = flag_extin[i]
    current_flag_ss = flag_ss[i]

    # Controllo per ciascun flag separatamente
    if current_flag_extin == 1:  # Se flag_extin è 1
        ax.scatter(dc, dx, dy, color='black', alpha=0)
        #print(current_flag_extin)
    elif current_flag_neg == 1:  # Se flag_neg è 1
        ax.scatter(dc, dx, dy, color='black', alpha=0)
    elif current_flag_ss == 1:   # Se flag_ss è 1
        ax.scatter(dc, dx, dy, color='green')
    elif current_flag_chaos == 1:   # Se flag_chaos è 1
        ax.scatter(dc, dx, dy, color='red')
        #print(current_flag_chaos)
    else:
        ax.scatter(dc, dx, dy, color='gold')

    # Stampa solo se flag_extin e flag_neg sono diversi da 0
    if current_flag_extin != 0 and current_flag_neg != 0:
        print(f"dx = {dx}, dc = {dc} -> chaos: {current_flag_chaos}, neg: {current_flag_neg}, extin: {current_flag_extin}, ss: {current_flag_ss}")

# Etichetta e titolo
ax.set_xlabel('dC1')
ax.set_ylabel('dC2')
ax.set_zlabel('dX1')
ax.set_xlim(0.05, 2.0)
ax.set_ylim(0.05, 2.0)
ax.set_zlim(0.05, 2.0)
ax.set_title('Flag-based Scatter Plot')
plt.legend()
    
# Mostra il grafico
plt.show()
