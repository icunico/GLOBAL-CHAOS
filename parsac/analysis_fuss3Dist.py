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

puntiplot=20000
flag_chaos = flag_chaos[-puntiplot:]
flag_neg = flag_neg[-puntiplot:]
flag_ss = flag_ss[-puntiplot:]
flag_extin = flag_extin[-puntiplot:]

#print("flag_chaos =", flag_chaos)
#print("flag_neg =", flag_neg)
#print("flag_ss =", flag_ss)
#print("flag_extin =", flag_extin)

# Creazione delle combinazioni con itertools.product
#combinations = x[:, 0] 
combinations = x[:, [0, 1, 2]]
print(combinations)


fig, ax = plt.subplots(figsize=(8, 6), subplot_kw={'projection': '3d'})
#fig, ax = plt.subplots(figsize=(8, 6))  # Puoi anche modificare la dimensione della figura
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
        ax.scatter(dc, dx, dy, color=(1.0, 100/255, 0/255))
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
ax.set_xlim(0.0, 2.0)
ax.set_ylim(0.0, 2.0)
ax.set_zlim(0.0, 2.0)
#ax.set_title('Flag-based Scatter Plot')
plt.legend()
    
# Mostra il grafico
plt.show()


import numpy as np

# Assume che combinations abbia 3 colonne: dx, dy, dz
dx = combinations[:, 0]
dy = combinations[:, 1]
dz = combinations[:, 2]

# Assegna i colori come prima
colori = []
for i in range(len(dx)):
    if flag_extin[i] == 1 or flag_neg[i] == 1:
        colori.append('black')
    elif flag_ss[i] == 1:
        colori.append('green')
    elif flag_chaos[i] == 1:
        colori.append('red')
    else:
        colori.append('gold')
colori = np.array(colori)

# Definizione dei colori da analizzare e bins
colori_da_usare = ['red', 'green', 'gold']
bins = 20  
x_edges = np.linspace(0.0, 0.6, bins + 1)
y_edges = np.linspace(0.0, 0.6, bins + 1)
z_edges = np.linspace(0.0, 0.6, bins + 1)

edges = [x_edges, y_edges, z_edges]

# Calcola histogramdd per ogni colore
conteggi = {}
for colore in colori_da_usare:
    mask = (colori == colore)
    coords = np.vstack((dx[mask], dy[mask], dz[mask])).T
    H, _ = np.histogramdd(coords, bins=edges)
    conteggi[colore] = H

# Stack dei conteggi con shape (n_colori, bins, bins, bins)
stack = np.stack([conteggi[colore] for colore in colori_da_usare], axis=0)

# Maschera dei voxel occupati (almeno un punto)
mask_voxel_occupati = np.sum(stack, axis=0) > 0

# Per ogni voxel occupato, trova il colore predominante
colore_predominante = np.nanargmax(stack, axis=0)

# Conta quanti voxel predominanti per ciascun colore SOLO tra quelli occupati
conteggi_voxel = {
    colore: np.sum((colore_predominante == i) & mask_voxel_occupati)
    for i, colore in enumerate(colori_da_usare)
}

totale_voxel = np.sum(list(conteggi_voxel.values()))
print(f"Totale voxel occupati: {totale_voxel}")

percentuali = {colore: 100 * n / totale_voxel for colore, n in conteggi_voxel.items()}

print("\n=== Percentuali volume (senza sovrapposizioni, solo voxel occupati, somma 100%) ===")
for colore in colori_da_usare:
    print(f"{colore}: {percentuali[colore]:.2f}%")

# Calcolo volume di un singolo voxel
volume_voxel = (x_edges[1] - x_edges[0]) * (y_edges[1] - y_edges[0]) * (z_edges[1] - z_edges[0])
print(f"Volume singolo voxel: {volume_voxel:.6f}")

# Calcola volume totale occupato
volume_totale = totale_voxel * volume_voxel
print(f"\nVolume totale occupato dai colori red, green e gold: {volume_totale:.6f}")

