import pickle
import numpy as np
import matplotlib.pyplot as plt

# -------------------------------
# Caricamento dati
# -------------------------------
pkname = 'fussmann.pickle'
with open(pkname, 'rb') as infile:
    new_dict = pickle.load(infile)

# Stampa informazioni sulle chiavi
for key, value in new_dict.items():
    print(f"Chiave: {key}, Tipo: {type(value)}")
    if isinstance(value, (list, np.ndarray)):
        print(f"Primi 10 valori: {value[:10]}")
    print("-" * 40)

# -------------------------------
# Estrazione flag e combinazioni
# -------------------------------
y = new_dict['Y']
x = new_dict['X']

flag_chaos = y[:, 0]
flag_neg = y[:, 1]
flag_ss = y[:, 2]
flag_extin = y[:, 3]

combinations = x[:, [0, 1, 2]]

# -------------------------------
# Taglio degli ultimi punti
# -------------------------------
puntiplot = 100000
N = min(puntiplot, len(flag_chaos), len(combinations))
flag_chaos = flag_chaos[-N:]
flag_neg = flag_neg[-N:]
flag_ss = flag_ss[-N:]
flag_extin = flag_extin[-N:]
combinations = combinations[-N:]

print("Shape combinazioni:", combinations.shape)

# -------------------------------
# Plot 3D dei flag
# -------------------------------
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

# Colori e label
color_map = {
    'extin_neg': 'none',
    'ss': 'green',
    'chaos': (1.0, 100/255, 0/255),
    'other': 'gold'
}

for i, (dc, dx, dy) in enumerate(combinations):
    current_flag_chaos = flag_chaos[i]
    current_flag_neg = flag_neg[i]
    current_flag_extin = flag_extin[i]
    current_flag_ss = flag_ss[i]

    # Controllo flag
    if current_flag_extin == 1 or current_flag_neg == 1:
        color = color_map['extin_neg']
    elif current_flag_ss == 1:
        color = color_map['ss']
    elif current_flag_chaos == 1:
        color = color_map['chaos']
    else:
        color = color_map['other']

    ax.scatter(dc, dx, dy, color=color)

    # Stampa solo se entrambi i flag negativi ed estinzione sono diversi da 0
    if current_flag_extin != 0 and current_flag_neg != 0:
        print(f"dx={dx}, dc={dc} -> chaos={current_flag_chaos}, neg={current_flag_neg}, extin={current_flag_extin}, ss={current_flag_ss}")

# Etichette assi
ax.set_xlabel('dC1')
ax.set_ylabel('dC2')
ax.set_zlabel('dX1')
ax.set_xlim(0.0, 2.0)
ax.set_ylim(0.0, 2.0)
ax.set_zlim(0.0, 2.0)

plt.show()

# -------------------------------
# Histogram 3D per colori
# -------------------------------
dx_vals = combinations[:, 0]
dy_vals = combinations[:, 1]
dz_vals = combinations[:, 2]

colori = np.array([
    'none' if flag_extin[i] == 1 or flag_neg[i] == 1 else
    'green' if flag_ss[i] == 1 else
    'red' if flag_chaos[i] == 1 else
    'gold'
    for i in range(N)
])

colori_da_usare = ['red', 'green', 'gold']
bins = 200
x_edges = np.linspace(0.0, 2.0, bins + 1)
y_edges = np.linspace(0.0, 2.0, bins + 1)
z_edges = np.linspace(0.0, 2.0, bins + 1)
edges = [x_edges, y_edges, z_edges]

conteggi = {}
for colore in colori_da_usare:
    mask = (colori == colore)
    coords = np.vstack((dx_vals[mask], dy_vals[mask], dz_vals[mask])).T
    if coords.shape[0] > 0:
        H, _ = np.histogramdd(coords, bins=edges)
    else:
        H = np.zeros((bins, bins, bins))
    conteggi[colore] = H

# Stack dei conteggi
stack = np.stack([conteggi[colore] for colore in colori_da_usare], axis=0)

# Maschera dei voxel occupati
mask_voxel_occupati = np.sum(stack, axis=0) > 0

# Colore predominante per voxel occupato
colore_predominante = np.nanargmax(stack, axis=0)

# Conteggio voxel per colore
conteggi_voxel = {
    colore: np.sum((colore_predominante == i) & mask_voxel_occupati)
    for i, colore in enumerate(colori_da_usare)
}

totale_voxel = np.sum(list(conteggi_voxel.values()))
percentuali = {colore: 100 * n / totale_voxel for colore, n in conteggi_voxel.items()}

print(f"\nTotale voxel occupati: {totale_voxel}")
print("\n=== Percentuali volume (senza sovrapposizioni, somma 100%) ===")
for colore in colori_da_usare:
    print(f"{colore}: {percentuali[colore]:.2f}%")

# Volume voxel
volume_voxel = (x_edges[1] - x_edges[0]) * (y_edges[1] - y_edges[0]) * (z_edges[1] - z_edges[0])
volume_totale = totale_voxel * volume_voxel
print(f"\nVolume singolo voxel: {volume_voxel:.6f}")
print(f"Volume totale occupato dai colori red, green e gold: {volume_totale:.6f}")

