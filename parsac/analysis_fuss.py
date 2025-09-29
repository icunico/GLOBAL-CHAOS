import pickle 
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import itertools

pkname = 'fussmann3G_legP1.pickle'
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
#combinations = x[:, 1:3] 
#combinations = x[:, 0] 
#combinations = x[:, [0, 1]]
combinations = x[-puntiplot:, [0, 1]]
#combinations = x[:, 0:3]

print(combinations)
#print(combinations)

#fig, ax = plt.subplots(figsize=(8, 6))  # Puoi anche modificare la dimensione della figura
fig, ax = plt.subplots(figsize=(6, 8))  # raddoppiato ma stesso rapporto
#fig, ax = plt.subplots(figsize=(6, 8))  # meno stretta, altezza simile

    # Prendi il valore per ogni flag associato alla combinazione corrente
for i, (dx, dc) in enumerate(combinations):
    current_flag_chaos = flag_chaos[i]
    current_flag_neg = flag_neg[i]
    current_flag_extin = flag_extin[i]
    current_flag_ss = flag_ss[i]

    # Controllo per ciascun flag separatamente
    if current_flag_extin == 1:  # Se flag_extin è 1
        ax.scatter(dx, dc, color='black', alpha=0, s=20)
        #print(current_flag_extin)
    elif current_flag_neg == 1:  # Se flag_neg è 1
        ax.scatter(dx, dc, color='black', alpha=0, s=20)
    elif current_flag_ss == 1:   # Se flag_ss è 1
        ax.scatter(dx, dc, color='green', s=20) 
    elif current_flag_chaos == 1:   # Se flag_chaos è 1
        ax.scatter(dx, dc, color=(1.0, 100/255, 0/255),  s=20)
        #print(current_flag_chaos)
    else:
        ax.scatter(dx, dc, color='gold', s=20)

    # Stampa dei valori finali associati a questa combinazione
#    print(f"dx = {dx}, dc = {dc} -> chaos: {current_flag_chaos}, neg: {current_flag_neg}, extin: {current_flag_extin}, ss: {current_flag_ss}")

# Etichetta e titolo
#ax.set_xlabel('dC1', fontsize=18)
#ax.set_ylabel('dX1', fontsize=18)
ax.set_xlim(0.05, 2.0)
ax.set_ylim(0.05, 2.0)
ax.tick_params(axis='both', which='major', labelsize=18)
#ax.set_title('Flag-based Scatter Plot')
#plt.legend()
    
# Mostra il grafico
plt.show()
plt.savefig("Result3B_paper.png")  
# Maschera combinata: flag_ss è 1 e flag_extin è 0
#mask = (flag_ss == 0) & (flag_chaos == 0) &  (flag_extin == 0)  

# Applichiamo la maschera a x e y
#x_selected = x[mask]
#y_selected = y[mask]

#sorted_x = np.sort(x_selected[:, 0])
#print(f"Numero di combinazioni con flag_chaos == 0 e flag_extin == 0: {np.sum(mask)}")
#print("Prime 10 combinazioni di parametri ordinate in ordine crescente:")
#print(sorted_x)


#CALCOLO PERCENTUALI
import numpy as np

dx = combinations[:, 0]
dc = combinations[:, 1]

colori = []
for i in range(len(dc)):
    if flag_extin[i] == 1 or flag_neg[i] == 1:
        colori.append('black')
    elif flag_ss[i] == 1:
        colori.append('green')
    elif flag_chaos[i] == 1:
        colori.append('red')
    else:
        colori.append('gold')
colori = np.array(colori)

import numpy as np

colori_da_usare = ['red', 'green', 'gold']
bins = 80
x_edges = np.linspace(0.01, 2.0, bins + 1)
y_edges = np.linspace(0.01, 2.0, bins + 1)

# Calcola histogram2d per ogni colore
conteggi = {}
for colore in colori_da_usare:
    mask = (colori == colore)
    H, _, _ = np.histogram2d(dc[mask], dx[mask], bins=[x_edges, y_edges])
    conteggi[colore] = H

# Stack per decidere per ogni bin quale colore è predominante
stack = np.stack([conteggi[colore] for colore in colori_da_usare], axis=-1)

# Crea una maschera per selezionare SOLO i bin OCCUPATI (almeno un punto in uno dei colori)
mask_bin_occupati = np.sum(stack, axis=-1) > 0

# Per quei bin occupati, trova quale colore è predominante
colore_predominante = np.nanargmax(stack, axis=-1)

# Conta quanti bin predominanti per ciascun colore SOLO tra quelli occupati
conteggi_bin = {
    colore: np.sum((colore_predominante == i) & mask_bin_occupati)
    for i, colore in enumerate(colori_da_usare)
}

totale_bin = np.sum(list(conteggi_bin.values()))
print(totale_bin)
percentuali = {colore: 100 * n / totale_bin for colore, n in conteggi_bin.items()}

print("\n=== Percentuali area (senza sovrapposizioni, solo bin occupati, somma 100%) ===")
for colore in colori_da_usare:
    print(f"{colore}: {percentuali[colore]:.2f}%")


# Calcolo area singolo bin
area_bin = (x_edges[1] - x_edges[0]) * (y_edges[1] - y_edges[0])
area_1D = (y_edges[1] - y_edges[0])
print(area_bin)
# Calcola area totale occupata
area_totale = totale_bin * area_bin
area_1D=totale_bin * area_1D
print(f"\nArea totale occupata dai colori red, green e gold: {area_totale:.4f}")
print(area_1D)

