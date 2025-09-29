import pickle
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt
import itertools

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
#combinations = x[:, 1:3] 
#combinations = x[:, 0] 
combinations = x[:, [1, 0]]
#combinations = x[:, 0:3]
print(combinations)
