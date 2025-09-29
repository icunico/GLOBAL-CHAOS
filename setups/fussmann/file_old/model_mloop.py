import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Parametri estratti dal documento
#P1, P2, P3
r = [2.5, 2.5, 2.5]
K = [1, 1, 1]


#C1
aP1C1 = 7.5
aP1C2 = 2.5
aP2C1 = 5.0
aP2C2 = 5.0
aP3C1 = 2.5
aP3C2 = 7.5
bP1C1 = 5.0
bP1C2 = 5.0
bP2C1 = 5.0
bP2C2 = 5.0
bP3C1 = 5.0
bP3C2 = 5.0
dC1 = 0.0022
dC2 = 0.0

# X1
aP1X1 = 0.0
aP2X1 = 0.0
aP3X1 = 0.0
aC1X1 = 1.0
aC2X1 = 0.0
bP1X1 = 0.0
bP2X1 = 0.0
bP3X1 = 0.0
bC1X1 = 2.0
bC2X1 = 0.0
dX1 = 0.223

# Y
aX1Y1 = 0.25
bX1Y1 = 0.50
dY1 = 1.0

# Z
aY1Z1 = 0.125
bY1Z1 = 0.25
dZ1 = 1.0


def ecological_network(t, vars):
    P1, P2, P3, C1, C2, X, Y, Z = vars
    
    # Equazioni per P1, P2, P3
    dP1dt = P1 * (r[0] * (1 - P1/K[0]) 
                 - (aP1C1 * C1) / (1 + bP1C1 * P1 + bP2C1 * P2 + bP3C1 * P3) 
                 - (aP1C2 * C2) / (1 + bP1C2 * P1 + bP2C2 * P2 + bP3C2 * P3) 
                 - (aP1X1 * X) / (1 + bP1X1 * P1 + bP2X1 * P2 + bP3X1 * P3 + bC1X1 * C1 + bC2X1 * C2))
    
    dP2dt = P2 * (r[1] * (1 - P2/K[1]) 
                 - (aP2C1 * C1) / (1 + bP1C1 * P1 + bP2C1 * P2 + bP3C1 * P3) 
                 - (aP2C2 * C2) / (1 + bP1C2 * P1 + bP2C2 * P2 + bP3C2 * P3) 
                 - (aP2X1 * X) / (1 + bP1X1 * P1 + bP2X1 * P2 + bP3X1 * P3 + bC1X1 * C1 + bC2X1 * C2))
    
    dP3dt = P3 * (r[2] * (1 - P3/K[2]) 
                 - (aP3C1 * C1) / (1 + bP1C1 * P1 + bP2C1 * P2 + bP3C1 * P3) 
                 - (aP3C2 * C2) / (1 + bP1C2 * P1 + bP2C2 * P2 + bP3C2 * P3) 
                 - (aP3X1 * X) / (1 + bP1X1 * P1 + bP2X1 * P2 + bP3X1 * P3 + bC1X1 * C1 + bC2X1 * C2))
    
    # Equazioni per C1 e C2
    dC1dt = C1 * ((aP1C1 * P1 + aP2C1 * P2 + aP3C1 * P3) / (1 + bP1C1 * P1 + bP2C1 * P2 + bP3C1 * P3) 
                  - (aC1X1 * X) / (1 + bP1X1 * P1 + bP2X1 * P2 + bP3X1 * P3 + bC1X1 * C1 + bC2X1 * C2) 
                  - dC1)
    
    dC2dt = C2 * ((aP1C2 * P1 + aP2C2 * P2 + aP3C2 * P3) / (1 + bP1C2 * P1 + bP2C2 * P2 + bP3C2 * P3)
                  - (aC2X1 * X) / (1 + bP1X1 * P1 + bP2X1 * P2 + bP3X1 * P3 + bC1X1 * C1 + bC2X1 * C2)
                  - dC2) 
    
    # Equazione per X
    dXdt = X * ((aP1X1 * P1 + aP2X1 * P2 + aP3X1 * P3 + aC1X1 * C1 + aC2X1 * C2) / 
                (1 + bP1X1 * P1 + bP2X1 * P2 + bP3X1 * P3 + bC1X1 * C1 + bC2X1 * C2) 
                - (aX1Y1* Y) / (1 + bX1Y1 * X) 
                - dX1)
    
    # Equazione per Y
    dYdt = Y * ((aX1Y1 * X) / (1 + bX1Y1 * X) 
               - (aY1Z1 * Z) / (1 + bY1Z1 * Y) 
               - dY1)
    
    # Equazione per Z
    dZdt = Z * ((aY1Z1 * Y) / (1 + bY1Z1 * Y) 
               - dZ1)
    
    return [dP1dt, dP2dt, dP3dt, dC1dt, dC2dt, dXdt, dYdt, dZdt]


# Grafico dei risultati
#plt.figure(figsize=(10, 6))
#labels = ["P1", "P2", "P3", "C1", "C2", "X", "Y", "Z"]
#for i in range(len(labels)):
#    plt.plot(dynamics.t, dynamics.y[i], label=labels[i])
#plt.xlabel("Tempo")
#plt.ylabel("Popolazione")
#plt.legend()
#plt.title("Dinamica della rete ecologica")
#plt.show()

#ANALISI STABILITA'
# Valori di dc e dx da esplorare
dx_values = np.linspace(0.10, 0.10, 1)
dc_values = np.linspace(0.65, 0.65, 1)


# Crea una figura per il grafico
fig, ax = plt.subplots(figsize=(10, 6))

# Ciclo sui valori di dc e dx
for dx in dx_values:
    for dc in dc_values:	
       # Aggiorna i parametri del modello
        dX1 = dx  # Assegna il valore di dx a dX1
        dC1 = dc  # Assegna il valore di dc a dC1 (se necessario)
        # Stampa i parametri in esecuzione
        print(f"Esecuzione con dX1={dx}, dC1={dc}")
        # Condizioni iniziali
        initial_conditions = [0.5, 0.5, 0.5, 0.5, 0.0, 0.5, 0.0, 0.0]

        # Intervallo di tempo
        time_span = (0, 6000)
        time_eval = np.linspace(*time_span, 6000)

        # Risoluzione numerica
        dynamics = solve_ivp(ecological_network, time_span, initial_conditions, t_eval=time_eval)

        # Estrai i risultati
        y = dynamics.y.T  # Trasponi per avere forma (n_time_steps, n_variables)
        specie_da_controllare = [0, 1, 2, 3, 5]  # Indici delle specie da controllare
        flag_neg = False
        for iv in specie_da_controllare:
            if np.any(dynamics.y[iv] < 0):  # Se ci sono valori negativi in qualsiasi punto temporale
               flag_neg = True
#        print(f"Attenzione: la specie {iv} ha valori negativi durante la simulazione!")

        # Prendi solo gli ultimi 300 punti temporali per l'analisi
        y = y[-500:]

        # Controllo sulla stabilità dello stato stazionario
        threshold = 0.01
        mean_y = np.mean(y, axis=0)
        std_y = np.std(y, axis=0)
        ratio = np.where(mean_y >= threshold, std_y / mean_y, 0)
        flag_ss = np.max(ratio) <= 0.01

        # Controllo se una specie si estingue
        flag_extin = False
        for iv in specie_da_controllare:
            if y[-1, iv] < 0.001:
#            if np.all(y[:, iv] < 0.001):
                flag_extin = True

        # Aggiungi un punto al grafico in base ai risultati
        if flag_extin:
            ax.scatter(dc, dx, c='k', label='Extinct' if flag_extin else '')
        elif flag_ss:
            ax.scatter(dc, dx, c='r', label='Stable Steady State' if flag_ss else '')
        elif flag_neg:
            ax.scatter(dc, dx, c='k', label='value neg' if flag_neg else '')        
        else:
            ax.scatter(dc, dx, c='y', label='Unstable' if not flag_ss else '')

# Configurazione del grafico
ax.set_xlabel('dc')
ax.set_ylabel('dx')
ax.set_title('Stability and Extinction Analysis')
plt.show()

# Salva il grafico
fig.savefig('bubblef2.png')
