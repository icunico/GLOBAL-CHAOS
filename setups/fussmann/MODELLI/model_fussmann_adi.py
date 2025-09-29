import ordpy as od
import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Parametri estratti dal documento
#P1, P2, P3
#r = [2.5, 2.5, 2.5]
#K = [1, 1, 1]

#C1
aP1C1 = 3.0
aP1C2 = 1.0
aP2C1 = 2.0
aP2C2 = 2.0
aP3C1 = 1.0
aP3C2 = 3.0
bP1C1 = 5.0
bP1C2 = 5.0
bP2C1 = 5.0
bP2C2 = 5.0
bP3C1 = 5.0
bP3C2 = 5.0
dC1 = 0.0022
dC2 = 0.3

# X1
aP1X1 = 0.0
aP2X1 = 0.0
aP3X1 = 0.0
aC1X1 = 0.4
aC2X1 = 0.0
bP1X1 = 0.0
bP2X1 = 0.0
bP3X1 = 0.0
bC1X1 = 2.0
bC2X1 = 0.0
dX1 = 0.2

# Y
aX1Y1 = 0.25
bX1Y1 = 0.50
dY1 = 1.0

# Z
aY1Z1 = 0.125
bY1Z1 = 0.25
dZ1 = 1.0


def ecological_network(t, vars):
    # Estrai tutte le variabili
    P1, P2, P3, C1, C2, X, Y, Z = vars

    # EQUAZIONI PER IL PLANCTON (P1, P2, P3) con dipendenza da N
    # ---------------------------------------------------------------
    dP1dt = P1 * ((1 - P1)  # Crescita dipendente da N
        - (aP1C1 * C1) / (1 + bP1C1 * P1 + bP2C1 * P2 + bP3C1 * P3)
        - (aP1C2 * C2) / (1 + bP1C2 * P1 + bP2C2 * P2 + bP3C2 * P3)
        - (aP1X1 * X) / (1 + bP1X1 * P1 + bP2X1 * P2 + bP3X1 * P3 + bC1X1 * C1 + bC2X1 * C2)
        )  
    # Equazione per P2
    dP2dt = P2 * ((1 - P2)
        - (aP2C1 * C1) / (1 + bP1C1 * P1 + bP2C1 * P2 + bP3C1 * P3)
        - (aP2C2 * C2) / (1 + bP1C2 * P1 + bP2C2 * P2 + bP3C2 * P3)
        - (aP2X1 * X) / (1 + bP1X1 * P1 + bP2X1 * P2 + bP3X1 * P3 + bC1X1 * C1 + bC2X1 * C2)
        )   
    # Equazione per P3
    dP3dt = P3 * ((1 - P3)
        - (aP3C1 * C1) / (1 + bP1C1 * P1 + bP2C1 * P2 + bP3C1 * P3)
        - (aP3C2 * C2) / (1 + bP1C2 * P1 + bP2C2 * P2 + bP3C2 * P3)
        - (aP3X1 * X) / (1 + bP1X1 * P1 + bP2X1 * P2 + bP3X1 * P3 + bC1X1 * C1 + bC2X1 * C2)
        )

    # ---------------------------------------------------------------
    # EQUAZIONI PER I CONSUMATORI (C1, C2)
    # ---------------------------------------------------------------

    dC1dt = C1 * ((aP1C1 * P1 + aP2C1 * P2 + aP3C1 * P3) / (1 + bP1C1 * P1 + bP2C1 * P2 + bP3C1 * P3)
                  - (aC1X1 * X) / (1 + bP1X1 * P1 + bP2X1 * P2 + bP3X1 * P3 + bC1X1 * C1 + bC2X1 * C2)
                  - dC1)

    dC2dt = C2 * ((aP1C2 * P1 + aP2C2 * P2 + aP3C2 * P3) / (1 + bP1C2 * P1 + bP2C2 * P2 + bP3C2 * P3)
                  - (aC2X1 * X) / (1 + bP1X1 * P1 + bP2X1 * P2 + bP3X1 * P3 + bC1X1 * C1 + bC2X1 * C2)
                  - dC2)
    # ---------------------------------------------------------------
    # EQUAZIONI PER I PREDATORI (X, Y, Z)
    # ---------------------------------------------------------------
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
        
    # ---------------------------------------------------------------
    # EQUAZIONI PER DETRITI E NUTRIENTI (BILANCIO DI MASSA)
    # ---------------------------------------------------------------
    
    return [dP1dt, dP2dt, dP3dt, dC1dt, dC2dt, dXdt, dYdt, dZdt]


#ANALISI STABILITA'
# Valori di dc e dx da esplorare
dtau_values = np.linspace(0.02, 0.18, 20)
dc_values = np.linspace(0.02, 0.48, 20)


# Crea una figura per il grafico
fig, ax = plt.subplots(figsize=(10, 6))

# Ciclo sui valori di dc e dx
for dtau in dtau_values:
    for dc in dc_values:
       # Aggiorna i parametri del modello
        dX1 = dtau  # Assegna il valore di dx a dX1
        dC1 = dc  # Assegna il valore di dc a dC1 (se necessario)
        # Stampa i parametri in esecuzione
        print(f"Esecuzione con tau={dtau}, dC1={dc}")
       # Condizioni iniziali
      
        initial_conditions = [0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 0.0]

        # Intervallo di tempo
        time_span = (0, 10000)
        time_eval = np.linspace(*time_span, 10000)

        # Risoluzione numerica
        dynamics = solve_ivp(ecological_network, time_span, initial_conditions, t_eval=time_eval, max_step=0.5)


        # Estrai i risultati
        y = dynamics.y.T  # Trasponi per avere forma (n_time_steps, n_variables)

        specie_da_controllare= [1, 3, 5] #'P1', 'P2', 'P3', 'C1', 'C2', 'X1', 'Y1', $
        flag_neg = False
        for iv in specie_da_controllare:
                if np.any(y[:, iv] < -0.01):
                    flag_neg = True

        y=y[-1000:]

                # Seleziona solo le specie da controllare
        y_controllo = y[:, specie_da_controllare]
        threshold = 0.01
        # Calcola la media e la deviazione standard SOLO per le specie selezionate
        mean_y = np.mean(y_controllo, axis=0)
        std_y = np.std(y_controllo, axis=0)
        
      # Controllo del rapporto std/mean solo per le specie selezionate

        ratio = np.where(mean_y >= threshold, std_y / mean_y, 0)
      # Controllo della stabilità
        flag_ss = np.max(ratio) <= 0.01

        #check if a specie gets extinct
        flag_extin=False
        for  iv in specie_da_controllare:
                   if np.max(y[:, iv] < 0.001):
                           flag_extin=True

        #chaos
        taux = np.logspace(np.log(1), np.log10(len(time_eval)//2), 200, dtype=int)
        #calcolo indici
 
        for iv in specie_da_controllare:
            T = y[:, iv]  # Usa i dati della specie da controllare
            H = np.zeros(len(taux))  # Inizializza array per H
            C = np.zeros(len(taux))  # Inizializza array per C
            # Calcola gli indici per ogni valore di taux
            flag_chaos=False

            for it, tx in enumerate(taux):
                try:
                    H[it], C[it] = od.complexity_entropy(T, taux=tx, dx=6)
        #            print(f"Specie {iv} - taux={tx}: H={H[it]}, C={C[it]}")
                except:
                    H[it], C[it] = np.nan, np.nan  # Gestisci eventuali errori
            cmax = 0.5 
#            print(cmax_interp)
            for it, tx in enumerate(taux):
                if not np.isnan(H[it]) and not np.isnan(C[it]):
                    if 0.45 <= H[it] <= 0.7 and np.abs(C[it] - cmax) < 0.1 * cmax:  # C vicino al massimo
                        flag_chaos = True
                        break  # Se trovi almeno una condizione che soddisfa il chaos, esci dal ciclo
                    if flag_chaos:
                        break  # Esci dal ciclo delle specie se la flag è stata settata
            #print(flag_chaos)

        if flag_extin:
                 ax.scatter(dc, dtau, c='k')
        elif flag_neg:
                 ax.scatter(dc, dtau, c='k')
        elif flag_ss:
                 ax.scatter(dc, dtau, c='g')
        elif flag_chaos:
                 ax.scatter(dc, dtau, c='r')
        else:
                 ax.scatter(dc, dtau, c="gold")
fig.savefig('bubble_ml_base.png')
plt.subplots(figsize=(15, 10)) 
