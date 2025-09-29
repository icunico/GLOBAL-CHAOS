import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Parametri 
#P1, P2, P3
r = [2.5, 2.5, 2.5]
K = [1, 1, 1]
tau=0.01

#Parametri per C1
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
dC1 = 0.5
dC2 = 0.0

# Parametri per X1
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
dX1 = 0.1

# Parametri Y
aX1Y1 = 0.25
bX1Y1 = 0.50
dY1 = 1.0

# Parametri Z
aY1Z1 = 0.125
bY1Z1 = 0.25
dZ1 = 1.0

# Parametri microbial loop
alpha1=0.3
alpha2=0.4
alpha3=0.3

def ecological_network(t, vars):
    # Estrai tutte le variabili 
    P1, P2, P3, C1, C2, X, Y, Z, D, N = vars

    # EQUAZIONI PER IL PLANCTON (P1, P2, P3) con dipendenza da N
    # ---------------------------------------------------------------
    dP1dt = P1 * ((r[0] *N *alpha1* (1 - P1/K[0]))  # Crescita dipendente da N
        - (aP1C1 * C1) / (1 + bP1C1 * P1 + bP2C1 * P2 + bP3C1 * P3)
        - (aP1C2 * C2) / (1 + bP1C2 * P1 + bP2C2 * P2 + bP3C2 * P3)
        - (aP1X1 * X) / (1 + bP1X1 * P1 + bP2X1 * P2 + bP3X1 * P3 + bC1X1 * C1 + bC2X1 * C2)   
        )	
    # Equazione per P2 
    dP2dt = P2 * ((r[1]* N * alpha2* (1 - P2/K[1]))
        - (aP2C1 * C1) / (1 + bP1C1 * P1 + bP2C1 * P2 + bP3C1 * P3)
        - (aP2C2 * C2) / (1 + bP1C2 * P1 + bP2C2 * P2 + bP3C2 * P3)
        - (aP2X1 * X) / (1 + bP1X1 * P1 + bP2X1 * P2 + bP3X1 * P3 + bC1X1 * C1 + bC2X1 * C2)
        )
    # Equazione per P3
    dP3dt = P3 * ((r[2]* N * alpha3*(1 - P3/K[2]))
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
    # Termini di morte totale (contributo di tutte le variabili)
    death_terms = (
        dC1 * C1 + dC2 * C2 +
        dX1 * X + dY1 * Y + dZ1 * Z
    )

    # Equazione per i Detriti
    dDdt = death_terms - (1/tau) * D

    # Equazione per i Nutrienti (N)
    dNdt = (1/tau) * D - N * (alpha1 * P1 *r[0]*(1 - P1/K[0]) + alpha2 * P2* r[1]*(1 - P2/K[1])  + alpha3 * P3*r[2]*(1 - P3/K[2]))

    # ---------------------------------------------------------------

    return [dP1dt, dP2dt, dP3dt, dC1dt, dC2dt, dXdt, dYdt, dZdt, dDdt, dNdt]



# Condizioni iniziali
initial_conditions = [0.5, 0.0, 0.0, 0.5, 0.0, 0.0, 0, 0, 0.5, 10]
# Intervallo di tempo
time_span = (0, 12000)
time_eval = np.linspace(*time_span, 6000)

# Risoluzione numerica
dynamics = solve_ivp(ecological_network, time_span, initial_conditions, t_eval=time_eval)
for iv in range(dynamics.y.shape[0]):
    dynamics.y[iv][dynamics.y[iv] < 0] = 0  # Imposta a zero tutti i valori negativi per ogni specie
# Grafico dei risultati

plt.figure(figsize=(10, 6))
labels = ["P1", "P2", "P3", "C1", "C2", "X", "Y", "Z", "D", "N"]
for i in range(len(labels)):
    plt.plot(dynamics.t, dynamics.y[i], label=labels[i])
plt.xlabel("Tempo")
plt.ylabel("Popolazione")
plt.legend()
plt.title("Dinamica della rete ecologica")
plt.show()

total_mass = np.sum(dynamics.y, axis=0)  # Somma delle variabili
# Plot della somma totale
plt.plot(dynamics.t, total_mass, color="black", linewidth=2.5, linestyle="--")
plt.xlabel("Tempo")
plt.ylabel("Massa tot")
plt.ylim(0, 10)
plt.title("Bilancio di Massa")
plt.grid(True)
plt.show()

