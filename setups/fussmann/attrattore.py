import numpy as np
import scipy.integrate
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Per il plot 3D

# Caricamento modello
import pyfabm
model = pyfabm.Model('fabm_netcom_tot.yaml')
model.cell_thickness = 1.
assert model.checkReady(), 'Model not ready'

# Estrai le specie da controllare
specie_da_controllare_nomi = ['P1/DW1','P1/DW3', 'C1/DWC1', 'C1/DWC2', 'X/DWX1']
variable_names = [var.path for var in model.state_variables]
specie_da_controllare = [i for i, name in enumerate(variable_names) if name in specie_da_controllare_nomi]

# Funzione derivata temporale
def dy(t0, y):
    model.state[:] = y
    return model.getRates()

# Integrazione temporale
t_ex = 3600 * 300000
t_eval = np.linspace(0, t_ex, 50000)
sol = scipy.integrate.solve_ivp(dy, [0., t_ex], model.state, t_eval=t_eval, max_step=5000)

# Estrai i dati
t = sol.t
y = sol.y.T
y = y[40000:]  # Utilizza solo gli ultimi valori, ad esempio dalla 40.000-esima iterazione in avanti
t_chaos = t[40000:]

# Seleziona le specie da plottare (3 specie in questo caso)
specie_selezionate = [0, 1, 2]  # Modifica questi indici per selezionare le specie che ti interessano
y_controllo = y[:, specie_selezionate]

# Plot 3D dell'attrattore
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# Traccia le tre specie selezionate
ax.plot(y_controllo[:, 0], y_controllo[:, 1], y_controllo[:, 2], label='Attrattore 3D')

# Etichette degli assi
ax.set_xlabel('Specie 1')
ax.set_ylabel('Specie 2')
ax.set_zlabel('Specie 3')

# Titolo e legenda
ax.set_title('Attrattore 3D')
ax.legend()

plt.show()

