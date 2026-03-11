import numpy as np
import matplotlib.pyplot as plt

# carica i dati ignorando la riga di intestazione
data = np.genfromtxt("lyapunov_vs_rho_direct_no_noise.dat", comments="#")

# colonne
rho = data[:,0]
lyap_trad = data[:,1]
lyap_tau = data[:,2]

# plot
plt.figure()

plt.scatter(rho, lyap_trad, label="lyap_trad")
plt.scatter(rho, lyap_tau, label="lyap_tau")

plt.xlabel("rho")
plt.ylabel("Lyapunov exponent")
plt.title("Lyapunov Exponents vs rho")

plt.legend()
plt.grid(True)

plt.show()
