import numpy as np
import matplotlib.pyplot as plt

# ===============================
# Parametri
# ===============================
num_steps = 100000
x0 = 0.2
epsilon = 0.0001
t_steps = 300
DELTA = 1e-2
delta0 = 1e-4

# ===============================
# Logistic map
# ===============================
def f_map(x, r, sigma=1e-7):
    return r * x * (1 - x) + sigma * np.random.uniform(-0.5, 0.5)

# ===============================
# Selezione y tramite ricorrenza spaziale
# ===============================
def select_y_by_recurrence(x_vals, epsilon, t_steps):
    for t1 in range(len(x_vals)-t_steps):
        for t2 in range(t1+t_steps, len(x_vals)):
            if abs(x_vals[t2]-x_vals[t1]) < epsilon:
                max_len = len(x_vals)-max(t1,t2)
                return x_vals[t1:t1+max_len], x_vals[t2:t2+max_len]
    return [], []

# ===============================
# Lyap_tau (versione continua)
# ===============================
def lyap_tau_continuous(x_vals, y_vals, DELTA, delta0):
    tau_list = []  # Lista dei tau misurati
    d0_list = []   # Lista delle d0 (distanze iniziali)

    counting = False
    tau_current = 0
    d0_current = 0

    # Loop per calcolare i riavvicinamenti
    for i in range(len(x_vals)):
        dist = abs(x_vals[i] - y_vals[i])

        if not counting:
            if dist < delta0:  # Cerchiamo un riavvicinamento
                counting = True
                tau_current = 0
                d0_current = dist
        else:
            tau_current += 1

            if dist >= DELTA:
                tau_list.append(tau_current)
                d0_list.append(d0_current)

                counting = False
                tau_current = 0

    if tau_list:
        mean_tau = np.mean(tau_list)
        mean_d0 = np.mean(d0_list)

        if mean_tau > 0 and mean_d0 > 0:
            lyap_val = (1.0 / mean_tau) * np.log(DELTA / mean_d0)
        else:
            lyap_val = np.nan
    else:
        mean_tau = np.nan
        lyap_val = np.nan

    return tau_list, mean_tau, lyap_val, d0_list  # Aggiungi d0_list qui

# ===============================
# Loop su valori di r
# ===============================
r_values = np.linspace(3.5, 4.0, 100)
lyap_trad_list = []
lyap_tau_list = []
mean_tau_list = []

# Apri file per scrivere SOLO i risultati
with open("lyapunov_results.txt", "w") as f:
    f.write("r_value | tau_medio | Lyap_trad | Lyap_tau\n")  # Apri file per scrivere SOLO i risultati

    for r in r_values:
        # Genera la traiettoria
        x_series = [x0]
        for _ in range(num_steps):
            x_series.append(f_map(x_series[-1], r))

        # Calcolo Lyap tradizionale
        lyap_sum = 0
        count = 0
        for x in x_series[:-1]:
            if 0 < x < 1:
                lyap_sum += np.log(abs(r * (1 - 2 * x)))  
                count += 1
        lyap_trad = lyap_sum / count if count > 0 else np.nan

        # Trova y per ricorrenza
        x_aligned, y_current = select_y_by_recurrence(x_series, epsilon, t_steps)

        if len(x_aligned) > 0:
            # Stampa la coppia scelta
            print(f"Coppia trovata: x = {x_aligned[0]}, y = {y_current[0]}")

            # Calcola Lyap tau
            tau_values, mean_tau, lyap_tau_val, d0_list = lyap_tau_continuous(x_aligned, y_current, DELTA, delta0)

            # Calcola la distanza media
            mean_d0 = np.mean(d0_list)
            print(f"Distanza media tra le ricorrenze = {mean_d0:.6f}")

            # Stampa le ricorrenze con i rispettivi tau
            for i, tau in enumerate(tau_values):
                print(f"  Ricorrenza {i+1}: Tau = {tau}, Distanza iniziale = {d0_list[i]}")
                # Stampa i punti x e y di ogni ricorrenza
                print(f"    Punti: x = {x_aligned[i]}, y = {y_current[i]}")

            # Stampa a schermo del tau medio e Lyap
            print(f"r={r:.3f} | tau medio={mean_tau:.3f} | Lyap trad={lyap_trad:.3f} | Lyap tau={lyap_tau_val:.3f} | Distanza media={mean_d0:.6f}")

            # Memorizza i valori nel list
            lyap_trad_list.append(lyap_trad)
            lyap_tau_list.append(lyap_tau_val)
            mean_tau_list.append(mean_tau)

            # Scrivi su file (SOLO i numeri) - nan diventa 0.000
            tau_for_file = 0.0 if np.isnan(mean_tau) else mean_tau
            lyap_for_file = 0.0 if np.isnan(lyap_tau_val) else lyap_tau_val

            f.write(f"{r:.3f} | {tau_for_file:.3f} | {lyap_trad:.3f} | {lyap_for_file:.3f}\n")
        else:
            print(f"r={r:.3f} | Nessuna ricorrenza trovata")
            lyap_trad_list.append(lyap_trad)
            lyap_tau_list.append(np.nan)            
            mean_tau_list.append(np.nan)

            # Scrivi su file (con 0.000 per valori mancanti)
            f.write(f"{r:.3f} | 0.000 | {lyap_trad:.3f} | 0.000\n")

# ===============================
# Plot risultati
# ===============================
plt.figure(figsize=(12, 8))

# Plot 1: Esponenti di Lyapunov
plt.subplot(2, 1, 1)
plt.plot(r_values, lyap_trad_list, 'b-', label='Lyap tradizionale', linewidth=2, marker='o', markersize=4)
plt.plot(r_values, lyap_tau_list, 'r--', label='Lyap tau', linewidth=2, marker='s', markersize=4)
plt.axhline(y=0, color='k', linestyle=':', alpha=0.5)
plt.ylabel('Esponente di Lyapunov')
plt.title('Confronto metodi Lyapunov')
plt.legend()
plt.grid(True, alpha=0.3)

# Plot 2: Tau medi
plt.subplot(2, 1, 2)
plt.plot(r_values, mean_tau_list, 'g-', label='Tau medio', linewidth=2, marker='^', markersize=4)
plt.xlabel('r')
plt.ylabel('Tau medio (passi)')
plt.title('Tau medio vs r')
plt.legend()
plt.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('lyapunov_comparison.png', dpi=150)
plt.show()

