import numpy as np
import matplotlib.pyplot as plt

# ===============================
# Parametri generali
# ===============================
num_steps = 100000
epsilon = 0.05
t_steps = 500
DELTA = 0.5
delta0 = 0.05
dt = 0.01
tol=0.2

sigma = 10.0
beta = 8.0 / 3.0

rho_values = np.linspace(34.2, 34.2, 1)

# ===============================
# Sistema di Lorenz
# ===============================
def f_map(x, rho):
    dx = sigma * (x[1] - x[0])
    dy = x[0] * (rho - x[2]) - x[1]
    dz = x[0] * x[1] - beta * x[2]
    return np.array([dx, dy, dz])

def jacobian(x, rho):
    return np.array([
        [-sigma, sigma, 0],
        [rho - x[2], -1, -x[0]],
        [x[1], x[0], -beta]
    ])

# ===============================
# Selezione ricorrenze (max 10)
# ===============================
def select_y_by_recurrence(x_vals, epsilon, t_steps, max_recurrences=500):
    recurrences_found = 0

    for t1 in range(len(x_vals) - t_steps):
        for t2 in range(t1 + t_steps, len(x_vals)):

            if np.linalg.norm(x_vals[t2] - x_vals[t1]) < epsilon:

                max_len = len(x_vals) - max(t1, t2)
                recurrences_found += 1

                if recurrences_found >= max_recurrences:
                    return (
                        x_vals[t1:t1 + max_len],
                        x_vals[t2:t2 + max_len]
                    )

                break

    return [], []

# ===============================
# Lyap_tau
# ===============================
def lyap_tau_continuous(x_vals, y_vals, DELTA, delta0):
    tau_list = []
    d0_list = []

    counting = False
    tau_current = 0
    d0_current = 0

    for i in range(len(x_vals)):
        dist = np.linalg.norm(x_vals[i] - y_vals[i])

        if not counting:
            if  delta0 * (1 - tol) <= dist <= delta0 :
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
        lyap_val = (1.0 / (mean_tau * dt)) * np.log(DELTA / mean_d0)
    else:
        lyap_val = np.nan

    return lyap_val

# ===============================
# Loop su rho
# ===============================
lyap_trad_list = []
lyap_tau_list = []

for rho in rho_values:

    # Genera traiettoria
    x_series = [np.array([1.0, 1.0, 1.0])]

    for _ in range(num_steps):
        x_new = x_series[-1] + dt * f_map(x_series[-1], rho)
        x_series.append(x_new)

    x_series = np.array(x_series)

    # -----------------
    # Lyapunov tradizionale (massimo)
    # -----------------
    v = np.array([1.0, 0.0, 0.0])
    lyap_sum = 0

    for x in x_series:
        J = jacobian(x, rho)
        v = v + dt * J @ v
        norm_v = np.linalg.norm(v)
        v = v / norm_v
        lyap_sum += np.log(norm_v)

    lyap_trad = lyap_sum / (num_steps * dt)
    lyap_trad_list.append(lyap_trad)

    # -----------------
    # Lyap_tau
    # -----------------
    x_aligned, y_current = select_y_by_recurrence(x_series, epsilon, t_steps)

    if len(x_aligned) > 0:
        lyap_tau_val = lyap_tau_continuous(
            x_aligned, y_current, DELTA, delta0
        )
    else:
        lyap_tau_val = np.nan

    lyap_tau_list.append(lyap_tau_val)

    print(f"rho={rho:.2f} | Lyap trad={lyap_trad:.3f} | Lyap tau={lyap_tau_val:.3f}")

# ===============================
# Plot finale
# ===============================
plt.figure(figsize=(8,6))
plt.scatter(rho_values, lyap_trad_list, color='black',  label='Lyap tradizionale')
plt.scatter(rho_values, lyap_tau_list, color= 'orange', label='Lyap tau')
plt.axhline(0, color='k', linestyle=':')
plt.xlabel("rho")
plt.ylim(-5, 5)
plt.ylabel("Lyapunov exponent")
plt.legend()
plt.grid(True)
plt.title("Lyapunov vs rho (Lorenz system)")
plt.savefig("lyaptau_LORENZ.png", dpi=300, bbox_inches='tight')
plt.show()
