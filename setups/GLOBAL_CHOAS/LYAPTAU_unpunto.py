import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import os
import shutil

# ===============================
# Time-delay embedding
# ===============================
def time_delay_embedding(x, m=3, tau=10):
    N = len(x)
    if N < (m - 1) * tau:
        return np.array([])
    return np.array([x[i : i + m * tau : tau] for i in range(N - (m - 1) * tau)])

# ===============================
# Selezione ricorrenze
# ===============================
def select_y_by_recurrence(x_vals, epsilon, t_steps, max_recurrences=50):
    recurrences_found = 0
    for t1 in range(len(x_vals) - t_steps):
        for t2 in range(t1 + t_steps, len(x_vals)):
            if np.linalg.norm(x_vals[t2] - x_vals[t1]) < epsilon:
                max_len = len(x_vals) - max(t1, t2)
                recurrences_found += 1
                if recurrences_found >= max_recurrences:
                    return (x_vals[t1:t1 + max_len], x_vals[t2:t2 + max_len])
                break
    return [], []

# ===============================
# Lyap_tau
# ===============================
def lyap_tau_continuous(x_vals, y_vals, DELTA, delta0, dt=1.0, tol=0.2):
    tau_list, d0_list = [], []
    counting = False
    tau_current = 0
    d0_current = 0
    for i in range(len(x_vals)):
        dist = np.linalg.norm(x_vals[i] - y_vals[i])
        if not counting:
            if delta0 * (1 - tol) <= dist <= delta0 * (1 + tol):
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
        return (1.0 / (mean_tau * dt)) * np.log(DELTA / mean_d0)
    else:
        return np.nan

# ===============================
# Carica dataset (link originale)
# ===============================
local_file = "regional_dataset_CHL_10_20_-30_-15.nc"
remote_file = "/g100_scratch/userexternal/gocchipi/GLOBAL/tmp_daily/regional_dataset_CHL_10_20_-30_-15.nc"

if not os.path.exists(local_file):
    shutil.copy2(remote_file, os.getcwd())
    print("copiato il file in working directory")

data = xr.open_dataset(local_file, decode_times=False)

# convert time
units, reference_date = data.time.attrs['units'].split('since')
data['time'] = pd.date_range(start=reference_date, periods=data.sizes['time'], freq='d')

# clean CHL
data['CHL'] = data['CHL'].where((data['CHL'] > 0) & (data['CHL'] < 100))

# ===============================
# Test su UN SOLO punto
# ===============================
lat_test = 17.5
lon_test = -17.5

point = data.sel(latitude=lat_test, longitude=lon_test, method="nearest")
xx = point['CHL'].values
xx = xx[~np.isnan(xx)]
print("Serie length:", len(xx))

# ===============================
# Parametri Lyap_tau
# ===============================
epsilon = 0.02
t_steps = 90
DELTA = 0.2
delta0 = 0.02

# ===============================
# Embedding
# ===============================
embedded = time_delay_embedding(xx, m=3, tau=10)
print("Embedded shape:", embedded.shape)

# ===============================
# Recurrence + Lyap_tau
# ===============================
if len(embedded) > 0:
    x_aligned, y_current = select_y_by_recurrence(embedded, epsilon, t_steps)
    if len(x_aligned) > 0:
        lyap_tau = lyap_tau_continuous(x_aligned, y_current, DELTA, delta0, dt=1.0, tol=0.2)
    else:
        lyap_tau = np.nan
else:
    lyap_tau = np.nan

print("Lyap_tau =", lyap_tau)

# ===============================
# Plot serie
# ===============================
plt.figure(figsize=(8,4))
plt.plot(xx)
plt.title(f"CHL time series ({lat_test:.2f}, {lon_test:.2f})")
plt.grid(True)
plt.savefig("single_point_series.png", dpi=300)
plt.close()
