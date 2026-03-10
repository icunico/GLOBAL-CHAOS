import numpy as np
import matplotlib.pyplot as plt

# ===============================
# Parametri generali
# ===============================
num_steps = 100000
t_steps = 300
DELTA = 1.0
delta0 = 0.1
dt = 0.01
tol = 0.4  # tolleranza per il calcolo Lyap_tau

sigma = 10.0
beta = 8.0 / 3.0
rho_values = np.linspace(0, 40, 100)

# Flag per scegliere il metodo
USE_EMBEDDING = False  # Se True usa embedding, se False usa output diretto del modello

print("=== PARAMETRI SIMULAZIONE ===")
print(f"Numero di passi: {num_steps}")
print(f"Passi minimi tra ricorrenze: {t_steps}")
print(f"DELTA (soglia massima): {DELTA}")
print(f"delta0 (soglia ricorrenze): {delta0}")
print(f"Tolleranza: {tol}")
print(f"dt: {dt}")
print(f"Valori di rho: {rho_values}")
print(f"Metodo: {'Embedding 3D' if USE_EMBEDDING else 'Output diretto del modello'}")
print()

# ===============================
# Sistema di Lorenz originale (full 3D)
# ===============================
def lorenz_system(x, rho):
    """
    Input:
        x: array numpy 1D di dimensione 3, stato del sistema [X, Y, Z]
        rho: float, parametro del sistema di Lorenz
    Output:
        array numpy 1D di dimensione 3, derivate del sistema di Lorenz [dX/dt, dY/dt, dZ/dt]
    """
    dX = sigma * (x[1] - x[0])  # dX/dt
    dY = x[0] * (rho - x[2]) - x[1]  # dY/dt
    dZ = x[0] * x[1] - beta * x[2]  # dZ/dt
    return np.array([dX, dY, dZ])

def create_3d_embedding(trajectory_1d, tau, embedding_dim=3):
    """
    Input:
        trajectory_1d: array numpy 1D, serie temporale di una singola componente
        tau: int, delay time per l'embedding
        embedding_dim: int, dimensione dello spazio di embedding (default=3)
    Output:
        array numpy 2D, traiettoria embedded di dimensione (n_points x embedding_dim)
        dove n_points = len(trajectory_1d) - (embedding_dim-1)*tau
    """
    n_points = len(trajectory_1d) - (embedding_dim - 1) * tau
    
    if n_points <= 0:
        raise ValueError(f"Tau troppo grande: tau={tau}, lunghezza={len(trajectory_1d)}, dim={embedding_dim}")
    
    embedded_traj = np.zeros((n_points, embedding_dim))
    
    for i in range(n_points):
        for j in range(embedding_dim):
            embedded_traj[i, j] = trajectory_1d[i + j * tau]
    
    print(f"  Embedding creato: {len(trajectory_1d)} -> {embedded_traj.shape[0]} punti")
    print(f"  Parametri: tau={tau}, dim={embedding_dim}")
    
    return embedded_traj

# ===============================
# Sistema di Lorenz per embedding (manteniamo per compatibilità)
# ===============================
def f_map(x, rho):
    """
    Input:
        x: array numpy 1D di dimensione 3, stato embedded del sistema
        rho: float, parametro del sistema di Lorenz
    Output:
        array numpy 1D di dimensione 3, derivate fittizie per il sistema embedded
    """
    # Per embedding, usiamo una dinamica fittizia basata sul sistema di Lorenz
    dx = sigma * (x[1] - x[0])
    dy = x[0] * (rho - x[2]) - x[1]
    dz = x[0] * x[1] - beta * x[2]
    return np.array([dx, dy, dz])

def jacobian(x, rho):
    """
    Input:
        x: array numpy 1D di dimensione 3, stato del sistema
        rho: float, parametro del sistema di Lorenz
    Output:
        array numpy 2D di dimensione 3x3, matrice Jacobiana del sistema di Lorenz
    """
    return np.array([
        [-sigma, sigma, 0],
        [rho - x[2], -1, -x[0]],
        [x[1], x[0], -beta]
    ])

# ===============================
# Selezione ricorrenze (max 200)
# ===============================
def select_y_by_recurrence(x_vals, delta0, tol, t_steps, max_recurrences=200):
    """
    Input:
        x_vals: array numpy 2D, traiettoria del sistema (num_points x 3)
        delta0: float, separazione iniziale di riferimento per identificare ricorrenze
        tol: float, tolleranza per identificare l'intorno di delta0
        t_steps: int, numero minimo di passi temporali tra ricorrenze
        max_recurrences: int, numero massimo di ricorrenze da considerare (default=50)
    Output:
        tupla (recurrence_pairs): lista di tuple (x_aligned, y_current) per ogni ricorrenza trovata,
        o [] se nessuna ricorrenza trovata
    """
    print(f"  Cercando ricorrenze nella traiettoria di {len(x_vals)} punti...")
    recurrences_found = 0
    soglia_min = delta0 * (1 - tol)
    soglia_max = delta0 * (1 + tol)
    print(f"  Soglia ricorrenze: {soglia_min:.4f} <= dist <= {soglia_max:.4f}")
    
    recurrence_pairs = []  # Lista per contenere tutte le ricorrenze
    last_t1_used = -t_steps  # Inizializza per permettere la prima ricorrenza
    
    for t1 in range(len(x_vals) - t_steps):
        # Verifica che t1 sia almeno t_steps passi dopo l'ultima ricorrenza trovata
        if t1 - last_t1_used < t_steps:
            continue
            
        for t2 in range(t1 + t_steps, len(x_vals)):
            dist = np.linalg.norm(x_vals[t2] - x_vals[t1])
            if soglia_min <= dist <= soglia_max:
                max_len = len(x_vals) - max(t1, t2)
                recurrences_found += 1
                last_t1_used = t1  # Aggiorna l'ultimo t1 utilizzato
                
                # Aggiungi la coppia di ricorrenze alla lista
                x_aligned = x_vals[t1:t1 + max_len]
                y_current = x_vals[t2:t2 + max_len]
                recurrence_pairs.append((x_aligned, y_current))
                
                print(f"  Ricorrenza {recurrences_found} trovata: t1={t1}, t2={t2}, dist={dist:.4f}, lunghezza={max_len}")
                
                if recurrences_found >= max_recurrences:
                    print(f"  Raggiunto limite massimo di ricorrenze ({max_recurrences})")
                    break
                break
        
        # Uscita anche dal loop esterno se raggiunto il massimo
        if recurrences_found >= max_recurrences:
            break
    
    if recurrences_found == 0:
        print("  Nessuna ricorrenza trovata!")
        return []
    else:
        print(f"  Totale ricorrenze trovate: {len(recurrence_pairs)}")
        return recurrence_pairs

def lyap_tau_continuous_multi(recurrence_pairs, DELTA, delta0, tol=0.2):
    """
    Input:
        recurrence_pairs: lista di tuple (x_vals, y_vals), coppie di traiettorie ricorrenti
        DELTA: float, soglia massima di separazione
        delta0: float, separazione iniziale di riferimento
        tol: float, tolleranza per identificare l'intorno di delta0 (default=0.2)
    Output:
        float, esponente di Lyapunov calcolato con il metodo tau su tutte le ricorrenze,
        o np.nan se non calcolabile
    """
    if not recurrence_pairs:
        print("  Nessuna ricorrenza fornita per il calcolo!")
        return np.nan
    
    print(f"  Calcolo Lyapunov tau su {len(recurrence_pairs)} ricorrenze...")
    
    all_tau_list = []
    all_d0_list = []
    soglia_min = delta0 * (1 - tol)
    soglia_max = delta0 * (1 + tol)
    
    # Processa ogni coppia di ricorrenze
    for pair_idx, (x_vals, y_vals) in enumerate(recurrence_pairs):
        print(f"    Processando ricorrenza {pair_idx + 1}/{len(recurrence_pairs)} ({len(x_vals)} punti)...")
        
        tau_list, d0_list = [], []
        counting, tau_current, d0_current = False, 0, 0
        
        for i in range(len(x_vals)):
            dist = np.linalg.norm(x_vals[i] - y_vals[i])
            
            if not counting:
                if soglia_min <= dist <= soglia_max:
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
        
        # Aggiungi i tau e d0 di questa ricorrenza al totale
        all_tau_list.extend(tau_list)
        all_d0_list.extend(d0_list)
        print(f"      Trovati {len(tau_list)} intervalli tau per questa ricorrenza")
    
    # Calcola il Lyapunov finale usando tutti i tau raccolti
    if all_tau_list:
        mean_tau = np.mean(all_tau_list)
        mean_d0 = np.mean(all_d0_list)
        lyap_val = (1.0 / (mean_tau * dt)) * np.log(DELTA / mean_d0)
        print(f"  Totale intervalli tau trovati: {len(all_tau_list)}")
        print(f"  Tau medio finale: {mean_tau:.2f}, d0 medio finale: {mean_d0:.4f}")
    else:
        lyap_val = np.nan
        print("  Nessun intervallo tau trovato in tutte le ricorrenze!")
    
    return lyap_val

# ===============================
# Loop su rho (con flag per scegliere il metodo)
# ===============================
print("=== INIZIO SIMULAZIONI ===")
lyap_trad_list = []
lyap_tau_list = []

# Parametri embedding (usati solo se USE_EMBEDDING = True)
embedding_tau = 10  # delay time per l'embedding
embedding_dim = 3   # dimensione embedding
component_to_embed = 0  # componente da embeddare (0=X, 1=Y, 2=Z)

if USE_EMBEDDING:
    print(f"Parametri embedding: tau={embedding_tau}, dim={embedding_dim}, componente={component_to_embed}")
else:
    print("Usando output diretto del sistema di Lorenz (X, Y, Z)")
print()

for i, rho in enumerate(rho_values):
    print(f"\n--- Simulazione {i+1}/{len(rho_values)}: rho = {rho:.2f} ---")
    
    # Genera traiettoria completa del sistema di Lorenz
    print("Generando traiettoria completa del sistema di Lorenz...")
    lorenz_series = [np.array([1.0, 1.0, 1.0])]  # condizioni iniziali [X, Y, Z]
    for step in range(num_steps):
        if step % 20000 == 0:
            print(f"  Passo {step}/{num_steps}")
        x_new = lorenz_series[-1] + dt * lorenz_system(lorenz_series[-1], rho)
        lorenz_series.append(x_new)
    lorenz_series = np.array(lorenz_series)
    print(f"Traiettoria completa generata: {len(lorenz_series)} punti")
    
    # Scegli il metodo in base al flag
    if USE_EMBEDDING:
        # Metodo embedding
        print(f"METODO EMBEDDING:")
        # Estrai la componente desiderata
        component_series = lorenz_series[:, component_to_embed]  # estrai X, Y o Z
        print(f"Componente {['X', 'Y', 'Z'][component_to_embed]} estratta: {len(component_series)} punti")
        
        # Crea embedding 3D dalla singola componente
        print(f"Creando embedding 3D dalla componente {['X', 'Y', 'Z'][component_to_embed]}...")
        try:
            x_series = create_3d_embedding(component_series, embedding_tau, embedding_dim)
            print(f"Embedding 3D creato: {x_series.shape}")
        except ValueError as e:
            print(f"Errore nell'embedding: {e}")
            lyap_trad_list.append(np.nan)
            lyap_tau_list.append(np.nan)
            continue
    else:
        # Metodo diretto
        print(f"METODO DIRETTO:")
        x_series = lorenz_series  # Usa direttamente la traiettoria 3D [X, Y, Z]
        print(f"Usando traiettoria diretta 3D: {x_series.shape}")

    # Lyapunov tradizionale (sempre usando il sistema originale)
    print("Calcolando Lyapunov tradizionale...")
    v = np.array([1.0, 0.0, 0.0])
    lyap_sum = 0
    # Usa sempre la traiettoria originale per il Jacobiano
    for x in lorenz_series:
        J = jacobian(x, rho)
        v = v + dt * J @ v
        norm_v = np.linalg.norm(v)
        v = v / norm_v
        lyap_sum += np.log(norm_v)
    lyap_trad = lyap_sum / (num_steps * dt)
    lyap_trad_list.append(lyap_trad)
    print(f"Lyapunov tradizionale: {lyap_trad:.4f}")

    # Lyap_tau (usando il metodo scelto)
    method_name = "embedding" if USE_EMBEDDING else "diretto"
    print(f"Calcolando Lyapunov tau su metodo {method_name}...")
    recurrence_pairs = select_y_by_recurrence(x_series, delta0, tol, t_steps)
    if recurrence_pairs:
        lyap_tau_val = lyap_tau_continuous_multi(recurrence_pairs, DELTA, delta0, tol=tol)
    else:
        lyap_tau_val = np.nan
        print("  Lyapunov tau: Non calcolabile (nessuna ricorrenza)")
    lyap_tau_list.append(lyap_tau_val)

    print(f"RISULTATO: rho={rho:.2f} | Lyap trad={lyap_trad:.3f} | Lyap tau={lyap_tau_val:.3f}")

# ===============================
# Plot finale
# ===============================
print("\n=== CREAZIONE GRAFICO ===")
plt.figure(figsize=(8,6))
plt.scatter(rho_values, lyap_trad_list, color='black', label='Lyap tradizionale')
plt.scatter(rho_values, lyap_tau_list, color='orange', label='Lyap tau')
plt.axhline(0, color='k', linestyle=':')
plt.xlabel("rho")
plt.ylabel("Lyapunov exponent")
plt.ylim(-5, 5)
plt.legend()
plt.grid(True)

# Titolo e nome file dipendono dal metodo usato
if USE_EMBEDDING:
    component_name = ['X', 'Y', 'Z'][component_to_embed]
    title = f"Lyapunov vs rho (Embedding 3D della componente {component_name})"
    filename = f"lyapunov_vs_rho_embedding_{component_name}.png"
else:
    title = "Lyapunov vs rho (Output diretto del sistema di Lorenz)"
    filename = "lyapunov_vs_rho_direct.png"

plt.title(title)
plt.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()
print(f"Grafico salvato come '{filename}'")
print("=== SIMULAZIONE COMPLETATA ===")
