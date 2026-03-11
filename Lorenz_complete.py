import numpy as np
import matplotlib.pyplot as plt

# ===============================
# Parametri generali
# ===============================
num_steps = 10000
t_steps = 300
DELTAMAX = 1.0  # Soglia massima di separazione (rinominata da DELTA)
delta0 = 0.1
dt = 0.01
tol = 0.4  # tolleranza per il calcolo Lyap_tau

sigma = 10.0
beta = 8.0 / 3.0
rho_values = np.linspace(40, 40, 1)

# Parametri rumore (processo di Wiener)
ADD_NOISE = False  # Flag per aggiungere rumore
noise_intensity = 0.01  # Intensità del rumore (sigma del processo di Wiener)

# Flag per scegliere il metodo
USE_EMBEDDING = False  # Se True usa embedding, se False usa output diretto del modello

# Flag per azzerare la parte deterministica
ZERO_DETERMINISTIC = False  # Se True azzera la derivata deterministica del sistema di Lorenz

print("=== PARAMETRI SIMULAZIONE ===")
print(f"Numero di passi: {num_steps}")
print(f"Passi minimi tra ricorrenze: {t_steps}")
print(f"DELTAMAX (soglia massima): {DELTAMAX}")
print(f"delta0 (soglia ricorrenze): {delta0}")
print(f"Tolleranza: {tol}")
print(f"dt: {dt}")
print(f"Valori di rho: {rho_values}")
print(f"Metodo: {'Embedding 3D' if USE_EMBEDDING else 'Output diretto del modello'}")
print(f"Rumore: {'Attivo' if ADD_NOISE else 'Disattivo'}")
if ADD_NOISE:
    print(f"Intensità rumore: {noise_intensity}")
print(f"Parte deterministica: {'AZZERATA' if ZERO_DETERMINISTIC else 'Attiva'}")
if ZERO_DETERMINISTIC:
    print("  ATTENZIONE: Sistema puramente stocastico (solo rumore di Wiener)")
print()

# ===============================
# Sistema di Lorenz originale (full 3D) con rumore
# ===============================
def lorenz_system(x, rho):
    """
    Input:
        x: array numpy 1D di dimensione 3, stato del sistema [X, Y, Z]
        rho: float, parametro del sistema di Lorenz
    Output:
        array numpy 1D di dimensione 3, derivate del sistema di Lorenz [dX/dt, dY/dt, dZ/dt]
    """
    if ZERO_DETERMINISTIC:
        # Azzera la parte deterministica - solo rumore
        return np.array([0.0, 0.0, 0.0])
    else:
        # Sistema di Lorenz normale
        dX = sigma * (x[1] - x[0])  # dX/dt
        dY = x[0] * (rho - x[2]) - x[1]  # dY/dt
        dZ = x[0] * x[1] - beta * x[2]  # dZ/dt
        return np.array([dX, dY, dZ])

def wiener_increment(dt, noise_intensity, dim=3):
    """
    Input:
        dt: float, passo temporale
        noise_intensity: float, intensità del rumore (sigma del processo di Wiener)
        dim: int, dimensionalità del rumore (default=3 per X, Y, Z)
    Output:
        array numpy 1D, incremento del processo di Wiener per ogni dimensione
    """
    # Incremento di Wiener: dW = sqrt(dt) * N(0,1)
    # Con intensità: noise_intensity * sqrt(dt) * N(0,1)
    return noise_intensity * np.sqrt(dt) * np.random.randn(dim)

def lorenz_system_sde(x, rho, dt, noise_intensity):
    """
    Input:
        x: array numpy 1D di dimensione 3, stato del sistema [X, Y, Z]
        rho: float, parametro del sistema di Lorenz
        dt: float, passo temporale
        noise_intensity: float, intensità del rumore
    Output:
        array numpy 1D di dimensione 3, nuovo stato dopo integrazione SDE
    """
    # Termine deterministico (drift) - può essere azzerato con il flag
    drift = lorenz_system(x, rho)
    
    # Termine stocastico (diffusione) - processo di Wiener
    if ADD_NOISE:
        diffusion = wiener_increment(dt, noise_intensity, dim=3)
    else:
        diffusion = np.zeros(3)
    
    # Integrazione secondo Itō: dx = f(x)dt + g(x)dW
    # Se ZERO_DETERMINISTIC=True: dx = 0*dt + g(x)dW = dW (moto Browniano puro)
    # Se ZERO_DETERMINISTIC=False: dx = f(x)dt + g(x)dW (sistema stocastico completo)
    x_new = x + drift * dt + diffusion
    
    return x_new

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
    if ZERO_DETERMINISTIC:
        return np.array([0.0, 0.0, 0.0])
    else:
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
    if ZERO_DETERMINISTIC:
        # Jacobiano nullo per sistema puramente stocastico
        return np.array([
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.0]
        ])
    else:
        # Jacobiano normale del sistema di Lorenz
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

def lyap_tau_continuous_multi(recurrence_pairs, DELTAMAX, delta0, tol=0.2):
    """
    Input:
        recurrence_pairs: lista di tuple (x_vals, y_vals), coppie di traiettorie ricorrenti
        DELTAMAX: float, soglia massima di separazione
        delta0: float, separazione iniziale di riferimento
        tol: float, tolleranza per identificare l'intorno di delta0 (default=0.2)
    Output:
        tuple: (tempo_predicibilita_tau, tempo_critico) 
               tempo_predicibilita_tau: float, tempo di predicibilità calcolato con il metodo tau,
               tempo_critico: float, tempo critico (mean_tau * dt)
               o (np.nan, np.nan) se non calcolabile
    """
    if not recurrence_pairs:
        print("  Nessuna ricorrenza fornita per il calcolo!")
        return np.nan, np.nan
    
    print(f"  Calcolo tempo di predicibilità tau su {len(recurrence_pairs)} ricorrenze...")
    
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
                if dist >= DELTAMAX:
                    tau_list.append(tau_current)
                    d0_list.append(d0_current)
                    counting = False
                    tau_current = 0
        
        # Aggiungi i tau e d0 di questa ricorrenza al totale
        all_tau_list.extend(tau_list)
        all_d0_list.extend(d0_list)
        print(f"      Trovati {len(tau_list)} intervalli tau per questa ricorrenza")
    
    # Calcola il tempo di predicibilità finale usando tutti i tau raccolti
    if all_tau_list:
        mean_tau = np.mean(all_tau_list)
        mean_d0 = np.mean(all_d0_list)
        
        # Tempo critico (tau medio in unità temporali)
        tempo_critico = mean_tau * dt
        
        # Tempo di predicibilità tau: T = (1/λ) * ln(DELTAMAX/delta0)
        # dove λ = 1/(mean_tau * dt) * ln(DELTAMAX/mean_d0)
        # Quindi: T = (mean_tau * dt) * ln(DELTAMAX/delta0) / ln(DELTAMAX/mean_d0)
        tempo_predicibilita_tau = tempo_critico * np.log(DELTAMAX / delta0) / np.log(DELTAMAX / mean_d0)
        
        print(f"  Totale intervalli tau trovati: {len(all_tau_list)}")
        print(f"  Tau medio finale: {mean_tau:.2f}, d0 medio finale: {mean_d0:.4f}")
        print(f"  Tempo critico: {tempo_critico:.3f}")
        print(f"  Tempo di predicibilità tau: {tempo_predicibilita_tau:.3f}")
    else:
        tempo_predicibilita_tau = np.nan
        tempo_critico = np.nan
        print("  Nessun intervallo tau trovato in tutte le ricorrenze!")
    
    return tempo_predicibilita_tau, tempo_critico

# ===============================
# Loop su rho (con flag per scegliere il metodo)
# ===============================
print("=== INIZIO SIMULAZIONI ===")
tempo_pred_trad_list = []
tempo_pred_tau_list = []
tempo_critico_list = []

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
    
    # Fissa il seed per riproducibilità (opzionale)
    np.random.seed(42 + i)  # seed diverso per ogni valore di rho
    
    # Genera traiettoria completa del sistema di Lorenz con rumore
    system_type = "puramente stocastico" if ZERO_DETERMINISTIC else ("con rumore" if ADD_NOISE else "deterministico")
    print(f"Generando traiettoria {system_type} del sistema di Lorenz...")
    lorenz_series = [np.array([1.0, 1.0, 1.0])]  # condizioni iniziali [X, Y, Z]
    for step in range(num_steps):
        if step % 20000 == 0:
            print(f"  Passo {step}/{num_steps}")
        
        # Usa la nuova funzione SDE che include il rumore
        x_new = lorenz_system_sde(lorenz_series[-1], rho, dt, noise_intensity)
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
            tempo_pred_trad_list.append(np.nan)
            tempo_pred_tau_list.append(np.nan)
            tempo_critico_list.append(np.nan)
            continue
    else:
        # Metodo diretto
        print(f"METODO DIRETTO:")
        x_series = lorenz_series  # Usa direttamente la traiettoria 3D [X, Y, Z]
        print(f"Usando traiettoria diretta 3D: {x_series.shape}")

    # Tempo di predicibilità tradizionale (basato su Lyapunov tradizionale)
    print("Calcolando tempo di predicibilità tradizionale...")
    v = np.array([1.0, 0.0, 0.0])
    lyap_sum = 0
    
    # Usa il sistema (deterministico o azzerato) per il Jacobiano
    for x in lorenz_series:
        J = jacobian(x, rho)
        v = v + dt * J @ v
        norm_v = np.linalg.norm(v)
        v = v / norm_v
        lyap_sum += np.log(norm_v)
    
    lyap_trad = lyap_sum / (num_steps * dt)
    
    # Converti in tempo di predicibilità: T = (1/λ) * ln(DELTAMAX/delta0)
    if ZERO_DETERMINISTIC:
        # Per sistema puramente stocastico, il Lyapunov è 0 (o molto piccolo)
        tempo_pred_trad = np.inf  # Tempo di predicibilità infinito per rumore puro
        print("  Sistema puramente stocastico: λ ≈ 0, tempo di predicibilità infinito")
    elif lyap_trad > 0:
        tempo_pred_trad = (1.0 / lyap_trad) * np.log(DELTAMAX / delta0)
    else:
        tempo_pred_trad = np.inf if lyap_trad < 0 else np.nan  # Sistema stabile o neutro
    
    tempo_pred_trad_list.append(tempo_pred_trad)
    print(f"Lyapunov tradizionale: {lyap_trad:.6f}")
    print(f"Tempo di predicibilità tradizionale: {tempo_pred_trad:.3f}")

    # Tempo di predicibilità tau (usando il metodo scelto)
    method_name = "embedding" if USE_EMBEDDING else "diretto"
    system_status = "puramente stocastico" if ZERO_DETERMINISTIC else ("con rumore" if ADD_NOISE else "deterministico")
    print(f"Calcolando tempo di predicibilità tau su metodo {method_name} {system_status}...")
    recurrence_pairs = select_y_by_recurrence(x_series, delta0, tol, t_steps)
    if recurrence_pairs:
        tempo_pred_tau_val, tempo_critico_val = lyap_tau_continuous_multi(recurrence_pairs, DELTAMAX, delta0, tol=tol)
    else:
        tempo_pred_tau_val = np.nan
        tempo_critico_val = np.nan
        print("  Tempo di predicibilità tau: Non calcolabile (nessuna ricorrenza)")
    
    tempo_pred_tau_list.append(tempo_pred_tau_val)
    tempo_critico_list.append(tempo_critico_val)

    print(f"RISULTATO: rho={rho:.2f} | T_pred trad={tempo_pred_trad:.3f} | T_pred tau={tempo_pred_tau_val:.3f} | T_critico={tempo_critico_val:.3f}")

# ===============================
# Plot finale
# ===============================
print("\n=== CREAZIONE GRAFICO ===")
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))

# Plot principale: Tempi di predicibilità
ax1.scatter(rho_values, tempo_pred_trad_list, color='black', label='Tempo pred. tradizionale', s=30, alpha=0.7)
ax1.scatter(rho_values, tempo_pred_tau_list, color='orange', label='Tempo pred. tau', s=30, alpha=0.7)
ax1.set_xlabel("ρ (rho)")
ax1.set_ylabel("Tempo di predicibilità")
ax1.set_ylim(0, 10)  # Limita per visualizzazione (ignora infiniti)
ax1.legend()
ax1.grid(True, alpha=0.3)
ax1.set_title("Tempi di Predicibilità vs ρ")

# Plot secondario: Tempo critico
valid_critico = ~np.isnan(tempo_critico_list)
if np.any(valid_critico):
    ax2.scatter(np.array(rho_values)[valid_critico], np.array(tempo_critico_list)[valid_critico], 
               color='red', label='Tempo critico (τ medio)', s=30, alpha=0.7)
ax2.set_xlabel("ρ (rho)")
ax2.set_ylabel("Tempo critico")
ax2.legend()
ax2.grid(True, alpha=0.3)
ax2.set_title("Tempo Critico vs ρ")

plt.tight_layout()

# Nome file dipende dal metodo usato
noise_suffix = f"_noise_{noise_intensity}" if ADD_NOISE else "_no_noise"
det_suffix = "_zero_det" if ZERO_DETERMINISTIC else ""

if USE_EMBEDDING:
    component_name = ['X', 'Y', 'Z'][component_to_embed]
    title_main = f"Tempi di Predicibilità vs rho (Embedding 3D della componente {component_name})"
    if ZERO_DETERMINISTIC:
        title_main += " - Sistema puramente stocastico"
    elif ADD_NOISE:
        title_main += f" - Rumore σ={noise_intensity}"
    filename = f"tempo_predicibilita_vs_rho_embedding_{component_name}{noise_suffix}{det_suffix}.png"
else:
    title_main = "Tempi di Predicibilità vs rho (Output diretto del sistema di Lorenz)"
    if ZERO_DETERMINISTIC:
        title_main += " - Sistema puramente stocastico"
    elif ADD_NOISE:
        title_main += f" - Rumore σ={noise_intensity}"
    filename = f"tempo_predicibilita_vs_rho_direct{noise_suffix}{det_suffix}.png"

fig.suptitle(title_main, fontsize=12)
plt.savefig(filename, dpi=300, bbox_inches='tight')
plt.close()
print(f"Grafico salvato come '{filename}'")

# Salva anche i dati in formato ASCII per l'analisi
data_filename = filename.replace('.png', '.dat')
with open(data_filename, 'w') as f:
    f.write("# rho tempo_pred_trad tempo_pred_tau tempo_critico\n")
    for i in range(len(rho_values)):
        f.write(f"{rho_values[i]:.6f} {tempo_pred_trad_list[i]:.6f} {tempo_pred_tau_list[i]:.6f} {tempo_critico_list[i]:.6f}\n")

print(f"Dati salvati in '{data_filename}'")
print("=== SIMULAZIONE COMPLETATA ===")
