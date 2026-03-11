import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# ===============================
# Parametri
# ===============================
num_steps = 50000
transient = 5000
dt = 0.01
epsilon = 0.02
t_steps = 500
DELTA = 0.2
delta0 = 0.02

# Parametri per l'embedding
embedding_dim = 3
embedding_delay = 3

# Parametri fissi di Lorenz
sigma = 10
beta = 8/3

# ===============================
# SISTEMA DI LORENZ
# ===============================
def lorenz_system(state, t, sigma, rho, beta):
    x, y, z = state
    dx = sigma * (y - x)
    dy = x * (rho - z) - y
    dz = x * y - beta * z
    return [dx, dy, dz]

def generate_lorenz(rho, sigma=10, beta=8/3, dt=0.01, n_points=20000):
    initial_state = [1.0, 1.0, 1.0]
    t = np.linspace(0, n_points*dt, n_points)
    trajectory = odeint(lorenz_system, initial_state, t, args=(sigma, rho, beta))
    return trajectory

# ===============================
# LYAPUNOV TRADIZIONALE
# ===============================
def lyapunov_lorenz_traditional(rho, sigma=10, beta=8/3, n_iter=2000, dt=0.01):
    state = np.array([1.0, 1.0, 1.0])
    pert = np.array([1.0, 0.0, 0.0])
    
    lyap_sum = 0
    
    # Transiente
    for _ in range(1000):
        state = state + dt * np.array(lorenz_system(state, 0, sigma, rho, beta))
    
    # Calcolo esponente
    for i in range(n_iter):
        state_next = state + dt * np.array(lorenz_system(state, 0, sigma, rho, beta))
        
        # Jacobiano approssimato
        J = np.zeros((3, 3))
        delta = 1e-6
        for j in range(3):
            e = np.zeros(3)
            e[j] = delta
            state_pert = state + e
            state_pert_next = state_pert + dt * np.array(lorenz_system(state_pert, 0, sigma, rho, beta))
            J[:, j] = (state_pert_next - state_next) / delta
        
        # Evolvi perturbazione
        pert = J @ pert
        norm = np.linalg.norm(pert)
        
        if norm > 0:
            lyap_sum += np.log(norm)
            pert = pert / norm
        
        state = state_next
    
    return lyap_sum / (n_iter * dt)

# ===============================
# FUNZIONI PER LYAPUNOV TAU (MULTIPLE RICORRENZE)
# ===============================
def create_embedding(series, dim, delay):
    if len(series) < (dim - 1) * delay + 1:
        return np.array([])
    
    n = len(series) - (dim - 1) * delay
    embedded = np.zeros((n, dim))
    for i in range(dim):
        embedded[:, i] = series[i * delay:i * delay + n]
    return embedded

def euclidean_distance(v1, v2):
    return np.sqrt(np.sum((v1 - v2)**2))

def find_multiple_recurrences(series, epsilon, min_sep, dim, delay, n_recurrences=20):
    """
    Trova multiple ricorrenze nella serie
    Restituisce una lista di (emb1, emb2, distanza_iniziale)
    """
    emb = create_embedding(series, dim, delay)
    
    if len(emb) < min_sep:
        return []
    
    n = len(emb)
    recurrences = []
    used_indices = set()  # Per evitare di usare gli stessi punti
    
    print(f"   Cerco {n_recurrences} ricorrenze in {n} punti embedded...")
    
    for i in range(0, n - min_sep, 5):
        if len(recurrences) >= n_recurrences:
            break
            
        for j in range(i + min_sep, n, 5):
            if j in used_indices or i in used_indices:
                continue
                
            dist = euclidean_distance(emb[i], emb[j])
            if dist < epsilon:
                # Calcola lunghezza massima disponibile
                max_len = min(n - i, n - j, 300)
                if max_len < 50:
                    continue
                
                # Estrai le traiettorie
                emb1 = emb[i:i + max_len]
                emb2 = emb[j:j + max_len]
                
                recurrences.append({
                    'emb1': emb1,
                    'emb2': emb2,
                    'dist0': dist,
                    'idx1': i,
                    'idx2': j
                })
                
                used_indices.add(i)
                used_indices.add(j)
                break
    
    print(f"   Trovate {len(recurrences)} ricorrenze valide")
    return recurrences

def lyap_tau_from_recurrences(recurrences, DELTA, delta0, dt=0.01):
    """
    Calcola Lyapunov tau usando multiple ricorrenze
    """
    all_taus = []
    all_d0s = []
    all_lyap_individual = []
    
    for idx, rec in enumerate(recurrences):
        emb1 = rec['emb1']
        emb2 = rec['emb2']
        
        tau_list = []
        d0_list = []
        
        counting = False
        tau_current = 0
        d0_current = 0
        
        for i in range(len(emb1)):
            dist = euclidean_distance(emb1[i], emb2[i])
            
            if not counting:
                if dist < delta0:
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
            mean_tau_rec = np.mean(tau_list)
            mean_d0_rec = np.mean(d0_list)
            if mean_tau_rec > 0 and mean_d0_rec > 0:
                lyap_rec = (1.0 / (mean_tau_rec * dt)) * np.log(DELTA / mean_d0_rec)
                all_lyap_individual.append(lyap_rec)
                all_taus.extend(tau_list)
                all_d0s.extend(d0_list)
    
    if all_taus:
        # Lyapunov dalla media di tutti i tau
        mean_tau_global = np.mean(all_taus)
        mean_d0_global = np.mean(all_d0s)
        lyap_global = (1.0 / (mean_tau_global * dt)) * np.log(DELTA / mean_d0_global)
        
        # Lyapunov come media dei Lyapunov individuali
        lyap_mean_individual = np.mean(all_lyap_individual) if all_lyap_individual else np.nan
        
        # Statistiche
        stats = {
            'n_recurrences': len(recurrences),
            'n_events': len(all_taus),
            'mean_tau': mean_tau_global,
            'std_tau': np.std(all_taus),
            'mean_d0': mean_d0_global,
            'lyap_global': lyap_global,
            'lyap_mean_individual': lyap_mean_individual,
            'lyap_std': np.std(all_lyap_individual) if len(all_lyap_individual) > 1 else 0,
            'all_taus': all_taus,
            'all_lyap': all_lyap_individual
        }
        
        return stats
    else:
        return None

# ===============================
# LOOP SU DIVERSI VALORI DI RHO
# ===============================
print("="*70)
print("CALCOLO LYAPUNOV CON MULTIPLE RICORRENZE (10+ punti)")
print("="*70)

rho_values = np.linspace(40, 40, 1)
lyap_trad_list = []
lyap_tau_list = []
lyap_tau_std_list = []
rho_valid = []
n_recurrences_found = []

print(f"\nParametri: σ={sigma}, β={beta}")
print(f"Embedding: dim={embedding_dim}, delay={embedding_delay}")
print(f"Cerco 10 ricorrenze per ogni ρ")

for i, rho in enumerate(rho_values):
    print(f"\n[{i+1}/{len(rho_values)}] ρ = {rho:.1f}")
    print("-" * 50)
    
    # Genera traiettoria
    traj = generate_lorenz(rho, sigma, beta, dt, num_steps)
    traj = traj[transient:]
    series_x = traj[:, 0]
    
    # Lyapunov tradizionale
    lyap_trad = lyapunov_lorenz_traditional(rho, sigma, beta, n_iter=2000, dt=dt)
    print(f"  Tradizionale: {lyap_trad:.4f}")
    
    # Trova multiple ricorrenze
    recurrences = find_multiple_recurrences(series_x, epsilon, t_steps, 
                                           embedding_dim, embedding_delay, 
                                           n_recurrences=10)
    
    if len(recurrences) >= 3:  # Almeno 3 ricorrenze
        stats = lyap_tau_from_recurrences(recurrences, DELTA, delta0, dt)
        
        if stats:
            print(f"  Trovate: {stats['n_recurrences']} ricorrenze")
            print(f"  Eventi di divergenza: {stats['n_events']}")
            print(f"  Tau medio: {stats['mean_tau']:.1f} ± {stats['std_tau']:.1f}")
            print(f"  Lyap tau (media eventi): {stats['lyap_global']:.4f}")
            print(f"  Lyap tau (media ricorrenze): {stats['lyap_mean_individual']:.4f} ± {stats['lyap_std']:.4f}")
            print(f"  Differenza da tradizionale: {abs(lyap_trad - stats['lyap_global']):.4f}")
            
            lyap_trad_list.append(lyap_trad)
            lyap_tau_list.append(stats['lyap_global'])
            lyap_tau_std_list.append(stats['lyap_std'])
            rho_valid.append(rho)
            n_recurrences_found.append(stats['n_recurrences'])
            
# ===============================
# PLOT FINALE
# ===============================
if len(rho_valid) > 0:
    plt.figure(figsize=(14, 8))
    
    # Plot con barre di errore
    plt.errorbar(rho_valid, lyap_tau_list, yerr=lyap_tau_std_list, 
                 fmt='r-s', linewidth=2, markersize=8, capsize=5, capthick=2,
                 label=f'Lyapunov tau (media su {int(np.mean(n_recurrences_found))} ricorrenze)', 
                 markeredgecolor='darkred', ecolor='gray', elinewidth=2)
    
    plt.plot(rho_valid, lyap_trad_list, 'b-o', linewidth=2.5, markersize=8, 
             label='Lyapunov tradizionale (Wolf)', markeredgecolor='darkblue')
    
    plt.axhline(y=0, color='k', linestyle='--', linewidth=1.5, alpha=0.7)
    plt.axvline(x=24.06, color='gray', linestyle=':', linewidth=2, alpha=0.8)
    
    plt.xlabel('ρ (parametro di Rayleigh)', fontsize=14, fontweight='bold')
    plt.ylabel('Esponente di Lyapunov massimo', fontsize=14, fontweight='bold')
    plt.title('Confronto Lyapunov tradizionale vs metodo dei tau (multiple ricorrenze)', 
              fontsize=16, fontweight='bold')
    
    # Aggiungi info sul numero di ricorrenze
    for i, (rho, n_rec) in enumerate(zip(rho_valid, n_recurrences_found)):
        plt.annotate(f'n={n_rec}', xy=(rho, lyap_tau_list[i]), 
                    xytext=(5, 5), textcoords='offset points', fontsize=9,
                    bbox=dict(boxstyle="round,pad=0.2", facecolor='white', alpha=0.7))
    
    plt.legend(loc='best', fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.xlim(20, 50)
    
    plt.tight_layout()
    plt.savefig('lyapunov_lorenz_multiple_recurrences.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    # Statistiche finali
    print("\n" + "="*70)
    print("STATISTICHE FINALI")
    print("="*70)
    
    lyap_trad_array = np.array(lyap_trad_list)
    lyap_tau_array = np.array(lyap_tau_list)
    
    print(f"\nPunti analizzati: {len(rho_valid)}/{len(rho_values)}")
    print(f"Media ricorrenze per punto: {np.mean(n_recurrences_found):.1f}")
    
    correlation = np.corrcoef(lyap_trad_array, lyap_tau_array)[0, 1]
    mae = np.mean(np.abs(lyap_trad_array - lyap_tau_array))
    
    print(f"\nCorrelazione: {correlation:.4f}")
    print(f"MAE: {mae:.4f}")
    
    print("\nTabella risultati:")
    print(f"{'ρ':>6} | {'Lyap trad':>10} | {'Lyap tau':>10} | {'Std dev':>8} | {'n_rec':>6}")
    print("-" * 50)
    
    for rho, ltrad, ltau, lstd, nrec in zip(rho_valid, lyap_trad_list, lyap_tau_list, lyap_tau_std_list, n_recurrences_found):
        print(f"{rho:6.1f} | {ltrad:10.4f} | {ltau:10.4f} | {lstd:8.4f} | {nrec:6d}")
