import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint

# ===============================
# Parametri
# ===============================
num_steps = 50000
transient = 50
dt = 0.1
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

def generate_lorenz(rho, sigma=10, beta=8/3, dt=0.01, n_points=50000):
    initial_state = [1.0, 1.0, 1.0]
    t = np.linspace(0, n_points*dt, n_points)
    trajectory = odeint(lorenz_system, initial_state, t, args=(sigma, rho, beta))
    return trajectory

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

def find_multiple_recurrences(series, epsilon, min_sep, dim, delay, n_recurrences=10):
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
    τ è direttamente il tempo in secondi
    Accetta solo eventi con d₀ entro ±20% di delta0
    """
    all_taus = []
    all_d0s = []
    all_lyap_individual = []
    
    # Liste per stampare tutti i valori
    tutti_i_tau_per_rho = []
    tutte_le_d0_per_rho = []
    
    # Definisci l'intervallo accettabile per d₀ (±20% di delta0)
    d0_min = delta0 * 0.8
    d0_max = delta0 * 1.2
    
    print(f"   Accetto solo d₀ comprese tra {d0_min:.6f} e {d0_max:.6f}")
    
    for idx, rec in enumerate(recurrences):
        emb1 = rec['emb1']
        emb2 = rec['emb2']
        
        tau_list = []
        d0_list = []
        
        counting = False
        tau_current = 0.0  # τ in secondi
        d0_current = 0
        
        for i in range(len(emb1)):
            dist = euclidean_distance(emb1[i], emb2[i])
            
            if not counting:
                if dist < delta0:
                    counting = True
                    tau_current = 0.0  # Resetta il tempo
                    d0_current = dist
            else:
                tau_current += dt  # Aggiunge il tempo (in secondi)
                if dist >= DELTA:
                    # Accetta solo se d₀ è vicino a δ₀ (±20%)
                    if d0_min <= d0_current <= d0_max:
                        tau_list.append(tau_current)  # τ già in secondi
                        d0_list.append(d0_current)
                    counting = False
                    tau_current = 0.0
        
        if tau_list:
            # Per questa ricorrenza, calcola λ usando i valori
            for t, d in zip(tau_list, d0_list):
                if t > 0 and d > 0:
                    lyap_evento = (1.0 / t) * np.log(DELTA / d)
                    all_lyap_individual.append(lyap_evento)
            
            all_taus.extend(tau_list)
            all_d0s.extend(d0_list)
            
            # Salva per stampa globale
            tutti_i_tau_per_rho.extend(tau_list)
            tutte_le_d0_per_rho.extend(d0_list)
            
            print(f"     Ricorrenza {idx+1}: {len(tau_list)} eventi validi")
            for k, (t, d) in enumerate(zip(tau_list, d0_list)):
                lyap = (1.0 / t) * np.log(DELTA / d)
                print(f"       Evento {k+1}: τ = {t:.3f}s, d₀ = {d:.6f}, λ = {lyap:.4f}")
    
    if all_taus:
        # Lyapunov dalla media di tutti i τ (in secondi)
        mean_tau_global = np.mean(all_taus)
        mean_d0_global = np.mean(all_d0s)
        lyap_global = (1.0 / mean_tau_global) * np.log(DELTA / mean_d0_global)
        
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
            'all_d0s': all_d0s,
            'all_lyap': all_lyap_individual,
            'tutti_i_tau': tutti_i_tau_per_rho,
            'tutte_le_d0': tutte_le_d0_per_rho
        }
        
        # STAMPA DETTAGLIATA DI TUTTI I TAU PER QUESTO ρ
        print(f"\n   📋 LISTA COMPLETA DI TUTTI I {len(all_taus)} τ TROVATI (in secondi):")
        for idx, t in enumerate(all_taus):
            if idx % 5 == 0 and idx > 0:
                print(f"\n   ", end="")
            print(f"{t:.3f}s ", end="")
        print(f"\n")
        
        # Statistiche base
        print(f"   📊 STATISTICHE τ:")
        print(f"      Min: {min(all_taus):.3f} s")
        print(f"      Max: {max(all_taus):.3f} s")
        print(f"      Media: {mean_tau_global:.3f} s")
        print(f"      Mediana: {np.median(all_taus):.3f} s")
        print(f"      Dev std: {np.std(all_taus):.3f} s")
        
        return stats
    else:
        print(f"   Nessun evento valido trovato (d₀ fuori dall'intervallo ±20%)")
        return None

# ===============================
# VALORI TEORICI NOTI (da letteratura)
# ===============================
def theoretical_lyapunov(rho):
    """Restituisce il valore teorico dell'esponente di Lyapunov per Lorenz"""
    if rho < 24.06:
        return 0.0  # Non caotico
    else:
        # Approssimazione lineare dai valori noti in letteratura
        if rho <= 28:
            return 0.9 * (rho - 24.06) / (28 - 24.06)
        elif rho <= 30:
            return 0.9 + 0.21 * (rho - 28) / 2
        elif rho <= 35:
            return 1.11 + 0.19 * (rho - 30) / 5
        elif rho <= 40:
            return 1.3 + 0.2 * (rho - 35) / 5
        else:
            return 1.5 + 0.1 * (rho - 40) / 5

# ===============================
# LOOP SU DIVERSI VALORI DI RHO
# ===============================
print("="*70)
print("CALCOLO LYAPUNOV CON MULTIPLE RICORRENZE - d₀ ±20% di δ₀")
print("="*70)

rho_values = np.linspace(40, 40, 1)
lyap_tau_list = []
lyap_tau_std_list = []
rho_valid = []
n_recurrences_found = []

# Valori teorici per il plot
rho_teorici = np.linspace(10, 45, 100)
lyap_teorici = [theoretical_lyapunov(r) for r in rho_teorici]

print(f"\nParametri: σ={sigma}, β={beta}")
print(f"Embedding: dim={embedding_dim}, delay={embedding_delay}")
print(f"δ₀={delta0}, Δ={DELTA}, ε={epsilon}, dt={dt}")
print(f"Accetto solo d₀ tra {delta0*0.8:.6f} e {delta0*1.2:.6f}")

# Dizionari per salvare tutti i tau e d₀ per ogni ρ
tutti_i_tau_per_rho = {}
tutte_le_d0_per_rho = {}

for i, rho in enumerate(rho_values):
    print(f"\n{'='*60}")
    print(f"[{i+1}/{len(rho_values)}] ρ = {rho:.1f}")
    print(f"{'='*60}")
    
    # Genera traiettoria
    traj = generate_lorenz(rho, sigma, beta, dt, num_steps)
    traj = traj[transient:]
    series_x = traj[:, 0]
    
    print(f"\n🔍 Cerco ricorrenze...")
    recurrences = find_multiple_recurrences(series_x, epsilon, t_steps, 
                                           embedding_dim, embedding_delay, 
                                           n_recurrences=10)
    
    if len(recurrences) >= 3:
        print(f"\n📋 DETTAGLIO EVENTI PER ρ = {rho:.1f}:")
        print(f"   {'-'*50}")
        
        stats = lyap_tau_from_recurrences(recurrences, DELTA, delta0, dt)
        
        if stats:
            # Salva tutti i tau e d₀ per questo ρ
            tutti_i_tau_per_rho[rho] = stats['all_taus']
            tutte_le_d0_per_rho[rho] = stats['all_d0s']
            
            teorico = theoretical_lyapunov(rho)
            
            print(f"\n📊 STATISTICHE FINALI PER ρ = {rho:.1f}:")
            print(f"   Numero ricorrenze trovate: {stats['n_recurrences']}")
            print(f"   Numero totale eventi validi: {stats['n_events']}")
            print(f"   τ medio: {stats['mean_tau']:.3f} s")
            print(f"   d₀ media: {stats['mean_d0']:.6f}")
            print(f"   Lyapunov τ calcolato: {stats['lyap_global']:.6f} ± {stats['lyap_std']:.6f}")
            print(f"   Valore teorico: {teorico:.6f}")
            print(f"   Differenza dal teorico: {abs(teorico - stats['lyap_global']):.6f}")
            
            lyap_tau_list.append(stats['lyap_global'])
            lyap_tau_std_list.append(stats['lyap_std'])
            rho_valid.append(rho)
            n_recurrences_found.append(stats['n_recurrences'])
    else:
        print(f"\n❌ Non trovate abbastanza ricorrenze (solo {len(recurrences)})")

# ===============================
# STAMPA FINALE DI TUTTI I TAU PER OGNI ρ (RIEPILOGO)
# ===============================
print("\n" + "="*70)
print("RIEPILOGO FINALE - TUTTI I τ E d₀ PER OGNI ρ")
print("="*70)

for rho in rho_valid:
    if rho in tutti_i_tau_per_rho:
        tau_list = tutti_i_tau_per_rho[rho]
        d0_list = tutte_le_d0_per_rho[rho]
        
        print(f"\n{'─'*60}")
        print(f"ρ = {rho:.1f} - {len(tau_list)} eventi validi totali")
        print(f"{'─'*60}")
        print(f"{'N°':<4} | {'τ (s)':<10} | {'d₀':<12} | {'λ individuale':<12} | {'λ teorico':<10}")
        print(f"{'─'*60}")
        
        teorico_rho = theoretical_lyapunov(rho)
        for idx, (t, d) in enumerate(zip(tau_list, d0_list)):
            lyap_ind = (1.0 / t) * np.log(DELTA / d)
            print(f"{idx+1:<4} | {t:<10.3f} | {d:<12.6f} | {lyap_ind:<12.4f} | {teorico_rho:<10.4f}")
        
        # Statistiche di riepilogo per questo ρ
        print(f"\n   Media τ: {np.mean(tau_list):.3f} s")
        print(f"   Mediana τ: {np.median(tau_list):.3f} s")
        print(f"   Min τ: {min(tau_list):.3f} s")
        print(f"   Max τ: {max(tau_list):.3f} s")
        print(f"   Dev std τ: {np.std(tau_list):.3f} s")

# ===============================
# PLOT FINALE
# ===============================
if len(rho_valid) > 0:
    plt.figure(figsize=(14, 8))
    
    # Plot dei valori teorici (linea tratteggiata)
    plt.plot(rho_teorici, lyap_teorici, 'g--', linewidth=2.5, 
             label='Valori teorici (da letteratura)', alpha=0.8)
    
    # Plot Lyapunov τ con barre di errore
    plt.errorbar(rho_valid, lyap_tau_list, yerr=lyap_tau_std_list, 
                 fmt='r-s', linewidth=2, markersize=8, capsize=5, capthick=2,
                 label=f'Lyapunov τ (solo d₀ ±20% di δ₀)', 
                 markeredgecolor='darkred', ecolor='gray', elinewidth=2)
    
    plt.axhline(y=0, color='k', linestyle='--', linewidth=1.5, alpha=0.7)
    plt.axvline(x=24.06, color='gray', linestyle=':', linewidth=2, alpha=0.8, label='ρ = 24.06 (biforcazione)')
    
    plt.xlabel('ρ', fontsize=14, fontweight='bold')
    plt.ylabel('Esponente di Lyapunov massimo', fontsize=14, fontweight='bold')
    plt.title('Confronto: Metodo dei τ vs valori teorici noti', 
              fontsize=16, fontweight='bold')
    
    plt.legend(loc='best', fontsize=12)
    plt.grid(True, alpha=0.3)
    plt.xlim(10, 45)
    plt.ylim(-0.2, 2.0)
    
    plt.tight_layout()
    plt.savefig('lyapunov_lorenz_vs_teorico.png', dpi=150, bbox_inches='tight')
    plt.show()
    
    # Statistiche finali
    print("\n" + "="*70)
    print("STATISTICHE FINALI - TUTTI I ρ")
    print("="*70)
    
    print("\nTabella risultati:")
    print(f"{'ρ':>6} | {'Lyapunov τ':>12} | {'Std dev':>10} | {'Teorico':>10} | {'Differenza':>12} | {'Eventi':>8}")
    print("-" * 75)
    
    for rho, ltau, lstd in zip(rho_valid, lyap_tau_list, lyap_tau_std_list):
        teorico = theoretical_lyapunov(rho)
        diff = abs(teorico - ltau)
        n_eventi = len(tutti_i_tau_per_rho[rho])
        print(f"{rho:6.1f} | {ltau:12.6f} | {lstd:10.6f} | {teorico:10.6f} | {diff:12.6f} | {n_eventi:8d}")
