import ordpy as od
import numpy as np
import scipy.integrate
import pyfabm

# =========================
# Lista dei valori di dC1 da testare
# =========================
dC1_values = np.linspace(0.25, 0.4, 20)

# =========================
# Apri file di output per i risultati
# =========================
output_file = 'risultati_dC1.txt'
with open(output_file, 'w') as f:
    f.write("# dC1\tN1/DWN1\tP1/DWP2\tC1/DWC1\tD/DWD1\n")
    f.write("#" + "-"*60 + "\n")

# =========================
# Ciclo sui parametri
# =========================
for dC1_val in dC1_values:
    print(f"\nTESTANDO dC1 = {dC1_val}")
    
    # =========================
    # Caricamento modello
    # =========================
    model = pyfabm.Model('/leonardo/home/userexternal/icunico0/seamless-notebooks/setups/fussmann/fabm.yaml')
    model.cell_thickness = 1.
    
    # =========================
    # IMPOSTA IL PARAMETRO dC1
    # =========================
    try:
        model.findParameter('fussmann/D/dC1').value = dC1_val
        print(f"  Parametro impostato a {dC1_val}")
    except:
        try:
            model.findParameter('D/dC1').value = dC1_val
            print(f"  Parametro impostato a {dC1_val}")
        except:
            print(f"  ERRORE: Impossibile impostare dC1 = {dC1_val}")
            continue
    
    # =========================
    # Inizializza il modello
    # =========================
    model.start()
    assert model.checkReady(), 'Model not ready'
    
    # =========================
    # Specie da controllare
    # =========================
    specie_da_controllare_nomi = ['P1/DW2', 'C1/DWC1', 'N/DWN1', 'D/DWD1']
    variable_names = [var.path for var in model.state_variables]
    specie_da_controllare = [i for i, name in enumerate(variable_names)
                             if name in specie_da_controllare_nomi]

    # =========================
    # Contatori evento estinzione
    # =========================
    event_counters = {iv: 0 for iv in specie_da_controllare}

    # =========================
    # Funzione derivata temporale
    # =========================
    def dy(t0, y):
        model.state[:] = y
        return model.getRates()

    # =========================
    # Evento estinzione
    # =========================
    def event(t, y):
        for iv in specie_da_controllare:
            if y[iv] < 0.001:
                event_counters[iv] += 1
                if event_counters[iv] >= 20000:
                    return 0
            else:
                event_counters[iv] = 0
        return 1

    event.terminal = True
    event.direction = -1

    # =========================
    # Integrazione temporale
    # =========================
    t_ex = 3600 * 100000
    t_eval = np.linspace(0, t_ex, 50000)

    try:
        sol = scipy.integrate.solve_ivp(
            dy,
            [0., t_ex],
            model.state,
            t_eval=t_eval,
            max_step=5000,
            events=event
        )
    except Exception as e:
        print(f"  ERRORE integrazione: {e}")
        continue

    # =========================
    # Analisi risultati
    # =========================
    if any(ev.size > 0 for ev in sol.t_events):
        print(f"  ESTINZIONE - nessun equilibrio")
        with open(output_file, 'a') as f:
            f.write(f"{dC1_val:.4f}\tNaN\tNaN\tNaN\tNaN\n")
    else:
        t = sol.t
        y = sol.y.T
        
        if len(y) < 10000:
            print(f"  SIMULAZIONE INCOMPLETA - nessun equilibrio")
            with open(output_file, 'a') as f:
                f.write(f"{dC1_val:.4f}\tNaN\tNaN\tNaN\tNaN\n")
            continue
            
        y = y[-10000:]  # prendi ultimi 10000 punti per l'analisi

        # Controllo negatività
        flag_neg = False
        for iv in specie_da_controllare:
            if np.any(y[:, iv] < -0.01):
                flag_neg = True
                break

        if flag_neg:
            print(f"  NEGATIVITÀ - nessun equilibrio")
            with open(output_file, 'a') as f:
                f.write(f"{dC1_val:.4f}\tNaN\tNaN\tNaN\tNaN\n")
            continue

        # Steady-state detection
        threshold = 0.01
        mean_y = np.mean(y[:, specie_da_controllare], axis=0)
        std_y = np.std(y[:, specie_da_controllare], axis=0)
        ratio = np.where(mean_y >= threshold, std_y / mean_y, 0)
        flag_ss = np.max(ratio) <= 0.01

        if flag_ss:
            # Salva i valori all'equilibrio
            p1_val = y[-1, specie_da_controllare[0]]
            c1_val = y[-1, specie_da_controllare[1]]
            n_val = y[-1, specie_da_controllare[2]]
            d_val = y[-1, specie_da_controllare[3]]
            
            print(f"  EQUILIBRIO: P1={p1_val:.8f}, C1={c1_val:.8f}, N={n_val:.8f}, D={d_val:.8f}")
            
            with open(output_file, 'a') as f:
                f.write(f"{dC1_val:.8f}\t{p1_val:.8f}\t{c1_val:.8f}\t{n_val:.8f}\t{d_val:.8f}\n")
        else:
            print(f"  NON EQUILIBRIO (caos/oscillazioni) - nessun equilibrio")
            with open(output_file, 'a') as f:
                f.write(f"{dC1_val:.4f}\tNaN\tNaN\tNaN\tNaN\n")

print(f"\n{'='*50}")
print(f"RISULTATI SALVATI IN: {output_file}")
print(f"{'='*50}")

# Mostra il contenuto del file
print("\nCONTENUTO DEL FILE:")
with open(output_file, 'r') as f:
    print(f.read())
