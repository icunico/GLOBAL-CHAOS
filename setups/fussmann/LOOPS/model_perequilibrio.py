import ordpy as od
import numpy as np
import scipy.integrate
import pyfabm

# =========================
# Lista dei valori di dC1 da testare
# =========================
dC1_values = np.linspace(0.0, 1.6, 10)

# =========================
# Apri file di output per i risultati
# =========================
output_file = 'risultati_dC1.txt'
with open(output_file, 'w') as f:
    # Intestazione con tutte le specie
    f.write("# dC1\tP1\tP2\tC1\tX\tN\tD\n")
    f.write("#" + "-"*70 + "\n")

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
    # Specie da controllare (TUTTE)
    # =========================
    specie_da_controllare_nomi = ['P1/DW1', 'P1/DW2', 'C1/DWC1', 'X/DWX1', 'N/DWN1', 'D/DWD1']
    variable_names = [var.path for var in model.state_variables]
    
    # Mappa degli indici per tutte le specie
    specie_indices = {}
    for nome in specie_da_controllare_nomi:
        try:
            idx = variable_names.index(nome)
            specie_indices[nome] = idx
            print(f"  Trovato {nome} all'indice {idx}")
        except ValueError:
            print(f"  ATTENZIONE: {nome} non trovato!")
            specie_indices[nome] = None
    
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
            # Scrivi NaN per tutte le 6 specie
            f.write(f"{dC1_val:.8f}\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\n")
    else:
        t = sol.t
        y = sol.y.T

        if len(y) < 10000:
            print(f"  SIMULAZIONE INCOMPLETA - nessun equilibrio")
            with open(output_file, 'a') as f:
                f.write(f"{dC1_val:.8f}\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\n")
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
                f.write(f"{dC1_val:.8f}\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\n")
            continue

        # Steady-state detection su TUTTE le specie
        threshold = 0.01
        mean_y = np.mean(y[:, specie_da_controllare], axis=0)
        std_y = np.std(y[:, specie_da_controllare], axis=0)
        
        # Calcola ratio solo dove mean_y >= threshold
        ratio = np.zeros_like(mean_y)
        mask = mean_y >= threshold
        ratio[mask] = std_y[mask] / mean_y[mask]
        
        flag_ss = np.max(ratio) <= 0.01

        if flag_ss:
            # Estrai i valori all'equilibrio per TUTTE le specie
            p1_val = y[-1, specie_indices['P1/DW1']] if specie_indices['P1/DW1'] is not None else np.nan
            p2_val = y[-1, specie_indices['P1/DW2']] if specie_indices['P1/DW2'] is not None else np.nan
            c1_val = y[-1, specie_indices['C1/DWC1']] if specie_indices['C1/DWC1'] is not None else np.nan
            x_val = y[-1, specie_indices['X/DWX1']] if specie_indices['X/DWX1'] is not None else np.nan
            n_val = y[-1, specie_indices['N/DWN1']] if specie_indices['N/DWN1'] is not None else np.nan
            d_val = y[-1, specie_indices['D/DWD1']] if specie_indices['D/DWD1'] is not None else np.nan

            print(f"  EQUILIBRIO: P1={p1_val:.8f}, P2={p2_val:.8f}, C1={c1_val:.8f}, X={x_val:.8f}, N={n_val:.8f}, D={d_val:.8f}")

            with open(output_file, 'a') as f:
                f.write(f"{dC1_val:.8f}\t{p1_val:.8f}\t{p2_val:.8f}\t{c1_val:.8f}\t{x_val:.8f}\t{n_val:.8f}\t{d_val:.8f}\n")
        else:
            print(f"  NON EQUILIBRIO (caos/oscillazioni) - nessun equilibrio")
            with open(output_file, 'a') as f:
                f.write(f"{dC1_val:.8f}\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\n")

print(f"\n{'='*50}")
print(f"RISULTATI SALVATI IN: {output_file}")
print(f"{'='*50}")

# Mostra il contenuto del file
print("\nCONTENUTO DEL FILE:")
with open(output_file, 'r') as f:
    print(f.read())
