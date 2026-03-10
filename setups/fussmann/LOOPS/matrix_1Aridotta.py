import sympy as sp
import numpy as np
import pandas as pd

# =========================
# Variabili indipendenti (D è espresso tramite bilancio di massa)
# =========================
N, P1, C1 = sp.symbols('N P1 C1')
vars_selected = [N, P1, C1]

# =========================
# Parametri (incluso dC1 come simbolo)
# =========================
r1 = 2.0
K1 = 1.0
aP1C1 = 2
bP1C1 = 5
tau = 50.0
M = 6.0  # Bilancio di massa totale
dC1 = sp.symbols('dC1')  # dC1 come simbolo

# =========================
# Bilancio di massa: esprimere D in funzione delle altre variabili
# =========================
# M = P1 + C1 + N + D  -->  D = M - P1 - C1 - N
D_expr = M - P1 - C1 - N
print(f"D in funzione delle altre variabili (M={M}):")
print(f"D = {D_expr}")
print()

# =========================
# Equazioni del sistema (sostituendo D)
# =========================
dNdt = D_expr/tau - N*(r1*P1/(N+K1))
dP1dt = P1*( r1*N/(N+K1) - aP1C1*C1/(1 + bP1C1*P1) )
dC1dt = C1*( (aP1C1*P1)/(1 + bP1C1*P1) - dC1 )

equations = [dNdt, dP1dt, dC1dt]

# =========================
# Jacobiana 3x3 simbolica
# =========================
J = sp.Matrix([[sp.diff(eq, var) for var in vars_selected] for eq in equations])
J_simplified = J.applyfunc(sp.simplify)

print("Jacobiana simbolica (sistema ridotto a 3 variabili):")
sp.pprint(J_simplified, use_unicode=True)

# =========================
# Funzione per analizzare un singolo punto di equilibrio
# =========================
def analyze_equilibrium(row, dC1_value):
    print(f"\n{'='*50}")
    print(f"ANALISI PER dC1 = {dC1_value:.6f}")
    print(f"{'='*50}")
    
    # Verifica se ci sono i valori necessari
    if pd.isna(row['P1']) or pd.isna(row['C1']):
        print("❌ Valori P1 o C1 NaN - impossibile calcolare")
        return None
    
    # Estrai i valori dal file
    P1_val = row['P1']
    C1_val = row['C1']
    
    # Per N, se non disponibile usiamo un valore tipico
    if pd.isna(row['N']):
        print("⚠️ N non disponibile nel file - uso valore tipico")
        N_val = 0.02  # Valore tipico dai dati
    else:
        N_val = row['N']
    
    # Calcola D dal bilancio di massa
    D_val = M - P1_val - C1_val - N_val
    print(f"P1 = {P1_val:.6f}")
    print(f"C1 = {C1_val:.6f}")
    print(f"N = {N_val:.6f}")
    print(f"D calcolato = {D_val:.6f}")
    
    # Verifica se D dal file è disponibile per confronto
    if not pd.isna(row['D']):
        print(f"D dal file = {row['D']:.6f}")
        print(f"Differenza: {abs(D_val - row['D']):.6f}")
    
    # Sostituisci i valori
    subs_dict = {
        N: N_val, 
        P1: P1_val, 
        C1: C1_val,
        dC1: dC1_value
    }
    
    J_numeric = J_simplified.subs(subs_dict)
    
    print("\nJacobiana numerica all'equilibrio (3x3):")
    # Arrotonda per una visualizzazione più chiara
    J_rounded = sp.Matrix([[round(float(J_numeric[i,j]), 6) for j in range(3)] for i in range(3)])
    sp.pprint(J_rounded, use_unicode=True)
    
    # =========================
    # Traccia < 0 (F1)
    # =========================
    trace = sum(J_numeric[i, i] for i in range(3))
    print(f"\nTraccia = {float(trace):.6f}")
    if trace < 0:
        print("Traccia < 0 ✅")
    else:
        print("Traccia >= 0 ❌")
    
    # =========================
    # Somma dei minori principali di ordine 2 > 0 (F2)
    # =========================
    # Calcolo i minori principali di ordine 2
    minor12 = J_numeric[0,0]*J_numeric[1,1] - J_numeric[0,1]*J_numeric[1,0]
    minor13 = J_numeric[0,0]*J_numeric[2,2] - J_numeric[0,2]*J_numeric[2,0]
    minor23 = J_numeric[1,1]*J_numeric[2,2] - J_numeric[1,2]*J_numeric[2,1]
    
    sum_minors = minor12 + minor13 + minor23
    print(f"\nSomma minori principali ordine 2 = {float(sum_minors):.6f}")
    if sum_minors > 0:
        print("Somma minori > 0 ✅")
    else:
        print("Somma minori <= 0 ❌")
    
    # =========================
    # Determinante < 0 (F3)
    # =========================
    det = float(J_numeric.det())
    print(f"\nDeterminante = {det:.6f}")
    if det < 0:
        print("Determinante < 0 ✅")
    else:
        print("Determinante >= 0 ❌")
    
    # =========================
    # Calcolo autovalori
    # =========================
    try:
        J_array = np.array(J_numeric).astype(np.float64)
        eigenvals = np.linalg.eigvals(J_array)
        print("\nAutovalori:")
        stable = True
        routh_hurwitz = (trace < 0 and sum_minors > 0 and det < 0)
        
        for i, ev in enumerate(eigenvals):
            if np.iscomplex(ev):
                print(f"  λ{i} = {ev.real:.6f} + {ev.imag:.6f}i")
                if ev.real >= 0:
                    stable = False
            else:
                print(f"  λ{i} = {ev.real:.6f}")
                if ev >= 0:
                    stable = False
        
        print(f"\nCriterio di Routh-Hurwitz: {'✅' if routh_hurwitz else '❌'}")
        print(f"Stabilità effettiva: {'✅' if stable else '❌'}")
        
    except Exception as e:
        print(f"\n⚠️ Impossibile calcolare gli autovalori: {e}")
        stable = None
        routh_hurwitz = (trace < 0 and sum_minors > 0 and det < 0)
        print(f"Criterio di Routh-Hurwitz: {'✅' if routh_hurwitz else '❌'}")
    
    return {
        'dC1': dC1_value,
        'P1': P1_val,
        'C1': C1_val,
        'N': N_val,
        'D_calc': D_val,
        'traccia': float(trace),
        'somma_minori': float(sum_minors),
        'det': det,
        'Routh_Hurwitz': routh_hurwitz,
        'stable': stable if stable is not None else routh_hurwitz
    }

# =========================
# Lettura del file e analisi per diversi valori di dC1
# =========================
print("\n" + "="*60)
print("ANALISI PER DIVERSI VALORI DI dC1 (SISTEMA RIDOTTO 3x3)")
print("="*60)

# Leggi il file
try:
    df = pd.read_csv('risultati_dC1.txt', comment='#', delim_whitespace=True, 
                     names=['dC1', 'P1', 'P2', 'C1', 'X', 'N', 'D'])
    
    print(f"\nColonne trovate: {list(df.columns)}")
    
    # Filtra righe dove P1 e C1 non sono NaN
    df_valid = df.dropna(subset=['P1', 'C1'])
    
    print(f"\nTrovate {len(df_valid)} righe con P1 e C1 validi su {len(df)} totali")
    
    if len(df_valid) > 0:
        print("\nValori di dC1 con dati validi:")
        for val in df_valid['dC1'].values:
            print(f"  {val:.6f}")
        
        # Analizza ogni punto di equilibrio
        results = []
        for idx, row in df_valid.iterrows():
            result = analyze_equilibrium(row, row['dC1'])
            if result:
                results.append(result)
        
        # =========================
        # Riepilogo finale
        # =========================
        print("\n" + "="*60)
        print("RIEPILOGO DEI RISULTATI (SISTEMA 3x3)")
        print("="*60)
        
        if results:
            summary_df = pd.DataFrame(results)
            pd.set_option('display.max_columns', None)
            pd.set_option('display.width', None)
            print("\n" + summary_df.to_string())
            
            # Salva i risultati
            summary_df.to_csv('analisi_stabilita_3x3.csv', index=False, float_format='%.6f')
            print("\n✅ Risultati salvati in 'analisi_stabilita_3x3.csv'")
            
            # Statistiche
            n_stabili = sum(summary_df['stable'])
            print(f"\nPunti stabili: {n_stabili}/{len(results)}")
        else:
            print("❌ Nessun risultato valido da analizzare")
    else:
        print("❌ Nessuna riga con P1 e C1 validi trovata nel file")
        
except FileNotFoundError:
    print("❌ File 'risultati_dC1.txt' non trovato!")
    print("Assicurati che il file sia nella stessa directory dello script")
except Exception as e:
    print(f"❌ Errore nella lettura del file: {e}")
    import traceback
    traceback.print_exc()
