import sympy as sp
import numpy as np
import pandas as pd
from itertools import combinations

# =========================
# Variabili di stato (6 variabili totali)
# =========================
D, N, P1, P2, C1, X = sp.symbols('D N P1 P2 C1 X')
state_vars = [D, N, P1, P2, C1, X]

# =========================
# Bilancio di massa totale
# =========================
M_tot = sp.symbols('M_tot')  # massa totale conservata

# Relazione di bilancio: M_tot = D + N + P1 + P2 + C1 + X
# Esprimiamo D in funzione delle altre
bilancio_massa = sp.Eq(M_tot, D + N + P1 + P2 + C1 + X)
D_expr = sp.solve(bilancio_massa, D)[0]  # D = M_tot - (N + P1 + P2 + C1 + X)
print(f"\nD espressa in funzione delle altre: D = {D_expr}")

# =========================
# Parametro di controllo
# =========================
dC1_sym = sp.symbols('dC1')

# =========================
# Parametri per i predatori (P1, P2)
# =========================
r1, K1 = 3.0, 1.0
r2, K2 = 2.0, 1.0

aP1C1, aP2C1 = 9.0, 2.0
bP1C1, bP2C1 = 5.0, 5.0

aP1X, aP2X = 0.25, 0.25
bP1X, bP2X = 0.5, 0.5

# =========================
# Parametri per C1
# =========================
aC1X = 3.0
bC1X = 2.0

# =========================
# Parametri per X
# =========================
dX = 0.5
tau = 50.0

# =========================
# Equazioni del sistema (senza dDdt, che è ridondante)
# =========================
# dDdt non è più un'equazione indipendente (D è determinato dal bilancio)
dNdt = D/tau - N*(r1*P1/(N+K1) + r2*P2/(N+K2))

dP1dt = P1*(
    r1*N/(N+K1)
    - aP1C1*C1/(1 + bP1C1*P1 + bP2C1*P2)
    - aP1X*X/(1 + bP1X*P1 + bP2X*P2 + bC1X*C1)
)

dP2dt = P2*(
    r2*N/(N+K2)
    - aP2C1*C1/(1 + bP1C1*P1 + bP2C1*P2)
    - aP2X*X/(1 + bP1X*P1 + bP2X*P2 + bC1X*C1)
)

dC1dt = C1*(
    (aP1C1*P1 + aP2C1*P2)/(1 + bP1C1*P1 + bP2C1*P2)
    - aC1X*X/(1 + bP1X*P1 + bP2X*P2 + bC1X*C1)
    - dC1_sym
)

dXdt = X*(
    (aP1X*P1 + aP2X*P2 + aC1X*C1) /
    (1 + bP1X*P1 + bP2X*P2 + bC1X*C1)
    - dX
)

# Lista delle 5 equazioni (senza dDdt)
equations = [dNdt, dP1dt, dP2dt, dC1dt, dXdt]
reduced_vars = [N, P1, P2, C1, X]  # Variabili indipendenti

# =========================
# Sostituiamo D con la sua espressione in tutte le equazioni
# =========================
equations_subs = [eq.subs(D, D_expr) for eq in equations]

# =========================
# Jacobiana 5x5 (sulle variabili indipendenti)
# =========================
J_reduced = sp.Matrix([[sp.diff(eq, var) for var in reduced_vars] for eq in equations_subs])
J_reduced_simplified = J_reduced.applyfunc(sp.simplify)

print("\n" + "="*80)
print("JACOBIANA RIDOTTA 5x5 (con D espresso dal bilancio di massa)")
print("="*80)
sp.pprint(J_reduced_simplified, use_unicode=True)

# =========================
# Lettura dati - CORRETTA
# =========================
df = pd.read_csv(
    'risultati_dC1.txt',
    sep='\s+',
    skipinitialspace=True,
    comment='#',
    header=None
)

# Assegnazione corretta delle colonne in base all'header mostrato
# dC1, P1, P2, C1, X, N, D
df.columns = ['dC1', 'P1', 'P2', 'C1', 'X', 'N', 'D']

print("Colonne nel dataframe:", df.columns.tolist())
print("\nPrime righe del file:")
print(df.head())

# =========================
# Analisi di stabilità (sistema ridotto)
# =========================
risultati = []
M_tot_val = 10.0  # Valore della massa totale conservata (da aggiustare)

for idx, row in df.iterrows():
    if row.isnull().any():
        continue

    try:
        dC1_val = float(row['dC1'])
        P1_val = float(row['P1'])
        P2_val = float(row['P2'])
        C1_val = float(row['C1'])
        X_val = float(row['X'])
        N_val = float(row['N'])
        # D_val non serve più dall'input, lo calcoliamo
        
        # Calcoliamo D dal bilancio di massa
        D_val = M_tot_val - (N_val + P1_val + P2_val + C1_val + X_val)
        
        # Verifica che D sia positivo (opzionale)
        if D_val < 0:
            print(f"  Attenzione: D negativo ({D_val:.4f}) per dC1={dC1_val:.4f}")
        
    except Exception as e:
        print(f"Errore alla riga {idx}: {e}")
        continue

    subs_dict = {
        N: N_val,
        P1: P1_val,
        P2: P2_val,
        C1: C1_val,
        X: X_val,
        dC1_sym: dC1_val,
        M_tot: M_tot_val  # M_tot serve anche per le derivate
    }

    try:
        J_numeric = J_reduced_simplified.subs(subs_dict).evalf()

        # Calcolo delle forze per sistema 5x5 (F1-F5)
        n_vars = 5
        
        F1 = sum(J_numeric[i,i] for i in range(n_vars))

        two_cycles = 0
        for i,j in combinations(range(n_vars),2):
            minor = J_numeric[i,j]*J_numeric[j,i]
            two_cycles += minor
        diag_products = 0
        for i,j in combinations(range(n_vars),2):
            diag_products += J_numeric[i,i]*J_numeric[j,j]
        F2 = (two_cycles - diag_products)

        F3 = 0
        for comb in combinations(range(n_vars),3):
            minor_det = J_numeric.extract(comb,comb).det()
            F3 += minor_det

        F4 = 0
        for comb in combinations(range(n_vars),4):
            minor_det = J_numeric.extract(comb,comb).det()
            F4 += minor_det
        F4 = -F4

        F5 = J_numeric.det()

        stabile = all(float(Fk) < 0 for Fk in [F1, F2, F3, F4, F5])

        risultati.append([
            dC1_val,
            'STABILE' if stabile else 'INSTABILE',
            float(F1), float(F2), float(F3), float(F4), float(F5)
        ])
    except Exception as e:
        print(f"Errore nel calcolo per dC1={dC1_val}: {e}")
        continue

# =========================
# Salvataggio file TXT
# =========================
with open('forze_stabilita_ridotto.txt', 'w') as f:
    f.write("# dC1\tSTATO\tF1\tF2\tF3\tF4\tF5\n")
    f.write("# Sistema ridotto 5x5 con D = M_tot - (N+P1+P2+C1+X)\n")
    f.write(f"# M_tot = {M_tot_val}\n")
    f.write("# Tutte le forze devono essere < 0 per stabilità\n")
    for r in risultati:
        f.write(f"{r[0]:.6f}\t{r[1]}\t{r[2]:.6e}\t{r[3]:.6e}\t{r[4]:.6e}\t{r[5]:.6e}\t{r[6]:.6e}\n")

print(f"\n✅ Salvati {len(risultati)} punti in forze_stabilita_ridotto.txt")
if risultati:
    print("\nPrime 5 righe:")
    print(f"{'dC1':>10} {'STATO':>12} {'F1':>12} {'F2':>12} {'F3':>12} {'F4':>12} {'F5':>12}")
    for r in risultati[:5]:
        print(f"{r[0]:10.6f} {r[1]:>12} {r[2]:12.6e} {r[3]:12.6e} {r[4]:12.6e} {r[5]:12.6e} {r[6]:12.6e}")
