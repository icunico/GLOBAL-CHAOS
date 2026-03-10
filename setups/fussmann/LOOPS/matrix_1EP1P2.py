import sympy as sp
import numpy as np
import pandas as pd
from itertools import combinations

# =========================
# Variabili di stato (6)
# =========================
D, N, P1, C1, P2, C2 = sp.symbols('D N P1 C1 P2 C2')
state_vars = [D, N, P1, C1, P2, C2]

# =========================
# Parametro di controllo
# =========================
dC1_sym = sp.symbols('dC1')

# =========================
# Parametri fissi
# =========================
r1, K1 = 3.0, 1.0
r2, K2 = 2.0, 1.0
aP1C1, aP2C1 = 9.0, 2.0
aP1C2, aP2C2 = 3.0, 6.0
bP1C1, bP2C1 = 5.0, 5.0
bP1C2, bP2C2 = 5.0, 5.0
dC2, tau = 0.6, 50.0

# =========================
# Equazioni del sistema
# =========================
dDdt = dC1_sym*C1 + dC2*C2 - D/tau
dNdt = D/tau - N*(r1*P1/(N+K1) + r2*P2/(N+K2))

dP1dt = P1*( 
    r1*N/(N+K1)
    - aP1C1*C1/(1 + bP1C1*P1 + bP2C1*P2)
    - aP1C2*C2/(1 + bP1C2*P1 + bP2C2*P2)
)

dP2dt = P2*( 
    r2*N/(N+K2)
    - aP2C1*C1/(1 + bP1C1*P1 + bP2C1*P2)
    - aP2C2*C2/(1 + bP1C2*P1 + bP2C2*P2)
)

dC1dt = C1*((aP1C1*P1 + aP2C1*P2)/(1 + bP1C1*P1 + bP2C1*P2) - dC1_sym)
dC2dt = C2*((aP1C2*P1 + aP2C2*P2)/(1 + bP1C2*P1 + bP2C2*P2) - dC2)

equations = [dDdt, dNdt, dP1dt, dC1dt, dP2dt, dC2dt]

# =========================
# Jacobiana 6x6
# =========================
J = sp.Matrix([[sp.diff(eq, var) for var in state_vars] for eq in equations])
J_simplified = J.applyfunc(sp.simplify)

print("\n" + "="*80)
print("JACOBIANA SIMBOLICA 6x6")
print("="*80)
sp.pprint(J_simplified, use_unicode=True)

# =========================
# Lettura dati
# =========================
df = pd.read_csv(
    'risultati_dC1.txt',
    sep='\s+',
    skipinitialspace=True,
    comment='#',
    header=None
)

if df.shape[1] == 8:
    df.columns = ['idx', 'dC1', 'P1', 'P2', 'C1', 'C2', 'N', 'D']
else:
    df.columns = ['dC1', 'P1', 'P2', 'C1', 'C2', 'N', 'D']

# =========================
# Analisi di stabilità
# =========================
risultati = []

for idx, row in df.iterrows():
    if row.isnull().any():
        continue

    try:
        dC1_val = float(row['dC1'])
        P1_val = float(row['P1'])
        P2_val = float(row['P2'])
        C1_val = float(row['C1'])
        C2_val = float(row['C2'])
        N_val  = float(row['N'])
        D_val  = float(row['D'])
    except:
        continue

    subs_dict = {
        D: D_val,
        N: N_val,
        P1: P1_val,
        C1: C1_val,
        P2: P2_val,
        C2: C2_val,
        dC1_sym: dC1_val
    }

    try:
        J_numeric = J_simplified.subs(subs_dict).evalf()

        # =========================
        # F1 = traccia
        # =========================
        F1 = sum(J_numeric[i,i] for i in range(6))

        # =========================
        # F2 = somma dei minori 2x2, segno invertito
        # =========================
        two_cycles = 0
        for i,j in combinations(range(6),2):
            minor = J_numeric[i,j]*J_numeric[j,i]
            two_cycles += minor
        diag_products = 0
        for i,j in combinations(range(6),2):
            diag_products += J_numeric[i,i]*J_numeric[j,j]
        F2 = (two_cycles - diag_products)

        # =========================
        # F3 = somma dei minori 3x3
        # =========================
        F3 = 0
        for comb in combinations(range(6),3):
            minor_det = J_numeric.extract(comb,comb).det()
            F3 += minor_det

        # =========================
        # F4 = somma dei minori 4x4, segno invertito
        # =========================
        F4 = 0
        for comb in combinations(range(6),4):
            minor_det = J_numeric.extract(comb,comb).det()
            F4 += minor_det
        F4 = -F4

        # =========================
        # F5 = somma dei minori 5x5
        # =========================
        F5 = 0
        for comb in combinations(range(6),5):
            minor_det = J_numeric.extract(comb,comb).det()
            F5 += minor_det

        # =========================
        # F6 = determinante 6x6, segno invertito
        # =========================
        F6 = -J_numeric.det()

        # =========================
        # Verifica stabilità: tutte le forze < 0
        # =========================
        stabile = all(float(Fk) < 0 for Fk in [F1,F2,F3,F4,F5,F6])

        risultati.append([
            dC1_val,
            'STABILE' if stabile else 'INSTABILE',
            float(F1), float(F2), float(F3), float(F4), float(F5), float(F6)
        ])

    except Exception as e:
        continue

# =========================
# Salvataggio file TXT
# =========================
with open('forze_stabilita.txt', 'w') as f:
    f.write("# dC1\tSTATO\tF1\tF2\tF3\tF4\tF5\tF6\n")
    f.write("# Tutte le forze devono essere < 0 per stabilità\n")
    for r in risultati:
        f.write(f"{r[0]:.6f}\t{r[1]}\t{r[2]:.6e}\t{r[3]:.6e}\t{r[4]:.6e}\t{r[5]:.6e}\t{r[6]:.6e}\t{r[7]:.6e}\n")

print(f"\n✅ Salvati {len(risultati)} punti in forze_stabilita.txt")
print("\nPrime 5 righe:")
print(f"{'dC1':>10} {'STATO':>12} {'F1':>12} {'F2':>12} {'F3':>12} {'F4':>12} {'F5':>12} {'F6':>12}")
for r in risultati[:5]:
    print(f"{r[0]:10.6f} {r[1]:>12} {r[2]:12.6e} {r[3]:12.6e} {r[4]:12.6e} {r[5]:12.6e} {r[6]:12.6e} {r[7]:12.6e}")
