import sympy as sp
import numpy as np

# =========================
# Leggi il file risultati_dC1.txt
# =========================
input_file = 'risultati_dC1.txt'
output_file = 'forze_dC1.txt'

with open(output_file, 'w') as f_out:
    f_out.write("# dC1\tF1\tF2\tF3\tF4\tF1<0\tF2<0\tF3<0\tF4<0\tstabilità\n")
    f_out.write("#" + "-"*90 + "\n")

with open(input_file, 'r') as f_in:
    lines = f_in.readlines()

data_lines = [line for line in lines if not line.startswith('#') and line.strip()]

print(f"Trovate {len(data_lines)} righe nel file.")
print("="*70)

# =========================
# Variabili selezionate - simboliche
# =========================
D, N, P1, C1 = sp.symbols('D N P1 C1')
vars_selected = [D, N, P1, C1]

# =========================
# Parametri simbolici
# =========================
r1 = 2.0
K1 = 1.0
aP1C1 = 2
bP1C1 = 5
dC1_sym = sp.symbols('dC1')
tau = 50.0

# =========================
# Equazioni del sistema simboliche
# =========================
dDdt = dC1_sym*C1 - D/tau
dNdt = D/tau - N*(r1*P1/(N+K1))
dP1dt = P1*( r1*N/(N+K1) - aP1C1*C1/(1 + bP1C1*P1) )
dC1dt = C1*( (aP1C1*P1)/(1 + bP1C1*P1) - dC1_sym )

equations = [dDdt, dNdt, dP1dt, dC1dt]

# =========================
# Jacobiana simbolica
# =========================
J = sp.Matrix([[sp.diff(eq, var) for var in vars_selected] for eq in equations])
J_simplified = J.applyfunc(sp.simplify)

# =========================
# STAMPA LA MATRICE SIMBOLICA SUL TERMINALE
# =========================
print("\n" + "="*70)
print("JACOBIANA SIMBOLICA (in funzione di D, N, P1, C1, dC1)")
print("="*70)
sp.pprint(J_simplified, use_unicode=True)
print("="*70 + "\n")

contatore_equilibri = 0

for line in data_lines:
    parts = line.strip().split()
    if len(parts) < 5:
        continue

    try:
        dC1_val = float(parts[0])  # dC1
        P1_val = float(parts[2])   # P1/DW2
        C1_val = float(parts[3])   # C1/DWC1
        N_val = float(parts[1])    # N/DWN1
        D_val = float(parts[4])    # D/DWD1
    except:
        continue

    if np.isnan(P1_val) or np.isnan(C1_val) or np.isnan(N_val) or np.isnan(D_val):
        continue

    contatore_equilibri += 1

    # =========================
    # Sostituzione valori di equilibrio
    # =========================
    equilibrium_values = {
        D: D_val,
        N: N_val,
        P1: P1_val,
        C1: C1_val,
        dC1_sym: dC1_val
    }

    J_numeric = J_simplified.subs(equilibrium_values).evalf()

    # =========================
    # Calcolo forze
    # =========================
    F1 = sum(J_numeric[i, i] for i in range(4))
    F1_cond = F1 < 0

    F2 = (
        J_numeric[1,0]*J_numeric[0,1] + J_numeric[2,0]*J_numeric[0,2] +
        J_numeric[2,1]*J_numeric[1,2] + J_numeric[3,0]*J_numeric[0,3] +
        J_numeric[3,1]*J_numeric[1,3] + J_numeric[3,2]*J_numeric[2,3]
        -
        (J_numeric[0,0]*J_numeric[1,1] + J_numeric[1,1]*J_numeric[2,2] +
         J_numeric[2,2]*J_numeric[3,3] + J_numeric[0,0]*J_numeric[2,2] +
         J_numeric[0,0]*J_numeric[3,3] + J_numeric[1,1]*J_numeric[3,3])
    )
    F2_cond = F2 < 0

    F3 = (
        J_numeric[1,0]*J_numeric[2,1]*J_numeric[0,2] +
        J_numeric[2,0]*J_numeric[1,2]*J_numeric[0,1] +
        J_numeric[3,0]*J_numeric[2,1]*J_numeric[0,3] +
        J_numeric[3,1]*J_numeric[2,3]*J_numeric[0,1] +
        J_numeric[3,2]*J_numeric[1,3]*J_numeric[0,2] +
        J_numeric[3,0]*J_numeric[1,2]*J_numeric[0,3]
        -
        (J_numeric[1,0]*J_numeric[0,1]*J_numeric[2,2] +
         J_numeric[2,0]*J_numeric[0,2]*J_numeric[1,1] +
         J_numeric[2,1]*J_numeric[1,2]*J_numeric[0,0] +
         J_numeric[3,0]*J_numeric[0,3]*J_numeric[1,1] +
         J_numeric[3,1]*J_numeric[1,3]*J_numeric[2,2] +
         J_numeric[3,2]*J_numeric[2,3]*J_numeric[0,0])
        +
        J_numeric[0,0]*J_numeric[1,1]*J_numeric[2,2] +
        J_numeric[0,0]*J_numeric[1,1]*J_numeric[3,3] +
        J_numeric[0,0]*J_numeric[2,2]*J_numeric[3,3] +
        J_numeric[1,1]*J_numeric[2,2]*J_numeric[3,3]
    )
    F3_cond = F3 < 0

    F4 = J_numeric.det()
    F4_cond = F4 < 0

    stabile = F1_cond and F2_cond and F3_cond and F4_cond

    # =========================
    # Stampa a schermo (solo dC1 e forze)
    # =========================
    print(f"\ndC1 = {dC1_val:.6f}")
    print(f"F1 = {float(F1):.6e}  {'✅' if F1_cond else '❌'}")
    print(f"F2 = {float(F2):.6e}  {'✅' if F2_cond else '❌'}")
    print(f"F3 = {float(F3):.6e}  {'✅' if F3_cond else '❌'}")
    print(f"F4 = {float(F4):.6e}  {'✅' if F4_cond else '❌'}")
    print(f"Stabilità: {'STABILE ✅' if stabile else 'INSTABILE ❌'}")

    # =========================
    # Scrivi su file (SOLO dC1 e forze)
    # =========================
    with open(output_file, 'a') as f_out:
        f_out.write(f"{dC1_val:.6f}\t")
        f_out.write(f"{float(F1):.6e}\t{float(F2):.6e}\t{float(F3):.6e}\t{float(F4):.6e}\t")
        f_out.write(f"{F1_cond}\t{F2_cond}\t{F3_cond}\t{F4_cond}\t")
        f_out.write(f"{'STABILE' if stabile else 'INSTABILE'}\n")

print(f"\n{'='*70}")
print(f"ANALISI COMPLETATA!")
print(f"Trovati {contatore_equilibri} equilibri nel file {input_file}")
print(f"Risultati salvati in: {output_file}")
print(f"{'='*70}")
