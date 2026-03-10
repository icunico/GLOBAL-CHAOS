import sympy as sp

# =========================
# Variabili selezionate
# =========================
D, N, P1, C1 = sp.symbols('P1 C1 N D')
vars_selected = [D, N, P1, C1]

# =========================
# Parametri
# =========================
r1 = 2.0
K1 = 1.0
aP1C1 = 2
bP1C1 = 5
dC1 = 0.35
tau = 50.0

# =========================
# Equazioni del sistema
# =========================
dDdt = dC1*C1 - D/tau
dNdt = D/tau - N*(r1*P1/(N+K1))
dP1dt = P1*( r1*N/(N+K1) - aP1C1*C1/(1 + bP1C1*P1) )
dC1dt = C1*( (aP1C1*P1)/(1 + bP1C1*P1) - dC1 )

equations = [dDdt, dNdt, dP1dt, dC1dt]

# =========================
# Jacobiana 4x4 simbolica
# =========================
J = sp.Matrix([[sp.diff(eq, var) for var in vars_selected] for eq in equations])
J_simplified = J.applyfunc(sp.simplify)

print("Jacobiana simbolica:")
sp.pprint(J_simplified, use_unicode=True)

# =========================
# Sostituzione valori di equilibrio
# =========================
equilibrium_values = {D: 3.38, N: 0.02, P1: 1.40, C1: 0.19}


J_numeric = J_simplified.subs(equilibrium_values)

print("\nJacobiana numerica all'equilibrio:")
sp.pprint(J_numeric, use_unicode=True)

# =========================
# F1 = somma dei diagonali < 0
# =========================
F1 = sum(J_numeric[i, i] for i in range(4))
print("\nF1 =", F1)
if F1 < 0:
    print("F1 < 0 ✅")
else:
    print("F1 >= 0 ❌")

# =========================
# F2 = combinazioni di prodotti a due < 0
# =========================

F2 = (
    # termini incrociati
    J_numeric[1,0]*J_numeric[0,1] + J_numeric[2,0]*J_numeric[0,2] + J_numeric[2,1]*J_numeric[1,2] +
    J_numeric[3,0]*J_numeric[0,3] + J_numeric[3,1]*J_numeric[1,3] + J_numeric[3,2]*J_numeric[2,3]
    -
    # termini diagonali
    (J_numeric[0,0]*J_numeric[1,1] + J_numeric[1,1]*J_numeric[2,2] + J_numeric[2,2]*J_numeric[3,3] +
     J_numeric[0,0]*J_numeric[2,2] + J_numeric[0,0]*J_numeric[3,3] + J_numeric[1,1]*J_numeric[3,3])
)
print("\nF2 =", F2)
if F2 < 0:
    print("F2 < 0 ✅")
else:
    print("F2 >= 0 ❌")

# =========================
# F3 = combinazioni triple < 0
# =========================

F3 = (
    # termini incrociati positivi (loop di 3 elementi)
    J_numeric[1,0]*J_numeric[2,1]*J_numeric[0,2] +
    J_numeric[2,0]*J_numeric[1,2]*J_numeric[0,1] +
    J_numeric[3,0]*J_numeric[2,1]*J_numeric[0,3] +
    J_numeric[3,1]*J_numeric[2,3]*J_numeric[0,1] +
    J_numeric[3,2]*J_numeric[1,3]*J_numeric[0,2] +
    J_numeric[3,0]*J_numeric[1,2]*J_numeric[0,3]
    -
    # termini combinati negativi (prodotti di 2 elementi * diagonale)
    (J_numeric[1,0]*J_numeric[0,1]*J_numeric[2,2] +
     J_numeric[2,0]*J_numeric[0,2]*J_numeric[1,1] +
     J_numeric[2,1]*J_numeric[1,2]*J_numeric[0,0] +
     J_numeric[3,0]*J_numeric[0,3]*J_numeric[1,1] +
     J_numeric[3,1]*J_numeric[1,3]*J_numeric[2,2] +
     J_numeric[3,2]*J_numeric[2,3]*J_numeric[0,0])
    +
    # termini completamente diagonali positivi
    J_numeric[0,0]*J_numeric[1,1]*J_numeric[2,2] +
    J_numeric[0,0]*J_numeric[1,1]*J_numeric[3,3] +
    J_numeric[0,0]*J_numeric[2,2]*J_numeric[3,3] +
    J_numeric[1,1]*J_numeric[2,2]*J_numeric[3,3]
)
print("\nF3 =", F3)
if F3 < 0:
    print("F3 < 0 ✅")
else:
    print("F3 >= 0 ❌")

# =========================
# F4 = -determinante < 0
# =========================
F4 = J_numeric.det()
print("\nF4 =", F4)
if F4 < 0:
    print("F4 < 0 ✅")
else:
    print("F4 >= 0 ❌")

