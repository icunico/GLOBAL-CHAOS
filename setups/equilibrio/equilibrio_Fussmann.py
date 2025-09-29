import sympy as sp

# Definiamo le variabili simboliche
P1, P2, P3, C, D, N, K1, K2, K3, a1, b1, d, r1, tau, a2, b2, r2, a3, b3, r3 = sp.symbols('P1 P2 P3 C D N K1 K2 K3 a1 b1 d r1 tau a2 b2 r2 a3 b3 r3')

# Equazioni del sistema
dP1_dt = 2.5 * P1 * N / (1+N) - P1 * 7.5 * C / (1 +  5*P1+ 5*P2 +5*P3) # Equazione per dx/dt = 0
dP2_dt = 2.5 * P2 * N / (1+N) - P2 * 2.5 * C / (1 +  5*P2+5*P1+ 5*P3) # Equazione per dx/dt = 0
dP3_dt = 2.5 * P3 * N / (1+N) - P3 * 5.0 * C / (1 +  5*P1+5*P2+5*P3) # Equazione per dx/dt = 0
dC_dt = C * ((7.5 * P1 + 2.5 * P2+ 5.0*P3) / (1 + 5*P1+ 5*P2 +5*P3)  - d)  # Equazione per dy/dt = 0
dD_dt = d*C - (1/50)*D
dN_dt = (1/50)*D - 2.5 * P1 * N/ (1+N)-  2.5 * P2 * N/ (1+N) -2.5 * P3 * N/ (1+N)

#dP1_dt = r1 * P1 * N / (K1+N) - P1 * a1 * C / (1 + b1 * P1+ b2*P2) # Equazione per dx/dt = 0
#dP2_dt = r2 * P2 * N / (K2+N) - P2 * a2 * C / (1 + b2 * P2+b1*P1) # Equazione per dx/dt = 0
#dC_dt = C * ((a1 * P1 + a2 * P2) / (1 + b1 * P1+ b2 * P2)  - d)  # Equazione per dy/dt = 0
#dD_dt = d*C - (1/tau)*D 
#dN_dt = (1/tau)*D - r1 * P1 * N/ (K1+N)-  r2 * P2 * N/ (K2+N)  
# Risolviamo il sistema per dx/dt = 0 e dy/dt = 0
equilibrio = sp.solve([dP1_dt, dP2_dt, dP3_dt, dC_dt, dD_dt, dN_dt], (P1, P2, P3, C, D, N))

# Mostriamo i punti di equilibrio
print("Punti di equilibrio:")
for punto in equilibrio:
    print(f"P = {punto[0]}, C = {punto[1]}, D = {punto[2]}, N = {punto[3]}")
