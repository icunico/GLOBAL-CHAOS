import sympy as sp

# Definiamo le variabili simboliche
x, y, r, K, a, b, d = sp.symbols('x y r K a b d')

# Equazioni del sistema
dx_dt = r * x * (1 - x / K) - y * (a * x / (1 + b * x))  # Equazione per dx/dt = 0
dy_dt = y * (a * x / (1 + b * x)- d)  # Equazione per dy/dt = 0

# Risolviamo il sistema per dx/dt = 0 e dy/dt = 0
equilibrio = sp.solve([dx_dt, dy_dt], (x, y))

# Mostriamo i punti di equilibrio
print("Punti di equilibrio:")
for punto in equilibrio:
    print(f"x = {punto[0]}, y = {punto[1]}")
