import jax.numpy as jnp
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt
from jax import random, jit
import networkx as nx
import scienceplots

plt.style.use('science')

def build_system(N, f, pars, balance_mass=True):
    """
    Builds a system of N equations in a callable format for ODE solvers.
    
    Parameters:
        N (int): Number of equations.
        f (function): Non-linearity function.
        pars (array-like): Fixed parameters for each equation.
        balance_mass (bool): Whether to use mass balance constraint.
    
    Returns:
        system (function): Callable system function that takes (t, x) for integration.
    """
    
    if balance_mass:
        @jit
        def system(t, x):
            dxdt = jnp.zeros(N)
            for i in range(N):
                dxdt = dxdt.at[i].set(pars[i-1] * f(x[i-1]) * x[i-1] - pars[i] * f(x[i]) * x[i])
            return dxdt
    else:
        @jit
        def system(t, x):
            dxdt = jnp.zeros(N)
            for i in range(N):
                dxdt = dxdt.at[i].set(pars[i-1] * f(x[i-1]) * x[i-1])
            return dxdt
    
    return system

def total_mass(N, f, pars, balance_mass=True):

    if balance_mass:
        @jit
        def system(t, x):
            dxdt = jnp.zeros(N)
            for i in range(N):
                dxdt = dxdt.at[i].set(pars[i-1] * f(x[i-1]) * x[i-1] - pars[i] * f(x[i]) * x[i])
            return dxdt
    else:
        @jit
        def system(t, x):
            dxdt = jnp.zeros(N)
            for i in range(N):
                dxdt = dxdt.at[i].set(pars[i-1] * f(x[i-1]) * x[i-1])
            return dxdt
    
    return system


def func(x):
    return (jnp.exp(-x)/(1+jnp.exp(-x))**2)

# Define N nodes
N = 6

# Create directed loop graph
G = nx.DiGraph()
for i in range(N):
    G.add_edge(i, (i + 1) % N)
plt.figure(figsize=(8, 6))
pos = nx.spring_layout(G)
nx.draw(G, pos, with_labels=True, node_color='lightblue', edge_color='gray', node_size=500, font_size=10, 
font_weight='bold')
plt.title('Directed Loop Graph')
plt.show()
plt.clf()
plt.close()

# Define seed and random parameters
seed = 420
key = random.key(seed)
key, subkey = random.split(key)
pars = random.uniform(subkey, (N,), minval=-15, maxval=5)
print(pars)

# Build system of equations
system = build_system(N, func, pars, balance_mass=True)

# Example integration
t_span = (0, 1000)
t_eval = jnp.linspace(t_span[0], t_span[1], int(1e5))
key, subkey = random.split(key)
x0 = random.uniform(subkey, (N,), minval=0, maxval=10)
print(jnp.sum(x0))
solution = solve_ivp(system, t_span, x0, method='RK45', t_eval=t_eval)

# Plot solutions
plt.figure(figsize=(10, 6))
for i in range(N):
    plt.plot(solution.t, solution.y[i], label=f'Variable {i+1}')
plt.xlabel('Time')
plt.ylabel('Solution Values')
plt.title('Time Evolution of the System')
plt.legend()
plt.grid(True)
plt.show()


plt.figure(figsize=(10, 6))
plt.plot(solution.t, jnp.sum(solution.y, axis=0))
plt.xlabel('Time')
plt.ylabel('Solution Values')
plt.title('Time Evolution of the System')
plt.grid(True)
plt.show()
