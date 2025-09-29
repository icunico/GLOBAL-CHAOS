import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def lyapunov_exponent(series, m=3, tau=2, dt=1, eps=1e-1, min_tsep=1):
    """
    Calcolo dell'esponente di Lyapunov massimo da una serie temporale scalare.
    
    Parametri:
    - series: array-like, serie temporale (es. biomassa)
    - m: embedding dimension
    - tau: time delay
    - dt: tempo tra campioni (default = 1)
    - eps: soglia massima per il vicino iniziale
    - min_tsep: separazione temporale minima tra punti (per evitare autocorrelazione)

    Ritorna:
    - lambda_max: esponente di Lyapunov stimato
    - num_pairs: numero di coppie usate
    """
    x = np.asarray(series)
    N = len(x) - (m - 1) * tau

    if N <= 0:
        raise ValueError("Serie troppo corta per embedding.")

    # Ricostruzione spazio delle fasi
    X = np.empty((N, m))
    for i in range(m):
        X[:, i] = x[i * tau : i * tau + N]

    sum_log = 0
    count = 0

    for i in range(N - 1):
        xi = X[i]
        # escludi vicini troppo vicini nel tempo
        dists = np.linalg.norm(X - xi, axis=1)
        dists[:i + min_tsep] = np.inf
        dmin = np.min(dists)
        j = np.argmin(dists)

        if dmin < eps and j + 1 < N and i + 1 < N:
            xi1 = X[i + 1]
            xj1 = X[j + 1]
            dist_next = np.linalg.norm(xi1 - xj1)
            if dist_next > 0 and dmin > 0:
                sum_log += np.log(dist_next / dmin)
                count += 1

    if count == 0:
        return np.nan, 0

    lambda_max = sum_log / (count * dt)
    return lambda_max, count


if __name__ == "__main__":
    # Carica dati da CSV
    file_path = 'chlorophyll_-65.0_0.0.csv'  # Modifica con il tuo file
    df = pd.read_csv(file_path, header=None)
    series = df.iloc[:, 0].values  # prima colonna

    # Calcola esponente
    lambda_max, n_pairs = lyapunov_exponent(series)

    print(f"Esponente di Lyapunov massimo: {lambda_max:.6f}")
    print(f"Numero di coppie usate: {n_pairs}")

    # Visualizza la serie
    plt.plot(series)
    plt.title("Serie temporale")
    plt.xlabel("Tempo")
    plt.ylabel("Biomassa")
    plt.grid(True)
    plt.show()

