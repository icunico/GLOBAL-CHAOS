#!/usr/bin/env python3
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

print("=== Plot Tempo di Predicibilità (τ) ===")

# Leggi i dati
try:
    data = pd.read_csv('lyapunov_paladin_mpi.dat', sep='\t', comment='#',
                       names=['lat', 'lon', 'tau'])
    print(f"Lette {len(data)} righe")
    print(f"Valori validi: {data['tau'].count()}")
    print(f"Range lat: [{data['lat'].min():.2f}, {data['lat'].max():.2f}]")
    print(f"Range lon: [{data['lon'].min():.2f}, {data['lon'].max():.2f}]")
    
    # Filtra valori validi (non NaN)
    valid_data = data.dropna()
    print(f"\nStatistiche τ (tempo di predicibilità):")
    print(f"  Min: {valid_data['tau'].min():.4f}")
    print(f"  Max: {valid_data['tau'].max():.4f}")
    print(f"  Media: {valid_data['tau'].mean():.4f}")
    print(f"  Mediana: {valid_data['tau'].median():.4f}")
    print(f"  Std: {valid_data['tau'].std():.4f}")
    
except Exception as e:
    print(f"Errore lettura file: {e}")
    exit(1)

# 1. Mappa scatter plot 2D
plt.figure(figsize=(12, 8))
sc = plt.scatter(valid_data['lon'], valid_data['lat'], 
                c=valid_data['tau'], cmap='plasma', 
                s=30, alpha=0.8, edgecolors='none')
cbar = plt.colorbar(sc, label='Tempo di predicibilità τ')
cbar.ax.tick_params(labelsize=10)
plt.xlabel('Longitudine', fontsize=12)
plt.ylabel('Latitudine', fontsize=12)
plt.title('Tempo di predicibilità (metodo Paladin)', fontsize=14)
plt.grid(True, alpha=0.3, linestyle='--')
plt.tight_layout()
plt.savefig('taupred_mappa.png', dpi=150, bbox_inches='tight')
print("\nMappa salvata: taupred_mappa.png")

# 2. Istogramma della distribuzione
plt.figure(figsize=(10, 6))
n, bins, patches = plt.hist(valid_data['tau'], bins=50, 
                            alpha=0.7, color='steelblue', 
                            edgecolor='black', linewidth=0.5)
plt.axvline(valid_data['tau'].mean(), color='red', 
            linestyle='--', linewidth=2, 
            label=f'Media: {valid_data["tau"].mean():.4f}')
plt.axvline(valid_data['tau'].median(), color='green', 
            linestyle='--', linewidth=2, 
            label=f'Mediana: {valid_data["tau"].median():.4f}')
plt.xlabel('Tempo di predicibilità τ', fontsize=12)
plt.ylabel('Frequenza', fontsize=12)
plt.title('Distribuzione dei tempi di predicibilità', fontsize=14)
plt.legend()
plt.grid(True, alpha=0.3, linestyle='--')
plt.tight_layout()
plt.savefig('taupred_istogramma.png', dpi=150, bbox_inches='tight')
print("Istogramma salvato: taupred_istogramma.png")

# 3. Mappa con proiezione geografica (se vuoi qualcosa di più professionale)
try:
    import cartopy.crs as ccrs
    import cartopy.feature as cfeature
    
    plt.figure(figsize=(14, 8))
    ax = plt.axes(projection=ccrs.PlateCarree())
    
    # Aggiungi elementi geografici
    ax.add_feature(cfeature.COASTLINE, linewidth=0.8)
    ax.add_feature(cfeature.BORDERS, linewidth=0.5, alpha=0.5)
    ax.add_feature(cfeature.LAND, facecolor='lightgray', alpha=0.3)
    ax.add_feature(cfeature.OCEAN, facecolor='lightblue', alpha=0.3)
    
    # Griglia
    gl = ax.gridlines(draw_labels=True, dms=True, 
                      x_inline=False, y_inline=False,
                      linewidth=0.5, color='gray', alpha=0.5)
    gl.top_labels = False
    gl.right_labels = False
    
    # Scatter plot
    sc = ax.scatter(valid_data['lon'], valid_data['lat'], 
                   c=valid_data['tau'], cmap='plasma', 
                   s=30, alpha=0.8, transform=ccrs.PlateCarree(),
                   edgecolors='none')
    
    cbar = plt.colorbar(sc, ax=ax, orientation='horizontal', 
                       pad=0.05, aspect=40, shrink=0.8)
    cbar.set_label('Tempo di predicibilità τ', fontsize=12)
    
    plt.title('Tempo di predicibilità - Mappa geografica', fontsize=14)
    plt.tight_layout()
    plt.savefig('taupred_mappa_geografica.png', dpi=150, bbox_inches='tight')
    print("Mappa geografica salvata: taupred_mappa_geografica.png")
    
except ImportError:
    print("\nNota: Installa cartopy per mappe geografiche più accurate:")
    print("  conda install -c conda-forge cartopy")
    print("  oppure")
    print("  pip install cartopy")

# 4. Heatmap 2D (se i dati sono su griglia regolare)
try:
    # Crea una griglia regolare
    lats_unique = np.sort(valid_data['lat'].unique())
    lons_unique = np.sort(valid_data['lon'].unique())
    
    nlat = len(lats_unique)
    nlon = len(lons_unique)
    
    if nlat * nlon == len(valid_data):
        print(f"\nDati su griglia regolare {nlat} x {nlon}")
        
        # Crea matrice 2D
        tau_grid = np.full((nlat, nlon), np.nan)
        
        for _, row in valid_data.iterrows():
            i = np.where(lats_unique == row['lat'])[0][0]
            j = np.where(lons_unique == row['lon'])[0][0]
            tau_grid[i, j] = row['tau']
        
        # Plot heatmap
        plt.figure(figsize=(12, 8))
        plt.imshow(tau_grid, extent=[lons_unique.min(), lons_unique.max(),
                                     lats_unique.min(), lats_unique.max()],
                  origin='lower', aspect='auto', cmap='plasma')
        plt.colorbar(label='Tempo di predicibilità τ')
        plt.xlabel('Longitudine')
        plt.ylabel('Latitudine')
        plt.title('Tempo di predicibilità - Heatmap')
        plt.tight_layout()
        plt.savefig('taupred_heatmap.png', dpi=150, bbox_inches='tight')
        print("Heatmap salvata: taupred_heatmap.png")
        
except Exception as e:
    print(f"Nota: Non è possibile creare heatmap (forse griglia irregolare): {e}")

print("\n=== Plot completati ===")
print("File generati:")
print("  - taupred_mappa.png (scatter plot base)")
print("  - taupred_istogramma.png (distribuzione)")
print("  - taupred_mappa_geografica.png (mappa con coste)")
print("  - taupred_heatmap.png (heatmap se griglia regolare)")
