import netCDF4 as nc
import matplotlib.pyplot as plt

# Apri il file NetCDF
fileoutput = 'result.nc'
f = nc.Dataset(fileoutput, mode='r')

# Verifica se la variabile esiste nel file
if 'P3/c' in f.variables:
    # Estrai i dati della variabile 'flPIR2c'
    P3c = f.variables['P3/c'][:]
    
    # Estrai anche la dimensione del tempo (o altre dimensioni se necessario)
    time = f.variables['time'][:]  # Modifica se il nome della variabile tempo è diverso

    # Crea il grafico
    plt.plot(time, P3c)
    plt.xlabel('Tempo (giorni)')  # Modifica a seconda delle unità
    plt.ylabel('P3/c')
    plt.title('Evoluzione di P3c nel tempo')
    plt.grid(True)
    plt.show()

else:
    print("La variabile 'P3c' non esiste nel file.")

# Chiudi il file
f.close()

