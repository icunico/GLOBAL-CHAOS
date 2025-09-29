import numpy as np
import scipy.integrate
import pyfabm
import netCDF4 as nc
import matplotlib.pyplot as plt

# Create model (loads fabm.yaml)
model = pyfabm.Model('fabm_nop.yaml')

# Configure the environment
# Note: the set of environmental dependencies depends on the loaded biogeochemical model.
model.dependencies['cell_thickness'].value = 1.
model.dependencies['temperature'].value = 15.
model.dependencies['practical_salinity'].value = 30.
model.dependencies['density'].value = 1000.
model.dependencies['depth'].value = 1.
model.dependencies['pressure'].value = 1.
model.dependencies['longitude'].value = 0.
model.dependencies['latitude'].value = 0.
model.dependencies['surface_downwelling_shortwave_flux'].value = 50.
model.dependencies['surface_air_pressure'].value = 1.
model.dependencies['wind_speed'].value = 5.
model.dependencies['mole_fraction_of_carbon_dioxide_in_air'].value = 390.
model.dependencies['number_of_days_since_start_of_the_year'].value = 1.

# Verify the model is ready to be used
model.cell_thickness=1.

assert model.checkReady(), 'One or more model dependencies have not been fulfilled.'

# Time derivative
def dy(t0, y):
    model.state[:] = y
    return model.getRates()

SEC_PER_DAY=86400
NYEARS=10
t_eval = np.linspace(0, SEC_PER_DAY*365*NYEARS, 100000)
#t_eval = np.linspace(0, SEC_PER_DAY*365*NYEARS, 300000)
sol = scipy.integrate.solve_ivp(dy, [0., SEC_PER_DAY*365*NYEARS], model.state,
t_eval=t_eval, max_step=1000)


# Plot results
#import pylab
#pylab.plot(t, y)
#pylab.legend([variable.path for variable in model.state_variables])
#pylab.show()

t = sol.t/86400
y = sol.y.T


#GRAFICI
import matplotlib.pyplot as plt  # Usa matplotlib.pyplot per il plotting
import numpy as np  # Assicurati che numpy sia importato

# Definisci le variabili da tracciare
variables_to_track = ['O3/c', 'O5/c', 'B1/c', 'R1/c','R2/c','R3/c', 'X1/c', 'X2/c', 'X3/c', 'R6/c','R8/c', 'P1/c', 'P2/c', 'P3/c', 'P4/c', 'Z3/c', 'Z4/c', 'Z5/c', 
'Z6/c']

# Crea un dizionario per memorizzare gli indici delle variabili
indices = {var: None for var in variables_to_track}

# Trova gli indici per le variabili
for i, variable in enumerate(model.state_variables):
    if variable.name in indices:  # Controlla se la variabile è nella lista
        indices[variable.name] = i  # Memorizza l'indice

# Verifica che tutte le variabili siano state trovate
if None in indices.values():
    print("Una o più variabili non sono state trovate!")
else:
    # Somma dei dati per tutte le variabili
    sum_data = np.zeros_like(t)  # Crea un array vuoto per la somma

    # Somma i dati per ogni variabile
    for var_name in variables_to_track:
        if indices[var_name] is not None:  # Verifica che l'indice non sia None
            var_data = y[:, indices[var_name]]  # Estrai i dati dalla matrice y
            sum_data += var_data  # Somma i dati

    # Crea una nuova figura per il grafico
    plt.figure(figsize=(10, 6))

    # Plot della somma dei dati
    plt.plot(t, sum_data, label='Carbon mass', color='black', linewidth=2)

    # Aggiungi etichette, titolo e legenda
    plt.xlabel('Time (days)')
    plt.ylabel('Sum of Variables')
    plt.title('Temporal Evolution of Sum of Selected Variables')
    plt.legend()

    # Migliora la disposizione del grafico
    plt.tight_layout()

    # Mostra il grafico
    plt.show()

#Nt=t.shape[0]
#deltaT=t[1]-t[0]
#laststeps = int(Nt/55) #compute the indicators just for this last timesteps
#freq = (1/deltaT) * np.linspace(0,laststeps/2,int(laststeps/2))  
#laststeps

# Save results

fileoutput = 'result.nc'
f = nc.Dataset(fileoutput, mode='w')

lat_dim = f.createDimension('lat', 1)
lon_dim = f.createDimension('lon', 1)
dep_dim = f.createDimension('z', 1)
time_dim = f.createDimension('time', Nt)

lat = f.createVariable('lat', np.float32, ('lat',))
lat.units = 'degrees_north'
lat.long_name = 'latitude'
f.variables['lat'][:]=0

lon = f.createVariable('lon', np.float32, ('lon',))
lon.units = 'degrees_east'
lon.long_name = 'longitude'
f.variables['lon'][:]=0

time = f.createVariable('time', np.float64, ('time',))
time.units = 'days'
time.long_name = 'time'
f.variables['time'][:]=t

depth = f.createVariable('z', np.float32, ('z',))
depth.units = 'meters'
depth.long_name = 'depth'
f.variables['z'][:]=1

for v,variable in enumerate(model.state_variables):
   ncvar = variable.name.replace("/","_")
   var = f.createVariable(ncvar, np.float64, ('time', 'z','lat','lon'))
   var.units = variable.units
   var.long_name = variable.long_name
   f.variables[ncvar][:]=y[:,v]


f.close()

