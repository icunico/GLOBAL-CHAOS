import numpy as np
import scipy.integrate
import pyfabm
import netCDF4 as nc
import matplotlib.pyplot as plt

# Create model (loads fabm.yaml)
model = pyfabm.Model('fabm_baczoo.yaml')

# Configure the environment
# Note: the set of environmental dependencies depends on the loaded biogeochemical model.
model.dependencies['cell_thickness'].value = 1.
model.dependencies['temperature'].value = 15.
#model.dependencies['temperature'].value = 15.
model.dependencies['practical_salinity'].value = 30.
model.dependencies['density'].value = 1000.
model.dependencies['depth'].value = 1.
model.dependencies['pressure'].value = 1.
model.dependencies['longitude'].value = 0.
model.dependencies['latitude'].value = 0.
model.dependencies['surface_downwelling_shortwave_flux'].value = 50.
#model.dependencies['surface_downwelling_shortwave_flux'].value = 50.
model.dependencies['surface_air_pressure'].value = 1.
model.dependencies['wind_speed'].value = 5.
model.dependencies['mole_fraction_of_carbon_dioxide_in_air'].value = 390.
model.dependencies['number_of_days_since_start_of_the_year'].value = 1.

# Verify the model is ready to be used
model.cell_thickness=1.

p_sum=[0.02, 0.2, 0.6]
for i,value in enumerate(p_sum):
	model.parameters['B1/p_sum'].value=value

	assert model.checkReady(), 'One or more model dependencies have not been fulfilled.'
	
	# Time derivative
	def dy(t0, y):
	    model.state[:] = y
	    return model.getRates()
	
	SEC_PER_DAY=86400
	NYEARS=100
	t_eval = np.linspace(0, SEC_PER_DAY*365*NYEARS, 100000)
	#t_eval = np.linspace(0, SEC_PER_DAY*365*NYEARS, 300000)
	#sol = scipy.integrate.solve_ivp(dy, [0., SEC_PER_DAY*365*NYEARS], model.state,
	#t_eval=t_eval)
	sol = scipy.integrate.solve_ivp(dy, [0., SEC_PER_DAY*365*NYEARS], model.state,
	t_eval=t_eval, max_step=50000)
	
	
	# Plot results
	#import pylab
	#pylab.plot(t, y)
	#pylab.legend([variable.path for variable in model.state_variables])
	#pylab.show()
	
	t = sol.t/86400
	y = sol.y.T
	
	
	#GRAFICI
	import matplotlib.pyplot as plt  # Usa matplotlib.pyplot invece di pylab
	
	# Definisci le variabili da tracciare
	
	#variables_to_track = [ 'R6/c', 'R8/c', 'R3/c', 'R1/c', 'R2/c']  # Variabili di interesse
	variables_to_track = ['B1/c', 'P1/c', 'P3/c']  # Variabili di interesse
	# Crea un dizionario per memorizzare gli indici delle variabili
	indices = {var: None for var in variables_to_track}
	
	# Trova gli indici per le variabili
	for i, variable in enumerate(model.state_variables):
	    print(f"Variabile {variable.name} trovata all'indice {i}")  # Stampa tutte le variabili
	    if variable.name in indices:  # Controlla se la variabile è nella lista
	        indices[variable.name] = i  # Memorizza l'indice
	
	# Verifica che tutte le variabili siano state trovate
	if None in indices.values():
	    print("Una o più variabili non sono state trovate!")
	else:
	    # Estrai i dati per tutte le variabili
	    linewidths = [7, 6, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1]  # Puoi cambiare gli spessori delle linee
	    
	
	    sum_data = np.zeros_like(t)  # Crea un array di zeri della stessa lunghezza di t
	
	    # Crea una nuova figura per il grafico
	    fig, ax1 = plt.subplots(figsize=(10, 6))
	    for idx, var_name in enumerate(variables_to_track):
	        var_data = y[:, indices[var_name]]  # Estrai i dati dalla matrice y
	        ax1.plot(t, var_data, label=var_name, linewidth=linewidths[idx])  # Plotta con linewidth personalizzato
	
	    # Aggiungi etichette, titolo e legenda
	    ax1.set_xlabel('Time (days)')
	    ax1.set_ylabel('Variables value')
	    ax1.set_title('Temporal evolution')
	    ax1.legend()
	  #  ax1.set_xlim(left=500)
	  #  ax1.set_yscale('log')
	    
	    # Migliora la disposizione del grafico
	    plt.tight_layout()
	
	    # Mostra il grafico
	    plt.show()
	
	    #GRAFICO CON SCALE DIVERSE
	
	# Plot della prima variabile sull'asse y primario
	#    var_data1 = y[:, indices['R1/c']]  # Estrai i dati per R1/c
	#    ax1.plot(t, var_data1, label='R1/c', linewidth=linewidths[0], color='tab:blue')
	#    ax1.set_xlabel('Time (t)')
	#    ax1.set_ylabel('R1/c', color='tab:blue')
	#    ax1.tick_params(axis='y', labelcolor='tab:blue')
	
	    # Crea il secondo asse y
	#    ax2 = ax1.twinx()
	#    var_data2 = y[:, indices['B1/p']]  # Estrai i dati per B1/p
	#    ax2.plot(t, var_data2, label='B1/p', linewidth=linewidths[1], color='tab:red')
	#    ax2.set_ylabel('B1/p', color='tab:red')
	#    ax2.tick_params(axis='y', labelcolor='tab:red')
	
	    # Aggiungi titolo e legenda
	#    ax1.set_title( 'R1/c e  B1/p')
	#    ax1.legend(loc='upper left')
	#    ax2.legend(loc='upper right')
	
	#    plt.tight_layout()  # Migliora la disposizione del grafico
	#    plt.show()
	
	Nt=t.shape[0]
	deltaT=t[1]-t[0]
	laststeps = int(Nt/55) #compute the indicators just for this last timesteps
	freq = (1/deltaT) * np.linspace(0,laststeps/2,int(laststeps/2))  
	#laststeps
	
	# Save results
	
	fileoutput = 'result_'+str(value)+'.nc'
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
	
