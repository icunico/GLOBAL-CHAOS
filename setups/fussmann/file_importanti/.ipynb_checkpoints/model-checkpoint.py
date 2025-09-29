import os,sys
import numpy as np
import scipy.integrate
import pyfabm
import netCDF4 as nc
from scipy.fft import fft, fftfreq
from scipy.signal import find_peaks
import math
import subprocess
#import statsmodels.api as sm
from scipy.stats import variation
import sys
import matplotlib.pyplot as plt
#ensamble_counter=int(np.loadtxt('ensamble_counter.txt'))

#command='/g100_scratch/userexternal/gocchipi/FUSSMAN_5c/' + str(ensamble_counter).zfill(6) + '/'

fig,ax=plt.subplots()

dxvalue=np.linspace(0.1, 0.45, 70) 
dcvalue=np.linspace(0.1, 1.2, 70)

for i in range(len(dxvalue)):
	for j in range(len(dcvalue)):
		print(f"Running with dxvalue[{i}] = {dxvalue[i]}, dcvalue[{j}] = {dcvalue[j]}")			
		#subprocess.run(["mkdir","-p",command])
		#subprocess.run(["cp","fabm.yaml",command])
		# Create model (loads fabm.yaml)
		model = pyfabm.Model('fabm_fussmann.yaml')
		
		# Configure the environment
		# Note: the set of environmental dependencies depends on the loaded biogeochemical model.
		#model.dependencies['cell_thickness'].value = 1.
		#model.dependencies['temperature'].value = 15.
		#model.dependencies['practical_salinity'].value = 30.
		#model.dependencies['density'].value = 1000.
		#model.dependencies['depth'].value = 1.
		#model.dependencies['pressure'].value = 1.
		#model.dependencies['isBen'].value = 1.
		#model.dependencies['longitude'].value = 0.
		#model.dependencies['latitude'].value = 0.
		#model.dependencies['surface_downwelling_shortwave_flux'].value = 50.
		#model.dependencies['surface_air_pressure'].value = 1.
		#model.dependencies['wind_speed'].value = 5.
		#model.dependencies['mole_fraction_of_carbon_dioxide_in_air'].value = 390.
		#model.dependencies['number_of_days_since_start_of_the_year'].value = 1.
		
		
		dx=dxvalue[i]
		dc=dcvalue[j]
		
		model.parameters['X/dX1'].value=dx
		model.parameters['C1/dC1'].value=dc
		
		
		# Verify the model is ready to be used
		model.cell_thickness=1.
		
		assert model.checkReady(), 'One or more model dependencies have not been fulfilled.'
		
		# Time derivative
		def dy(t0, y):
		    model.state[:] = y
		    return model.getRates()

		t_ex= 3600*50000
		# Time-integrate over 1000 days (note: FABM's internal time unit is seconds!)
		t_eval = np.linspace(0, t_ex, 50000)
		sol = scipy.integrate.solve_ivp(dy, [0., t_ex], model.state, t_eval=t_eval)
		
		
		# Plot results
		#import pylab
		#pylab.plot(t, y)
		#pylab.legend([variable.path for variable in model.state_variables])
		#pylab.show()
		
		t = sol.t
		y = sol.y.T
		#t = t_eval
		
		specie_da_controllare= [0, 1, 2, 3, 5] #['P1', 'P2', 'P3', 'C1', 'C2', 'X1', 'Y1', 'Z1']
		flag_neg = False
		for iv in specie_da_controllare:
		    if np.any(y[:, iv] < 0.0):  
		        flag_neg = True

		y=y[-1500:]

		variable_names = [var.path for var in model.state_variables]

		# Seleziona solo le specie da controllare
		y_controllo = y[:, specie_da_controllare]
		threshold = 0.01
		# Calcola la media e la deviazione standard SOLO per le specie selezionate
		mean_y = np.mean(y_controllo, axis=0)
		std_y = np.std(y_controllo, axis=0)

		# Controllo del rapporto std/mean solo per le specie selezionate
		ratio = np.where(mean_y >= threshold, std_y / mean_y, 0)
		# Controllo della stabilità
		flag_ss = np.max(ratio) <= 0.01

		print("Ratio (std/mean):", mean_y)
		print("Ratio (std/std):", std_y)
		print("Ratio (std/ratio):", ratio)
		
		#check if a specie gets extinct
		flag_extin=False
		for  iv in specie_da_controllare: 
			if np.max(y[:, iv] < 0.001):
				flag_extin=True

		if flag_extin:
			ax.scatter(dc, dx, c='k')
		elif flag_neg:
			ax.scatter(dc, dx, c='k')
		elif flag_ss:
			ax.scatter(dc, dx, c='r')  
		else:
			ax.scatter(dc, dx, c='y')

fig.savefig('bubblef3.png')
#plt.savefig('bubblef.png')
plt.subplots(figsize=(15, 10))

Nt=t.shape[0]
deltaT=t[1]-t[0]
laststeps =1000# int(Nt/10) #compute the indicators just for this last timesteps
freq = (1/deltaT) * np.linspace(0,laststeps/2,int(laststeps/2)) / laststeps

## Save results


fileoutput = 'result.nc'
f = nc.Dataset(fileoutput, mode='w')

##lat_dim = f.createDimension('lat', 1)
##lon_dim = f.createDimension('lon', 1)
##dep_dim = f.createDimension('z', 1)
##time_dim = f.createDimension('time', Nt)

##lat = f.createVariable('lat', np.float32, ('lat',))
##lat.units = 'degrees_north'
##lat.long_name = 'latitude'
##f.variables['lat'][:]=0

##lon = f.createVariable('lon', np.float32, ('lon',))
##lon.units = 'degrees_east'
##lon.long_name = 'longitude'
##f.variables['lon'][:]=0

##time = f.createVariable('time', np.float64, ('time',))
##time.units = 'days'
##time.long_name = 'time'
##f.variables['time'][:]=t

##depth = f.createVariable('z', np.float32, ('z',))
##depth.units = 'meters'
##depth.long_name = 'depth'
##f.variables['z'][:]=1

##for v,variable in enumerate(model.state_variables):
##   ncvar = variable.name.replace("/","_")
##   var = f.createVariable(ncvar, np.float64, ('time', 'z','lat','lon'))
##   var.units = variable.units
##   var.long_name = variable.long_name
##   f.variables[ncvar][:]=y[:,v]

#GRAFICI
import matplotlib.pyplot as plt  # Usa matplotlib.pyplot invece di pylab

# Definisci le variabili da tracciare
#variables_to_track = ['Y/DWY1']  # Variabili di interesse
#variables_to_track = ['C1/DWC1', 'C1/DWC2', 'P1/DW3', 'P1/DW2', 'P1/DW1', 'X/DWX1', 'Y/DWY1', 'Z/DWZ1']  # Variabili di interesse
# Crea un dizionario per memorizzare gli indici delle variabili
##indices = {var: None for var in variables_to_track}

# Trova gli indici per le variabili
##for i, variable in enumerate(model.state_variables):
##    print(f"Variabile {variable.name} trovata all'indice {i}")  # Stampa tutte le variabili
##    if variable.name in indices:  # Controlla se la variabile è nella lista
##        indices[variable.name] = i  # Memorizza l'indice

# Verifica che tutte le variabili siano state trovate
##if None in indices.values():
##    print("Una o più variabili non sono state trovate!")
##else:
##    # Estrai i dati per tutte le variabili
##    linewidths = [7, 6, 5, 4, 3, 2, 1, 1, 1, 1, 1, 1]  # Puoi cambiare gli spessori delle linee

    # Crea una nuova figura per il grafico
##    fig, ax1 = plt.subplots(figsize=(10, 6))
##    for idx, var_name in enumerate(variables_to_track):
##        var_data = y[:, indices[var_name]]  # Estrai i dati dalla matrice y
##        ax1.plot(t, var_data, label=var_name, linewidth=linewidths[idx])  # Plotta con linewidth personalizzato

    # Aggiungi etichette, titolo e legenda
##    ax1.set_xlabel('Time (days)')
##    ax1.set_ylabel('Variables value')
##    ax1.set_title('Temporal evolution')
##    ax1.legend()
##   # ax1.set_xlim(left=500)
##    #ax1.set_yscale('log')

    # Migliora la disposizione del grafico
##    plt.tight_layout()

    # Mostra il grafico#
##    plt.show()

f.close()
