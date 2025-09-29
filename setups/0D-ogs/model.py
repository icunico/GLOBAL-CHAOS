import numpy as np
import scipy.integrate
import pyfabm
import netCDF4 as nc

# Create model (loads fabm.yaml)
model = pyfabm.Model('fabm.yaml')

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

# Time-integrate over 1000 days (note: FABM's internal time unit is seconds!)
t_eval = np.linspace(0, 3650.*86400, 10000) 
#t_eval = np.linspace(0, 3650.*86400, 300000) 
sol = scipy.integrate.solve_ivp(dy, [0., 3650.*86400], model.state, t_eval=t_eval)
#y = scipy.integrate.odeint(dy, model.state, t*86400)

# Plot results
#import pylab
#pylab.plot(t, y)
#pylab.legend([variable.path for variable in model.state_variables])
#pylab.show()

t = sol.t/86400
y = sol.y.T


Nt=t.shape[0]
deltaT=t[1]-t[0]
laststeps = int(Nt/10) #compute the indicators just for this last timesteps
freq = (1/deltaT) * np.linspace(0,laststeps/2,int(laststeps/2)) / laststeps

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
