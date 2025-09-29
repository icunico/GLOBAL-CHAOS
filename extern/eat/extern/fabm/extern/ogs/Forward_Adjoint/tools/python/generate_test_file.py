import numpy as np

nw=33

wavelenght=[ 250, 325, 350, 375, 400, 425, 450, 475, 500, 525, 550, 575, 600, 
             625, 650, 675, 700, 725, 775, 850, 950, 1050, 1150, 1250, 1350,
            1450, 1550, 1650, 1750, 1900, 2200, 2900, 3700]

dt = np.dtype([('wl', np.int), ('Rrs', np.float), ('Ed0m', np.float), ('Es0m', np.float)])

indata=np.zeros((nw),dtype=dt)

for i in range(nw):

    indata[i]['wl']=wavelenght[i] #wavelength

    indata[i]['Rrs'] =0.05        #Rrs0p_sat

    indata[i]['Ed0m']=0.35        #Ed0mOASIM

    indata[i]['Es0m']=0.35        #Es0mOASIM

fileout="surfdata.txt"

fmt = '%d', '%1.4f', '%1.4f', '%1.4f'

np.savetxt(fileout,indata,fmt=fmt)
