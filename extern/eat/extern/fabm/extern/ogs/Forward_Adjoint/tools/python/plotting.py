import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt


#iter, z(lvl), z(lvl+1), a(lvl), b(lvl), bb(lvl), vd(lvl), vs, vu, rd, rs, ru
dt = np.dtype([('iter', np.int), ('z_up', np.float), ('z_bot', np.float), ('a', np.float),
               ('b', np.float), ('bb', np.float), ('vd', np.float), ('vs', np.float),
               ('vu', np.float), ('rd', np.float), ('rs', np.float), ('ru', np.float)])

infile="/galileo/home/userexternal/plazzari/Forward_Adjoint/bin/sol_0250.txt"
indata=np.loadtxt(infile,dtype=dt)

Nite=np.max(indata['iter'])
Nlev=len(indata['iter'])/Nite
z_up  = np.reshape(indata['z_up'],(Nite,Nlev))
z_bot = np.reshape(indata['z_bot'],(Nite,Nlev))

fig,axs = plt.subplots(3,3, gridspec_kw = {'wspace':.5, 'hspace':.5})
fig.set_size_inches(10,10)
listavar=['a','b','bb']
z=0.5*(z_bot[-1,:]+z_up[-1,:])

for j,var in enumerate(listavar):

    var2plot = np.reshape(indata[var],(Nite,Nlev))

    index_row=int(j / 3)
    index_col=int(j % 3)
    ax=axs[index_row, index_col]
    
    ax.scatter(var2plot[-1,:],-z,marker='*')

fileout='inversion.png'
fig.savefig(fileout, format='png',dpi=150)

fileout='inversion.eps'
fig.savefig(fileout, format='eps')
   
