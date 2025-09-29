import os,sys
import glob
from datetime import datetime
import matplotlib.dates as mdates
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def create_init_params(lam):
    df_abw   = pd.read_csv("bcs/abw25_boundaries.dat", delimiter="\s+",skiprows = 5,header=None, engine='python')
    lam_bin  = df_abw.iloc[:,0]
    lam_low  = df_abw.iloc[:,1]
    lam_hig  = df_abw.iloc[:,2]
    idx=10000
    for i in range(len(lam_low)):
#       print(lam)
#       print(lam_low[i])
#       print(lam_hig[i])
        if (lam >= lam_low[i]) and (lam < lam_hig[i]):
            idx=i
#           print(idx)

    aw       = df_abw.iloc[idx,3]
    bw       = df_abw.iloc[idx,4]
    bbw      = df_abw.iloc[idx,5]
    return aw, bw, bbw

#########################
### MAIN CODE 

#### Create output file
result_file="inversion_result.txt"
df_res  = pd.read_csv(result_file, delimiter="\t",skiprows = 0, engine='python')
#profile = data.loc[(data['time'] == hhmm) & (data[var] > -998.0)]

WL_list=[412.5,442.5,490.,510.,555.,560.,665.,670.,681.25]

fig,axs = plt.subplots(9,3, gridspec_kw = {'wspace':.5, 'hspace':.5})
fig.set_size_inches(10,10)
hivernal_lines=[]
estival_lines=[]
for yyyy in range(2004,2013):
    hivernal_lines.append(datetime(yyyy, 1, 1))
    estival_lines.append(datetime(yyyy, 7, 1))
for i,WL in enumerate(WL_list):

    ag,bg,bbg = create_init_params(WL)

    df_filter_WL=df_res.loc[(df_res['WL'] == WL)]
    df_filter_depth=df_filter_WL.loc[(df_filter_WL['Depth_Up'] == 0.)]
    yyyymmdd = df_filter_depth['yyyymmdd'].values
    HHMMSS   = df_filter_depth['HHMMSS'].values
    date_list=[]
    for tt,mydate in enumerate(yyyymmdd):
        datestr= mydate + ' ' + HHMMSS[tt]
        dateobj=datetime.strptime(datestr, "%Y-%m-%d %H-%M-%S")
        date_list.append(dateobj)

    a        = df_filter_depth['a'].values
    b        = df_filter_depth['b'].values
    bb       = df_filter_depth['bb'].values
    rrs      = df_filter_depth['rrs'].values

    index_row=i
### First Column "a"
    ax=axs[index_row, 0]
    axT=ax.twinx()
    ax.scatter(date_list,a,s=0.8,marker='*',color='k')
#   axT.scatter(date_list,rrs,s=0.8,marker='*',color='r',alpha=0.5)
    ax.axhline(y=ag,color='g')
    for vline in hivernal_lines:
        ax.axvline(x=vline,color='c',linestyle='--',alpha=0.5)
    for vline in estival_lines:
        ax.axvline(x=vline,color='m',linestyle='--',alpha=0.5)
    ax.plot(date_list,a,marker='*',color='k')
    if index_row < 8:
        ax.set_xticklabels([])
    else:
        myFmt = mdates.DateFormatter('%Y')
        ax.xaxis.set_major_formatter(myFmt)
        ax.xaxis.set_tick_params(rotation=45)
    ax.set_ylim([0, None])
    ax.set_ylabel(r'$a_{%.2f}$' %WL)
### Second Column "b"
    ax=axs[index_row, 1]
    axT=ax.twinx()
    ax.scatter(date_list,b,s=0.8,marker='*',color='k')
#   axT.scatter(date_list,rrs,s=0.8,marker='*',color='r',alpha=0.5)
    ax.axhline(y=bg,color='g')
    for vline in hivernal_lines:
        ax.axvline(x=vline,color='c',linestyle='--',alpha=0.5)
    for vline in estival_lines:
        ax.axvline(x=vline,color='m',linestyle='--',alpha=0.5)
    if index_row < 8:
        ax.set_xticklabels([])
    else:
        myFmt = mdates.DateFormatter('%Y')
        ax.xaxis.set_major_formatter(myFmt)
        ax.xaxis.set_tick_params(rotation=45)
    ax.set_ylim([0, None])
    ax.set_ylabel(r'$b_{%.2f}$' %WL)
### Third Column "bb"
    ax=axs[index_row, 2]
    axT=ax.twinx()
    ax.scatter(date_list,bb, s=0.8, marker='*',color='k')
#   axT.scatter(date_list,rrs,s=0.8,marker='*',color='r',alpha=0.5)
    ax.axhline(y=bbg,color='g')
    for vline in hivernal_lines:
        ax.axvline(x=vline,color='c',linestyle='--',alpha=0.5)
    for vline in estival_lines:
        ax.axvline(x=vline,color='m',linestyle='--',alpha=0.5)
    if index_row < 8:
        ax.set_xticklabels([])
    else:
        myFmt = mdates.DateFormatter('%Y')
        ax.xaxis.set_major_formatter(myFmt)
        ax.xaxis.set_tick_params(rotation=45)
    ax.set_ylim([0, None])
    ax.set_ylabel(r'$bb_{%.2f}$' %WL)

fileout='inversion.png'
fig.savefig(fileout, format='png',dpi=150)

#fileout='inversion.eps'
#fig.savefig(fileout, format='eps')

print("EOB")
