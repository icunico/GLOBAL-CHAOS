import os,sys
import glob
import datetime
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


DIRPATH="/m100_work/OGS21_PRACE_P_0/plazzari/Forward_Adjoint/BOUSSOLE_DATA/2POINTFIT/"
fig,axs = plt.subplots(9,1, gridspec_kw = {'wspace':.5, 'hspace':.8})
fig.set_size_inches(10,10)

kd_list=['412.5', '442.5', '490', '510','555','560', '665','670', '681.25']
kd_list_m= ['m412.50','m442.50','m490.00','m510.00','m555.00','m560.00','m665.00','m670.00','m681.25']
for kd_idx in range(9):
# BOUSSOLE DATA
    file_boussole=DIRPATH + "Kd."+ kd_list[kd_idx] + "BOUSSOLEFit_AntoineMethod.txt"
    dfB = pd.read_csv(file_boussole, delimiter="\t",skiprows = 0, engine='python')
    Kd_b   = dfB['Kd'].values
    yearsB = dfB['yyyy'].values
    monthB = dfB['mm'].values
    dayB   = dfB['dd'].values

    date_listB=[]
    for i in range(len(Kd_b)):
        datestr= str(int(2000)).zfill(4) + str(int(monthB[i])).zfill(2) + str(int(dayB[i])).zfill(2)
        dateobj=datetime.datetime.strptime(datestr, "%Y%m%d")
        date_listB.append(dateobj)
  
    df = pd.read_csv("inversion_result_kd_L-BFGS-B.txt", delimiter="\s+",skiprows = 0, engine='python')

    df_1= df.loc[df['WL'] > 300.] # filter 2

    yyyy = df_1['yyyy'].values
    mm = df_1['mm'].values
    dd = df_1['dd'].values
    HH   = df_1['HH'].values
    MM   = df_1['MM'].values
    SS   = df_1['SS'].values
    date_list=[]
    for tt,myyear in enumerate(yyyy):
        datestr= '2000' + '-' + str(mm[tt]) + '-' + str(dd[tt]) + ' ' + str(HH[tt]) + '-' + str(MM[tt]) + '-' + str(SS[tt])
        dateobj=datetime.datetime.strptime(datestr, "%Y-%m-%d %H-%M-%S")
        date_list.append(dateobj)

    y1=df_1[kd_list_m[kd_idx]].values 

    axs[kd_idx].scatter(date_listB,Kd_b,c='r',s=0.1,label='Kd - Data',alpha=0.1)
    axs[kd_idx].scatter(date_list,y1,c='k',s=0.1,label='Kd - Model')
    annotation=r'$\lambda=$ ' + kd_list[kd_idx]
    axs[kd_idx].text(0.8, 0.8, annotation, horizontalalignment='center',
     verticalalignment='center', transform=axs[kd_idx].transAxes)

months = mdates.MonthLocator(interval=1)
months_fmt = mdates.DateFormatter('%b')
for p in range(9):
    axs[p].xaxis.set_major_locator(months)
    axs[p].xaxis.set_major_formatter(months_fmt)
    axs[p].xaxis.set_tick_params(rotation=45)
    axs[p].set_xlim([datetime.date(2000, 1, 1), datetime.date(2001, 1, 1)])
axs[0].set_ylim([0., 0.3])
#plt.show()
plt.savefig("inversion_kd_seasonal_L-BFGS-B.png")

print("EOB")
