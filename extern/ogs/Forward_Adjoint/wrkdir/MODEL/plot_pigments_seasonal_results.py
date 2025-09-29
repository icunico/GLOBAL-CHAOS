import os,sys
import glob
import datetime
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


fig,axs = plt.subplots(4,2, gridspec_kw = {'wspace':.5, 'hspace':.5})
fig.set_size_inches(10,10)
# BOUSSOLE DATA
DATADIR='.'
DATAFILE=DATADIR+'/'+'boussole_pigments_jul08.csv'
dfB = pd.read_csv(DATAFILE, delimiter=",",skiprows = 41, engine='python')
#Photoprotective carotenoids (PPC) are Alloxanthin, Lutein, Violaxanthin, Zeaxanthin, diadinoxanthin(DD), diatoxanthin(DT) and alpha-beta carotenes.
allo            = dfB['allo'].values
lut             = dfB['lut'].values
viola           = dfB['viola'].values
zea             = dfB['zea'].values
diadino         = dfB['diadino'].values
diato           = dfB['diato'].values
alpha_beta_car  = dfB['alpha-beta-car'].values
yearsB = dfB['year'].values
monthB = dfB['month'].values
dayB   = dfB['day'].values

date_listB=[]
for i in range(len(allo)):
    datestr= str(int(2000)).zfill(4) + str(int(monthB[i])).zfill(2) + str(int(dayB[i])).zfill(2)
    dateobj=datetime.datetime.strptime(datestr, "%Y%m%d")
    date_listB.append(dateobj)
  
#df = pd.read_csv("inversion_result_bio_L-BFGS-B.txt", delimiter="\s+",skiprows = 0, engine='python')

#df_0= df.loc[df['NAP'] < 0.9] # filter 1
#df_1= df_0.loc[df_0['WL'] > 300.] # filter 2

#yyyy = df_1['yyyy'].values
#mm = df_1['mm'].values
#dd = df_1['dd'].values
#HH   = df_1['HH'].values
#MM   = df_1['MM'].values
#SS   = df_1['SS'].values
#date_list=[]
#for tt,myyear in enumerate(yyyy):
#    datestr= '2000' + '-' + str(mm[tt]) + '-' + str(dd[tt]) + ' ' + str(HH[tt]) + '-' + str(MM[tt]) + '-' + str(SS[tt])
#    dateobj=datetime.datetime.strptime(datestr, "%Y-%m-%d %H-%M-%S")
#    date_list.append(dateobj)

#p1=df_1["P1"].values 
#p2=df_1["P2"].values
#p3=df_1["P3"].values
#p4=df_1["P4"].values
#nap=df_1["NAP"].values 

#axs[0,0].scatter(date_list,y1,c='k',s=0.1,label='Total Chla - Model')
allo            = dfB['allo'].values
lut             = dfB['lut'].values
viola           = dfB['viola'].values
zea             = dfB['zea'].values
diadino         = dfB['diadino'].values
diato           = dfB['diato'].values
alpha_beta_car  = dfB['alpha-beta-car'].values

axs[0,0].scatter(date_listB,allo,c='r',s=0.1,label='Alloxanthin')
axs[1,0].scatter(date_listB,lut,c='r',s=0.1,label='Lutein')
axs[2,0].scatter(date_listB,viola,c='r',s=0.1,label='Violaxanthin')
axs[3,0].scatter(date_listB,zea,c='r',s=0.1,label='Zeaxanthin')
axs[0,1].scatter(date_listB,diadino,c='r',s=0.1,label='Diadinoxanthin(DD)')
axs[1,1].scatter(date_listB,diato,c='r',s=0.1,label='Diatoxanthin(DT)')
axs[2,1].scatter(date_listB,alpha_beta_car,c='r',s=0.1,label='Alpha-beta carotenes')
#axs[0,0].scatter(date_list,bbw_m_442,c='k',s=0.1,label='Model')
#axs[0,0].set_ylabel('$m^{-1}$')
#axs[0,0].text(0.3, 1.1, '$Data=bbw_{442}, model=bbw_{442}$', horizontalalignment='center',
#    verticalalignment='center', transform=axs[0,0].transAxes)
#axs[0,0].legend(ncol=2, markerscale=12.)

#axs[1,0].scatter(date_listB,bbw_488,c='r',s=0.1,label='Data')
#axs[1,0].scatter(date_list,bbw_m_490,c='k',s=0.1,label='Model')
#axs[1,0].text(0.3, 1.1, '$Data=bbw_{488}, model=bbw_{490}$', horizontalalignment='center',
#    verticalalignment='center', transform=axs[1,0].transAxes)
#axs[1,0].set_ylabel('$m^{-1}$')
#axs[1,0].legend(ncol=2, markerscale=12.)

months = mdates.MonthLocator(interval=1)
months_fmt = mdates.DateFormatter('%b')
for j in range(2):
    for i in range(4):
        axs[i,j].xaxis.set_major_locator(months)
        axs[i,j].xaxis.set_major_formatter(months_fmt)
        axs[i,j].xaxis.set_tick_params(rotation=45)
        axs[i,j].set_xlim([datetime.date(2000, 1, 1), datetime.date(2001, 1, 1)])
        axs[i,j].set_ylabel('$mg.m^{-3}$')
        axs[i,j].legend(ncol=1, markerscale=12.)
#axs[0].set_ylim([0., 5.])
#plt.show()
plt.savefig("pigments_seasonal.png")

print("EOB")
