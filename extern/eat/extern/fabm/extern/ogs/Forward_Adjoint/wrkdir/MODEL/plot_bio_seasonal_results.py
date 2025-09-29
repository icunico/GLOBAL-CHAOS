import os,sys
import glob
import datetime
import matplotlib.dates as mdates
import matplotlib.ticker as ticker
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


fig,axs = plt.subplots(4,1, gridspec_kw = {'wspace':.5, 'hspace':.5})
fig.set_size_inches(10,10)
# BOUSSOLE DATA
dfB = pd.read_csv("BOUSSOLE_PIGMENTS.csv", delimiter=",",skiprows = 0, engine='python')
dfB_surf = dfB.loc[dfB['Depth'] < 10]
Tchla  = dfB_surf['Tchla'].values
yearsB = dfB_surf['year'].values
monthB = dfB_surf['month'].values
dayB   = dfB_surf['day'].values

date_listB=[]
for i in range(len(Tchla)):
    datestr= str(int(2000)).zfill(4) + str(int(monthB[i])).zfill(2) + str(int(dayB[i])).zfill(2)
    dateobj=datetime.datetime.strptime(datestr, "%Y%m%d")
    date_listB.append(dateobj)
  
df = pd.read_csv("inversion_result_bio_L-BFGS-B.txt", delimiter="\s+",skiprows = 0, engine='python')
#df_1= df.loc[df['NAP'] < 1.9] # filter 1

df_0= df.loc[df['NAP'] < 0.9] # filter 1
df_1= df_0.loc[df_0['WL'] > 300.] # filter 2

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

y1=df_1["P1"].values + df_1["P2"].values + df_1["P3"].values +df_1["P4"].values 
ywl=df_1["WL"].values
y1_1=df_1["P1"].values 
y1_2=df_1["P2"].values 
y1_3=df_1["P3"].values 
y2=df_1["NAP"].values 
y3=df_1["CDOM"].values 

axs[0].scatter(date_list,y1,c='k',s=0.1,label='Total Chla - Model')
axs[0].scatter(date_listB,Tchla,c='r',s=0.1,label='Total Chla - Data')
axsT=axs[0].twinx()
axsT.scatter(date_list,ywl,c='g',s=0.1,label='WL')
axs[0].legend(ncol=2,markerscale=12.)
axs[1].scatter(date_list,y1_1,c='k',s=0.1,label='diatoms')
axs[1].scatter(date_list,y1_2,c='m',s=0.1,label='flagellates')
axs[1].scatter(date_list,y1_3,c='c',s=0.1,label='picophytoplankton')
axs[1].scatter(date_listB,Tchla,c='r',s=0.1,label='Total Chla - Data')
axs[1].legend(ncol=4, markerscale=12.)
#axs[1].scatter(date_list,y1_2,c='g',s=0.1)
#axs[1].scatter(date_list,y1_3,c='c',s=0.1)
axs[2].scatter(date_list,y2,c='k',s=0.1,label='NAP')
axs[2].legend(markerscale=12.)
axs[3].scatter(date_list,y3,c='k',s=0.1,label='CDOM')
axs[3].legend(markerscale=12.)

months = mdates.MonthLocator(interval=1)
months_fmt = mdates.DateFormatter('%b')
for p in range(4):
    axs[p].xaxis.set_major_locator(months)
    axs[p].xaxis.set_major_formatter(months_fmt)
    axs[p].xaxis.set_tick_params(rotation=45)
    axs[p].set_xlim([datetime.date(2000, 1, 1), datetime.date(2001, 1, 1)])
axs[0].set_ylim([0., 2.5])
#plt.show()
plt.savefig("inversion_bio_seasonal_L-BFGS-B.png")

print("EOB")
