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
DATADIR='/m100_work/OGS21_PRACE_P_0/plazzari/Forward_Adjoint/BOUSSOLE_DATA/orig'
DATAFILE=DATADIR+'/'+'boussole_multi_rrs_bbpw_T10_IES20_2006-2012.csv'
dfB = pd.read_csv(DATAFILE, delimiter=";",skiprows = 0, engine='python')
bbp_442  = dfB['bbp_442'].values
bbw_442  = dfB['bbw_442'].values
bbp_488  = dfB['bbp_488'].values
bbw_488  = dfB['bbw_488'].values
bbp_550  = dfB['bbp_550'].values
bbw_550  = dfB['bbw_550'].values
bbp_620  = dfB['bbp_620'].values
bbw_620  = dfB['bbw_620'].values
yearsB = dfB['YEAR'].values
monthB = dfB['MONTH'].values
dayB   = dfB['DAY'].values

date_listB=[]
for i in range(len(bbp_442)):
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

abso  = pd.read_csv("abs_coeff.txt", sep="\s+",skiprows = 0, engine='python')
scat  = pd.read_csv("scat_coeff.txt", sep="\s+",skiprows = 0, engine='python')
bscat = pd.read_csv("back_scat_coeff.txt", sep="\s+",skiprows = 0, engine='python')
Nrows = abso.shape[0]
Nlambdas = 9 # maximum number of lambdas at BOUSSOLE

a_c =np.zeros((Nrows,7))
b_c =np.zeros((Nrows,7))
bb_c=np.zeros((Nrows,7))
count=0
for a_it,b_it,bb_it in zip(abso.iterrows(),scat.iterrows(),bscat.iterrows()):
    a_c[count,:]= [a_it[1][1]  , a_it[1][2]  , a_it[1][3]  , a_it[1][4]  , a_it[1][5]  , a_it[1][6]  , a_it[1][7]  ]
    b_c[count,:]= [b_it[1][1]  , b_it[1][2]  , b_it[1][3]  , b_it[1][4]  , b_it[1][5]  , b_it[1][6]  , b_it[1][7]  ]
    bb_c[count,:]=[bb_it[1][1] , bb_it[1][2] , bb_it[1][3] , bb_it[1][4] , bb_it[1][5] , bb_it[1][6] , bb_it[1][7] ]
    count +=1

bbw_m_442=np.zeros(len(date_list)) + bb_c[1,0]
bbw_m_490=np.zeros(len(date_list)) + bb_c[2,0]
bbw_m_555=np.zeros(len(date_list)) + bb_c[4,0]
bbw_m_665=np.zeros(len(date_list)) + bb_c[6,0]
p1=df_1["P1"].values 
p2=df_1["P2"].values
p3=df_1["P3"].values
p4=df_1["P4"].values
nap=df_1["NAP"].values 
bbp_m_442=3.0/0.4*( p1*bb_c[1,1]+p2*bb_c[1,2]+p3*bb_c[1,3]+p4*bb_c[1,4]+nap*bb_c[1,6])
bbp_m_490=3.0/0.4*( p1*bb_c[2,1]+p2*bb_c[2,2]+p3*bb_c[2,3]+p4*bb_c[2,4]+nap*bb_c[2,6])
bbp_m_555=3.0/0.4*( p1*bb_c[4,1]+p2*bb_c[4,2]+p3*bb_c[4,3]+p4*bb_c[4,4]+nap*bb_c[4,6])
bbp_m_665=3.0/0.4*( p1*bb_c[6,1]+p2*bb_c[6,2]+p3*bb_c[6,3]+p4*bb_c[6,4]+nap*bb_c[6,6])

#axs[0,0].scatter(date_list,y1,c='k',s=0.1,label='Total Chla - Model')

axs[0,0].scatter(date_listB,bbw_442,c='r',s=0.1,label='Data')
axs[0,0].scatter(date_list,bbw_m_442,c='k',s=0.1,label='Model')
axs[0,0].set_ylabel('$m^{-1}$')
axs[0,0].text(0.3, 1.1, '$Data=bbw_{442}, model=bbw_{442}$', horizontalalignment='center',
    verticalalignment='center', transform=axs[0,0].transAxes)
axs[0,0].legend(ncol=2, markerscale=12.)

axs[1,0].scatter(date_listB,bbw_488,c='r',s=0.1,label='Data')
axs[1,0].scatter(date_list,bbw_m_490,c='k',s=0.1,label='Model')
axs[1,0].text(0.3, 1.1, '$Data=bbw_{488}, model=bbw_{490}$', horizontalalignment='center',
    verticalalignment='center', transform=axs[1,0].transAxes)
axs[1,0].set_ylabel('$m^{-1}$')
#axs[1,0].legend(ncol=2, markerscale=12.)

axs[2,0].scatter(date_listB,bbw_550,c='r',s=0.1,label='Data')
axs[2,0].scatter(date_list,bbw_m_555,c='k',s=0.1,label='Model')
axs[2,0].text(0.3, 1.1, '$Data=bbw_{550}, model=bbw_{555}$', horizontalalignment='center',
    verticalalignment='center', transform=axs[2,0].transAxes)
axs[2,0].set_ylabel('$m^{-1}$')
#axs[2,0].legend(ncol=2, markerscale=12.)

axs[3,0].scatter(date_listB,bbw_620,c='r',s=0.1,label='Data$')
axs[3,0].scatter(date_list,bbw_m_665,c='k',s=0.1,label='Model')
axs[3,0].text(0.3, 1.1, '$Data=bbw_{620}, model=bbw_{665}$', horizontalalignment='center',
    verticalalignment='center', transform=axs[3,0].transAxes)
axs[3,0].set_ylabel('$m^{-1}$')
#axs[3,0].legend(ncol=2, markerscale=12.)

axs[0,1].scatter(date_listB,bbp_442,c='r',s=0.1,label='Data')
axs[0,1].scatter(date_list,bbp_m_442,c='k',s=0.1,label='Model')
axs[0,1].text(0.3, 1.1, '$Data=bbp_{442}, model=bbp_{442}$', horizontalalignment='center',
    verticalalignment='center', transform=axs[0,1].transAxes)
axs[0,1].set_ylabel('$m^{-1}$')
#axs[0,1].legend(ncol=2, markerscale=12.)

axs[1,1].scatter(date_listB,bbp_488,c='r',s=0.1,label='Data')
axs[1,1].scatter(date_list,bbp_m_490,c='k',s=0.1,label='Model')
axs[1,1].text(0.3, 1.1, '$Data=bbp_{488}, model=bbp_{490}$', horizontalalignment='center',
    verticalalignment='center', transform=axs[1,1].transAxes)
axs[1,1].set_ylabel('$m^{-1}$')
#axs[1,1].legend(ncol=2, markerscale=12.)

axs[2,1].scatter(date_listB,bbp_550,c='r',s=0.1,label='Data')
axs[2,1].scatter(date_list,bbp_m_555,c='k',s=0.1,label='Model')
axs[2,1].text(0.3, 1.1, '$Data=bbp_{550}, model=bbp_{555}$', horizontalalignment='center',
    verticalalignment='center', transform=axs[2,1].transAxes)
axs[2,1].set_ylabel('$m^{-1}$')
#axs[2,1].legend(ncol=2, markerscale=12.)

axs[3,1].scatter(date_listB,bbp_620,c='r',s=0.1,label='Data')
axs[3,1].scatter(date_list,bbp_m_665,c='k',s=0.1,label='Model')
axs[3,1].text(0.3, 1.1, '$Data=bbp_{620}, model=bbp_{665}$', horizontalalignment='center',
    verticalalignment='center', transform=axs[3,1].transAxes)
axs[3,1].set_ylabel('$m^{-1}$')
#axs[3,1].legend(ncol=2, markerscale=12.)

#axsT.scatter(date_list,ywl,c='g',s=0.1,label='WL')
#axs[0].legend(ncol=2,markerscale=12.)
#axs[1].scatter(date_list,y1_1,c='k',s=0.1,label='diatoms')
#axs[1].scatter(date_list,y1_2,c='m',s=0.1,label='flagellates')
#axs[1].scatter(date_list,y1_3,c='c',s=0.1,label='picophytoplankton')
#axs[1].scatter(date_listB,Tchla,c='r',s=0.1,label='Total Chla - Data')
#axs[1].legend(ncol=4, markerscale=12.)
#axs[1].scatter(date_list,y1_2,c='g',s=0.1)
#axs[1].scatter(date_list,y1_3,c='c',s=0.1)
#axs[2].scatter(date_list,y2,c='k',s=0.1,label='NAP')
#axs[2].legend(markerscale=12.)
#axs[3].scatter(date_list,y3,c='k',s=0.1,label='CDOM')
#axs[3].legend(markerscale=12.)

months = mdates.MonthLocator(interval=1)
months_fmt = mdates.DateFormatter('%b')
for j in range(2):
    for i in range(4):
        axs[i,j].xaxis.set_major_locator(months)
        axs[i,j].xaxis.set_major_formatter(months_fmt)
        axs[i,j].xaxis.set_tick_params(rotation=45)
        axs[i,j].set_xlim([datetime.date(2000, 1, 1), datetime.date(2001, 1, 1)])
#axs[0].set_ylim([0., 5.])
#plt.show()
plt.savefig("inversion_bb_seasonal_L-BFGS-B.png")

print("EOB")
