import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generic averager for sat files.
    It works with one mesh, without performing interpolations.
    Files with dates used for the average provided (dirdates).
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    parser.add_argument(   '--opt_method', '-m',
                                type = str,
                                required = False,
                                default='SLSQP',
                                choices = ['L-BFGS-B', 'TNC', 'SLSQP', 'Powell', 'trust-constr'],
                                help = ''' Optimization method choose between ['L-BFGS-B', 'TNC', 'SLSQP', 'Powell', 'trust-constr']'''

                                )
    return parser.parse_args()

args = argument()

import os,sys
import glob
import datetime
import matplotlib.dates as mdates
import subprocess
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


fig,axs = plt.subplots(12,1, gridspec_kw = {'wspace':.5, 'hspace':.5})
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
    datestr= str(int(yearsB[i])).zfill(4) + str(int(monthB[i])).zfill(2) + str(int(dayB[i])).zfill(2)
    dateobj=datetime.datetime.strptime(datestr, "%Y%m%d")
    date_listB.append(dateobj)
  
datavar  = ["d412.50","d442.50","d490.00","d510.00","d555.00","d560.00","d665.00","d670.00","d681.25"]
modelvar = ["m412.50","m442.50","m490.00","m510.00","m555.00","m560.00","m665.00","m670.00","m681.25"]
xlabels  = ["412.50","442.50","490.00","510.00","555.00","560.00","665.00","670.00","681.25"]
months   = ["JAN", "FEB", "MAR", "APR", "MAY", "JUN", "JUL", "AUG", "SEP", "OCT", "NOV", "DEC"]

#d412.50 d442.50 d490.0	d510.00	d555.00	d560.00	d665.00	d670.00	d681.25	m412.50	m442.50	m490.00	m510.00	m555.00	m560.00	m665.00	m670.00	m681.25
method=args.opt_method
infile="inversion_result_rrs_" + method + ".txt"
df = pd.read_csv(infile, delimiter="\s+",skiprows = 0, engine='python')

for mm in range(12):

    df_month = df.loc[df['mm'] ==  mm+1]

    rrsd_ave =[]
    rrsm_ave =[]
    rrsd_std =[]
    rrsm_std =[]

    for v in datavar:
        rrsd_ave.append(np.nanmean(df_month[v].values)) 
        rrsd_std.append(np.nanstd(df_month[v].values)) 

    for v in modelvar:
        rrsm_ave.append(np.nanmean(df_month[v].values)) 
        rrsm_std.append(np.nanstd(df_month[v].values)) 

    x=range(9)
    axs[mm].errorbar(x,rrsd_ave,yerr=rrsd_std,c='r',label='Data')
    axs[mm].errorbar(x,rrsm_ave,yerr=rrsm_std,c='k',label='Model')
    axs[mm].text(0.8, 0.6, months[mm], horizontalalignment='center',
     verticalalignment='center', transform=axs[mm].transAxes)
    if mm == 0:
        axs[mm].legend(ncol=2)
    if mm < 11:
        axs[mm].set_xticks(range(len(xlabels)))
        empty_string_labels = ['']*len(xlabels)
        axs[mm].set_xticklabels(empty_string_labels)
    if mm == 5:
        axs[mm].set_ylabel(r'$R_{rs} [st^{-1}]$')
    if mm == 11:
        axs[mm].set_xticks(range(len(xlabels)))
        axs[mm].set_xticklabels(xlabels)
        axs[mm].set_xlabel(r'$\lambda [nm]$')

plt.suptitle(r'$R_{rs}(0-)$')

outfile="inversion_rrs_"+ method + ".png"
plt.savefig(outfile)

print("EOB")
