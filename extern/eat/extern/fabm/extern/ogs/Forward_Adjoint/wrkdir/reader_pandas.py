import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import datetime
from datetime import timedelta

def myint(xa,xb,ya,yb,xc):
    yc = ya*(xb-xc)/(xb-xa) + yb*(xa-xc)/(xa-xb)
    return yc

def put_on_array(date_start,df,look_up_table,N,var):

    result=np.zeros(N)

    data=df[var].values

    yyyy_boussole=df['YEAR'].values
    mm_boussole=df['MONTH'].values
    dd_boussole=df['DAY'].values
    HH_boussole=df['GMTHOUR'].values
    MM_boussole=df['MINUTE'].values

    for tt,dat in enumerate(data):
        yyyy=yyyy_boussole[tt]
        mm=mm_boussole[tt]-1
        dd=dd_boussole[tt]-1
        HH=HH_boussole[tt]-1
        MM=MM_boussole[tt]
        date_boussole=datetime.datetime(yyyy,mm+1,dd+1,HH+1,MM,0)
        if date_boussole > date_start:
          yyyy_rel = yyyy-date_start.year
          idx=look_up_table[yyyy_rel,mm,dd,HH,MM]
          result[idx]=data[tt]

    return result

def put_on_array_OASIM(date_start,df,look_up_table,N,var):

    #yyyymmdd	HHMMSS	L0250	L0325
    #2004-01-01	00:07:30

    result=np.zeros(N)

    data=df[var].values

    yyyymmdd=df['yyyymmdd'].values
    HHMMSS=df['HHMMSS'].values

    for tt,dat in enumerate(data):
        yyyy=int(yyyymmdd[tt][0:4])
        mm  =int(yyyymmdd[tt][5:7])-1
        dd  =int(yyyymmdd[tt][8:10])-1
        HH  =int(HHMMSS[tt][0:2])-1
        MM  =int(HHMMSS[tt][3:5])
        date_boussole=datetime.datetime(yyyy,mm+1,dd+1,HH+1,MM,0)

        if date_boussole > date_start:
          yyyy_rel = yyyy-date_start.year
          idx=look_up_table[yyyy_rel,mm,dd,HH,MM]
          result[idx]=data[tt]

    return result

def write_line(fid,wl ,value_rrs, value_Ed, value_Es,solar_zenith):
    if ~ np.isnan(value_rrs):
        fid.write(wl)
        fid.write(" ")
        fid.write(str(value_rrs))
        fid.write(" ")
        fid.write(str(value_Ed))
        fid.write(" ")
        fid.write(str(value_Es))
        fid.write(" ")
        fid.write(str(solar_zenith))
        fid.write("\n")


###MAIN CODE

## input file header
#YEAR    MONTH   DAY     GMTHOUR      MINUTE      SS      date    sunzen  tangent Wt      sal     gamma   cond    depth.ctd       At      SST     pressure..atm   hygrometry      winddir windspeed       windspeed.dailymean     gust    WaveHeight      WavePeriod      ePAR    ed.0p.412.5     ed.0p.442.5     ed.0p.490       ed.0p.510       ed.0p.555       ed.0p.560       ed.0p.665       ed.0p.670       ed.0p.681.25    ed.412.5        ed.442.5        ed.490  ed.510  ed.555  ed.560  ed.665  ed.670  ed.681.25       i.es.412.5      i.es.442.5      i.es.490        i.es.510        i.es.555        i.es.560        i.es.665        i.es.670        i.es.681.25     depth.down4     ed4.412.5       ed4.442.5       ed4.490 ed4.510 ed4.555 ed4.560 ed4.665 ed4.670 ed4.681.25      depth.down9     ed9.412.5       ed9.442.5       ed9.490 ed9.510 ed9.555 ed9.560 ed9.665 ed9.670 ed9.681.25      lw.412.5        lw.442.5        lw.490  lw.510  lw.555  lw.560  lw.665  lw.670  lw.681.25       chl

df_ALL = pd.read_csv("../BOUSSOLE_DATA/orig/buoy.DPFF.2003-09-06_2012-12-31.dat", sep="\t",skiprows = 0, engine='python')
df_RRS = pd.read_csv("../BOUSSOLE_DATA/orig/boussole_multi_rrs_bbpw_T10_IES20_2006-2012.csv", sep=";",skiprows = 0, engine='python')
df_Ed_OASIM = pd.read_csv("../OASIM_DATA/Ed.txt", sep="\t",skiprows = 0, engine='python')
df_Es_OASIM = pd.read_csv("../OASIM_DATA/Es.txt", sep="\t",skiprows = 0, engine='python')

date_start=datetime.datetime(2004, 1, 1, 0, 0, 0)
date___end=datetime.datetime(2013, 1, 1, 0, 0, 0)

lista_date=[]
current_date = date_start
yyyy0        = date_start.year
print(yyyy0)
yyyye        = date___end.year
print(yyyye)
N_years      = yyyye - yyyy0 
print(N_years)
idx_counter  = 0

look_up_table=np.zeros((N_years,12,31,24,60),dtype=int)
look_up_table[:,:,:,:,:] = 0

step = datetime.timedelta(seconds=int(15*60))

while current_date < date___end:

    lista_date.append(current_date)
    print(current_date)

    yyyy=int(current_date.year)
    mm=int(current_date.month)-1
    dd=int(current_date.day)-1
    HH=int(current_date.hour)-1
    MM=int(current_date.minute)

    yyyy_rel=yyyy-yyyy0
    look_up_table[yyyy_rel,mm,dd,HH,MM]=idx_counter

    current_date=current_date + step
    idx_counter=idx_counter+1

N_rows=idx_counter

## Data Header
WAVELENGHTS=['412.5','442.5','490','510','555','560','665','670','681.25']
# csv var names rrs_412.5;rrs_442.5;rrs_490;rrs_510;rrs_555;rrs_560;rrs_665;rrs_670;rrs_681.25
r412_5   = put_on_array(date_start,df_RRS,look_up_table,N_rows,'rrs_412.5')
r442_5   = put_on_array(date_start,df_RRS,look_up_table,N_rows,'rrs_442.5')
r490     = put_on_array(date_start,df_RRS,look_up_table,N_rows,'rrs_490')
r510     = put_on_array(date_start,df_RRS,look_up_table,N_rows,'rrs_510')
r555     = put_on_array(date_start,df_RRS,look_up_table,N_rows,'rrs_555')
r560     = put_on_array(date_start,df_RRS,look_up_table,N_rows,'rrs_560')
r665     = put_on_array(date_start,df_RRS,look_up_table,N_rows,'rrs_665')
r670     = put_on_array(date_start,df_RRS,look_up_table,N_rows,'rrs_670')
r681_25  = put_on_array(date_start,df_RRS,look_up_table,N_rows,'rrs_681.25')

#OASIM_DATA
Ed_400   = put_on_array_OASIM(date_start,df_Ed_OASIM,look_up_table,N_rows,'L0400')
Ed_425   = put_on_array_OASIM(date_start,df_Ed_OASIM,look_up_table,N_rows,'L0425')
Ed_450   = put_on_array_OASIM(date_start,df_Ed_OASIM,look_up_table,N_rows,'L0450')
Ed_475   = put_on_array_OASIM(date_start,df_Ed_OASIM,look_up_table,N_rows,'L0475')
Ed_500   = put_on_array_OASIM(date_start,df_Ed_OASIM,look_up_table,N_rows,'L0500')
Ed_525   = put_on_array_OASIM(date_start,df_Ed_OASIM,look_up_table,N_rows,'L0525')
Ed_550   = put_on_array_OASIM(date_start,df_Ed_OASIM,look_up_table,N_rows,'L0550')
Ed_575   = put_on_array_OASIM(date_start,df_Ed_OASIM,look_up_table,N_rows,'L0575')
Ed_600   = put_on_array_OASIM(date_start,df_Ed_OASIM,look_up_table,N_rows,'L0600')
Ed_625   = put_on_array_OASIM(date_start,df_Ed_OASIM,look_up_table,N_rows,'L0625')
Ed_650   = put_on_array_OASIM(date_start,df_Ed_OASIM,look_up_table,N_rows,'L0650')
Ed_675   = put_on_array_OASIM(date_start,df_Ed_OASIM,look_up_table,N_rows,'L0675')
Ed_700   = put_on_array_OASIM(date_start,df_Ed_OASIM,look_up_table,N_rows,'L0700')

Es_400   = put_on_array_OASIM(date_start,df_Es_OASIM,look_up_table,N_rows,'L0400')
Es_425   = put_on_array_OASIM(date_start,df_Es_OASIM,look_up_table,N_rows,'L0425')
Es_450   = put_on_array_OASIM(date_start,df_Es_OASIM,look_up_table,N_rows,'L0450')
Es_475   = put_on_array_OASIM(date_start,df_Es_OASIM,look_up_table,N_rows,'L0475')
Es_500   = put_on_array_OASIM(date_start,df_Es_OASIM,look_up_table,N_rows,'L0500')
Es_525   = put_on_array_OASIM(date_start,df_Es_OASIM,look_up_table,N_rows,'L0525')
Es_550   = put_on_array_OASIM(date_start,df_Es_OASIM,look_up_table,N_rows,'L0550')
Es_575   = put_on_array_OASIM(date_start,df_Es_OASIM,look_up_table,N_rows,'L0575')
Es_600   = put_on_array_OASIM(date_start,df_Es_OASIM,look_up_table,N_rows,'L0600')
Es_625   = put_on_array_OASIM(date_start,df_Es_OASIM,look_up_table,N_rows,'L0625')
Es_650   = put_on_array_OASIM(date_start,df_Es_OASIM,look_up_table,N_rows,'L0650')
Es_675   = put_on_array_OASIM(date_start,df_Es_OASIM,look_up_table,N_rows,'L0675')
Es_700   = put_on_array_OASIM(date_start,df_Es_OASIM,look_up_table,N_rows,'L0700')

sunzen   = put_on_array(date_start,df_ALL,look_up_table,N_rows,'sunzen')

#INTERPOLATE OASIM
#'681.25']
#412.5
Ed_412_5  =  myint(400.,425., Ed_400 , Ed_425, 412.5)
Es_412_5  =  myint(400.,425., Es_400 , Es_425, 412.5)
#442.5
Ed_442_5  =  myint(425.,450., Ed_425 , Ed_450, 442.5)
Es_442_5  =  myint(425.,450., Es_425 , Es_450, 442.5)
#490
Ed_490  =  myint(475.,500., Ed_475 , Ed_500, 490.)
Es_490  =  myint(475.,500., Es_475 , Es_500, 490.)
#510
Ed_510  =  myint(500.,525., Ed_500 , Ed_525, 510.)
Es_510  =  myint(500.,525., Es_500 , Es_525, 510.)
#555
Ed_555  =  myint(550.,575., Ed_550 , Ed_575, 555.)
Es_555  =  myint(550.,575., Es_550 , Es_575, 555.)
#560
Ed_560  =  myint(550.,575., Ed_550 , Ed_575, 560.)
Es_560  =  myint(550.,575., Es_550 , Es_575, 560.)
#665
Ed_665  =  myint(650.,675., Ed_650 , Ed_675, 665.)
Es_665  =  myint(650.,675., Es_650 , Es_675, 665.)
#670
Ed_670  =  myint(650.,675., Ed_650 , Ed_675, 670.)
Es_670  =  myint(650.,675., Es_650 , Es_675, 670.)
#681.25
Ed_681_25  =  myint(675., 700.,Ed_675 , Ed_700, 681.25)
Es_681_25  =  myint(675., 700.,Es_675 , Es_700, 681.25)

#WRITE OUTPUT DATA

dir_output='SURFACE_DATA'

for n in range(N_rows):

    sec=int(lista_date[n].strftime("%H"))*3600 + int(lista_date[n].strftime("%M"))*60

    if ( sec > 10*3600 ) and ( sec < 14*3600 ):  

        if np.nanmax((r412_5[n],r442_5[n],r490[n],r510[n],r555[n],r560[n],r665[n],r670[n],r681_25[n])) > 0.0:

            mydate=lista_date[n].strftime("%Y-%m-%d_%H-%M-%S")
            file_out = dir_output + '/surface.' + mydate + '.txt'
        
            fid = open(file_out,'w')
        
            write_line(fid, WAVELENGHTS[0], r412_5[n], Ed_412_5[n], Es_412_5[n],  sunzen[n])
            write_line(fid, WAVELENGHTS[1], r442_5[n], Ed_442_5[n], Es_442_5[n],  sunzen[n])
            write_line(fid, WAVELENGHTS[2],   r490[n], Ed_490[n],   Es_490[n],    sunzen[n])
            write_line(fid, WAVELENGHTS[3],   r510[n], Ed_510[n],   Es_510[n],    sunzen[n])
            write_line(fid, WAVELENGHTS[4],   r555[n], Ed_555[n],   Es_555[n],    sunzen[n])
            write_line(fid, WAVELENGHTS[5],   r560[n], Ed_560[n],   Es_560[n],    sunzen[n])
            write_line(fid, WAVELENGHTS[6],   r665[n], Ed_665[n],   Es_665[n],    sunzen[n])
            write_line(fid, WAVELENGHTS[7],   r670[n], Ed_670[n],   Es_670[n],    sunzen[n])
            write_line(fid, WAVELENGHTS[8],r681_25[n], Ed_681_25[n],Es_681_25[n], sunzen[n])
    
            fid.close()

