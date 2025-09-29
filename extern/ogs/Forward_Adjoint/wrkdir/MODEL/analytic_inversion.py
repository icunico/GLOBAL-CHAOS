import argparse
def argument():
    parser = argparse.ArgumentParser(description = '''
    Generic averager for sat files.
    It works with one mesh, without performing interpolations.
    Files with dates used for the average provided (dirdates).
    ''',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    parser.add_argument(   '--file_bio', '-fb',
                                type = str,
                                required = False,
                                default='inversion_result_bio_simple.txt',
                                help = ''' tab separeted file containing biological results'''

                                )

    parser.add_argument(   '--file_rrs', '-fr',
                                type = str,
                                required = False,
                                default='inversion_result_rrs_simple.txt',
                                help = ''' tab separeted file containing reflectance results'''

                                )

    parser.add_argument(   '--file_kd', '-fk',
                                type = str,
                                required = False,
                                default='inversion_result_kd_simple.txt',
                                help = ''' tab separeted file containing downward attenuation results'''

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
from datetime import datetime
import numpy as np
from scipy.optimize import minimize
import pandas as pd
import matplotlib.pyplot as plt
from argparse import ArgumentParser

def getSurface(flnm):
    dt = pd.read_csv(flnm, delimiter="\s+",skiprows = 0,header=None, engine='python')
    wl_list=[412.5,442.5, 490., 510., 555., 560., 665.,  670.,  681.25]
    RF_a   =[0.003,0.004,0.011,0.013,0.017,0.017,0.018,0.018,  0.018]
    RF_b1  =[0.014,0.015,0.010,0.010,0.010,0.010,0.0010,0.010,  0.010]
    RF_b2  =[-0.022,-0.023,-0.051,-0.060,-0.080,-0.080,-0.081,-0.081, -0.081]
    Nlambdas=len(wl_list)
    SURF_DATA=np.zeros((Nlambdas,5))

    lambdas=[False, False, False, False, False, False, False, False, False]

    for j,ww in enumerate(wl_list):
        for i in range(dt.shape[0]):
            if dt[0][i] == ww:
                lambdas[j] = True
    c=0
    for j,ww in enumerate(wl_list):
        if lambdas[j]:
            SURF_DATA[j,:]=[dt[0][c],dt[1][c],dt[2][c],dt[3][c],dt[4][c]]
            c +=1

    return c, lambdas,SURF_DATA

def getmud(sunz):
    refrac_idx=1.341
    rad    = 180.0/np.arccos(-1.0) # radians
    rsza = sunz/rad
    sinszaw = np.sin(rsza)/refrac_idx
    szaw = np.arcsin(sinszaw)
    avgcos=np.cos(szaw) 
    return avgcos

def G(IOP):
    Eu_m=Eu_model(IOP)
    error = 0.
    for i in range(Nrows):
      if lambdas[i]:
        error += (Eu_m[i]- Eu_0sat[i])**2.0
    return error

def Eu_model(IOP):

    res=np.zeros((Nrows))

    for i in range(Nrows):
      if lambdas[i]:

#        FACT=[0.3,0.3,0.3,0.3,1.0,1.0]
         FACT=[3.,3.,3.,3.,1.,1.]
         a=a_c[i,0]   + np.dot(a_c[i,1:7],FACT*IOP)
         b=b_c[i,0]   + np.dot(b_c[i,1:7],FACT*IOP)
         bb=bb_c[i,0] + np.dot(bb_c[i,1:7],FACT*IOP)
         Ad = (a+b)/vd
         Fd = (b-rd*bb)/vd
         Bd = rd*bb/vd
         Cs = (a+rs*bb)/vs
         Cu = (a+ru*bb)/vu
         Bs = rs*bb/vs
         Bu = ru*bb/vu
#        !eigenvalues:  Ad = Ad,
         DD = 0.5*(Cs+Cu+np.sqrt((Cs+Cu)*(Cs+Cu)-4.0*Bs*Bu))
         kp = DD-Cu
         km = DD-Cs
#        !eigenvectors
         rkp = Bs/DD
         rkm = Bu/DD
#        Inhomogenous solution
         Mi=np.array([[Ad+Cu,-Bu],[-Bs,Ad-Cs]])
         Vi=np.array([-Fd,Bd])
         Ni = ((Ad-Cs)*(Ad+Cu) +Bs*Bu )
         if abs(Ni) > eps :
             Ni=((Ad-Cs)*(Ad+Cu) +Bs*Bu ) # normalization
             [x,y]=1.0/Ni*np.matmul(Mi,Vi)
             Eu_0=(Es_0[i] - x*Ed_0[i])*rkp + y*Ed_0[i]
         else :
             [x,y]=np.array([Fd,-Bd])
             Eu_0=(Es_0[i])*rkp 
         res[i] = Eu_0

    return res
def EdEsEu_model(z,IOP):

    Ed_z=np.zeros((Nrows))
    Es_z=np.zeros((Nrows))
    Eu_z=np.zeros((Nrows))

    for i in range(Nrows):
      if lambdas[i]:

         FACT=[3.,3.,3.,3.,1.,1.]
         a=a_c[i,0]   + np.dot(a_c[i,1:7],FACT*IOP)
         b=b_c[i,0]   + np.dot(b_c[i,1:7],FACT*IOP)
         bb=bb_c[i,0] + np.dot(bb_c[i,1:7],FACT*IOP)
         Ad = (a+b)/vd
         Fd = (b-rd*bb)/vd
         Bd = rd*bb/vd
         Cs = (a+rs*bb)/vs
         Cu = (a+ru*bb)/vu
         Bs = rs*bb/vs
         Bu = ru*bb/vu
#        !eigenvalues:  Ad = Ad,
         DD = 0.5*(Cs+Cu+np.sqrt((Cs+Cu)*(Cs+Cu)-4.0*Bs*Bu))
         kp = DD-Cu
         km = DD-Cs
#        !eigenvectors
         rkp = Bs/DD
         rkm = Bu/DD
#        Inhomogenous solution
         Mi=np.array([[Ad+Cu,-Bu],[-Bs,Ad-Cs]])
         Vi=np.array([-Fd,Bd])
         Ni = ((Ad-Cs)*(Ad+Cu) +Bs*Bu )
         Ed_z[i]=Ed_0[i]*np.exp(-Ad*z)
         if abs(Ni) > eps :
             Ni=((Ad-Cs)*(Ad+Cu) +Bs*Bu ) # normalization
             [x,y]=1.0/Ni*np.matmul(Mi,Vi)
             Es_z[i]=(Es_0[i] - x*Ed_0[i])    *np.exp(-kp*z) + x*Ed_z[i]
             Eu_z[i]=(Es_0[i] - x*Ed_0[i])*rkp*np.exp(-kp*z) + y*Ed_z[i]
         else :
             [x,y]=np.array([Fd,-Bd])
             Es_z[i]= Es_0[i]     *np.exp(-kp*z) + x*z*Ed_z[i]
             Eu_z[i]=(Es_0[i])*rkp*np.exp(-kp*z) + y*z*Ed_z[i]

    return Ed_z,Es_z,Eu_z
def Kd_model(z,IOP):
    Kd_z=np.zeros((Nrows))
    Ed_z,Es_z,Eu_z=EdEsEu_model(z,IOP)
    for i in range(Nrows):
        if lambdas[i]:
            Kd_z[i] = -np.log( (Ed_z[i]+Es_z[i]) / (Ed_0[i]+Es_0[i]) ) / z
        else:
            Kd_z[i] = np.nan
    return Kd_z

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

#b_c[:,1:]=2.0*b_c[:,1:]
#bb_c[:,1:]=0.2*bb_c[:,1:]
#==> ../SURFACE_DATA/surface.2006-06-26_13-15-00.txt <==
Ed_0=np.zeros((Nlambdas))
Es_0=np.zeros((Nlambdas))


rd = 1.0
rs = 1.5
ru = 3.0
vs = 0.83
vu = 0.4

eps= 0.001

W   =1.0
P1  =0.1
P2  =0.1
P3  =0.1
P4  =0.05
CDOM=0.5
NAP =0.05

state=np.array((P1,P2,P3,P4,CDOM,NAP))
state_w=np.array((0.,0.,0.,0.,0.,0.))


#########################
### MAIN CODE


file_list=[]

for input_file in glob.iglob('../SURFACE_DATA/surface.*.txt'):
#for input_file in glob.iglob('../SURFACE_DATA/surface.2008-07-25_11-30-00.txt'):
    file_list.append(input_file)

file_list.sort()

fid_bio = open(args.file_bio,'w')
header='yyyy\tmm\tdd\tHH\tMM\tSS\tWL\tP1\tP2\tP3\tP4\tCDOM\tNAP'
fid_bio.write(header)
fid_bio.write("\n")

fid_rrs = open(args.file_rrs,'w')
header='yyyy\tmm\tdd\tHH\tMM\tSS\tWL\td412.50\td442.50\td490.00\td510.00\td555.00\td560.00\td665.00\td670.00\td681.25\tm412.50\tm442.50\tm490.00\tm510.00\tm555.00\tm560.00\tm665.00\tm670.00\tm681.25'
fid_rrs.write(header)
fid_rrs.write("\n")

fid_kd = open(args.file_kd,'w')
header='yyyy\tmm\tdd\tHH\tMM\tSS\tWL\tm412.50\tm442.50\tm490.00\tm510.00\tm555.00\tm560.00\tm665.00\tm670.00\tm681.25'
fid_kd.write(header)
fid_kd.write("\n")

for input_file in file_list:
    print("input_file -->" + input_file)
    date_time_str=os.path.basename(input_file)[8:27]
#   0000000000111111111
#   0123456789012345678   
#   yyyy-mm-dd_HH-MM-SS    
    HH=int(date_time_str[11:13])
    if HH<10 and HH>14:
        continue
    date_time_obj = datetime.strptime(date_time_str, '%Y-%m-%d_%H-%M-%S')


    c,lambdas,SURF_DATA =  getSurface(input_file)

    b=np.asarray(list(map(int, lambdas))) # map presence of lambda to int [0 or 1]
    WL=b.dot(2**np.arange(b.size)[::-1]) # convert bynary to integer

    if c < 7:
        continue

#   SURF_DATA[0,:]=[412.5 ,0.00550517   ,7.544209500000001  ,12.9893255         ,29.28330459]
#   SURF_DATA[1,:]=[442.5 ,0.00558697   ,10.012499400000001 ,15.458903          ,29.28330459]
#   SURF_DATA[2,:]=[490   ,0.00514099   ,11.576611799999998 ,15.5211058         ,29.28330459]
#   SURF_DATA[3,:]=[510   ,0.00388067   ,11.661843399999999 ,14.8509742         ,29.28330459]
#   SURF_DATA[5,:]=[560   ,0.00193718   ,11.911050199999998 ,13.812709599999998 ,29.28330459]
#   SURF_DATA[6,:]=[665   ,0.000162026  ,10.3260218         ,10.641849          ,29.28330459]
#   SURF_DATA[8,:]=[681.25,0.000149838  ,8.7379275          ,8.89418125         ,29.28330459]

    Rrs0p =SURF_DATA[:,1]

    wl_list=[412.5,442.5, 490., 510., 555., 560., 665.,  670.,  681.25]
    RF_a   =[0.003,0.004,0.011,0.013,0.017,0.017,0.018,0.018,  0.018]
    RF_b1  =[0.014,0.015,0.010,0.010,0.010,0.010,0.0010,0.010,  0.010]
    RF_b2  =[-0.022,-0.023,-0.051,-0.060,-0.080,-0.080,-0.081,-0.081, -0.081]
    RF=np.zeros(Nlambdas)

    if lambdas[1] and lambdas[4]:
        for ii in range(Nlambdas):
            RF[ii] = RF_a[ii]*(Rrs0p[1]/Rrs0p[4]) + RF_b1[ii] * np.power(Rrs0p[4],RF_b2[ii])
    else:
            RF[:]  = 0.0
    print(RF)

    Rrs0p_c = Rrs0p/(1.0 + RF) # Correction for Raman Scattering effect Lee et al., 2013

#derive Rrs0m using the correction by Lee et al. 2002
    T=0.52
    GammaQ=1.7
#   Rrs0m = Rrs0p
    Rrs0m = Rrs0p_c/(T+GammaQ*Rrs0p_c) # assume input is already Rrs0m
    Ed_0=SURF_DATA[:,2]
    Es_0=SURF_DATA[:,3]
    Zenith=SURF_DATA[0,4]
# Formula to retrieve QRrs Aas and Hojerslev 1999
# 90 -solz = solar elevation
    QRrs    = 5.33*np.exp(-0.45*np.sin(np.pi/180.*(90.0-Zenith)))
    Eu_0sat = QRrs * (Ed_0+Es_0)*Rrs0m
    vd = getmud(Zenith)

#   state=np.array((P1,P2,P3,P4,CDOM,NAP))
    bnds = ((0, 10), (0, 10), (0, 10), (0, 0.1),(0, 10), (0, 1))

    res = minimize(G, state, method=args.opt_method, bounds=bnds,tol=1.e-10) 

#   res = minimize(G, state, method='L-BFGS-B', bounds=bnds,tol=1.e-10) #interesting results
#   res = minimize(G, state, method='Nelder-Mead', bounds=bnds,tol=1.e-10) # gives negative results
#   res = minimize(G, state, method='SLSQP', bounds=bnds,tol=1.e-10)

    fid_bio.write(date_time_obj.strftime("%Y\t%m\t%d"))
    fid_bio.write("\t")
    fid_bio.write(date_time_obj.strftime("%H\t%M\t%S"))
    fid_bio.write("\t")
    fid_bio.write(str(WL))
    fid_bio.write("\t")
    fid_bio.write('{0: >#016.5f}'. format(float(res.x[0])))
    fid_bio.write("\t")
    fid_bio.write('{0: >#016.5f}'. format(float(res.x[1])))
    fid_bio.write("\t")
    fid_bio.write('{0: >#016.5f}'. format(float(res.x[2])))
    fid_bio.write("\t")
    fid_bio.write('{0: >#016.5f}'. format(float(res.x[3])))
    fid_bio.write("\t")
    fid_bio.write('{0: >#016.5f}'. format(float(res.x[4])))
    fid_bio.write("\t")
    fid_bio.write('{0: >#016.5f}'. format(float(res.x[5])))
    fid_bio.write("\n")

    fid_rrs.write(date_time_obj.strftime("%Y\t%m\t%d"))
    fid_rrs.write("\t")
    fid_rrs.write(date_time_obj.strftime("%H\t%M\t%S"))
    fid_rrs.write("\t")
    fid_rrs.write(str(WL))
    for kk in range(9):
        fid_rrs.write("\t")
        if lambdas[kk]:
           fid_rrs.write('{0: >#012.5f}'. format(float(Rrs0m[kk])))
        else:
           fid_rrs.write('{0: >#012.5f}'. format(float(np.nan)))
#   rrs_model=Eu_model(state_w)/(Ed_0+Es_0)/4.0
    rrs_model=Eu_model(res.x)/(Ed_0+Es_0)/4.0
    for kk in range(9):
        fid_rrs.write("\t")
        if lambdas[kk]:
           fid_rrs.write('{0: >#012.5f}'. format(float(rrs_model[kk])))
        else:
           fid_rrs.write('{0: >#012.5f}'. format(float(np.nan)))
    fid_rrs.write("\n")

    fid_kd.write(date_time_obj.strftime("%Y\t%m\t%d"))
    fid_kd.write("\t")
    fid_kd.write(date_time_obj.strftime("%H\t%M\t%S"))
    fid_kd.write("\t")
    fid_kd.write(str(WL))
    kd_m=Kd_model(9.,res.x)
    for kk in range(9):
        fid_kd.write("\t")
        if lambdas[kk]:
           fid_kd.write('{0: >#012.5f}'. format(float(kd_m[kk])))
        else:
           fid_kd.write('{0: >#012.5f}'. format(float(np.nan)))
    fid_kd.write("\n")

fid_bio.close()
fid_rrs.close()
fid_kd.close()
print("EOB")
