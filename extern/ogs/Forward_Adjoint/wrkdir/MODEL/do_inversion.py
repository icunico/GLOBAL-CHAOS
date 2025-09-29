import os,sys
import glob
from datetime import datetime
import subprocess
import numpy as np
import pandas as pd

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

    print("Writing init param file")
    param_file="init_params.txt"
    fid_p = open(param_file,'w')
    fid_p.write(str(aw))
    fid_p.write(" ")
    fid_p.write(str(bw))
    fid_p.write("  ")
    fid_p.write(str(bbw))
    fid_p.close()

def create_init_opt_const():
    # 7 optically active constituents
    #Water, P1, P2, P3, P4, CDOM, NAP
    W   =1.0
    P1  =0.5
    P2  =0.5
    P3  =0.5
    P4  =0.5
    CDOM=0.5
    NAP =0.05
    param_file="init_opt_const.txt"
    fid_p = open(param_file,'w')
    header="W       P1      P2      P3      P4      CDOM    NAP    "
    data  ="1.00000 0.50000 0.50000 0.50000 0.10000 0.50000 0.05000"
    data  ="{:1.5f}".format(W)    + ' ' + \
           "{:1.5f}".format(P1)   + ' ' + \
           "{:1.5f}".format(P2)   + ' ' + \
           "{:1.5f}".format(P3)   + ' ' + \
           "{:1.5f}".format(P4)   + ' ' + \
           "{:1.5f}".format(CDOM) + ' ' + \
           "{:1.5f}".format(NAP ) 
           
    fid_p.write(header)
    fid_p.write("\n")
    fid_p.write(data)
    fid_p.close()


def write_row(fid,input_file):
    #                           012345678901234567890123456789
#    input_file="../SURFACE_DATA/surface.2006-02-03_13-30-00.txt"

    local_file="surfdata.txt"
    subprocess.run(["rm -f " + local_file], shell=True)
    subprocess.run(["cp  " + input_file + " " + local_file], shell=True)
    
    date_time_str=os.path.basename(input_file)[8:27]
    
    date_time_obj = datetime.strptime(date_time_str, '%Y-%m-%d_%H-%M-%S')
    
    with open(input_file) as f:
        lines = [line.rstrip() for line in f]
    
    nwl = 0 # number of wavelengths
    for line in lines:
        nwl += 1
#       input_str=line.split(" ")
    
    nwl_s   = str(nwl)

    ## clean existing files
    try:
        res_bio_file = "res_bio.txt"
        subprocess.run(['rm ' + res_bio_file], shell=True)
        res_opt_file = "res_opt*.txt"
        subprocess.run(['rm ' + res_opt_file], shell=True)
    except:
        print("File already clean!")

    
#    command=["./adj.xx 7"]
    command=["./adj.xx", nwl_s]
    
    print(command)
    subprocess.Popen(command) # creates results file

    return
    
#       if os.path.isfile(result_file): 

#           with open(result_file) as g:
#              rlines = [rline.rstrip() for rline in g]
        
#           WFUNC=rlines[0].split()
#           depths=rlines[1].split()
#           NDEPTH=len(depths)
#           NLAYERS=NDEPTH-1
#           for k in range(NLAYERS):
#               coeff=rlines[k+2].split()
#               print(coeff)
#               a =coeff[0]
#               b =coeff[1]
#               bb=coeff[2]
#               fid.write(date_time_obj.strftime("%Y-%m-%d")) 
#               fid.write("\t")
#               fid.write(date_time_obj.strftime("%H-%M-%S")) 
#               fid.write("\t")
#               fid.write(depths[k])
#               fid.write("\t")
#               fid.write(depths[k+1])
#               fid.write("\t")
#               fid.write(wl)
#               fid.write("\t")
#               fid.write('{0: >#016.5f}'. format(float(rrs)))
#               fid.write("\t")
#               fid.write('{0: >#016.2f}'. format(float(ed)))
#               fid.write("\t")
#               fid.write('{0: >#016.2f}'. format(float(es)))
#               fid.write("\t")
#               fid.write('{0: >#016.1f}'. format(float(sunz)))
#               fid.write("\t")
#               fid.write(a)
#               fid.write("\t")
#               fid.write(b)
#               fid.write("\t")
#               fid.write(bb)
#               fid.write("\t")
#               fid.write('{:.2e}'. format(float(WFUNC[0])))
#               fid.write("\n")


#########################
### MAIN CODE 

create_init_opt_const()

file_list=[]

for input_file in glob.iglob('../SURFACE_DATA/surface.*.txt'):
#for input_file in glob.iglob('../SURFACE_DATA/surface.2008-07-25_11-30-00.txt'):
    file_list.append(input_file)

file_list.sort()


#### Create output file
output_file="inversion_result.txt"

fid = open(output_file,'w')
header='yyyymmdd\tHHMMSS\tDepth_Up\tDepth_Down\tWL\trrs\ted\tes\tsunz\ta\tb\tbb\tERROR_NORM'
fid.write(header)
fid.write("\n")


#for input_file in file_list[0:2]:
for input_file in file_list:
    print(input_file)
    write_row(fid,input_file)
    sys.exit()

fid.close()    
print("EOB")
