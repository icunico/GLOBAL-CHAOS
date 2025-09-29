import os
import numpy as np
import re

nombrito = "fabm_"

# Set input and output directories
indir = os.getcwd() + "/"
outdir = os.getcwd() + "/"

#########################################################################
## COMPOSE & SAVE fabm.yaml with several size-classes per functional type
#########################################################################

# Define size-classes in an octave scale (2^n)
# Based on: "Dealing with size-spectra: some conceptual and mathematical problems" 
# JM Blanco, F Echevarría, CM García - Scientia Marina 1994
clases = np.arange(-8, 33)   # n
volumes = np.power(2.,clases)          # biovolume in um^3
diameters = ((volumes * 6) / np.pi)**(1/3)  # Equivalent spherical diameter (ESD)

# Define size range (bins) for each functional type:
clasesB1 = [-5]
clasesP6 = np.arange(-4, -2)  # Equivalent to -4:-3 in R
clasesP9 = np.arange(-2, 1)   # Equivalent to -2:0 in R
clasesP3 = np.arange(0, 3)    # Equivalent to 0:2 in R
clasesP7 = np.arange(2, 8)    # Equivalent to 2:7 in R
clasesP8 = np.arange(2, 8)    # Same as P7
clasesP2 = np.arange(4, 14)   # Equivalent to 4:13 in R
clasesP5 = np.arange(6, 16)   # Equivalent to 6:15 in R
clasesP1 = np.arange(11, 21)  # Equivalent to 11:20 in R
clasesP4 = np.arange(11, 21)  # Same as P1
clasesZ6 = np.arange(5, 12)   # Equivalent to 5:11 in R
clasesZ5 = np.arange(9, 19)   # Equivalent to 9:18 in R
clasesZ4 = np.arange(15, 25)  # Equivalent to 15:24 in R
clasesZ3 = np.arange(22, 32)  # Equivalent to 22:31 in R

# Combine all groups
all_groups = np.concatenate([
    clasesB1,
    clasesP1, clasesP2, clasesP3, clasesP4, clasesP5, clasesP6, clasesP7, clasesP8, clasesP9,
    clasesZ6, clasesZ5, clasesZ4, clasesZ3
])

# Histogram: Count how many types (all, phyto and zoo) fall in each bin
histo_all, bins_all = np.histogram(all_groups, bins=clases)
clases_init = bins_all[1:]
counts_init = histo_all

# Phytoplankton groups
phyto_groups = np.concatenate([
    clasesP1, clasesP2, clasesP3, clasesP4, clasesP5, clasesP6, clasesP7, clasesP8, clasesP9
])
histo_phyto, bins_phyto = np.histogram(phyto_groups, bins=clases)
clases_init_phyto = bins_phyto[1:]
counts_init_phyto = histo_phyto

# Zooplankton groups
zoo_groups = np.concatenate([clasesZ6, clasesZ5, clasesZ4, clasesZ3])
histo_zoo, bins_zoo = np.histogram(zoo_groups, bins=clases)
clases_init_zoo = bins_zoo[1:]
counts_init_zoo = histo_zoo

# To display the combined counts (for debugging or logging)
counts_combined = np.column_stack((counts_init, counts_init_phyto, counts_init_zoo))
print(counts_combined)


###############################################################################################
## Initialization biomass (fabm) and reference state (gotm)
## If reference state is given (in gotm.yaml for relaxation), the value in fabm.yaml is not used
## Several options tested below
###############################################################################################

# Initialize each size-class to the same biomass c=2 mg m^-3 (all phyto and zoo together)
# Note: B1 does not get the initialization value
c_init = 2

# Initialize each type to the same biomass
c_init_phyto = 0.5
c_init_zoo = 0.1

# Initialize each size-class (bin) to the same biomass, phyto and zoo separately
# Realistic initialization values (similar to Serra-Pompei et al. 2020, c=5 mg m^-3)
c_init_phyto = 1
c_init_zoo = 0.5

# Very small initialization values (Bruggeman & Kooijman 2007, c=0.0096 mg m^-3 per type)
c_init_phyto = 0.01
c_init_zoo = 0.005

# Intermediate initialization values (null net flow, on average, with migration in and out)
c_init_phyto = 0.1
c_init_zoo = 0.08

# Choose one
equal_initialization = False
separate_phytozoo = True

###############################################

## Read the fabm.yaml template that contains the 9PFTs and 4Zoo
## Several text tags have been included to locate the different "lego bricks"
template_file = os.path.join(indir, "fabm_multispectral_9PFTs_template.yaml")

# Read the lines from the template file
with open(template_file, 'r') as file:
    modulo_productor = file.readlines()

## Locate the common part: lightspectral, bethicLayer, initialization for pelagic_base, PelagicCSYS, PelOxygen, B1, PelChem, and CalciteDissolution
inicio = next(i for i, line in enumerate(modulo_productor) if "start_base" in line)
final = next(i for i, line in enumerate(modulo_productor) if "end_base" in line)

# START the fabm.yaml
output_file = os.path.join(outdir, f"{nombrito}.yaml")

# Write the common part
with open(output_file, 'w') as out_file:
    out_file.writelines(modulo_productor[inicio+1:final])

###########################################
### Create BRICKS for PRODUCERS
###########################################
    
#names and parameters for the 9 phytoplankton groups
old_names = ["P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P9"]
old_long_names = ["diatoms", "prymnesiophyta", "smalleukariotes", "dinoflagellates", "coccolithophores",
              "prochlorococcus", "chlorophyceae", "prasinophyceae", "synechococcus"]
clases_names = [clasesP1, clasesP2, clasesP3, clasesP4, clasesP5, clasesP6, clasesP7, clasesP8, clasesP9]
# Parameters to modify
parameters = ["p_sum", "p_srs", "p_res", "p_qlcPPY", "p_qup", "p_qun", "p_qplc", 
              "p_qpcPPY", "p_qnlc", "p_qncPPY", "p_quantum_yield"]
# Old parameter values
old_values = [
             [3.7, 0.06, 5.0, 0.039, 0.0046, 0.025, 0.00057, 0.000786, 0.00687, 0.0126, 0.504e-3],
             [2.6, 0.09, 0.0, 0.045, 0.0033, 0.025, 0.000352, 0.000556, 0.00687, 0.0126, 0.563e-3],
             [3.5, 0.1, 0.0, 0.010, 0.0036, 0.25, 0.0004288, 0.00072,  0.00687, 0.0126, 1.043e-3],
             [1.5, 0.1, 2.5, 0.015, 0.0038, 0.025, 0.0004288, 0.000786, 0.00687, 0.0126, 0.617e-3],
             [2.6, 0.09, 0.0, 0.045, 0.0034, 0.025, 0.000352, 0.000556,  0.00687, 0.0126, 0.447e-3],
             [3.5, 0.1, 0.0, 0.010, 0.0027, 0.25, 0.0004288, 0.00072, 0.00687, 0.0126, 1.489e-3],
             [2.6, 0.09, 0.0, 0.045, 0.0034, 0.025, 0.000352, 0.000556, 0.00687, 0.0126, 0.447e-3],
             [2.6, 0.09, 0.0, 0.040, 0.0038,  0.025,  0.000352, 0.000556, 0.00687, 0.0126, 0.591e-3],
             [3.5, 0.1, 0.0, 0.014, 0.0034, 0.25, 0.0004288, 0.00072, 0.00687, 0.0126, 1.570e-3]
             ]
# Number of digits for formatting
digito = [3, 3, 2, 3, 4, 4, 6, 6, 6, 6, 5]
#digito = [4, 4, 3, 4, 5, 5, 7, 7, 7, 7, 6]

#match Eva's order
#old_names_eva = ["P1","P4","P2","P5","P7","P8","P3","P9","P6"]
id_eva         = [0,3,1,4,6,7,2,8,5]
old_names      = np.array(old_names)[id_eva]
old_long_names = np.array(old_long_names)[id_eva]
#clases_names   = np.array(clases_names)[id_eva]
old_values     = np.array(old_values)[id_eva]
clases_names   = [clases_names[i] for i in id_eva]


for iPFT, old_name in enumerate(old_names):

    old_long_name = old_long_names[iPFT]

    clasesP = clases_names[iPFT]  # ESD = (((2^clasesP)*6)/pi)^(1/3)

    # Replace "-" with "m" in clasesP, the name of the instance cannot have "-"
    clasesN = [re.sub("-", "m", str(c)) for c in clasesP]

    old_value = old_values[iPFT]

    # Extract module block (template) for Px
    inicio = next(i for i, line in enumerate(modulo_productor) if f"start_{old_name}" in line)
    final = next(i for i, line in enumerate(modulo_productor) if f"end_{old_name}" in line)    

    # Calculate lower and upper limits for the size class
    lower_limit = 2. ** np.array(clasesP)
    upper_limit = 2. ** (np.array(clasesP) + 1)

    # Calculate diameters
    low_diameter = ((lower_limit * 6) / np.pi) ** (1/3)

    # Define initial values
    p_sum = 3.8 * np.power(lower_limit,(-0.17))  # Tang 1995 & Irwin 2006
    p_srs = 0.063 - 0.008 * np.log10(lower_limit)  # Shimoda 2016
    carbon = (10. ** (-0.933)) * np.power(lower_limit,0.881)  # Menden-Deuer and Lessard, 2000 (diatoms > 3000)
    QPmin = (10. ** (-1.04)) * np.power(lower_limit,0.714)   # Minimum P:C quota (Grover 1989), fmol cell-1
    p_qplc = (QPmin / carbon) * 1e9 * 1e-12  # mmolP mgC-1
    QPmax = (10. ** (-0.29)) * np.power(lower_limit,0.767)  # Maximum P:C quota (Grover 1989), fmol cell-1
    p_qpcPPY = ((QPmax / carbon) * 1e9 * 1e-12) / 2  # mmolP mgC-1
    QNmin = (1.36e-9) * np.power(lower_limit,0.77)    # Minimum N:C quota (Lichtman et al 2007), umol cell-1
    p_qnlc = (QNmin / carbon) * 1e9 * 1e-3  # mmolN mgC-1
    QNmax = (4.64e-9) * np.power(lower_limit, 0.81)    # Maximal N:C quota (Montagnes & Franklin 2001), umol cell-1, divided by luxury storage to get optimum quota
    p_qncPPY = ((QNmax / carbon) * 1e9 * 1e-3) / 2  # mmolN mgC-1
    p_res = 0.024 * np.power(lower_limit,0.37)  # Durante et al, 2019 (cylinders)
    aP = (10. ** (-8.1)) * np.power(lower_limit,0.73)  # P affinity (Edwards 2012)
    p_qup = (aP / carbon) * 1e-3 * 1e9 / 4  # m3 mgC-1 d-1
    aN = (10. ** (-8.2)) * np.power(lower_limit, 0.75)  # N affinity (Edwards 2012)
    p_qun = (aN / carbon) * 1e-3 * 1e9  # m3 mgC-1 d-1

    # Adjust P:C and N:C quotas
    p_qplc = p_qpcPPY * 0.55
    p_qncPPY = (p_qpcPPY * 0.0126) / 0.000786
    p_qnlc = p_qncPPY * 0.55

    # Default values for Chl-a:C and photochemical efficiency
    p_qlcPPY = np.full(len(lower_limit), 0.02)
    p_quantum_yield = np.full(len(lower_limit), 0.480e-3)

    # Loop over each size class (j)
    for j in range(len(clasesP)):
        modulo_modificado = np.copy(modulo_productor)

        # Calculate minimum diameter
        diameter_min = np.power((lower_limit[j] * 6) / np.pi,1/3)

        # Replace old name and long name in the module block
        modulo_modificado[inicio:final][2] =  modulo_modificado[inicio:final][2].split(':\n')[0]+f"_{clasesN[j]}:\n"
        modulo_modificado[inicio:final][3] =  modulo_modificado[inicio:final][3].split('\n')[0]+f"_{round(diameter_min, 1)}\n"
        
        # Replace parameter values
        for i, param in enumerate(parameters):
            # Locate the parameter in the block
            filas = [re.search(f"{param}:", line) for line in modulo_modificado[inicio:final]]
            fila = next((line for match, line in zip(filas, modulo_modificado[inicio:final]) if match), None)
            if fila:
                #get the index of the line
                cual = np.where(modulo_modificado[inicio:final] == fila)[0][0]
                # Replace old value with the new calculated one
                parametro = eval(param)[j]
                modulo_modificado[inicio:final][cual] = re.sub(str(old_value[i]), str(round(parametro, digito[i])), fila)

        # Compute initialization values
        if equal_initialization:
            carbon_value = c_init_phyto
        elif separate_phytozoo:
            carbon_value = c_init_phyto / counts_init_phyto[np.where(clases_init == clasesP[j])[0][0]]
        else:
            carbon_value = c_init / counts_init[np.where(clases_init == clasesP[j])[0][0]]

        nitro = carbon_value * (0.1008 / 8)
        phosp = carbon_value * (0.006288 / 8)
        silica = carbon_value * (0.08 / 8)
        chloro = carbon_value * (0.16 / 8)

        # Locate and replace initialization chunk
        location = [re.search("initialization:", line) for line in modulo_modificado[inicio:final]]
        cual = next((i for i, match in enumerate(location) if match), None)
        if cual is not None:
            modulo_modificado[inicio:final][cual+1] = re.sub("8.0", str(round(carbon_value, 11)), modulo_modificado[inicio:final][cual+1])
            modulo_modificado[inicio:final][cual+2] = re.sub("0.1008", str(round(nitro, 11)), modulo_modificado[inicio:final][cual+2])
            modulo_modificado[inicio:final][cual+3] = re.sub("0.006288", str(round(phosp, 11)), modulo_modificado[inicio:final][cual+3])
            modulo_modificado[inicio:final][cual+4] = re.sub("0.08", str(round(silica, 11)), modulo_modificado[inicio:final][cual+4])
            modulo_modificado[inicio:final][cual+5] = re.sub("0.16", str(round(chloro, 11)), modulo_modificado[inicio:final][cual+5])

        # Append the modified module to the fabm.yaml
        with open(f"{outdir}{nombrito}.yaml", "a") as outfile:
            outfile.writelines(modulo_modificado[inicio+1:final])


###########################################
### Create BRICKS for PREDATORS
###########################################
            
# Z6
# Variables setup
old_name = "Z6"
old_long_name = "Heterotrophic Nanoflagellates (HNAN)"
clasesP = clasesZ6  # ESD = (((2**clasesP)*6)/pi)**(1/3)

# Parameters to modify
parameters = ["p_sum", "p_chuc"]
old_value = ["3.88", "90"]
digito = [2, 0]
#digito = [3, 0]

old_nprey = 12
# Define all potential prey choices
preys = np.concatenate([clasesB1, clasesP2, clasesP3, clasesP5, clasesP6, clasesP7, clasesP8, clasesP9, clasesZ6])
names_preys = ["B1"] + [f"P2_{str(cl).replace('-', 'm')}" for cl in clasesP2] + \
              [f"P3_{str(cl).replace('-', 'm')}" for cl in clasesP3] + \
              [f"P5_{str(cl).replace('-', 'm')}" for cl in clasesP5] + \
              [f"P6_{str(cl).replace('-', 'm')}" for cl in clasesP6] + \
              [f"P7_{str(cl).replace('-', 'm')}" for cl in clasesP7] + \
              [f"P8_{str(cl).replace('-', 'm')}" for cl in clasesP8] + \
              [f"P9_{str(cl).replace('-', 'm')}" for cl in clasesP9] + \
              [f"Z6_{str(cl).replace('-', 'm')}" for cl in clasesZ6]
names_preys = np.array(names_preys)
# Extract module block for Zx
inicio = [i for i, line in enumerate(modulo_productor) if re.search(f"start_{old_name}", line)]
final = [i for i, line in enumerate(modulo_productor) if re.search(f"end_{old_name}", line)]

lower_limit = np.power(2.,clasesP)
upper_limit = np.power(2.,(clasesP + 1))
diameter = ((lower_limit * 6) / np.pi)**(1/3)

# Compute parameter values (allometry for taxon-specific)
p_sum = 26 * (diameter ** -0.4)  # Maximum ingestion rate (Hansen 1997)
p_chuc = np.full(len(diameter), 120)  # Banas 2011, default value for this model

# Optimum prey size
diameter_preys = ((np.power(2.,preys) * 6) / np.pi)**(1/3)
opt_prey_size = diameter / 10  # Kiorboe 2008, Ward 2012

# Loop over size classes of predators
for j in range(len(clasesP)):
    modulo_modificado = np.copy(modulo_productor)
    diameter_min = ((lower_limit[j] * 6) / np.pi) ** (1 / 3)
    # Replace module names
    modulo_modificado[inicio[0]+2] = modulo_modificado[inicio[0]+2].split(':\n')[0]+f"_{clasesP[j]}:\n"
    modulo_modificado[inicio[0]+4] = modulo_modificado[inicio[0]+4].split('\n')[0]+f"_{round(diameter_min, 1)}\n"
    
    #modulo_modificado[inicio[0]+2] = re.sub(old_name, f"{old_name}_{clasesP[j]}", modulo_modificado[inicio[0]+2])
    #modulo_modificado[inicio[0]+4] = re.sub(old_long_name, f"{old_long_name}_{round(diameter_min, 1)}", modulo_modificado[inicio[0]+4])
    
    # Replace parameter values
    for i, param in enumerate(parameters):
        param_value = locals()[param]  # dynamically retrieve variable
        filas = [idx for idx, line in enumerate(modulo_modificado[inicio[0]:final[0]]) if re.search(f"{param}:", line)]
        if filas:
            fila_idx = inicio[0] + filas[0]
            modulo_modificado[fila_idx] = re.sub(old_value[i], str(round(param_value[j], digito[i])), modulo_modificado[fila_idx])

    # Compute preferences, SD of gaussian = 0.20
    prefs = np.exp(-((np.log10(diameter_preys) - np.log10(opt_prey_size[j])) / 0.20)**2)
    prefs[prefs < 1e-2] = 0
    sorted_indices = np.argsort(np.power(2.,preys))
    
    # Replace number of preys (nprey)
    nprey_fila_idx = next(i for i in range(inicio[0], final[0]) if re.search("nprey", modulo_modificado[i]))
    modulo_modificado[nprey_fila_idx] = re.sub(str(old_nprey), str(len(names_preys[prefs != 0])), modulo_modificado[nprey_fila_idx])
    
    # Block of preferences
    pref_lines = [f"      suprey{idx+1}:   {round(pref, 2)}     # [-] Availability of {names_preys[idx]} to Z6_{clasesP[j]}"
                  for idx, pref in enumerate(prefs) if pref != 0]
    modulo_modificado[nprey_fila_idx+1:nprey_fila_idx+1+len(pref_lines)] = pref_lines
    
    # Block of calcifier
    boolean = np.zeros(len(names_preys[prefs != 0]), dtype=int)
    boolean[[i for i, name in enumerate(names_preys[prefs != 0]) if name.startswith("P5")]] = 1 #why P5?
    calcifier_lines = [f"      isP2{idx+1}:   {boolean[idx]}           # [-] identify P2 among the preys [1 for P2 otherwise 0]"
                       for idx in range(len(boolean))]
    calcifier_fila_idx = next(i for i in range(inicio[0], final[0]) if re.search("isP21:", modulo_modificado[i]))
    modulo_modificado[calcifier_fila_idx:calcifier_fila_idx+len(calcifier_lines)] = calcifier_lines
    
    # Block of coupling
    coupling_lines = [f"      prey{idx+1}:   {names_preys[idx]}" for idx, pref in enumerate(prefs) if pref != 0]
    prey_fila_idx = next(i for i in range(inicio[0], final[0]) if re.search("      prey1:", modulo_modificado[i]))
    modulo_modificado[prey_fila_idx:prey_fila_idx+len(coupling_lines)] = coupling_lines
    
    # Compute initialization values
    if equal_initialization:
        carbon = c_init_zoo
    elif separate_phytozoo:
        carbon = c_init_zoo / counts_init_zoo[np.where(clases_init == clasesP[j])[0][0]]
    else:
        carbon = c_init / counts_init[np.where(clases_init == clasesP[j])[0][0]]
        
    nitro = carbon * (0.0167 / 1.0)
    phosp = carbon * (0.00185 / 1.0)
    
    # Replace initialization values
    init_fila_idx = next(i for i in range(inicio[0], final[0]) if re.search("initialization:", modulo_modificado[i]))
    modulo_modificado[init_fila_idx + 1] = re.sub("1.0", str(round(carbon, 10)), modulo_modificado[init_fila_idx + 1])
    modulo_modificado[init_fila_idx + 2] = re.sub("0.0167", str(round(nitro, 10)), modulo_modificado[init_fila_idx + 2])
    modulo_modificado[init_fila_idx + 3] = re.sub("0.00185", str(round(phosp, 10)), modulo_modificado[init_fila_idx + 3])

    # Append to fabm.yaml
    with open(f"{outdir}/{nombrito}.yaml", 'a') as f:
        for line in modulo_modificado[inicio[0]+1:final[0]-1]:
            f.write(line + "\n")


# Z5
# Variables setup
old_name = "Z5"
old_long_name = "Microzooplankton"
clasesP = clasesZ5  # ESD = (((2**clasesP)*6)/pi)**(1/3)

# Parameters to modify
parameters = ["p_sum", "p_chuc"]
old_value = ["2.71",   "20"]
digito = [2, 0]
#digito = [3, 0]

old_nprey = 12
# Define all potential prey choices
preys = np.concatenate([clasesB1, clasesP2, clasesP3, clasesP4, clasesP5, clasesP6, clasesP7, clasesP8, clasesP9, clasesZ6, clasesZ5])
names_preys = ["B1"] + [f"P2_{str(cl).replace('-', 'm')}" for cl in clasesP2] + \
              [f"P3_{str(cl).replace('-', 'm')}" for cl in clasesP3] + \
              [f"P4_{str(cl).replace('-', 'm')}" for cl in clasesP4] + \
              [f"P5_{str(cl).replace('-', 'm')}" for cl in clasesP5] + \
              [f"P6_{str(cl).replace('-', 'm')}" for cl in clasesP6] + \
              [f"P7_{str(cl).replace('-', 'm')}" for cl in clasesP7] + \
              [f"P8_{str(cl).replace('-', 'm')}" for cl in clasesP8] + \
              [f"P9_{str(cl).replace('-', 'm')}" for cl in clasesP9] + \
              [f"Z6_{str(cl).replace('-', 'm')}" for cl in clasesZ6] + \
              [f"Z5_{str(cl).replace('-', 'm')}" for cl in clasesZ5]
names_preys = np.array(names_preys)
# Extract module block for Zx
inicio = [i for i, line in enumerate(modulo_productor) if re.search(f"start_{old_name}", line)]
final = [i for i, line in enumerate(modulo_productor) if re.search(f"end_{old_name}", line)]

lower_limit = np.power(2.,clasesP)
upper_limit = np.power(2.,(clasesP + 1))
diameter = np.power((lower_limit * 6) / np.pi, 1/3)

# Compute parameter values (allometry for taxon-specific)
p_sum = 26 * (np.power(diameter, -0.4))  # Maximum ingestion rate (Hansen 1997)
p_chuc = np.full(len(diameter), 120)  # Banas 2011, default value for this model

# Optimum prey size
diameter_preys = np.power((np.power(2.,preys) * 6) / np.pi, 1/3)

opt_prey_size = diameter / 10  # Kiorboe 2008, Ward 2012

# Loop over size classes of predators
for j in range(len(clasesP)):
    modulo_modificado = np.copy(modulo_productor)
    diameter_min = np.power((lower_limit[j] * 6) / np.pi, 1/3)
    
    # Replace module names
    modulo_modificado[inicio[0]+2] = re.sub(old_name, f"{old_name}_{clasesP[j]}", modulo_modificado[inicio[0]+2])
    modulo_modificado[inicio[0]+4] = re.sub(old_long_name, f"{old_long_name}_{round(diameter_min, 1)}", modulo_modificado[inicio[0]+4])
    
    # Replace parameter values
    for i, param in enumerate(parameters):
        param_value = locals()[param]  # dynamically retrieve variable
        filas = [idx for idx, line in enumerate(modulo_modificado[inicio[0]:final[0]]) if re.search(f"{param}:", line)]
        if filas:
            fila_idx = inicio[0] + filas[0]
            modulo_modificado[fila_idx] = re.sub(old_value[i], str(round(param_value[j], digito[i])), modulo_modificado[fila_idx])

    # Compute preferences, SD of gaussian = 0.20
    prefs = np.exp(-((np.log10(diameter_preys) - np.log10(opt_prey_size[j])) / 0.20)**2)
    prefs[prefs < 1e-2] = 0
    sorted_indices = np.argsort(2.**preys)
    
    # Replace number of preys (nprey)
    nprey_fila_idx = next(i for i in range(inicio[0], final[0]) if re.search("nprey", modulo_modificado[i]))
    modulo_modificado[nprey_fila_idx] = re.sub(str(old_nprey), str(len(names_preys[prefs != 0])), modulo_modificado[nprey_fila_idx])
    
    # Block of preferences
    pref_lines = [f"      suprey{idx+1}:   {round(pref, 2)}     # [-] Availability of {names_preys[idx]} to Z5_{clasesP[j]}"
                  for idx, pref in enumerate(prefs) if pref != 0]
    modulo_modificado[nprey_fila_idx+1:nprey_fila_idx+1+len(pref_lines)] = pref_lines


    # Block of calcifier
    boolean = np.zeros(len(names_preys[prefs != 0]), dtype=int)
    boolean[[i for i, name in enumerate(names_preys[prefs != 0]) if name.startswith("P5")]] = 1 #why P5?
    calcifier_lines = [f"      isP2{idx+1}:   {boolean[idx]}           # [-] identify P2 among the preys [1 for P2 otherwise 0]"
                       for idx in range(len(boolean))]
    
    calcifier_fila_idx = next(i for i in range(inicio[0], final[0]) if re.search("isP21:", modulo_modificado[i]))
    modulo_modificado[calcifier_fila_idx:calcifier_fila_idx+len(calcifier_lines)] = calcifier_lines

    # Block of coupling
    coupling_lines = [f"      prey{idx+1}:   {names_preys[idx]}" for idx, pref in enumerate(prefs) if pref != 0]
    prey_fila_idx = next(i for i in range(inicio[0], final[0]) if re.search("      prey1:", modulo_modificado[i]))
    modulo_modificado[prey_fila_idx:prey_fila_idx+len(coupling_lines)] = coupling_lines

    # Compute initialization values

    if equal_initialization:
        carbon = c_init_zoo
    elif separate_phytozoo:
        carbon = c_init_zoo / counts_init_zoo[np.where(clases_init == clasesP[j])[0][0]]
    else:
        carbon = c_init / counts_init[np.where(clases_init == clasesP[j])[0][0]]
    
    nitro = carbon * (0.0167 / 1.0)
    phosp = carbon * (0.00185 / 1.0)

    # Replace initialization values
    init_fila_idx = next(i for i in range(inicio[0], final[0]) if re.search("initialization:", modulo_modificado[i]))
    modulo_modificado[init_fila_idx + 1] = re.sub("1.0", str(round(carbon, 10)), modulo_modificado[init_fila_idx + 1])
    modulo_modificado[init_fila_idx + 2] = re.sub("0.0167", str(round(nitro, 10)), modulo_modificado[init_fila_idx + 2])
    modulo_modificado[init_fila_idx + 3] = re.sub("0.00185", str(round(phosp, 10)), modulo_modificado[init_fila_idx + 3])

    # Append to fabm.yaml
    with open(f"{outdir}/{nombrito}.yaml", 'a') as f:
        for line in modulo_modificado[inicio[0]+1:final[0]-1]:
            f.write(line + "\n")

# Z4

# Variables setup
old_name = "Z4"
old_long_name = "Omnivorous Mesozooplankton"
clasesP = clasesZ4  # ESD = (((2**clasesP)*6)/pi)**(1/3)

# Parameters to modify
parameters = ["p_sum"]
old_value = ["2.0"]
digito = [2]
#digito = [3]

old_nprey = 8
# Define all potential prey choices
preys = np.concatenate([clasesP1, clasesP2, clasesP4, clasesP5, clasesP7, clasesP8, clasesZ6, clasesZ5, clasesZ4])
names_preys = [f"P1_{str(cl).replace('-', 'm')}" for cl in clasesP1] + \
              [f"P2_{str(cl).replace('-', 'm')}" for cl in clasesP2] + \
              [f"P4_{str(cl).replace('-', 'm')}" for cl in clasesP4] + \
              [f"P5_{str(cl).replace('-', 'm')}" for cl in clasesP5] + \
              [f"P7_{str(cl).replace('-', 'm')}" for cl in clasesP7] + \
              [f"P8_{str(cl).replace('-', 'm')}" for cl in clasesP8] + \
              [f"Z6_{str(cl).replace('-', 'm')}" for cl in clasesZ6] + \
              [f"Z5_{str(cl).replace('-', 'm')}" for cl in clasesZ5] + \
              [f"Z4_{str(cl).replace('-', 'm')}" for cl in clasesZ4]
names_preys = np.array(names_preys)

# Extract module block for Zx
inicio = [i for i, line in enumerate(modulo_productor) if re.search(f"start_{old_name}", line)]
final = [i for i, line in enumerate(modulo_productor) if re.search(f"end_{old_name}", line)]

lower_limit = np.power(2.,clasesP)

upper_limit = np.power(2.,(clasesP + 1))
diameter = np.power((lower_limit * 6) / np.pi, 1/3)

# Compute parameter values (allometry for taxon-specific)
p_sum = 26 * (np.power(diameter, -0.4))  # Maximum ingestion rate (Hansen 1997)

# Optimum prey size
diameter_preys = np.power((np.power(2.,preys) * 6) / np.pi,1/3)
opt_prey_size = diameter / 10  # Kiorboe 2008, Ward 2012

# Loop over size classes of predators
for j in range(len(clasesP)):
    modulo_modificado = np.copy(modulo_productor)
    diameter_min = np.power((lower_limit[j] * 6) / np.pi, 1/3)
    
    # Replace module names
    modulo_modificado[inicio[0]+2] = re.sub(old_name, f"{old_name}_{clasesP[j]}", modulo_modificado[inicio[0]+2])
    modulo_modificado[inicio[0]+4] = re.sub(old_long_name, f"{old_long_name}_{round(diameter_min, 1)}", modulo_modificado[inicio[0]+4])
    
    # Replace parameter values
    for i, param in enumerate(parameters):
        param_value = locals()[param]  # dynamically retrieve variable
        filas = [idx for idx, line in enumerate(modulo_modificado[inicio[0]:final[0]]) if re.search(f"{param}:", line)]
        if filas:
            fila_idx = inicio[0] + filas[0]
            modulo_modificado[fila_idx] = re.sub(old_value[i], str(round(param_value[j], digito[i])), modulo_modificado[fila_idx])

    # Compute preferences, SD of gaussian = 0.20
    prefs = np.exp(-((np.log10(diameter_preys) - np.log10(opt_prey_size[j])) / 0.20)**2)
    prefs[prefs < 1e-2] = 0
    sorted_indices = np.argsort(np.power(2.,preys))
    
    # Replace number of preys (nprey)
    nprey_fila_idx = next(i for i in range(inicio[0], final[0]) if re.search("nprey", modulo_modificado[i]))
    modulo_modificado[nprey_fila_idx] = re.sub(str(old_nprey), str(len(names_preys[prefs != 0])), modulo_modificado[nprey_fila_idx])
    
    # Block of preferences
    pref_lines = [f"      suprey{idx+1}:   {round(pref, 2)}     # [-] Availability of {names_preys[idx]} to Z4_{clasesP[j]}"
                  for idx, pref in enumerate(prefs) if pref != 0]
    modulo_modificado[nprey_fila_idx+1:nprey_fila_idx+1+len(pref_lines)] = pref_lines

    # Block of calcifier
    boolean = np.zeros(len(names_preys[prefs != 0]), dtype=int)
    boolean[[i for i, name in enumerate(names_preys[prefs != 0]) if name.startswith("P5")]] = 1 #why P5?

    calcifier_lines = [f"      isP2{idx+1}:   {boolean[idx]}           # [-] identify P2 among the preys [1 for P2 otherwise 0]"
                          for idx in range(len(boolean))]
    calcifier_fila_idx = next(i for i in range(inicio[0], final[0]) if re.search("isP21:", modulo_modificado[i]))
    modulo_modificado[calcifier_fila_idx:calcifier_fila_idx+len(calcifier_lines)] = calcifier_lines

    # Block of coupling
    coupling_lines = [f"      prey{idx+1}:   {names_preys[idx]}" for idx, pref in enumerate(prefs) if pref != 0]
    prey_fila_idx = next(i for i in range(inicio[0], final[0]) if re.search("      prey1:", modulo_modificado[i]))
    modulo_modificado[prey_fila_idx:prey_fila_idx+len(coupling_lines)] = coupling_lines

    # Compute initialization values
    if equal_initialization:
        carbon = c_init_zoo
    elif separate_phytozoo:
        carbon = c_init_zoo / counts_init_zoo[np.where(clases_init == clasesP[j])[0][0]]
    else:
        carbon = c_init / counts_init[np.where(clases_init == clasesP[j])[0][0]]

    nitro = carbon * (0.0167 / 1.0)
    phosp = carbon * (0.00185 / 1.0)

    # Replace initialization values
    init_fila_idx = next(i for i in range(inicio[0], final[0]) if re.search("initialization:", modulo_modificado[i]))
    modulo_modificado[init_fila_idx + 1] = re.sub("1.0", str(round(carbon, 10)), modulo_modificado[init_fila_idx + 1])
    modulo_modificado[init_fila_idx + 2] = re.sub("0.0167", str(round(nitro, 10)), modulo_modificado[init_fila_idx + 2])
    modulo_modificado[init_fila_idx + 3] = re.sub("0.00185", str(round(phosp, 10)), modulo_modificado[init_fila_idx + 3])

    # Append to fabm.yaml
    with open(f"{outdir}/{nombrito}.yaml", 'a') as f:
        for line in modulo_modificado[inicio[0]+1:final[0]-1]:
            f.write(line + "\n")


# Z3

old_name = "Z3"
old_long_name = "Carnivorous Mesozooplankton"
clasesP = clasesZ3  # ESD = (((2**clasesP)*6)/pi)**(1/3)

# Parameters to modify
parameters = ["p_sum"]
old_value = ["2.0"]
digito = [2]
#digito = [3]

old_nprey = 3
# Define all potential prey choices
preys = np.concatenate([clasesP1, clasesP3, clasesZ5, clasesZ4, clasesZ3])
names_preys = [f"P1_{str(cl).replace('-', 'm')}" for cl in clasesP1] + \
              [f"P3_{str(cl).replace('-', 'm')}" for cl in clasesP3] + \
              [f"Z5_{str(cl).replace('-', 'm')}" for cl in clasesZ5] + \
              [f"Z4_{str(cl).replace('-', 'm')}" for cl in clasesZ4] + \
              [f"Z3_{str(cl).replace('-', 'm')}" for cl in clasesZ3]
names_preys = np.array(names_preys)

# Extract module block for Zx
inicio = [i for i, line in enumerate(modulo_productor) if re.search(f"start_{old_name}", line)]
final = [i for i, line in enumerate(modulo_productor) if re.search(f"end_{old_name}", line)]

lower_limit = np.power(2.,clasesP)
upper_limit = np.power(2.,(clasesP + 1))
diameter = np.power((lower_limit * 6) / np.pi, 1/3)

# Compute parameter values (allometry for taxon-specific)
p_sum = 26 * np.power(diameter, -0.4)  # Maximum ingestion rate (Hansen 1997)


# Optimum prey size
diameter_preys = np.power((np.power(2.,preys) * 6) / np.pi,1/3)
opt_prey_size = diameter / 10  # Kiorboe 2008, Ward 2012

# Loop over size classes of predators
for j in range(len(clasesP)):
    modulo_modificado = np.copy(modulo_productor)
    diameter_min = np.power((lower_limit[j] * 6) / np.pi, 1/3)
    
    # Replace module names
    modulo_modificado[inicio[0]+2] = re.sub(old_name, f"{old_name}_{clasesP[j]}", modulo_modificado[inicio[0]+2])
    modulo_modificado[inicio[0]+4] = re.sub(old_long_name, f"{old_long_name}_{round(diameter_min, 1)}", modulo_modificado[inicio[0]+4])
    
    # Replace parameter values
    for i, param in enumerate(parameters):
        param_value = locals()[param]  # dynamically retrieve variable
        filas = [idx for idx, line in enumerate(modulo_modificado[inicio[0]:final[0]]) if re.search(f"{param}:", line)]
        if filas:
            fila_idx = inicio[0] + filas[0]
            modulo_modificado[fila_idx] = re.sub(old_value[i], str(round(param_value[j], digito[i])), modulo_modificado[fila_idx])

    # Compute preferences, SD of gaussian = 0.20
    prefs = np.exp(-((np.log10(diameter_preys) - np.log10(opt_prey_size[j])) / 0.20)**2)
    prefs[prefs < 1e-2] = 0
    sorted_indices = np.argsort(2.**preys)
    
    # Replace number of preys (nprey)
    nprey_fila_idx = next(i for i in range(inicio[0], final[0]) if re.search("nprey", modulo_modificado[i]))
    modulo_modificado[nprey_fila_idx] = re.sub(str(old_nprey), str(len(names_preys[prefs != 0])), modulo_modificado[nprey_fila_idx])
    
    # Block of preferences
    pref_lines = [f"      suprey{idx+1}:   {round(pref, 2)}     # [-] Availability of {names_preys[idx]} to Z3_{clasesP[j]}"
                  for idx, pref in enumerate(prefs) if pref != 0]
    modulo_modificado[nprey_fila_idx+1:nprey_fila_idx+1+len(pref_lines)] = pref_lines

    # Block of calcifier

    boolean = np.zeros(len(names_preys[prefs != 0]), dtype=int)
    boolean[[i for i, name in enumerate(names_preys[prefs != 0]) if name.startswith("P5")]] = 1 #why P5?
    calcifier_lines = [f"      isP2{idx+1}:   {boolean[idx]}           # [-] identify P2 among the preys [1 for P2 otherwise 0]"
                       for idx in range(len(boolean))]
    
    calcifier_fila_idx = next(i for i in range(inicio[0], final[0]) if re.search("isP21:", modulo_modificado[i]))
    modulo_modificado[calcifier_fila_idx:calcifier_fila_idx+len(calcifier_lines)] = calcifier_lines

    # Block of coupling
    coupling_lines = [f"      prey{idx+1}:   {names_preys[idx]}" for idx, pref in enumerate(prefs) if pref != 0]
    prey_fila_idx = next(i for i in range(inicio[0], final[0]) if re.search("      prey1:", modulo_modificado[i]))
    modulo_modificado[prey_fila_idx:prey_fila_idx+len(coupling_lines)] = coupling_lines

    # Compute initialization values
    if equal_initialization:
        carbon = c_init_zoo
    elif separate_phytozoo:
        carbon = c_init_zoo / counts_init_zoo[np.where(clases_init == clasesP[j])[0][0]]
    else:
        carbon = c_init / counts_init[np.where(clases_init == clasesP[j])[0][0]]

    nitro = carbon * (0.0167 / 1.0)
    phosp = carbon * (0.00185 / 1.0)

    # Replace initialization values
    init_fila_idx = next(i for i in range(inicio[0], final[0]) if re.search("initialization:", modulo_modificado[i]))
    modulo_modificado[init_fila_idx + 1] = re.sub("1.0", str(round(carbon, 10)), modulo_modificado[init_fila_idx + 1])
    modulo_modificado[init_fila_idx + 2] = re.sub("0.0167", str(round(nitro, 10)), modulo_modificado[init_fila_idx + 2])
    modulo_modificado[init_fila_idx + 3] = re.sub("0.00185", str(round(phosp, 10)), modulo_modificado[init_fila_idx + 3])

    # Append to fabm.yaml
    with open(f"{outdir}/{nombrito}.yaml", 'a') as f:
        for line in modulo_modificado[inicio[0]+1:final[0]-1]:
            f.write(line + "\n")




                      







