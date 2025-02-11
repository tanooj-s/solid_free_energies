# calculate spring constants from MSD sims and also plot MSD as a function of time

import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
s = 4
temps = [1000,1040,1080,1120]
# first pull log info into csvs

for T in temps: os.system(f"python parseLammps.py -i MSD_T-{T}.log")

# create lists to plot spring constant as a function of time
K_nas = []
K_cls = []

# note that we want to write out lattice constant as a function of temperature to this file as well for later TI sims
with open(f"spring_constants.csv",'w') as f:
    f.write(f"T,lx,K_na,K_cl\n") # adjust as per number of elements
    for T in Ts:
        df = pd.read_csv(f"MSD_T-{T}.log.csv")
        nUse = int(0.5*df.shape[0])
        kT = 8.61733e-5 * T # kBT in eV 
        MSD_na = np.mean(np.array(df['c_msdNa[4]'])[-nUse:])
        MSD_cl = np.mean(np.array(df['c_msdCl[4]'])[-nUse:]) # note this is assuming rough constancy of MSD
        lx = (1/s) * np.mean(np.array(df['Lx'][-nUse:])) 
        K_na = 3*kT/MSD_na 
        K_cl = 3*kT/MSD_cl 
        f.write(f"{T},{lx},{K_na},{K_cl}\n")  
        K_nas.append(K_na)
        K_cls.append(K_cl)

# plot spring constant as a function of temperature, only two elements here so both on the same panel
plt.rcParams['figure.figsize'] = (6,4)
plt.scatter(x=Ts,y=K_nas,label="Na")
plt.scatter(x=Ts,y=K_cls,label="Cl")
plt.xlabel("Temperature (K)")
plt.ylabel("Spring constant (eV/Å^2)")
plt.ylim(0,1.2*np.max(np.concatenate((np.array(K_nas),np.array(K_cls)),axis=0)))
plt.legend(loc="lower left")
plt.grid()
plt.savefig(f"spring_constants.png")


