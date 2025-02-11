import numpy as np
import pandas as pd
import os

temps = [1000,1040,1080,1120]
s = 4 # this is an external parameter that needs to be consistent across scripts
with open("lattice_constants.csv",'w') as f:
	f.write("T,lx\n")
	for T in Ts:
		os.system(f"python parseLammps.py -i eq_solid_T-{T}.log")
		df = pd.read_csv(f"eq_solid_T-{T}.log.csv")
		nUse = int(0.5*df.shape[0]) # only average over last half of each sim 
		lx = (1/s) * np.mean(np.array(df['Lx'][-nUse:]))
		line = f"{T},{lx}\n"
		print(line)
		f.write(line)

