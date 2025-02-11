# run simulations to compute mean squared displacements for every (temperature, lattice constant) tuple
import numpy as np
import pandas as pd
import os
s = 4
with open(f"lattice_constants.csv",'r') as f:
    lines = f.readlines()
    for l in lines[1:]:
        tokens = l.strip('\n').split(',')
        T, lx = tokens[0], tokens[1]
        print(f"T={T}, lx={lx}")
        os.system(f"sbatch MSD.job {lx} {T} {s}")
