import os
import numpy as np
import argparse
parser = argparse.ArgumentParser(description="run harmonic crystal pathway simulations")
parser.add_argument('-s',action="store",dest="s")
args = parser.parse_args()

s = int(args.s) # supercell size
lambdas = np.arange(0.00,1.01,0.05)

# note that a separate LAMMPS input file is needed because fix spring/self errors out when you pass in a spring constant of 0
# or alternatively one could edit fix_spring_self.cpp to not error out and recompile LAMMPS
with open("spring_constants.csv","r") as f:
    lines = f.readlines()
    for line in lines[1:]:
        tokens = line.strip('\n').split(',')
        T, lc, K_na, K_cl = int(tokens[0]), float(tokens[1]), float(tokens[2]), float(tokens[3])
        for l in lambdas:
            l = round(l,2)
            print(f"Running TI sim at T={T}, lambda={l}, Kna={K_na*l}, Kcl={K_cl*l}")
            if l == 0.0:
                os.system(f"sbatch spring0.job {T} {lc} {s}")
            else:
                os.system(f"sbatch spring.job {T} {lc} {K_na} {K_cl} {l} {s}")
                # note that scaling of each spring constant is done within LAMMPS input script
