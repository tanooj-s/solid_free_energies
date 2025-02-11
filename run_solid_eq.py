import os

# list of temperatures needs to be consistent across all python scripts
temps = [1000,1040,1080,1120] 
s = 4

for T in temps:
    print(f"Running solid NPT equilibration at {T}K")
    os.system(f"sbatch eq_solid.job {T} {s}")
