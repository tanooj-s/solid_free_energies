# assuming TI sims have been run and free energies computed at multiple box sizes, within a different folder for each box
# pull free energy data and fit to logarithmic function to extrapolate to macroscopic values
# (may need to fit to a different type of function for different classes of materials)

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

# logarithmic function to fit G(N)
def log_correction(x,a,b): 
	return a + b*np.log(x)/x

temps = np.array([1000,1040,1080,1120])
numbers = np.array([512,1000,2744,5832])

# first read in all data to a single dataframe
dfs = []
for N in numbers:
    dfs.append(pd.read_csv(f"N{N}/solid_free_energies.csv"))
    # note naming convention used here for box sizes
    # also note that free energies computed in these csvs are per particle
df = pd.concat(dfs,axis=0,ignore_index=True)

plt.rcParams['figure.figsize'] = (6,2.5*len(temps))
fig, axes = plt.subplots(4,1,sharex=True)

Ginfs = [] # large N limit for each Gsol(N), as a function of temperature

for i in range(len(temps)):
    T = temps[i]
    print(f"T={T}K")
    Gsols = np.array((df[df['T']==T]['Gsol'])) 
    # Initial guesses for parameters of fit function
    a_guess = Gsols[-1]
    b_guess = np.mean((Gsols - a_guess) * numbers / np.log(numbers))
    popt, pcov = curve_fit(log_correction, numbers, Gsols, p0=[a_guess, b_guess])
    a_fit, b_fit = popt
    print(f"  Gsol(N): {Gsols}")
    print(f"  b_fit: {b_fit}") # "decay" parameter
    print(f"  a_fit: {a_fit}") # asymptotic value of G at large N

    # dense set of points to plot fitted trend
    Ndense = np.arange(numbers[0],numbers[-1]+1,5)
    axes[i].scatter(x=numbers,y=Gsols,label="collected MD data")
    axes[i].plot(Ndense,log_correction(Ndense,a_fit,b_fit),label="fitted lnN/N",color='k')
    axes[i].axhline(a_fit,label=f"G_inf={round(a_fit,3)}",color='r',linestyle='dashed')
    axes[i].set_title(f"T={T}K ${{G_{{sol}}(N)}}$")
    axes[i].grid()
    axes[i].legend(loc="upper right")
    axes[i].set_ylabel("kBT / particle")
    Ginfs.append(a_fit)
    
axes[-1].set_xlabel("Number of particles")
plt.tight_layout()
plt.savefig("large_N_extrapolation.png")

print("-------------")
print("Gsol(T) at large N")
print(temps)
print(Ginfs)