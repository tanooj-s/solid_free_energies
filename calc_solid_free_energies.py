# script to pull relevant LAMMPS log data and compute solid absolute free energies from reversible pathways
# using NaCl as template material to validate methodology
# pathway is run such that λ=0 is "actual material" and λ=1 is the harmonic crystal reference state
# LAMMPS simulations are set up in a way which should be as general as possible for whichever force field used
# pair_style hybrid/scaled to scale down interaction potential, and fix spring/self with scaled spring constants to tether particles to lattice sites 
# this is complete pending an interpretable insertion of the finite size correction to the center of mass constraint oscillator term, but that can also be sidestepped by running these sims at multiple box sizes

# note that there are two components of the "thermodynamic driving force" for this pathway 
# one, which is from deviations of particles from exact lattice sites
# and then there's an additional "large" driving force on top of that which is the change in the internal energy of the crystal when particles are at exact lattice sites
# at λ=0 this is the potential energy for particles frozen at lattice sites, what I call the "static crystal potential energy"
# at λ=1 there is no internal energy for particles at exact lattice sites
# (or really that λ=1 term would be infinite as per K(x-x0)^2, since springs would be "maximally compressed" in that case, but that is another rabbit hole irrelevant for this application)

# results validated against https://doi.org/10.1063/1.1522375
# target number: Gsol = -97.75 ± 0.02 kBT per ion pair at 1074K, N=512
# numbers obtained here for N=512
# T     | Gsol (kBT per ion pair)
# 1000K | -103.669
# 1040K | -100.297
# 1080K | -97.207
# 1120K | -94.400

import os
import json
import glob
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
parser = argparse.ArgumentParser(description="compute solid free energies")
parser.add_argument('-s', action="store", dest="s") # supercell size, external parameter because you should confirm size dependence of these calcs
args = parser.parse_args()
pwd = os.getcwd()
s = int(args.s)  


# ---- universal physical constants ----
kB = 1.380649e-23 # in joules = 8.61733e-5 eV
h = 6.62607015e-34  # in J.s
hbar = h / (2*np.pi)

# ---- unit conversions (set here for LAMMPS metal units) ----
# LAMMPS metal units - energy in eV, presure in bars, distance in angstroms, mass in AMUs
angstrom_to_meter = 1e-10
bar_to_pascal = 1e5
joule_to_eV = 6.241509e18
amu_to_kg = 1.6605402e-27

# ---- crystal specific parameters (using NaCl as template material, edit as appropriate) ----
m_Na = 22.9898
m_Cl = 35.453 
particles_per_unit_cell = 8 # 8 atoms in a rocksalt lattice unit cell as per introductory materials science courses
N = particles_per_unit_cell * (s ** 3) # total number of particles in simbox

# ---- pathway parameters (these should correspond to the actual simulations that were run) -----
temps = [1000,1040,1080,1120]
lambdas = np.arange(0.0,1.01,0.05)


# ---- analytical function for free energy of a harmonic crystal -----

def harmonicHelmholtz(K,T,m):
    '''
    return analytical free energy of single classical harmonic oscillator
    assume K in eV/A^2, m in AMUs, T in kelvin
    '''
    # convert everything to SI units to obtain effective oscillator frequency
    K *= 1/joule_to_eV # eV/A^2 => joules/A^2
    K /= (angstrom_to_meter ** 2) # joules/A^2 => joules/m^2
    m *= amu_to_kg # AMU => kg
    w = np.sqrt(K/m)  # angular frequency in Hz
    kT = kB * T # joules
    A = 3 * kT * np.log(hbar*w/kT) 
    A *= joule_to_eV # convert back to eV
    return [w,A] # return frequency and free energy
# TODO later - write a function to generalize above for the center of mass term as well, which is computed explicitly below


# ---- read in computed Keff(T) and assign to data structure ----
# these spring constants should be determined from prior NPT/NVT runs  
# edit below either as per specific material
spring_constants = { 'Na': {}, 'Cl': {} }
with open("spring_constants.csv","r") as f:
    lines = f.readlines()
    for l in lines[1:]:
        tokens = l.strip('\n').split(',')
        T, K_Na, K_Cl = int(tokens[0]), float(tokens[2]), float(tokens[3])
        spring_constants['Na'][T] = K_Na
        spring_constants['Cl'][T] = K_Cl

print("Effective spring constant as a function of temperature for each element")
print(f"Keff(T) (eV/Å^2)") # edit energy units per LAMMPS 
print(json.dumps(spring_constants, indent=2))
print("-------------------------------")


# ---- convert every TI log file to a csv ----
for fname in glob.iglob(f"{pwd}/T-*spring.log"):
    print(fname)
    os.system(f"python parseLammps.py -i {fname}")


# ---- now actually parse simulation logs and compute free energies ----
# plot dU/dλ(λ) for each temperature 
plt.rcParams['figure.figsize'] = (5*len(temps),4) 
plt.rcParams['font.size'] = 12
fig, axes = plt.subplots(1,len(temps),sharey=True) 
axes[0].set_ylabel('eV')
# note axes[i].plot commands will complain if only run at one temperature, so don't use plt.subplot in that case

# iterate over temperatures and pull relevant values from logs 
with open(f"{pwd}/solid_free_energies.csv",'w') as outf:
    outf.write(f"T,kT,dA,VdP,PV,wosc,wcom,Aosc,Acom,Aref,Gsol\n")
    # units for each column of output csv file (implicit)
    # kelvin,eV,kT,kT,kT,Hz,Hz,kT,kT,kT,kT
    # free energy terms in kBT per particle (edit these intensive numbers for "per molecular unit" values later if necessary)
    for i in range(len(temps)):
        T = temps[i]
        kT = kB * T * joule_to_eV 
        print(f"T={T} K | kT = {kT} eV")
        # create a separate list/array (later convert to np arrays) for each relevant energy component of dU/dλ
        # 1. spring energy for each element, scalar value output by fix spring/self 
        # 2. potential energy of scaled interaction potential, directly from LAMMPS PE
        # 3. static crystal potential energy at each value of λ - the interaction potential term for that value of λ when all particles are at exact lattice sites
        # 4. the static crystal potential energy at λ=0, internal energy of the crystal at fully scaled interaction potential
        # averaged energy terms as functions of λ
        spring1, spring2 = [], [] # extend list of lists as per number of elements, spring energy for every element at each λ
        pe = [] # LAMMPS computed PotEng at each λ
        stat = [] # static crystal energy at each λ 

        # pull static energy from start of pathway
        df0 = pd.read_csv(f"{pwd}/T-{T}_lam-0.0_spring.log.csv") 
        Ustat0 = df0['PotEng'][0] # static crystal potential energy at λ=0
        
        # now iterate over every simulated λ along pathway
        # note that all energy values are thermodynamically extensive, so that we deal with the logic of normalization ourselves
        # ensure this by confirming thermo_modify norm no in LAMMPS input scripts  
        for l in lambdas:
            df = pd.read_csv(f"{pwd}/T-{T}_lam-{round(l,2)}_spring.log.csv")
            nUse = int(0.5*df.shape[0])
            if l < 0.00001: 
                spring1.append(0)
                spring2.append(0)
            else:
                spring1.append(np.mean(np.array(df['f_0'])[-nUse:]))
                spring2.append(np.mean(np.array(df['f_1'])[-nUse:])) # different spring/self fix for each element
            pe.append(np.mean(np.array(df['PotEng'])[-nUse:]))
            stat.append(np.array(df['PotEng'])[0])

        spring1, spring2, pe, stat = np.array(spring1), np.array(spring2), np.array(pe), np.array(stat)
        # rename these arrays for physical clarity of computation below
        Uharmonic = spring1 + spring2 # + spring3 + ...
        Uinteraction = pe
        Ustatic = stat
        #print(f"Ustat0: {Ustat0}")
        #print(f"Ustat(λ): {Ustat}")
        #print(f"Uint(λ): {Uint}")
        #print(f"Uref(λ): {Uref}") # print for debugging

        # plot dU/dλ(λ), integrate under curve for transformation free energy
        tdfs = -Ustat0 + (Uharmonic - (Uinteraction - Ustatic)) # dU/dλ(λ)
        # as mentioned above, there would also be a + Ureference1 term here in theory, but it is 0
        dA_ext = np.trapz(x=lambdas,y=tdfs) # dA for entire system in whichever energy units 
        axes[i].scatter(x=lambdas,y=tdfs,label=f"${{ΔA_{{TI}}}}={round(dA_ext,3)}eV$") # using inline TeX formatting
        axes[i].set_title(f"$T={T}K$, $N={{{N}}}$ ∂U/∂λ(λ), extensive")
        axes[i].legend(loc="lower right")
        axes[i].set_xlabel('λ')
        axes[i].grid()
        dA_ext_kT = dA_ext / kT
        dA_int_kT = dA_ext_kT / N # this is the number we want to report, free energies per particle in terms of kBT
        # should assert df['Atoms'] == N earlier in this script
        
        # compute PV terms 
        df0 = pd.read_csv(f"{pwd}/T-{T}_lam-0.0_spring.log.csv")
        df1 = pd.read_csv(f"{pwd}/T-{T}_lam-1.0_spring.log.csv")
        V = float(df1['Volume'][0]) # should assert this is the same 
        P0, P1 = np.mean(np.array(df0['Press'])[-nUse:]), np.mean(np.array(df1['Press'])[-nUse:]) 
        # assuming every sim over the pathway is the same number of timesteps 
        dP = P1 - P0
        VdP_ext = (V * (angstrom_to_meter ** 3) * dP * (bar_to_pascal)) * joule_to_eV # change in simbox pressure due to transformation
        PVref_ext = (P1 * bar_to_pascal * V * (angstrom_to_meter ** 3)) * joule_to_eV # PV of the harmonic crystal
        VdP_ext_kT = VdP_ext / kT
        PVref_ext_kT = PVref_ext / kT
        VdP_int_kT = VdP_ext_kT / N
        PVref_int_kT = PVref_ext_kT / N


        # --- free energy of the harmonic crystal reference state ----
        # separate this out as center of mass constraint term and sum of oscillators term

        # ---- individual oscillators ----
        K_Na = spring_constants["Na"][T]
        K_Cl = spring_constants["Cl"][T]
        w_osc = 0.5 * (harmonicHelmholtz(K_Na,T,m_Na)[0] + harmonicHelmholtz(K_Cl,T,m_Cl)[0])
        A_osc = 0.5 * (harmonicHelmholtz(K_Na,T,m_Na)[1] + harmonicHelmholtz(K_Cl,T,m_Cl)[1])
        # note this oscillator term is a mean to obtain an intensive number

        # ---- center of mass oscillator ----
        # recompute whatever is done within harmonicHelmholtz with the appropriate frequency
        # this vibrates in the opposite "net" direction of the individual oscillators to keep the center of mass of the simulated crystal stationary
        M = (N/2)*m_Na + (N/2)*m_Cl # total mass 
        mu_Na = m_Na / (N/2)
        mu_Cl = m_Cl / (N/2) # fractional mass of each element
        # use these and spring constants for effective frequency of center of mass

        # wbar_com is extremely stoichiometry dependent, later write a function to compute this systematically 
        wbar_com = 1/(mu_Na*mu_Cl) * (((K_Na*(1/joule_to_eV)/(angstrom_to_meter**2))*(K_Cl*(1/joule_to_eV)/(angstrom_to_meter**2)) / (2*np.pi*M*amu_to_kg*(N/2)*(N/2))) ** 0.5)
        molecular_units = N/2 # this is also intensive now
        A_com = 3*kT/molecular_units * np.log(hbar*wbar_com/kT) 
        # no need to reconvert to eV since kT computed at beginning of temperature loop is already in eV 

        # convert to kT
        A_osc_kT = A_osc / kT
        A_com_kT = A_com / kT
        A_ref_kT = A_osc_kT - A_com_kT 

        # no finite size correction term needed - instead run at multiple box sizes and fit observed G(N) curve to either an lnN/N or power law decay curve for large N limit
        

        # ----------------

        Gsolid = A_ref_kT + PVref_int_kT - dA_int_kT - VdP_int_kT

        # output oscillator frequencies interpretably
        w_osc = "{:e}".format(w_osc)
        wbar_com = "{:e}".format(wbar_com) 

        # sanity check prints
        print(f"  V (cubic Å): {V}")
        print(f"  P0 (bars): {P0}")
        print(f"  P1 (bars): {P1}")
        print(f"  VdP (kT, intensive): {VdP_int_kT}")
        print(f"  PV_ref (kT, intensive): {PVref_int_kT}")
        print(f"  dA (kT, intensive): {dA_int_kT}")
        print(f"  w_osc (Hz): {w_osc}")
        print(f"  wbar_com (Hz): {wbar_com}")
        print(f"  A_osc (kT, intensive): {A_osc_kT}")
        print(f"  A_com (kT, intensive): {A_com_kT}")
        print(f"  A_ref (kT, intensive): {A_ref_kT}") 
        print(f"  Gsolid (kT, intensive): {Gsolid}")
        print("-------------------------------")
        # edit for extensive and direct energy unit values as appropriate in case of debugging
        csv_line = f"{T},{kT},{dA_int_kT},{VdP_int_kT},{PVref_int_kT},{w_osc},{wbar_com},{A_osc_kT},{A_com_kT},{A_ref_kT},{Gsolid}\n"
        outf.write(csv_line)
plt.tight_layout()
plt.savefig(f"{pwd}/tdf_curves.png")