Workflow:

1. run_solid_eq.py to run NPT sims for equilbrium lattice constants
   - calc_lattice_constants.py to compute above
2. run_MSD.py to run NVT sims to estimate spring constants
   - calc_spring_constants.py to compute above
3. run_tether_sims.py to run solid -> harmonic crystal pathway sims
4. calc_solid_free_energies.py to pull relevant data from logs and compute absolute free energies
5. extrapolate_large_N.py to extrapolate computed free energies to macroscopic length scales

Things to note:  

- this entire workflow (except step 5) should be run at multiple system sizes (at least 4 points for an observable trend)
- each system size should be run in a different folder, so for example I would mkdir N512 mkdir N1000 mkdir N2744 mkdir N5832, and then copy the relevant scripts in this directory (except extrapolate_large_N.py) to each of those subdirectories
- (update: this pdf needs to be updated) there is a pdf included in this directory that contains somewhat of a rederivation of the analytical free energy of the finite sized harmonic crystal, in a form that seems more easily interpretable to me than the expression typically used in the literature
- I have also included multiple other relevant background papers from the literature in the ref/ directory for the interested reader 
