## Study of the Ising model critical temperature and phase transitions

We will be simulating the Ising model using the Metropolis method, a Markov Chain Monte Carlo (MCMC) method. 
Our goal is to find the critical temperature of the Ising model in order to study its phase transitions and magneitzation.

We will be implementing the 2D lattice version of the Ising model, using a 2x2 case as our analytical solution. 
We increase the lattice size eventually and run for multiple Monte Carlo cycles.




Compile: 
```
make compile
```
OBS: make sure you have \emph{openmp} installed


Run example:
```
make run L=2 mc_cycles=100000 burn_pct=10 lower_temp=1 upper_temp=2 temp_step=0.5 align=0
```
- Parameters 
  - L: size of the grid of a total of N = LXL electrons
  - mc_cycles: number of N atemped randomly selected spin flips
  - burn_pct: percentage of mc_cycles to burn
  - lower_temp: temperature to begin the loop of temperatures (each temperature is saved in a different file)
  - upper_temp: temperature to stop the loop of temperatures
  - temp_step: step to make beging from lower_temp, inclusive
  - align: 0 is initialize spin configuration randomly, 1 is initialize all up, 2 is initialize all down.

