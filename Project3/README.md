# Project 3:

## Compiling and running
### inputs: 

method: Use "RK4" or "Euler"

n_particles: number of particles (int)

total_time: Total simulation time (int)

intercation: Use "on" or "off"

### Compiling one particle:

```
make compile_one
```
### Run:

```
./main method n_timesteps total_time 
```
### Compiling two particles:

```
make compile_two
```

### Run:

```
./main method n_timesteps total_time intercation 
```

### Compiling multiple particles:
```
make compile_multi
```

### Run

```
./main method n_timesteps total_time intercation n_particles
```

### Compiling with perturbation:

```
make compile_perturb
```

### Run

```
./main method n_timesteps total_time intercation n_particles
```

## Plotting
All plotting scripts are in '/plotting'-folder. Run from said folder after producing the data.
