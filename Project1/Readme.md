# Project1


## Reproducing results

### Problem 7
For compiling, run the following command in terminal

```
g++ P7.cpp -o P7
```

Then run the following command for producing the data files with 10, 100 and 1000 steps 

```
for nsteps in 10 100 1000 ; do ./P7 $nsteps; done
```

and finally for plotting the results run 

```
python3 Plot_P7.py 10 100 1000
```

### Problem 9

For compiling, run the following command in terminal

```
g++ P9.cpp -o P9
```

Then run the following command for producing the data file with whatever nsteps desired. The example for P9 only contains nsteps = 1000, but multiple can be used similarly to what was done in P7

```
for nsteps in 1000; do ./P9 $nsteps; done
```


### Problem 10
OBS: the iteration 
```
for (std::string mode : modes_vector)
```
is a C++11 range-based for loop. So in case g++ or Clang is not using C++11 as default, the compiling command needs to be 

```
g++ --std=c++11 problem10.cpp -o problem10
```

Then, run the following command for producing the data files with steps up to 10^6

```
for nsteps in 10 100 1000 10000 100000 1000000 ; do ./P10 $nsteps; done
```
