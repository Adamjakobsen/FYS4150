# Project1


## Reproducing results

### Problem 7
For compiling, run the following command in terminal

```
g++ P2.cpp -o P2
```

Then run the following command for producing the data file

```
./P2
```

and finally for plotting the results run 

```
python3 Plot_P2.py
```

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

### Problem 8

OBS: Before compiling problem 8, we need to re-execute P7 with more input steps as below

```
for nsteps in 10 100 1000 10000 100000 1000000 10000000 ; do ./P7 $nsteps; done
```

now we can properly plot the data from problem 8 as follows

```
python3 Plot_P8.py
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
OBS: the iteration in line 16
```
for (std::string mode : modes_vector)
```
is a C++11 range-based for loop. So in case g++ or Clang is not using C++11 as default, the compiling command needs to be 

```
g++ --std=c++11 P10.cpp -o P10
```

In case C++ is already default, for compliling and linking simply run the default

```
g++ P10.cpp -o P10
```

Then, run the following command for producing the data files with steps up to 10^6

```
for nsteps in 10 50 100 500 1000 5000 10000 50000 100000 500000 1000000 ; do ./P10 $nsteps; done
```

This command will generate 1 .txt file for each nstep, for each algorithm used.
For plotting the result graph in a pdf file then run 

```
python3 Plot_P10.py
```