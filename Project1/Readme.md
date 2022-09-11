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

