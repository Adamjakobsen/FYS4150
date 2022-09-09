# Project1


## Reproducing results

### Problem 7
For compiling, run the following command in terminal

```
g++ P7.cpp -o P7
```

Then run the following command for producing the data files 

```
for nsteps in 10 100 1000 10000; do ./P7 $nsteps; done
```

and finally for plotting the results run 


