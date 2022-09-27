# Project 2

All problems are solved in the the "main.cpp" file. To ask for the specific output and txt files you want generated, simply compile and link by running

```
g++ main.cpp utils.cpp -o main -std=c++11 -larmadillo  
```
and execute the code by running 
```
./main PN
```
Where N is the number of the desired problem. An example for executing problem 2 should then be  "./main P2".

**Exception:** 
The only exception for this process is problem 6, which requires 2 arguments. The second argument should be the dimension of the desired square matrix. An example for the execution of problem 6 would then be:

```
./main P6 10
```

## Ploting results
Given that all problems were properly executed and the .txt files are in the local directory, the plots for problem 5 and 6 can be saved by running

```
python plot.py PN
```

where, again, N is the desired problem. For the plots, only the arguments 5 and 6 are possible.


