## Simulating time-dependent two dimensional Schrödinger equation

We will be using our simulation to study the double-slit experiment. 


# compile

 ```
 make compile
 ```
 
 # Run
 Produce data for one problem at a time:
 ```
 make run_"problem"
 ```
where problem can be `P7_1`, `P7_2`, `P8`, `P9_single`, `P9_double` or `P9_triple`.

or simply produce all the data at one:

```
 make run_all
 ```
 
 # plot
 
 For plotting run:
 
 ```
 python3 plot.py "problem"
 ```
