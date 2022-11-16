compile: 
	g++ toymodel.cpp -o main -std=c++11 -larmadillo -I/opt/homebrew/opt/libomp/include -L/opt/homebrew/opt/libomp/lib -lomp -O2
	
run:
	./main $(L) $(mc_cycles) $(burn_pct) $(lower_temp) $(upper_temp) $(temp_step) $(align)
