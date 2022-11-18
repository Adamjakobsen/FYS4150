compile: 
	g++ ising_model.cpp -o main -std=c++11 -larmadillo -I/opt/homebrew/opt/libomp/include -L/opt/homebrew/opt/libomp/lib -lomp -O2
compileomp: 
	CPATH=/opt/homebrew/include LIBRARY_PATH=/opt/homebrew/lib g++-12 ising_model.cpp -std=c++11 -larmadillo  -fopenmp -o main -O2 

compilegrid:
	CPATH=/opt/homebrew/include LIBRARY_PATH=/opt/homebrew/lib g++-12 output_big_grid_state.cpp -std=c++11 -larmadillo  -fopenmp -o main -O2 

run:
	./main $(L) $(mc_cycles) $(burn_pct) $(lower_temp) $(upper_temp) $(temp_step) $(align) 
