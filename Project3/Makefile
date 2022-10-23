<<<<<<< HEAD
compile_two: 
	g++ main_two_particles.cpp ./src/particle.cpp ./include/Particle.hpp ./src/penningtrap.cpp ./include/PenningTrap.hpp -larmadillo -o main

compile_one:
	g++ main_one_particle.cpp ./src/particle.cpp ./include/Particle.hpp ./src/penningtrap.cpp ./include/PenningTrap.hpp -larmadillo -o main
=======
all: 
	g++ main.cpp ./src/particle.cpp ./src/penningtrap.cpp -std=c++11 -larmadillo -o main
>>>>>>> 84745ab42bfbcda10086fca4a7f2c131e3bd7286

run:
	./main $(method) $(n) $(t) $(interaction)

all: compile run
