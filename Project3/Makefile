compile_two: 
	g++ main_two_particles.cpp ./src/particle.cpp ./include/Particle.hpp ./src/penningtrap.cpp ./include/PenningTrap.hpp -larmadillo -o main

compile_one:
	g++ main_one_particle.cpp ./src/particle.cpp ./include/Particle.hpp ./src/penningtrap.cpp ./include/PenningTrap.hpp -larmadillo -o main

run:
	./main $(method) $(n) $(t) $(interaction)

all: compile run
