# plot a 2d grid of the L x L spins
from time import sleep 
import matplotlib.pyplot as plt
import random
import numpy as np

def init_random_config(L):
    """
    L is the length of the square lattice
    """
    config = np.zeros((L, L))
    for i in range(L):
        for j in range(L):
            config[i, j] = random.choice([+1, -1])
    return config

def find_interacting_pairs(index, raveled_config, L):
    side_neighbours = [(index - 1) % L + L*(index//L), (index + 1) % L+ L*(index//L)]
    vert_neighbours = [(index - L)%(L*L), (index + L)%(L*L)]
    neibours = side_neighbours + vert_neighbours
    interacting_pairs = [(raveled_config[index], raveled_config[x]) for x in neibours]
    return interacting_pairs


# evolve the system given the initial configuration and boltzmann probability
def evolve(config, beta):
    L = len(config)
    raveled_config = config.ravel()
    for i in range(L*L):
        interacting_pairs = find_interacting_pairs(i, raveled_config, L)
        energy = 0
        for a,b in interacting_pairs:
            energy += -(a*b)

        delta_E = -2*energy # because delta_E = E_new - E_old = -((-a)*b) - (-(a*b)) = 2*(a*b) = -2*energy
        if delta_E < 0: # if the energy is lower, accept the new configuration
            raveled_config[i] *= -1
        else:
            if random.random() < np.exp(-beta*delta_E): # if the energy is higher, accept the new configuration with a probability
                raveled_config[i] *= -1
    return raveled_config.reshape(L,L)

def find_energy_given_array_of_config(config, L):
    """
    notice we are taking into account the periodic boundary conditions
    """
    energy = 0
    raveled_config = config.ravel()
    ## this is dumb now but will be useful later when we select only one term to flip

    for i in range(len(raveled_config)):
        ## i will be the index of the term. k,l will be the coordinates of the term

        interacting_pairs = find_interacting_pairs(i, raveled_config, L)
        for a, b in interacting_pairs:
            energy += a*b
    return (-1/2)*energy ## notice the minus sign and that we count each pair twice




def plot_config(L, init=False, epochs=None):
    """
    config is a 2d array of the L x L spins
    """
    if init: # else we are evolving the system
        config = init_random_config(L)
        energy = find_energy_given_array_of_config(config, L)
        plt.imshow(config, cmap = 'Purples', interpolation = 'nearest')
        plt.title(f'Energy = {energy}')

        plt.show()
    else:
        config = init_random_config(L)

        kb = 1 ## only for testing
        T = 1
        beta = 1/(kb*T)
        plt.ion()
        plt.figure()
        plt.show()
        for i in range(epochs):
            config = evolve(config, beta)
            energy = find_energy_given_array_of_config(config, L)
            plt.imshow(config, cmap = 'Purples', interpolation = 'nearest')
            plt.title(f'Energy = {energy}')
            # animate the plot wit
            plt.draw()

            plt.pause(0.1)

if __name__ == '__main__':
    L = 100
    epochs = 100
    plot_config(L, init=False, epochs=epochs)