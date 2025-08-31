import json
import argparse
import os
import numpy as np

from animation import Animation
from wavefunction import GausWavefunction
from simulation import Simulation

# Set up argument parser
parser = argparse.ArgumentParser(description='Wavefunction Animation')
parser.add_argument('-s','--save', dest='save', action='store_true', default=False, help='Set this flag to not save the animation')
parser.add_argument('-f','--fname', dest='fname', type=str, default='wavefunction', help='Name of the animation file')
parser.add_argument('-c','--config', dest='config', type=str, default='test1.json', help='Name of the configuration file')

args = parser.parse_args()

# Load data from JSON file
with open(args.config, 'r') as f:
    data = json.load(f)

# Physics - FIXME: Use these
hbar = data["constants"]["hbar"]
m = data["wavefunction"]["m"]

# Wavefunction
x0 = data["wavefunction"]["x0"]
k0 = data["wavefunction"]["k0"]
sigma = data["wavefunction"]["sigma"]
epsilon = data["wavefunction"]["epsilon"]

# Simulation
dt = data["simulation"]["dt"]
dx = data["simulation"]["dx"]
N_t = data["simulation"]["N_t"]
L = data["simulation"]["L"]
evolve_method = data["simulation"]["evolve_method"]
grid = np.arange(-L / 2, L / 2, dx)

# Animation
interval = data["animation"]["interval"]
N_steps = data["animation"]["N_steps"]
fname = args.fname
folder = 'animations/'
if os.path.exists(folder) == False:
    os.mkdir(folder)


# Set up potential energy
V = np.zeros_like(grid)
vmode = data["potential_energy"]["mode"]
if vmode == "delta":
    V[int(len(grid) / 2)] = 1000.0  # Delta function potential
else:
    print(f'Potential energy mode "{vmode}" not recognized. Using default V(x) = 0.')

if __name__ == "__main__":
    # Initialize the wavefunction
    psi = GausWavefunction(grid, x0, k0, sigma, epsilon)
    
    # Initialize the simulation
    simulation = Simulation(dt, [psi], V, L, grid, evolve_method)
    
    # Initialize the animation
    ani = Animation(interval, N_t, simulation)

    
    if args.save:
        ani.save(folder + fname + '.mp4')
    else:
        ani.run()