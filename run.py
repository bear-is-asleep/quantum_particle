import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from functools import partial
import json
import argparse

# Set up argument parser
parser = argparse.ArgumentParser(description='Wavefunction Animation')
parser.add_argument('--save', dest='save', default=False, help='Set this flag to not save the animation')

args = parser.parse_args()

# Load data from JSON file
with open('config.json', 'r') as f:
    data = json.load(f)

# Physics
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
grid = np.arange(-L / 2, L / 2, dx)

# Animation
interval = data["animation"]["interval"]
N_steps = data["animation"]["N_steps"]
fname = data["animation"]["fname"]
folder = 'animations/'


# Set up potential energy
V = np.zeros_like(grid)
if data["potential_energy"]["mode"] == "delta":
    V[int(len(grid) / 2)] = 1000.0  # Delta function potential

# Prepare the figure
fig, ax = plt.subplots()
ax.set_xlim(-L/2, L/2)
ax.set_xlabel('x')
ax.plot(grid, V, label='V(x)', color = 'red')

def initialize_gaussian_wavepacket(x, x0, k0, sigma):
    """
    Returns a Gaussian wavepacket
    ----------------------
    Parameters:
    x: float
        Spatial coordinate
    x0: float
        Initial position of the wavepacket
    k0: float
        Initial wavevector of the wavepacket
    sigma: float
        Width of the wavepacket
    ----------------------
    Returns:
    psi: complex
        Value of the wavefunction at x
    """
    return np.exp(-0.5 * ((x - x0) / sigma) ** 2 + 1j * k0 * x)

def initialize_wavepacket(grid, x0, k0, sigma):
    """
    Initializes the wavefunction for a Gaussian wavepacket
    ----------------------
    Parameters:
    grid: np.array
        1D array of spatial coordinates
    x0: float
        Initial position of the wavepacket
    k0: float
        Initial wavevector of the wavepacket
    sigma: float
        Width of the wavepacket
    ----------------------
    Returns:
    psi: np.array
        1D array of the wavefunction
    """
    # Gaussian wavepacket
    psi = initialize_gaussian_wavepacket(grid, x0, k0, sigma)
    
    # Normalize the wavefunction
    psi = normalize(psi, dx)
    
    #Set end points to zero (assuming periodic boundary conditions)
    if psi[0] < epsilon and psi[-1] < epsilon:
        psi[0] = 0.0
        psi[-1] = 0.0
    else:
        print(f'Wavefunction boundaries are larger than epsilon: {epsilon}. Reducing the width until the boundaries are zero.')
        iterations = 0
        while psi[0] > epsilon or psi[-1] > epsilon:
            sigma -= 0.01 * np.abs(grid[0] - grid[-1])  # Reduce the width by 1% of the box size
            psi = initialize_gaussian_wavepacket(grid, x0, k0, sigma)
            psi = normalize(psi, dx)  # Normalize the wavefunction
            iterations += 1
        print(f'Wavefunction boundaries are now zero with sigma = {sigma} ({iterations} iterations).')
        psi[0] = psi[-1] = 0.0
    
    return psi

def normalize(psi, dx):
    """
    Normalizes the wavefunction
    ----------------------
    Parameters:
    psi: np.array
        1D array of the wavefunction
    dx: float
        Spacing of the spatial grid
    ----------------------
    Returns:
    psi: np.array
        1D array of the normalized wavefunction
    """

    # Compute the norm
    norm = np.sqrt(np.sum(compute_probability(psi) * dx))

    # Normalize the wavefunction
    psi /= norm

    return psi

def compute_probability(psi):
    """
    Computes the probability density from the wavefunction
    ----------------------
    Parameters:
    psi: np.array
        1D array of the wavefunction
    ----------------------
    Returns:
    probability: np.array
        1D array of the probability density
    """
    return np.real(psi * np.conj(psi))

#Replace with Crank-Nicolson method
#https://people.sc.fsu.edu/~jpeterson/5-CrankNicolson.pdf
def compute_laplacian(psi, dx, pad=True, order=1):
    """
    Computes the Laplacian of the wavefunction
    ----------------------
    Parameters:
    psi: np.array
        1D array of the wavefunction
    dx: float
        Spacing of the spatial grid
    order: int
        Order of the finite difference scheme
    ----------------------
    Returns:
    laplacian: np.array
        1D array of the Laplacian of the wavefunction
    """

    laplacian = np.zeros_like(psi)
    # Compute the Laplacian using central differences
    if order == 1:
        laplacian[1:-1] = (psi[:-2] - 2 * psi[1:-1] + psi[2:]) / dx ** 2
        laplacian[0:1] = laplacian[-1:-2] = 0.0 # Set the boundaries to zero
    elif order == 2:
        laplacian[2:-2] = (-1/12 * psi[:-4] + 4/3 * psi[1:-3] - 5/2 * psi[2:-2] + 4/3 * psi[3:-1] - 1/12 * psi[4:]) / dx ** 2
        laplacian[0:2] = laplacian[-2:] = 0.0 # Set the boundaries to zero
    return laplacian

def evolve(psi, dt, dx, V):
    """
    Evolves the wavefunction in time using Euler's method
    ----------------------
    Parameters:
    psi: np.array
        1D array of the wavefunction
    dt: float
        Time step
    dx: float
        Spacing of the spatial grid
    V: np.array
        1D array of the potential energy
    ----------------------
    Returns:
    psi: np.array
        1D array of the evolved wavefunction
    """
    
    # Compute the Laplacian
    laplacian = compute_laplacian(psi, dx)

    # Evolve the wavefunction
    psi += -1j * dt * (hbar / (2 * m) * laplacian + V * psi)

    # Normalize the wavefunction
    psi = normalize(psi, dx)

    #for i,x in enumerate(grid):
    #   print(f'grid: {x:.2e}, psi: {np.real(psi[i]):.2e}, prob: {compute_probability(psi)[i]:.2e}')
    #print(f'psi: {psi}')

    return psi

def run_simulation(psi, dt, dx, n, V):
    """
    Runs the simulation for a given number of time steps
    ----------------------
    Parameters:
    psi: np.array
        1D array of the wavefunction
    dt: float
        Time step
    dx: float
        Spacing of the spatial grid
    n: int
        Number of time steps
    V: np.array
        1D array of the potential energy
    ----------------------
    Returns:
    psi: np.array
        1D array of the evolved wavefunction
    """

    # Normalize the wavefunction
    psi = normalize(psi, dx)

    # Compute the initial probability density
    probability = compute_probability(psi)

    # Evolve the wavefunction
    for i in range(n):
        psi = evolve(psi, dt, dx, V)

        # # Check the normalization
        # if i % 1 == 0:
        #     psi = normalize(psi, dx)

    # Compute the final probability density
    probability = compute_probability(psi)

    return psi, probability

def animate(i, nsteps):
    """
    Animates the wavefunction
    ----------------------
    Parameters:
    i: int
        Frame number
    nsteps: int
        Number of time steps to evolve the wavefunction
    ----------------------
    Returns:
    line: matplotlib.lines.Line2D
        Line object to animate
    """
    global psi, probability
    psi, probability = run_simulation(psi, dt, dx, nsteps, V)
    line.set_ydata(probability)
    return line,

if __name__ == "__main__":
    # Initialize the wavefunction
    psi = initialize_wavepacket(grid, x0, k0, sigma)
    ax.set_ylim(0, 1.1 * np.max(compute_probability(psi)))
    
    # Animate
    make_animation = partial(animate, nsteps=N_t)
    line, = ax.plot(grid, compute_probability(psi))
    ani = FuncAnimation(fig, make_animation, frames=N_steps, interval=interval)
    if args.save:
        ani.save(folder+fname+'.mp4', writer='ffmpeg', fps=30)
    else:
        plt.show()
    
 