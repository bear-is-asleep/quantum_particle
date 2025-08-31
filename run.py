import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from functools import partial
import json
import argparse
import os
import time

# Set up argument parser
parser = argparse.ArgumentParser(description='Wavefunction Animation')
parser.add_argument('-s','--save', dest='save', action='store_true', default=False, help='Set this flag to not save the animation')
parser.add_argument('-f','--fname', dest='fname', type=str, default='wavefunction', help='Name of the animation file')
parser.add_argument('-c','--config', dest='config', type=str, default='test1.json', help='Name of the configuration file')

args = parser.parse_args()

# Load data from JSON file
with open(args.config, 'r') as f:
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

#Works
def evolve_crank_nicholson(psi, dt, dx, V):
    """
    Evolves the wavefunction in time using the Crank-Nicholson method.

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
    psi_new: np.array
        1D array of the evolved wavefunction
    """
    global steps
    t0 = time.time()
    # Number of spatial points
    N = len(psi)

    # Calculate the pre-factor
    lambda_coeff = hbar* dt / (2 * m * dx**2)

    # Construct matrix A for L * psi_next = R * psi_current
    A = np.zeros((N, N), dtype=complex)
    B = np.zeros((N, N), dtype=complex)

    # Populate matrix A and matrix B
    for i in range(N):
        # Construct the main diagonal of A
        A[i, i] = 2 + 2j * lambda_coeff + 2j * lambda_coeff * V[i] * dx**2
        B[i, i] = 2 - 2j * lambda_coeff - 2j * lambda_coeff * V[i] * dx**2

        # Construct the sub-diagonal (-i lambda) for A
        if i > 0:
            A[i, i-1] = -1j * lambda_coeff
            B[i, i-1] = 1j * lambda_coeff

        # Construct the super-diagonal (-i lambda) for A
        if i < N - 1:
            A[i, i+1] = -1j * lambda_coeff
            B[i, i+1] = 1j * lambda_coeff

    t1 = time.time()

    # Calculate the right hand side B @ psi
    b = B @ psi

    t2 = time.time()

    # Solve A @ psi_new = b
    psi_new = np.linalg.solve(A, b)
    
    t3 = time.time()
    
    if steps % 100 == 0:
        print(f'-cn step {steps} psi - psi_new: {np.sum(np.abs(psi - psi_new)):.2e}')
        print(f'-Construction of A and B: {t1-t0:.2e} s')
        print(f'-Calculation of b: {t2-t1:.2e} s')
        print(f'-Solution of A @ psi_new = b: {t3-t2:.2e} s')

    return psi_new

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

steps = 0
def evolve(psi, dt, dx, V, method='crank_nicholson'):
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
    method: str
        Method to use for time evolution ("euler" or "crank_nicholson")
    ----------------------
    Returns:
    psi: np.array
        1D array of the evolved wavefunction
    """
    global steps
    _psi = np.copy(psi)
    if method == 'euler':
        # Compute the Laplacian
        laplacian = compute_laplacian(psi, dx)

        # Evolve the wavefunction
        psi += -1j * dt * (hbar / (2 * m) * laplacian + V * psi)
    elif method == 'crank_nicholson':
        psi = evolve_crank_nicholson(psi, dt, dx, V)
    else:
        raise ValueError(f'Method "{method}" not recognized. Use "euler" or "crank_nicholson".')

    # Normalize the wavefunction
    psi = normalize(psi, dx)

    #for i,x in enumerate(grid):
    #   print(f'grid: {x:.2e}, psi: {np.real(psi[i]):.2e}, prob: {compute_probability(psi)[i]:.2e}')
    #print(f'psi: {psi}')
    
    if steps % 100 == 0:
        qqq = 0
        #print(f'- step {steps} psi - _psi: {np.sum(np.abs(psi - _psi)):.2e}')

    steps += 1
    return psi

def run_simulation(psi, dt, dx, n, V, method='crank_nicholson'):
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
    method: str
        Method to use for time evolution ("euler" or "crank_nicholson")
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
        psi = evolve(psi, dt, dx, V, method)

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
    t0 = time.time()
    psi, probability = run_simulation(psi, dt, dx, nsteps, V, method=evolve_method)
    t1 = time.time()
    line.set_ydata(probability)
    t2 = time.time()
    
    #print(f'Evolution of the wavefunction: {t1-t0:.2e} s')
    #print(f'Update of the plot: {t2-t1:.2e} s')
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