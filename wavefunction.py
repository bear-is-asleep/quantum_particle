import numpy as np

#FIXME: Set physics constants globally
hbar = 1.0
m = 1.0

class Wavefunction:
    def __init__(self, grid):
        """
        Initializes the wavefunction as a Gaussian wavepacket.
        ----------
        Parameters:
        grid : ndarray
            Grid points.
        """
        self.dx = grid[1] - grid[0] # Grid spacing
        self.psi = self.init_wavefunction(grid)
        
        #Bookkeeping
        self.steps = 0
        #assert (self.compute_norm() > 0.0) & (self.compute_norm() < 1e3), "Wavefunction is not normalized."
    
    def init_wavefunction(self, grid):
        raise NotImplementedError("Wavefunction initialization method not implemented.")

    def compute_norm(self):
        """
        Computes the norm of the wavefunction.
        ----------
        Returns:
        norm : float
            Norm of the wavefunction.
        """
        return np.sqrt(np.sum(self.compute_probability(self.psi) * self.dx))

    def normalize(self):
        """
        Normalizes the wavefunction.
        """
        self.psi = self.psi/self.compute_norm()
        return self.psi

    @staticmethod
    def compute_probability(psi):
        """
        Computes the probability density of the wavefunction.
        """
        return np.real(psi * np.conj(psi))
    
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
    
    def evolve(dt, V, method='crank_nicholson'):
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
        _psi = np.copy(self.psi)
        if method == 'euler':
            # Compute the Laplacian
            laplacian = compute_laplacian(self.psi, dx)

            # Evolve the wavefunction
            self.psi += -1j * dt * (hbar / (2 * m) * laplacian + V * self.psi)
        elif method == 'crank_nicholson':
            self.psi = evolve_crank_nicholson(self.psi, dt, dx, V)
        else:
            raise ValueError(f'Method "{method}" not recognized. Use "euler" or "crank_nicholson".')

        # Normalize the wavefunction
        self.normalize()

        #for i,x in enumerate(grid):
        #   print(f'grid: {x:.2e}, psi: {np.real(psi[i]):.2e}, prob: {compute_probability(psi)[i]:.2e}')
        #print(f'psi: {psi}')
        
        if self.steps % 100 == 0:
            qqq = 0
            #print(f'- step {self.steps} psi - _psi: {np.sum(np.abs(psi - _psi)):.2e}')

        self.steps += 1
        return psi

class GausWavefunction(Wavefunction):
    def __init__(self, grid, x0, k0, sigma, epsilon):
        """
        Initializes the wavefunction as a Gaussian wavepacket.
        ----------
        Parameters:
        grid : ndarray
            Grid points.
        x0 : float
            Initial position of the wavepacket.
        k0 : float
            Initial wave number of the wavepacket.
        sigma : float
            Width of the wavepacket.
        epsilon : float
            Threshold for the wavepacket amplitude at the boundaries.
        """
        #Initalize variables first
        self.x0 = x0
        self.k0 = k0
        self.sigma = sigma
        self.epsilon = epsilon
        super().__init__(grid)

    def init_gaus_packet(self, x, x0, k0, sigma, epsilon):
        """Initializes a Gaussian wavepacket.
        ----------
        Parameters:
        x : ndarray
            Grid points.
        x0 : float
            Initial position of the wavepacket.
        k0 : float
            Initial wave number of the wavepacket.
        sigma : float
            Width of the wavepacket.
        epsilon : float
            Threshold for the wavepacket amplitude at the boundaries
        ----------  
        Returns:
        psi : ndarray
            Wavefunction.
        """
        self.psi = np.exp(-0.5 * ((x - x0) / sigma)**2 + 1j * k0 * x)
        self.psi = self.normalize()
        #Set end points to zero (assuming periodic boundary conditions)
        if self.psi[0] < epsilon and self.psi[-1] < epsilon:
            self.psi[0] = 0.0
            self.psi[-1] = 0.0
        else:
            print(f'Wavefunction boundaries are larger than epsilon: {epsilon}. Reducing the width until the boundaries are zero.')
            iterations = 0
            while psi[0] > epsilon or self.psi[-1] > epsilon:
                sigma -= 0.01 * self.dx
                self.psi = np.exp(-0.5 * ((x - x0) / sigma)**2 + 1j * k0 * x)
                self.psi = self.normalize(x)
                iterations += 1
            print(f'Wavefunction boundaries are now zero with sigma = {sigma} ({iterations} iterations).')
            self.psi[0] = self.psi[-1] = 0.0
        return psi

    def init_wavefunction(self, grid):
        return self.init_gaus_packet(grid, self.x0, self.k0, self.sigma)