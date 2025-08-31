

class Simulation:
    def __init__(self, dt, wavefunctions, V, L, grid, method):
        """
        Initializes the simulation.
        ----------
        Parameters:
        dt : float
            Time step.
        wavefunctions : list
            List of wavefunctions.
        V : ndarray
            Potential energy.
        L : float
            Length of the infinite square well.
        grid : ndarray
            Grid points.
        method : str
            Method to evolve the wavefunction in time ("euler" or "crank_nicholson").
        """
        self.dt = dt
        self.wavefunctions = wavefunctions
        self.V = V
        self.L = L
        self.grid = grid
        self.method = method
        # Bookkeeping
        self.steps = 0
        
        def evolve(self,n):
            """
            Evolves the wavefunction in time.
            ----------
            Parameters:
            n : int
                Number of time steps to evolve.
            """
            for i in range(n):
                for wavefunction in self.wavefunctions:
                    wavefunction.evolve(self.dt, self.V, self.method)
                self.steps += 1
            