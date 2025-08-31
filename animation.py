import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from functools import partial

class Animation:
    def __init__(self, interval, nsteps, simulation, fps=30):
        """
        Initializes the animation.
        ----------
        Parameters:
        interval : int
            Time interval between frames in milliseconds.
        nsteps : int
            Number of time steps to evolve.
        simulation : Simulation
            Simulation object.
        fps : int
            Frames per second.
        """
        self.interval = interval
        self.nsteps = nsteps
        self.simulation = simulation
        self.init_figure()
    
    def init_figure(self):
        """
        Initializes the figure.
        """
        # Prepare the figure
        fig, self.ax = plt.subplots()
        self.ax.set_xlim(-simulation.L/2, simulation.L/2)
        self.ax.set_xlabel('x')
        self.ax.plot(simulation.grid,simulation.V, label='V(x)', color = 'red')
        
        max_y = max([psi.compute_probability() for psi in simulation.wavefunctions])
        self.ax.set_ylim(0, 1.1*max_y)
        
        self.lines = []
        for psi in self.simulation.wavefunctions:
            line, = self.ax.plot(self.simulation.grid, psi.compute_probability(), label=f'psi_{psi.id}')
            self.lines.append(line)
    
    def animate(self, i):
        """
        Animates the wavefunction.
        ----------
        Parameters:
        i : int
            Frame number.
        """
        # Evolve the wavefunction nsteps
        self.simulation.evolve(self.nsteps)
        
        # Set the data for each line
        for i, psi in enumerate(self.simulation.wavefunctions):
            self.lines[i].set_ydata(psi.compute_probability())
        
        return self.lines
    
    def make_animation(self):
        """
        Makes the animation.
        """
        self.ani = FuncAnimation(self.ax.figure, self.animate, frames=range(0, self.nsteps), interval=self.interval, blit=True)
    
    def run(self):
        """
        Runs the animation.
        """
        self.make_animation()
        plt.show()
    
    def save(self, filename):
        """
        Saves the animation to a file.
        ----------
        Parameters:
        filename : str
            Name of the file.
        """
        self.make_animation()
        self.ani.save(filename, writer='ffmpeg', fps=self.fps)
        
        