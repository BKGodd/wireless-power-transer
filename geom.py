"""
Created by: Brandon Goddard
Adapted from: Antonio Franco
"""
import numpy as np
from biot_savart import BiotSolver
import scipy as sc
import pylab as p
import mpl_toolkits.mplot3d.axes3d as p3


class Coil:
    """Creates an arbitrary shaped and orientated conductive coil."""

    def __init__(self, center, radius, num_points, orientation='yz'):
        """
        Constructor for the coil.

        Args:
            center (list): The center of the circle, given as a 3-element vector
                            containing the x, y, and z coordinates.
            radius (float): The radius of the circle.
            num_points (int): The number of points on the circle.
            Orientation (str, optional): The orientation of the circle, given as
                                        'xy', 'xz', or 'yz'. Default: 'yz'.

        Returns:
            None
        """
        # Initialize the complex current of the coil.
        self.I = complex(1, 0)
        # Initialize the parameters of the coil.
        self.center = center
        self.radius = radius
        self.num_points = num_points
        self.orientation = orientation
        # Initialize the coordinates that define the coil.
        self.init_coil()
        
        return
    

    def init_coil(self):
        """
        Create a coil with the specified center, radius, and orientation.

        Returns:
            None
        """
        # Parameterize the points on the circle.
        t = np.linspace(0, 2*np.pi, self.num_points)

        # Check the orientation of the circle in 3D space.
        if self.orientation == 'xy':
            X = self.center[0] + self.radius * np.sin(t)
            Y = self.center[1] + self.radius * np.cos(t)
            Z = np.zeros(self.num_points)
        elif self.orientation == 'xz':
            X = self.center[0] + self.radius * np.sin(t)
            Z = self.center[1] + self.radius * np.cos(t)
            Y = np.zeros(self.num_points)
        elif self.orientation == 'yz':
            X = np.zeros(self.num_points)
            Y = self.center[1] + self.radius * np.sin(t)
            Z = self.center[2] + self.radius * np.cos(t)

        self.coordz = np.array(X, Y, Z)

        return


    def calc_inductance(self, mutual_coil, disp=False):
        """
        Calculates the total magnetic flux and inductance of a coil.

        Args:
            mutual_coil (Coil): Coil object to calculate the mutual inductance.
            disp (bool): Whether to display/plot the results of integration or not.

        Returns:
            inductance (float): The inductance of the loop.
        """
        # Create points that form a line from the center to the radius.
        points = np.linspace(0, mutual_coil.radius, mutual_coil.num_points)
        x = np.ones(len(points)) * mutual_coil.center[0]
        y = np.ones(len(points)) * mutual_coil.center[1]
        z = (np.ones(len(points)) * mutual_coil.center[2]) + points
        integral_points = np.array([x,y,z])

        # Find the total magnetic field solution for all points (norm of the vector).
        Solver = BiotSolver()
        b = Solver.Solve(self, integral_points, 0.1)
        b_norm = np.linalg.norm(b, axis=1)
        
        # Account for error of width of wire, only integrate up to the max bm value (only for self inductance).
        if self == mutual_coil:
            max_index = np.where(b_norm==max(b_norm))[0][0]
            b_norm = b_norm[:max_index+1]
            points = points[:len(b_norm)]

        # Execute the intergration to find the flux and inductance.
        current = np.linalg.norm(self.I)
        area_circle = 2*np.pi*points
        integrand = area_circle * b_norm
        flux = sc.integrate.simps(integrand, points)
        inductance = flux/current

        # Display results, if desired
        if disp:
            fig = p.figure(None)
            ax = p3.Axes3D(fig)
            # Draw the two coils in question (only one if self inductance).
            ax.plot(mutual_coil.coordz[0], mutual_coil.coordz[1], mutual_coil.coordz[2], 
                    label='Coil', linewidth=2, c='k')
            if self == mutual_coil:
                ax.plot(self.coordz[0], self.coordz[1], self.coordz[2], 
                        label='Coil', linewidth=2, c='k')
            
            # Draw the magnetic vector field solution (B).
            ax.quiver(x, y, z, b[:,0], b[:,1], b[:,2], color='b', length=0.3, normalize=True)
            # Set axis labels
            ax.set_xlabel('X')
            ax.set_ylabel('Y')
            ax.set_zlabel('Z')
            p.draw()
            p.show()

        return inductance
