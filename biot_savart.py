"""
Created by: Brandon Goddard
Adapted from: Antonio Franco
"""
import numpy as np


class BiotSolver(object):
    """
    Calculates the magnetic field vector solution for an arbitrary geometry.
    """

    def solve(self, wire, points, dL):
        """
        Solves the B field at points points for a given wire wire discretized uniformly
        at a length dL

        Args:
            wire (Coil): A coil object containing the coordinates of the wire vertices.
            points (ndarray): All points at which to calculate the magnetic field.
            dL (float): The length of the discretized segments (m).
        
        Returns:
            (ndarray): The magnetic field solution at all points.
        """
        # Initialize magnetic field.
        B = np.zeros((len(points[0]), 3), dtype=np.complex_)
        
        # Pre-compute discretized segments for all wire segments.
        num_segments = int(len(wire.coordz[0])) - 1
        segments = np.zeros((num_segments, len(points[0]), 3), dtype=np.complex_)
        for i in range(num_segments):
            segments[i] = self.__discretize(wire, dL, i)
            # Calculate magnetic field vector at each point.
            for p in range (0, int(len(points[0]))):
                B += self.__solve_segment(segments[i], points[:,p], wire.I)

        return B


    def __solve_segment(self, d_segment, point, complex_current):
        """
        Solves the magnetic field at a given point for a given discretized segment.

        Args:
        d_segment (ndarray): Discretized segment for which to calculate the magnetic field.
        point (ndarray): Point at which to calculate the magnetic field.
        complex_current (complex): Complex current running through the segment.

        Returns:
        (ndarray): The calculated magnetic field at the given point due to the given discretized segment.
        """
        B = np.zeros((1,3))
        for i in range (0, int(len(d_segment[0]))-1):
            dL = [np.array([d_segment[0][i],d_segment[1][i],d_segment[2][i]]),
                  np.array([d_segment[0][i+1],d_segment[1][i+1],d_segment[2][i+1]])]
            B += self.__db_scale(dL, point, complex_current)

        return B
    

    def __db_scale(self, dL, point, complex_current):
        """
        Compute the magnetic field (dB) produced by a wire segment carrying a complex current.

        Args:
            dL (float): A 2D array containing the coordinates of the wire segment.
            point (ndarray): The coordinates of the point where the magnetic field is to be calculated.
            complex_current (complex): The complex current flowing through the wire segment.

        Returns:
            (ndarray) The magnetic field vector at the specified calculation point.
        """
        # Magnetic permittivity of the vacuum
        mu = 4 * np.pi * 1e-7

        # Compute the sub-segment vector and the reference point in the middle of the sub-segment
        vdL = dL[1] - dL[0]
        rp = (dL[1] + dL[0]) / 2

        # Compute the displacement vector
        r = point - rp

        # Compute the contribution of the sub-segment to the magnetic field using the Biot-Savart law
        r3 = np.linalg.norm(r)**3
        integrand = np.cross(vdL, r) / r3
        dB = mu / (4 * np.pi) * complex_current * integrand

        return dB
    
        
    def __discretize(self, wire, dL, n):
        """
        Discretize the given N-segment wire into smaller segments of length dL.

        Args:
            wire (Coil): A coil object containing the coordinates of the wire vertices.
            dL (float): The length of the discretized segments (m).
            n (int): The index of the wire segment to be discretized.

        Returns:
            (ndarray): An array containing the coordinates of the discretized wire segment.
        """
        # Check if n is valid.
        if len(wire.coordz[0]) < n + 1:
            return []

        # Get the segment coordinates.
        x1, y1, z1 = wire.coordz[0][n], wire.coordz[1][n], wire.coordz[2][n]
        x2, y2, z2 = wire.coordz[0][n+1], wire.coordz[1][n+1], wire.coordz[2][n+1]
        segment = np.array([[x1, x2], [y1, y2], [z1, z2]])

        # Calculate the segment vector and length.
        u = segment[:,1] - segment[:,0]
        segment_length = np.linalg.norm(u)

        # Return the segment if it is shorter than dL.
        if segment_length < dL:
            return segment

        # Calculate the number of segments needed.
        n_segments = int(np.floor(segment_length / dL))

        # Calculate the length vector.
        v = dL * u / segment_length

        # Discretize the segment.
        d_segment = segment[:,0][:,None] + v * np.arange(1, n_segments+1)

        # Check if the final segment needs to be added.
        if n_segments != segment_length / dL:
            d_segment = np.append(d_segment, segment[:,1][:,None], axis=1)

        return d_segment
