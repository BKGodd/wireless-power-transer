"""
Created by: Brandon Goddard
"""
from csv import DictReader
import mpl_toolkits.mplot3d.axes3d as p3
import pylab as p
import numpy as np


def calc_coefficients(csv, disp=False):
    """
    Calculate the coupling coefficients and output voltage of the coil WPT system.

    Args:
        csv (File): File object for reading the data in.
        disp (bool): Whether to plot and show the results.

    Returns:
        None
    """
    # Create a dictionary to store the data for each column.
    data = {}
    # Pre-computed self inductance of coils 1 (Tx) and 4 (Rx).
    s11 = 1.3738219953677048e-06
    # Set frequency of 10MHz
    w = 2*np.pi*10.0e6

    # Gather all data from the CSV file.
    reader = DictReader(csv, delimiter=',', skipinitialspace=True)
    for row in reader:
        # Collect each parameter in the new dictionary.
        for key, value in row.items():
            if key not in data:
                data[key] = np.array([float(value)])
            else:
                data[key] = np.append(data[key], float(value))

    # Calculate all coupling coefficients (K).
    k21 = k43 = data['m21']/(np.sqrt(s11 * data['s22']))
    k32 = data['m32']/(np.sqrt(data['s22'] * data['s22']))

    # Calculate the output voltage.
    factor1 = ((50*w) * (data['m21']*data['m32'])) * data['m21']
    factor2 = (data['m21']**2) * (data['m21']**2)
    factor3 = factor2 * w**2
    factor4 = 2500 * data['m32']**2
    vout = 2 * (factor1/(factor3 + factor4))

    # Plot results, if desired.
    if disp:
        fig = p.figure(None)
        ax = p3.Axes3D(fig)
        ax.scatter(data['coil3_x']-data['coil2_x'], data['coil23_rad'], vout, color='b')
        ax.set_xlabel('Distance (m)')
        ax.set_ylabel('Radius (m)')
        ax.set_zlabel('Vout (V)')
        p.draw()
        p.show()

    return k21, k32
    

def main(csv_filepath):
    """
    Entry point for the program. Calculates the coupling coefficients
    between all coils in the system.

    Args:
        csv_filepath (str): The file path of the CSV file where 
                            the results are saved to read from.

    Returns:
        None.
    """
    # Execute all coil permutations, saving results to file.
    with open(csv_filepath, 'r') as csv:
        # Use coupling coefficients however you need them.
        k21, k32 = calc_coefficients(csv)