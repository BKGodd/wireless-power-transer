'''
Created by: Brandon Goddard
'''
import shapes
import numpy as np
import multiprocessing as mp


def get_inductance(coil_one, coil_two):
    """
    Calculate the mutual inductance between two coils.
    
    Args:
        coil_one (Coil): The first coil object.
        coil_two (Coil): The second coil object.

    Returns:
        (float): The mutual inductance between the two coils.
    """
    return coil_one.calc_inductance(coil_two)


def write_results(csv, all_inductances, coil2, coil3):
    """
    Write all inductance calculation results to a CSV file.

    Args:
        csv (File): A CSV file object
        all_inductances (list(tuple)): A list of three pairs of inductances.
        coil2 (Coil): Third coil in the WPT system.
        coil3 (Coil): Second coil in the WPT system.

    Returns:
        None
    """
    # Save results to the csv file.
    # Note: M43 = M21
    m21, m32, s22 = round(all_inductances[0], 5), \
                    round(all_inductances[1], 5), \
                    round(all_inductances[2], 5)
    coil2_dist = round(coil2.center[0], 5)
    coil3_dist = round(coil3.center[0], 5)
    coil23_rad = round(coil2.radius, 5)
    csv.write(f'{coil2_dist}, {coil3_dist}, {coil23_rad}, {m21}, {m32}, {s22}\n')


def coil_purmutations(csv):
    """
    Generate a mesh of possible radii and positions for two intermediate 
    coils in a system of four coils, and calculate the mutual inductances 
    between the four coils for each combination of intermediate coil radii and positions.

    Args:
        csv (File): A CSV file object where the results of the calculations will be saved.

    Returns:
        None.
    """
    num_points = 500
    # Create transmit (Tx) coil.
    coil1 = shapes.Coil([0,0,0], 0.2, num_points)
    # create receive (Rx) coil (not needed due to symmetry).
    #coil4 = shapes.Coil([2,0,0], 0.2, num_points)

    # Create a mesh of all possible radii and positions of coils.
    loop_radii = np.linspace(0.2167, 0.55, 7)
    loop_pos = np.linspace(0.05, 0.975, 19)
    radii_grid, pos_grid = np.meshgrid(loop_radii, loop_pos)

    # Iterate over all possible combinations of varying radii and distances.
    # Maintaining symmetry reduces comp. complexity in this case.
    with mp.Pool() as pool:
        for i in range(radii_grid.shape[0]):
            for j in range(radii_grid.shape[1]):
                # Create intermediate coils for lilypad technique of magnetic coupling.
                coil2 = shapes.Coil([pos_grid[i, j], 0, 0], radii_grid[i, j], num_points)
                coil3 = shapes.Coil([2-pos_grid[i, j], 0, 0], radii_grid[i, j], num_points)

                # Calculate all mutual inductances in parallel (m12=m21 , m13=m31 , etc.).
                coils = [(coil1, coil2), (coil2, coil3), (coil2, coil2)]
                all_inductances = pool.map(get_inductance, coils)

                # Save results to CSV file for reference.
                write_results(csv, all_inductances, coils)


def main(csv_filepath):
    """
    Entry point for the program. Executes the coil permutations 
    and saves the results to a CSV file.

    Args:
        csv_filepath (str): The file path of the CSV file where 
                            the results will be saved.

    Returns:
        None.
    """
    # Execute all coil permutations, saving results to file.
    with open(csv_filepath, 'a+') as csv:
        csv.write('coil2_x, coil3_x, coil23_rad, m21, m32, s22\n')
        coil_purmutations(csv)
