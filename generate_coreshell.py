import os
import sys
import numpy as np
from ase import Atoms
from ase.io import read
from ase.data import atomic_masses, atomic_numbers

"""
This script reads a VASP structure file, generates a LAMMPS data file for the core-shell model.

Usage:
1. Run the script with:
   python generate_coreshell.py <input.vasp> <output.lammps>
2. Follow prompts to select core-shell elements and define atom order.

Email: Hao YANG (yanghao2023@shanghaitech.edu.cn) if any questions
Modified by Zihan YAN (yanzihan@westlake.edu.cn)
"""

def read_vasp_to_dict(filename):
    """
    Reads a VASP structure file and extracts information into a dictionary.
    Args:
        filename (str): Path to the VASP input structure file.
    Returns:
        dict: A dictionary containing box information, element counts, positions, and symbols.
    """
    atoms = read(filename, format='vasp')

    # Extract box dimensions, elements, and positions
    box_info = atoms.get_cell()  # Box information (cell matrix)
    symbols = atoms.get_chemical_symbols()  # List of element symbols
    unique_elements, counts = np.unique(symbols, return_counts=True)  # Count occurrences of each element
    positions = atoms.get_positions()  # Atomic positions

    # Store the extracted information in a dictionary
    structure_dict = {
        "box": box_info,
        "elements": dict(zip(unique_elements, counts)),
        "positions": positions,
        "symbols": symbols
    }

    return structure_dict

def add_coreshell_elements(structure_dict, elements):
    """
    Adds core and shell elements for specified elements in the structure.
    Core elements are a modification of the original element, and shell elements
    represent a surrounding shell of atoms.
    Args:
        structure_dict (dict): Dictionary containing atomic structure information.
        elements (list): List of elements to which core-shell modifications should be added.
    """
    symbols = structure_dict["symbols"]
    positions = structure_dict["positions"]

    for element in elements:
        # Find the indices of atoms corresponding to the selected element
        indices = [i for i, symbol in enumerate(symbols) if symbol == element]

        for index in indices:
            # Add core and shell atoms for each atom of the selected element
            symbols.append(f"{element}_core")
            symbols.append(f"{element}_shell")
            positions = np.vstack([positions, positions[index], positions[index]])

        # Update the element counts for core and shell atoms
        structure_dict["elements"][f"{element}_core"] = len(indices)
        structure_dict["elements"][f"{element}_shell"] = len(indices)

    # Update symbols and positions in the dictionary
    structure_dict["symbols"] = symbols
    structure_dict["positions"] = positions

if __name__ == "__main__":
    # Command-line arguments
    vasp_file = sys.argv[1]  # Input VASP file path
    output_file = sys.argv[2]  # Output file path for the LAMMPS input

    # Read and process the VASP structure
    structure_info = read_vasp_to_dict(vasp_file)

    # Print out the element counts
    print(", ".join([f"{element}: {count}" for element, count in structure_info["elements"].items()]))

    # Ask user to select elements for core-shell processing
    print("------------------->")
    selected_elements = input("Enter the coreshell species (eg. Ba Zr O): \n").split()

    add_coreshell_elements(structure_info, selected_elements)
    # Print updated element counts after adding core-shell atoms
    print("------------------->")
    print("Current elements in the dictionary:")
    print(", ".join([f"{element}: {count}" for element, count in structure_info["elements"].items()]))


    # Ask user for the output order of elements
    print("------------------->")
    output_order = input("Output order (eg. Li La Zr O_core O_shell): \n").split()

    positions = structure_info["positions"]
    symbols = structure_info["symbols"]
    output_matrix = []
    element_indices = {elem: i + 1 for i, elem in enumerate(output_order)}  # Map elements to atom type indices
    element_counts = {elem: 1 for elem in output_order}  # Count atoms per element type

    row_number = 1
    for element in output_order:
        for i, symbol in enumerate(symbols):
            if symbol == element:
                # Create a row for the LAMMPS atom list
                row = [int(row_number), 1, int(element_indices[symbol]), 0] + list(positions[i])
                output_matrix.append(row)
                element_counts[element] += 1
                row_number += 1

    # Create core-shell pairs (bonds between core and shell elements)
    core_shell_pairs = []
    group_number = 1
    for element in selected_elements:
        element_core = f"{element}_core"
        element_shell = f"{element}_shell"
        if element_core in output_order and element_shell in output_order:
            # Find core and shell atom indices for bonding
            core_indices = [int(row[0]) for row in output_matrix if row[2] == element_indices[element_core]]
            shell_indices = [int(row[0]) for row in output_matrix if row[2] == element_indices[element_shell]]
            for core_idx, shell_idx in zip(core_indices, shell_indices):
                core_shell_pairs.append([group_number, core_idx, shell_idx])
            group_number += 1

    core_shell_pairs = np.array(core_shell_pairs)

    # Add numbered core-shell pair IDs
    numbered_core_shell_pairs = []
    for idx, pair in enumerate(core_shell_pairs, start=1):
        numbered_core_shell_pairs.append([int(idx)] + pair.tolist())

    numbered_core_shell_pairs = np.array(numbered_core_shell_pairs)

    # Write the LAMMPS input file
    with open(output_file, 'w') as f:
        # Header and system information
        f.write("# Core-shell model generated by Hao YANG et al.\n\n")
        f.write(f"{len(output_matrix)} atoms\n")
        f.write(f"{len(output_order)} atom types\n")
        f.write(f"{len(numbered_core_shell_pairs)} bonds\n")
        f.write(f"{group_number-1} bond types\n\n")
        
        # Box dimensions
        box_info = structure_info["box"]
        f.write(f"0.0 {box_info[0, 0]:.6f} xlo xhi\n")
        f.write(f"0.0 {box_info[1, 1]:.6f} ylo yhi\n")
        f.write(f"0.0 {box_info[2, 2]:.6f} zlo zhi\n")
        f.write(f"{box_info[0, 1]:.6f} {box_info[0, 2]:.6f} {box_info[1, 2]:.6f} xy xz yz\n\n")
        
        # Masses section: Calculate masses for core-shell elements
        f.write("Masses\n\n")
        for idx, element in enumerate(output_order, start=1):
            if "_core" in element:
                base_element = element.replace("_core", "")
                mass = atomic_masses[atomic_numbers[base_element.capitalize()]] * 0.9  # Core atom mass (90%)
            elif "_shell" in element:
                base_element = element.replace("_shell", "")
                mass = atomic_masses[atomic_numbers[base_element.capitalize()]] * 0.1  # Shell atom mass (10%)
            else:
                mass = atomic_masses[atomic_numbers[element.capitalize()]]
            f.write(f"{idx} {mass:.3f} # Mass of {element}\n")
        f.write("\n")

        # Atoms section: List atom information
        f.write("Atoms  # full\n\n")
        for row in output_matrix:
            row[:4] = list(map(int, row[:4]))  # Convert row data to integers
            f.write(" ".join(map(str, row)) + "\n")
        f.write("\n")

        # Bonds section: List core-shell bonds
        f.write("Bonds\n\n")
        for row in numbered_core_shell_pairs:
            row = list(map(int, row))
            f.write(" ".join(map(str, row)) + "\n")

    print(f"LAMMPS input file written to: {output_file}")
