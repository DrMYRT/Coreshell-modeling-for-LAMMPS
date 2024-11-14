import numpy as np
from ase import Atoms
from ase.io import read
from ase.data import atomic_masses, atomic_numbers
import os

def read_vasp_to_dict(filename):
    # 使用ase库读取VASP结构文件
    atoms = read(filename, format='vasp')

    # 提取盒子信息、元素种类、数目和坐标
    box_info = atoms.get_cell()  # 盒子信息
    symbols = atoms.get_chemical_symbols()  # 元素符号列表
    unique_elements, counts = np.unique(symbols, return_counts=True)  # 元素种类和对应数目
    positions = atoms.get_positions()  # 坐标集合

    # 创建字典存储这些信息
    structure_dict = {
        "box": box_info,
        "elements": dict(zip(unique_elements, counts)),
        "positions": positions,
        "symbols": symbols
    }

    return structure_dict

def add_coreshell_elements(structure_dict, elements):

    # 获取原子符号和坐标
    symbols = structure_dict["symbols"]
    positions = structure_dict["positions"]

    # 添加core和shell元素
    for element in elements:
        # 找到选定元素的索引
        indices = [i for i, symbol in enumerate(symbols) if symbol == element]

        for index in indices:
            symbols.append(f"{element}_core")
            symbols.append(f"{element}_shell")
            positions = np.vstack([positions, positions[index], positions[index]])

        # 更新字典
        structure_dict["elements"][f"{element}_core"] = len(indices)
        structure_dict["elements"][f"{element}_shell"] = len(indices)

    structure_dict["symbols"] = symbols
    structure_dict["positions"] = positions

if __name__ == "__main__":
    vasp_file = '/path/to/your/POSCAR.vasp'  # 修改为你的VASP结构文件路径
    output_file = '/path/to/your/output.vasp'
    structure_info = read_vasp_to_dict(vasp_file)

    for element, count in structure_info["elements"].items():
        print(f"Element: {element}, Count: {count}")

    selected_elements = input("Enter the elements you want to process for coreshell (separated by space): ").split()
    add_coreshell_elements(structure_info, selected_elements)

    # 输出当前字典中的元素种类
    print("Current elements in the dictionary:")
    for element, count in structure_info["elements"].items():
        print(f"Element: {element}, Count: {count}")

    output_order = input("Enter the output order of elements (separated by space): ").split()

    positions = structure_info["positions"]
    symbols = structure_info["symbols"]
    output_matrix = []
    element_indices = {elem: i + 1 for i, elem in enumerate(output_order)}
    element_counts = {elem: 1 for elem in output_order} 

    row_number = 1
    for element in output_order:
        for i, symbol in enumerate(symbols):
            if symbol == element:
                row = [int(row_number), 1, int(element_indices[symbol]), 0] + list(positions[i])
                output_matrix.append(row)
                element_counts[element] += 1
                row_number += 1

    #output_matrix = np.array(output_matrix)

    #print("First 20 rows of the coordinate matrix:")
    #print(output_matrix[:20])

    core_shell_pairs = []
    group_number = 1
    for element in selected_elements:
        element_core = f"{element}_core"
        element_shell = f"{element}_shell"
        if element_core in output_order and element_shell in output_order:
            core_indices = [int(row[0]) for row in output_matrix if row[2] == element_indices[element_core]]
            shell_indices = [int(row[0]) for row in output_matrix if row[2] == element_indices[element_shell]]
            for core_idx, shell_idx in zip(core_indices, shell_indices):
                core_shell_pairs.append([group_number, core_idx, shell_idx])
            group_number += 1

    core_shell_pairs = np.array(core_shell_pairs)

    # 为每行core-shell对添加从1开始的编号
    numbered_core_shell_pairs = []
    for idx, pair in enumerate(core_shell_pairs, start=1):
        numbered_core_shell_pairs.append([int(idx)] + pair.tolist())

    numbered_core_shell_pairs = np.array(numbered_core_shell_pairs)

    #print("Core and shell index pairs with group numbers and row numbers:")
    #print(numbered_core_shell_pairs)

    with open(output_file, 'w') as f:
        f.write("#Lammps full yanghao2023\n\n")
        f.write(f"{len(output_matrix)} atoms\n")
        f.write(f"{len(output_order)} atom types\n")
        f.write(f"{len(numbered_core_shell_pairs)} bonds\n")
        f.write(f"{group_number-1} bond types\n\n")
        
        box_info = structure_info["box"]
        f.write(f"0.0 {box_info[0, 0]:.6f} xlo xhi\n")
        f.write(f"0.0 {box_info[1, 1]:.6f} ylo yhi\n")
        f.write(f"0.0 {box_info[2, 2]:.6f} zlo zhi\n")
        f.write(f"{box_info[0, 1]:.6f} {box_info[0, 2]:.6f} {box_info[1, 2]:.6f} xy xz yz\n\n")
        
        f.write("Masses\n\n")
        for idx, element in enumerate(output_order, start=1):
            if "_core" in element:
                base_element = element.replace("_core", "")
                mass = atomic_masses[atomic_numbers[base_element.capitalize()]] * 0.9
            elif "_shell" in element:
                base_element = element.replace("_shell", "")
                mass = atomic_masses[atomic_numbers[base_element.capitalize()]] * 0.1
            else:
                mass = atomic_masses[atomic_numbers[element.capitalize()]]
            f.write(f"{idx} {mass:.3f} # Mass of {element}\n")
        f.write("\n")

        f.write("Atoms  # full\n\n")
        for row in output_matrix:
            row[:4] = list(map(int, row[:4]))
            f.write(" ".join(map(str, row)) + "\n")
        f.write("\n")

        f.write("Bonds\n\n")
        for row in numbered_core_shell_pairs:
            row = list(map(int, row))
            f.write(" ".join(map(str, row)) + "\n")
        #f.write("\n")
