import os
from spglib import get_symmetry_dataset

class CrystalStructure:
    def __init__(self, lattice_matrix, atom_positions, atom_symbols, atom_types):
        self.lattice_matrix = lattice_matrix
        self.atom_positions = atom_positions
        self.atom_symbols = atom_symbols
        self.atom_types = atom_types

    @classmethod
    def from_POSCAR(cls, poscar_path):
        with open(poscar_path, 'r') as file:
            # 读取POSCAR文件的内容
            lines = file.readlines()

            lattice_matrix = [[float(lines[2].split()[i]) for i in range(3)],
                              [float(lines[3].split()[i]) for i in range(3)],
                              [float(lines[4].split()[i]) for i in range(3)]]

            # num_atoms = sum([int(x) for x in lines[6].split()])
            # atom_positions = []
            # atom_symbols = []
            # atom_types = []
            # atom_symbols.append(lines[5].split())
            # atom_types.extend([int(element) for element in lines[6].split()])
            # for i in range(8, 8 + num_atoms):
            #     data = lines[i].split()
            #     atom_positions.append([float(data[j]) for j in range(3)])

            num_atoms = sum([int(x) for x in lines[6].split()])
            atom_positions = []
            atom_symbols = lines[5].split()
            atom_types = []

            atom_type_counts = [int(x) for x in lines[6].split()]
            current_type_index = 0

            for i in range(8, 8 + num_atoms):
                data = lines[i].split()
                atom_positions.append([float(data[j]) for j in range(3)])

            # Update the atom_types list based on the current type index
                atom_types.append(current_type_index + 1)

            # Check if we have reached the end of the current type
                if len(atom_positions) == sum(atom_type_counts[:current_type_index + 1]):
                    current_type_index += 1


        return cls(lattice_matrix, atom_positions, atom_symbols, atom_types)

poscar_file_path = "POSCAR"
crystal = CrystalStructure.from_POSCAR(poscar_file_path)

print("Lattice Matrix:")
for row in crystal.lattice_matrix:
    print(row)
print(crystal.lattice_matrix)
print("Atom Positions:")
for position in crystal.atom_positions:
    print(position)

print("Atom Symbols:", crystal.atom_symbols)
print("Atom Types:", crystal.atom_types)



lattice = crystal.lattice_matrix
position = crystal.atom_positions
types = crystal.atom_types
symprec = 1e-3

cell = (lattice, position, types)
dataset = get_symmetry_dataset(cell, symprec=symprec)

print("Lattice:", lattice)
print("Position:", position)
print("Types:", types)
print("Cell:", cell)
print("dataset:", dataset)
print(f'International symbol: {dataset["international"]} ({dataset["number"]})')
print(f'Hall symbol: {dataset["hall"]}')
print("Wyckoff letters: ", end="")
print(" ".join([f"{w}" for w in dataset["wyckoffs"]]))
print("Equivalent atoms:")
for i, equiv_atom in enumerate(dataset["equivalent_atoms"]):
    print(f"{i} -> {equiv_atom}")
print("Space group operations:")
for i, (r, t) in enumerate(zip(dataset["rotations"], dataset["translations"])):
    print(f"--- {i + 1} ---")
    for vec in r:
        print(f"{vec[0]:2d} {vec[1]:2d} {vec[2]:2d}")
    print(f"{t[0]:.5f} {t[1]:.5f} {t[2]:.5f}")