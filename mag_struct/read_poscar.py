import os

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

            num_atoms = sum([int(x) for x in lines[6].split()])
            atom_positions = []
            atom_symbols = []
            atom_types = []
            atom_symbols.append(lines[5].split())
            atom_types.append(lines[6].split())
            for i in range(8, 8 + num_atoms):
                data = lines[i].split()
                atom_positions.append([float(data[j]) for j in range(3)])

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
