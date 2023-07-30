class Atom:
    def __init__(self, symbol, x, y, z):
        self.symbol = symbol
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.connections = []

class Carbon(Atom):
    def __init__(self, x, y, z):
        super().__init__("C", x, y, z)
        self.max_bonds = 4

class Hydrogen(Atom):
    def __init__(self, x, y, z):
        super().__init__("H", x, y, z)
        self.max_bonds = 1

class Oxygen(Atom):
    def __init__(self, x, y, z):
        super().__init__("O", x, y, z)
        self.max_bonds = 2

class Nitrogen(Atom):
    def __init__(self, x, y, z):
        super().__init__("N", x, y, z)
        self.max_bonds = 3

class Sulfur(Atom):
    def __init__(self, x, y, z):
        super().__init__("S", x, y, z)
        self.max_bonds = 2

class Fe(Atom):
    def __init__(self, x, y, z):
        super().__init__("Fe", x, y, z)
        self.max_bonds = 0

class Na(Atom):
    def __init__(self, x, y, z):
        super().__init__("Na", x, y, z)
        self.max_bonds = 0

def distance(atom1, atom2):
    dx = atom1.x - atom2.x
    dy = atom1.y - atom2.y
    dz = atom1.z - atom2.z
    return (dx**2 + dy**2 + dz**2)**0.5

def tree_search(atoms, atom_index):
    if atom_index >= len(atoms):
        return True

    atom = atoms[atom_index]

    # Prioritize satisfying all H- connections first
    if atom.symbol == "H":
        for i in range(len(atoms)):
            if i == atom_index:
                continue

            candidate_atom = atoms[i]
            dist = distance(atom, candidate_atom)
            if dist <= 2.0 and len(atom.connections) < atom.max_bonds and len(candidate_atom.connections) < candidate_atom.max_bonds:
                # Form a connection between atom and candidate_atom
                atom.connections.append(i)
                candidate_atom.connections.append(atom_index)

                # Move on to the next atom
                if tree_search(atoms, atom_index + 1):
                    return True

                # If the current connection didn't lead to a valid solution, backtrack and try another connection
                atom.connections.pop()
                candidate_atom.connections.pop()

    # If the atom is not H or all H- connections are satisfied, try connecting to other atoms
    else:
        for i in range(len(atoms)):
            if i == atom_index:
                continue

            candidate_atom = atoms[i]
            dist = distance(atom, candidate_atom)
            if dist <= 2.0 and len(atom.connections) < atom.max_bonds and len(candidate_atom.connections) < candidate_atom.max_bonds:
                # Form a connection between atom and candidate_atom
                atom.connections.append(i)
                candidate_atom.connections.append(atom_index)

                # Move on to the next atom
                if tree_search(atoms, atom_index + 1):
                    return True

                # If the current connection didn't lead to a valid solution, backtrack and try another connection
                atom.connections.pop()
                candidate_atom.connections.pop()

    # If no valid connection was found for the current atom, return False
    return False

def infer_connectivity(atoms):
    # Start tree search from the first atom
    tree_search(atoms, 0)

    return atoms

def parse_xyz_file(file_path):
    atoms = []
    with open(file_path, 'r') as f:
        lines = f.readlines()
        for i, line in enumerate(lines[2:]):  # Skip the first two lines (number of atoms and comment)
            symbol, x, y, z = line.split()
            if symbol == "C":
                atom = Carbon(x, y, z)
            elif symbol == "H":
                atom = Hydrogen(x, y, z)
            elif symbol == "O":
                atom = Oxygen(x, y, z)
            elif symbol == "N":
                atom = Nitrogen(x, y, z)
            elif symbol == "S":
                atom = Sulfur(x, y, z)
            elif symbol == "Fe":
                atom = Fe(x, y, z)
            elif symbol == "Na":
                atom = Na(x, y, z)
            elif symbol == "Y":
                continue
            elif symbol == "CF":
                continue
            else:
                raise ValueError(f"Unknown atom symbol: {symbol}")
            atoms.append(atom)
    return atoms

def main():
    xyz_file_path = r"C:\Users\lisap\Downloads\NPT-pos-1-trimmed.xyz"  # Replace this with the actual file path
    atoms = parse_xyz_file(xyz_file_path)
    atoms = infer_connectivity(atoms)

    # Print connectivity records as objects
    for atom in atoms:
        connected_atoms = [atoms[i] for i in atom.connections]
        connected_atoms_symbols = [connected_atom.symbol for connected_atom in connected_atoms]
        print(f"{atom.symbol} {atom.x:.3f} {atom.y:.3f} {atom.z:.3f} - Connected to: {', '.join(connected_atoms_symbols)}")

if __name__ == "__main__":
    main()
