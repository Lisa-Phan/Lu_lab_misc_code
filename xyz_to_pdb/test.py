from itertools import islice
import numpy as np

####### CLASS #######

class XYZ_atoms:
    """
    Atoms in XYZ file, contains name and coordinate of atoms
    """
    def __init__(self, name: str, coordinates: tuple, connected_atoms: list = []):
        self.name = name
        self.coordinates = coordinates
        self.connected_atoms = connected_atoms

    def __str__(self):
        return f'{self.name} {self.coordinates}'

    def __repr__(self):
        return f'{self.name} {self.coordinates}'
    
    

class XYZ_object: 
    """
    One xyz object contains the coordinates that makes up a molecule
    """
    def __init__(self, name: str, atoms: list):
        self.name = name
        self.atoms = atoms

    def __str__(self):
        return f'{self.name} {self.coordinates}'

    def __repr__(self):
        return f'{self.name} {self.coordinates}'
    
    def remove_atom(self, atom: XYZ_atoms):
        self.atoms.remove(atom)

    def get_xyz(self):
        return self.atoms
    
    def to_file(self, file):
        """
        Write xyz object to file
        """
        with open(file, 'w') as out_file:
            out_file.write(f'{len(self.atoms)}\n')
            out_file.write(f'{self.name}\n')
            for atom in self.atoms:
                out_file.write(f'{atom.name} {atom.coordinates[0]} {atom.coordinates[1]} {atom.coordinates[2]}\n')
            out_file.close()

####### FUNCTIONS #######
def create_coordinate(coord_list):
    """
    Create a coordinate object from a string
    """
    return tuple([float(value) for value in coord_list])

def get_distance(atom1: XYZ_atoms, atom2: XYZ_atoms):
    """
    Calculate distance between two atoms
    """
    return np.linalg.norm(np.array(atom1.coordinates) - np.array(atom2.coordinates))

####### HARD CODE FUNCTIONS #######

def create_XYZ_objs(file: str):
    """
    assume file format is 
    repeat:
        1 line of atom number
        1 line of i = comment
        atom_number lines of atom name and coordinates

            
    Read an xyz file
    Get number of atoms
    Create XYZ objects
    return list of xyz objects
    """
    object_list = []
    with open(file, 'r') as f:
        while True:
            try:
                atom_number = int(f.readline().strip())
            except ValueError: #end of line
                break
            number = f.readline().split()[5].strip(',')
            #print(number)
            atom_line = list(islice(f, int(atom_number)))
            assert(len(atom_line) == atom_number), "Number of atoms in XYZ file does not match the number of atoms in the file"
            XYZ_atom_list = [XYZ_atoms(line.split()[0], create_coordinate(line.split()[1:])) for line in atom_line]
            object_list.append(XYZ_object(number, XYZ_atom_list))
    return object_list


def get_sum_of_covalent_radii(atom1: XYZ_atoms, atom2: XYZ_atoms):
    """
    Get the sum of covalent radii of two atoms
    """
    covalent_radii = {
        'C': 0.762,
        'N': 0.71,
        'O': 0.661,
        'H': 0.31
    }
    if atom1.name not in covalent_radii or atom2.name not in covalent_radii:
        return None
    return covalent_radii[atom1.name] + covalent_radii[atom2.name]

def get_neighbors_within(atom: XYZ_atoms, obj: XYZ_object, distance: float):
    """
    Get neighbor atoms within a specified radii
    """
    neighbors = []
    for atom2 in obj.get_xyz():
        if get_distance(atom, atom2) < distance:
            neighbors.append(atom2)
    return neighbors

def infer_connectivity(obj: XYZ_object, search_distance: float):
    """
    C =	0.762
    N = 0.71
    O = 0.661(19)
    H = 0.31
    
    return a list of conected atoms based on the rule that
    atoms are connected if their distances are smaller than the sum of their covalent radii

    """
    for atoms in obj.get_xyz():
        neighbors = get_neighbors_within(atoms, obj, search_distance)
        for neighbor in neighbors:
            if get_distance(atoms, neighbor) < get_sum_of_covalent_radii(atoms, neighbor):
                atoms.connected_atoms.append(neighbor)
        print(atoms.connected_atoms)
        


    #get atoms around some parameters of iron


def trim_water(obj: XYZ_object, name = 'default_name'):
    """
    Trim water molecules from the object
    Assume that nothing exists after the stretch of water molecules
    """
    xyz_atoms = obj.get_xyz() 
    for count, atom in enumerate(xyz_atoms):
        if atom.name == 'O':
            #look 9 atoms forward
            nine_forward = xyz_atoms[count:count+9]
            nine_atom_names = [atom.name for atom in nine_forward]
            nine_atom_names = ''.join(nine_atom_names)
            if nine_atom_names == 'OHHOHHOHH':
                #delete everything that follows
                del xyz_atoms[count+3:]
    return  XYZ_object(name, xyz_atoms)    



if __name__ == '__main__':
    objs = create_XYZ_objs(r"C:\Users\lisap\Downloads\NPT-pos-1.xyz")
    infer_connectivity(objs[0], 1.7)
    #print(get_distance(obj[0].get_xyz()[0].coordinates, obj[0].get_xyz()[1].coordinates))
    

    # atoms_string_list = []
    # for obj in objs: 
    #     atoms_string = [atom.name for atom in obj.get_xyz()]
    #     #join the string together
    #     atoms_string = ''.join(atoms_string)
    #     atoms_string_list.append(atoms_string)
    
    # #check if all atoms positioning remains the same
    # assert(len(set(atoms_string_list)) == 1), "Atoms positioning is not the same for all objects"
    # print(set(atoms_string_list))
    #print(len(objs[0].get_xyz()))
    #trimmed = trim_water(objs[0])
    #trimmed.to_file(r"C:\Users\lisap\Downloads\NPT-pos-1-trimmed.xyz")

