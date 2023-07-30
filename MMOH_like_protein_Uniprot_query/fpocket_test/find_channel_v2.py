#implementing adaptive search for pocket finding
import pandas as pd
import numpy as np
import argparse
import os
import re
from operator import itemgetter

######## FILE FORMAT CORRECTION FUNCTIONS ########
def line_formatting(line):
    # Add a space between the 'HETATM' prefix and the number
    line = re.sub(r'(HETATM)(\d+)', r'\1 \2', line)

    # Replace the minus sign (-) with a space
    line = re.sub(r'(\S)-(\S)', r'\1 -\2', line)

    # Return the corrected line
    return line


########## NORMAL FUNCTIONS #########

FILE = None
POCKET_DICT = None

def get_coordinate_from_line(line):
    x = float(line.split()[6])
    y = float(line.split()[7])
    z = float(line.split()[8])
    coord = np.array([x,y,z])
    return coord



def file_to_pocket_dict(file):
    """
    read a file and returns pocket -> {pocket_number: {sphere: [sphere.coords]}}
    """
    pocket_dict = {}
    with open(file, 'r') as f:
        for line in f:
            if line.startswith('HETATM'):
                line = line_formatting(line)
                print(line)
                print(line.split())
                pocket_mark, sphere_number, pocket_num = itemgetter(2,1,5)(line.split())
                pocket_num = int(pocket_num)
                if pocket_mark == 'APOL':
                    coord = get_coordinate_from_line(line)
                    if pocket_num not in pocket_dict:
                        pocket_dict[pocket_num] = [{'sphere_number': int(sphere_number), 'coord': coord}]
                        #print(pocket_dict)
                    else:
                        if len(pocket_dict[pocket_num]) < int(sphere_number):
                            pocket_dict[pocket_num].append({'sphere_number': int(sphere_number), 'coord': coord})
    return pocket_dict

def run_fpocket(file, min, max):
    """
    run fpocket on command line
    assuming that file input contains only one chain
    file, min, max: string values for file name, min and max sphere radius
    """
    os.system(f'fpocket -f {file} -m {min} -M {max}')

def evaluate_fpocket_output():
    """
    Takes a directory of fpocket output files
    evaluate the pocket detected at active site based on specified conditions

    """
    pass

def get_distance(coord1, coord2):
    return np.linalg.norm(coord1-coord2)


def pocket_to_pocket(pocket_number: int, radius: float, search_mode = 'sphere', pocket_dict = POCKET_DICT ):
    """
    starting search state
    if pocket is found, search of nearby pockets
    iteratively search all atoms
    
    #returns a dictionary

    Not a great implementation, more of a binary report yes/no whether there exists a pocket nearby
    TODO later: search mode lining atoms
    """
    if search_mode == 'sphere':
        #get sphere coordinates of starting
        sphere_coords = [sphere['coord'] for sphere in pocket_dict[pocket_number]]

        #which pocket are located nearby ? 
        nearby_pockets = {}
        for pocket_key in pocket_dict:
            if pocket_key != pocket_number:
                for sphere in pocket_dict[pocket_key]:
                    for coord in sphere_coords:
                        distance = get_distance(coord, sphere['coord'])
                        if distance <= radius:
                            nearby_pockets[pocket_key] = {'query_coord':coord, 'target_coord':sphere['coord'] ,'distance':distance}
        return nearby_pockets
    

def get_min_distance(some_coordinates, one_coordinate):
    """
    return the minimum distance between a pocket and end coordinate
    """
    min_distance = 100000
    for coord in some_coordinates:
        distance = get_distance(coord, one_coordinate)
        if distance <= min_distance:
            min_distance = distance
    return min_distance
   
def find_onepocket(pocket_dict, one_coordinate):
    """
    Determine which pocket to use as starting and ending point
    """
    set_min = 100000
    for pocket_key in pocket_dict:
        some_coordinates = [sphere['coord'] for sphere in pocket_dict[pocket_key]]
        min_distance = get_min_distance(some_coordinates, one_coordinate)
        if min_distance <= set_min:
            set_min = min_distance
            pocket = pocket_key
    return pocket

def pocket_graph_search(start_coord, end_coord, radius, pocket_dict = POCKET_DICT):
    """
    search for a path from start to end
    
    enpoint = find_endpoint()
    start = parser()
    
    if distance(new pocket, destination) < distance(old pocket, destination), reject move
    else:
        move to new pocket
        generate graphs to test
        repeat until find a path from start to end 
    """
    start = find_onepocket(pocket_dict, start_coord)
    end = find_onepocket(pocket_dict, end_coord)
    initial_distance = get_min_distance([sphere['coord'] for sphere in pocket_dict[start]], end_coord)

    search_stack = [key for key in pocket_to_pocket(pocket_dict, start, radius)]
    
    #sort by distance
    search_stack.sort(key=lambda x: get_min_distance([sphere['coord'] for sphere in pocket_dict[x]], end_coord))
    
    visited = []
    while search_stack:
        pocket = search_stack.pop()
        if pocket == end:
            print(pocket, 'ending pocket found')
            return True
        else:
            visited.append(pocket)
            coords = [sphere['coord'] for sphere in pocket_dict[pocket]]
            #get closest distance between pocket and end
            min_distance = get_min_distance(coords, end_coord)
            if min_distance < initial_distance:
                #accept 
                initial_distance = min_distance
                search_stack = [key for key in pocket_to_pocket(pocket_dict, pocket, radius) if key not in visited]
                

#pocket dict is a global variable

class Node:
    def __init__(self, state, state_number, parent=None):
        """
        self.state is list of coords
        self.parent is also list of coords representing object
        self.state_number is pocket_number
        """
        self.state = state
        self.parent = parent
        self.state_number = state_number
    
    def calculate_depth(self):
        """
        Calculate how deep the node is in the graph
        """
        if self.parent is None:
            return 0
        else:
            return self.parent.calculate_depth() + 1 


    def is_goal(self, pocket_dict, ending):
        #find pocket closest to ending point
        one_pocket = find_onepocket(pocket_dict, ending)
        #check if all atoms are in the same pocket
        if one_pocket == self.state_number:
            return True
        return False



def search(node, ending, pocket_dict = POCKET_DICT, radius = 6):
    visited = []
    frontier = [node]
    results = []
    while frontier:
        node = frontier.pop(0) 
        if node.is_goal(pocket_dict, ending):
            results.append(node)
        neighbor_pocket = pocket_to_pocket(node.state_number, radius, 'sphere', pocket_dict)
        children = [Node(state = pocket_dict[p], state_number = p, parent = node) for p in neighbor_pocket if p not in visited]
        frontier.extend(children)
        visited.append(node.state_number)
    results.sort(key=lambda x: x.calculate_depth())
    return results

def path_tracing(node):
    """
    trace the path from the node to the root
    """
    path = []
    while node.parent is not None:
        print('parent node', node.state_number)
        path.append(node.state_number)
        node = node.parent
    print('root node', node.state_number)
    path.append(node.state_number)
    return path 


def get_info_write_info(pocket_dict, final_node, pathstring, outfile, min_sphere, max_sphere):
    """
    Get info and write
    outFile.write('total_pockets shortest_path_length pockets_in_path pdb sphere_params_min sphere_params_max \n')
    """
    with open(outfile, 'w') as outFile:
        total_pock = len(pocket_dict)
        path_length = final_node.calculate_depth() + 1
        pockets_in_path = ''.join([str(pocket) + ',' for pocket in path_tracing(final_node)])
        pdb = os.path.basename(pathstring).split('_')[0]
        print(total_pock, path_length, pockets_in_path, pdb, min_sphere, max_sphere)
        outFile.write('{} {} {} {} {} {}\n'.format(total_pock, path_length, pockets_in_path, pdb, min_sphere, max_sphere))
        
def main():
    parser = argparse.ArgumentParser(description='Pocket path search via specified coordinates')
    parser.add_argument('-f', '--FILE', required = False, type=str, help='path to fpocket output pdb file')
    parser.add_argument('-t', '--TEXTFILE', required = False, type=str, help='path to text file containing params and pocket file to analyze')
    parser.add_argument('-r', '--RADIUS', type=float, help='radius of search')
    parser.add_argument('-o', '--OUTFILE', required = False, type=str, help='path to output file')
    
    args = parser.parse_args()
    
    assert args.RADIUS > 0, 'radius must be greater than 0'
    assert args.FILE or args.TEXTFILE, 'must specify either a pdb file or a text file'
    assert not (args.FILE and args.TEXTFILE), 'must specify either a pdb file or a text file, not both'

    init_coord = np.array([146.95, -121.46, -46.27])
    end_coord = np.array([163.84, -138.45, -56.90])
    
    if args.TEXTFILE:
        textfile = args.TEXTFILE
        pocket_tab = pd.read_csv(textfile, sep='\t')
        for index, row in pocket_tab.iterrows():
            file = row['file_path']
            min_sphere = row['min']
            max_sphere = row['max']
            POCKET_DICT = file_to_pocket_dict(file)
            init_pocket_num = find_onepocket(POCKET_DICT, init_coord)
            end_pocket_num = find_onepocket(POCKET_DICT, end_coord)
            print('initial pocket number:', init_pocket_num)
            print('ending pocket number:', end_pocket_num)

            init_node = Node(state = POCKET_DICT[init_pocket_num], state_number = init_pocket_num)
            end_nodelist = search(init_node, end_coord, POCKET_DICT, args.RADIUS)
            for end_node in end_nodelist:
                get_info_write_info(POCKET_DICT, end_node, file, args.OUTFILE, min_sphere, max_sphere)


    else: 
        file = args.FILE
        POCKET_DICT = file_to_pocket_dict(file)
        
        init_pocket_num = find_onepocket(POCKET_DICT, init_coord)
        end_pocket_num = find_onepocket(POCKET_DICT, end_coord)
        print('initial pocket number:', init_pocket_num)
        print('ending pocket number:', end_pocket_num)

        init_node = Node(state = POCKET_DICT[init_pocket_num], state_number = init_pocket_num)
        end_node_list = search(init_node, end_coord, POCKET_DICT, args.RADIUS)
        for node in end_node_list:
            print('depth of node:', node.calculate_depth())
            path_tracing(node)
            print('----------------------')

    
            

def unit_test(file):
    init_coord = np.array([37.70, -93.78, 12.56])
    end_coord = np.array([27,-93,12])
    pocket_dict = file_to_pocket_dict(file)
    init_pocket_num = find_onepocket(pocket_dict, init_coord)
    init_node = Node(state = pocket_dict[init_pocket_num], state_number = init_pocket_num)
    end_node = search(init_node, end_coord, pocket_dict, 6)
    return end_node


if __name__ == '__main__':
    main()

