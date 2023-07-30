import os
import pandas as pd
import sys
import urllib
import urllib.request as request


def make_directory(path, namelist):
    
    """Create a series of directories in given path
        param: path: absolute path as a string
        param: namelist: list of strings of directory names
        return: a list of absolute paths of made directories
    """
    dirlist = []
    assert path[-1] == '/', 'Make sure path ends in backslash symbol'
    for name_index in range(len(namelist)):
        folder = str(path + namelist[name_index])
        try:
            if not os.path.exists(folder):
                os.makedirs(folder)
        except:
            print(namelist[name_index])
        
        dirlist.append(folder)
    return dirlist

#borrowed code
def download_alphafold_pdb(uniprot_id, datadir, version = 'v4'):
    """
    Downloads a PDB file from the Internet and saves it in a data directory.
    :param uniprot_id: UniprotID of structures to get from alphafold
    :param datadir: The directory where the downloaded file will be saved
    :return: the full path to the downloaded PDB file or None if something went wrong
    """
    pdbfn = uniprot_id + "_AlphaFold.pdb"
    url = 'https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_{}.pdb'.format(uniprot_id, version)
    outfnm = os.path.join(datadir, pdbfn)
    try:
        urllib.urlretrieve(url, outfnm)
        return outfnm
    except Exception as err:
        print(str(err))
        return None

    
def modify(tab):
    """
    Download experimental PDB files into specified directory
    param: tab: table to update
    return: updated table with new column for pdbfile, making it rectangular
    """
    pdb = tab['PDB'].to_list()
    uni = tab['Entry'].to_list()
    
    #create a dictionary for pdb id: uniprot id
    pdb_uniprot_dict = {}
    for i in range(len(uni)):
        for pdbcode in pdb[i].split(';'):
            if pdbcode != '':
                pdb_uniprot_dict[pdbcode] = uni[i]
    
    #make the dataset rectangular:
    rect = pd.DataFrame({'Single_pdb': list(pdb_uniprot_dict.keys()), 
                         'Entry': list(pdb_uniprot_dict.values())})
    fin = pd.merge(rect, tab)
    
    return fin
    
def download_pdb(pdbcode, datadir, downloadurl="https://files.rcsb.org/download/"):
    """
    Downloads a PDB file from the Internet and saves it in a data directory.
    :param pdbcode: The standard PDB ID e.g. '3ICB' or '3icb'
    :param datadir: The directory where the downloaded file will be saved
    :param downloadurl: The base PDB download URL, cf.
        `https://www.rcsb.org/pages/download/http#structures` for details
    :return: the full path to the downloaded PDB file or None if something went wrong
    """
    pdbfn = pdbcode + ".pdb"
    url = downloadurl + pdbfn
    outfnm = os.path.join(datadir, pdbfn)
    try:
        urllib.request.urlretrieve(url, outfnm)
        return outfnm
    except Exception as err:
        print(str(err))
        return None
    
def main():
    
    #read in file
    tab = pd.read_csv(intsv, sep = '\t')
    
    assert all(col in tab.columns for col in ['Entry', 'AlphaFoldDB']), 'Check column names'
    
    #uniprotid formating
    Alphafold = [str(i).strip(';') for i in tab['AlphaFoldDB']]
    Uniprotid = tab['Entry']
    
    ######## make the data rectangular, w.r.t pdb entries
    tab = modify(tab)
    
    
    ######## create directory if one does not already exist
    pdbdirectory = os.getcwd() + '/PDB_dump/'
    
    if not os.path.exists(pdbdirectory):
        os.makedirs(pdbdirectory)
    
    pdb_paths = []
    for pdbcode in tab['Single_pdb']:
        pdb_paths.append(download_pdb(pdbcode, pdbdirectory))
    
    ######## add new column
    tab['PDB_paths'] = pdb_paths   
    tab.to_csv('Uniprot_with_structure_cofactor_single_pdb.tsv', sep = '\t', index = False)
    
intsv = sys.argv[1]

if __name__ == '__main__':
    main()
