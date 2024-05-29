import pubchempy as pcp
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdmolops
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_useSVG=True
import numpy as np
import builtins
import copy
import matplotlib.pyplot as plt
from collections import defaultdict
from rdkit import RDLogger
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

def get_smiles(molecule_name : str):
    """
    Fetch the SMILES encoding for a molecule by its name using PubChem.

    Parameters:
        molecule_name (str): The name of the molecule to search for.

    Returns:
        str: The SMILES encoding of the molecule or an error message if not found.
    """
    try:
        # Search for compounds in PubChem by name
        results = pcp.get_compounds(molecule_name, 'name')

        # Check if we got any results
        if results:
            # Return the SMILES encoding of the first result
            return results[0].isomeric_smiles
        else:
            return "No results found for the given molecule name."
    except Exception as e:
        return f"An error occurred: {str(e)}"


def canonicalize_smiles(smiles : str):
    mol = Chem.MolFromSmiles(smiles)
    canonicalized_smiles = Chem.MolToSmiles(mol)
    return canonicalized_smiles



def multiplicity(smiles : str):
    """
    Calculates the multiplicity of hydrogen atoms for each carbon atom in a given molecule.

    This function takes a SMILES string representing a molecule, parses it into an RDKit molecule,
    and calculates the multiplicity of hydrogen atoms for each carbon atom that has hydrogens. 
    It also includes other heavy atoms with hydrogens that are not carbon.

    The multiplicity of a carbon atom is defined as the total number of hydrogen atoms on adjacent
    carbon atoms with hydrogens.

    Parameters:
    smiles (str): A valid SMILES string representing the molecule.

    Returns:
    dict_mult_Heavy (dict): A dictionary where keys are the indices of carbon atoms with hydrogens,
                            and values are the multiplicity of hydrogen atoms on adjacent carbon atoms.
    F_dict_Hs (dict): A dictionary where keys are the indices of atoms with hydrogens,
                      and values are the number of hydrogen atoms bonded to those atoms.

    Raises:
    ValueError: If the input SMILES string is invalid.

    Example:
    >>> multiplicity('CCO')
    ({0: 2, 1: 3, 2: 0}, {0: 3, 1: 2, 2: 1})
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError("Invalid SMILES string")
    
    # Dictionary of hydrogens bonded to each atom
    dict_Hs = {atom.GetIdx(): atom.GetTotalNumHs(includeNeighbors=True) for atom in mol.GetAtoms()}
    F_dict_Hs = {idx: count for idx, count in dict_Hs.items() if count != 0}
    
    dict_mult_Heavy = {atom.GetIdx(): 0 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C' and dict_Hs[atom.GetIdx()] != 0}
    
    # Calculate multiplicity for carbons with hydrogens next to carbons with hydrogens
    for bond in mol.GetBonds():
        begin_idx, end_idx = bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()
        if begin_idx in dict_mult_Heavy and end_idx in dict_mult_Heavy:
            dict_mult_Heavy[begin_idx] += dict_Hs[end_idx]
            dict_mult_Heavy[end_idx] += dict_Hs[begin_idx]
    
    # Include heavy atoms with hydrogens that are not carbon
    for atom in mol.GetAtoms():
        idx = atom.GetIdx()
        if dict_Hs[idx] != 0 and idx not in dict_mult_Heavy:
            dict_mult_Heavy[idx] = 0
    
    return dict_mult_Heavy, F_dict_Hs



def from_mol_to_shift_figure(smiles):
    """
    Generate a dictionary representing the type of carbon atoms in a molecule.
    
    Parameters:
        smiles (str): The SMILES string representation of the molecule.
    
    Returns:
        dict: A dictionary where keys are the indices of atoms with hydrogens
              and values represent their chemical environment:
              0 for non-carbon atoms,
              1 for aliphatic carbons,
              2 for aromatic carbons,
              3 for carbons double-bonded to other carbons.
    """
    mol = Chem.MolFromSmiles(smiles)
    
    atom_symbols = {atom.GetIdx(): atom.GetSymbol() for atom in mol.GetAtoms()}
    multiplicities, hydrogens_dict = multiplicity(smiles)
    
    aromatic_carbons = {}
    double_bonded_carbons = {}
    
    for bond in mol.GetBonds():
        begin_idx = bond.GetBeginAtomIdx()
        end_idx = bond.GetEndAtomIdx()
        
        if bond.GetIsAromatic():
            aromatic_carbons[begin_idx] = 1
            aromatic_carbons[end_idx] = 1
        elif bond.GetBondType() == Chem.BondType.DOUBLE and atom_symbols[begin_idx] == 'C' and atom_symbols[end_idx] == 'C':
            double_bonded_carbons[begin_idx] = 1
            double_bonded_carbons[end_idx] = 1
    
    shift_figure = {}
    
    for atom_idx in hydrogens_dict.keys():
        if atom_symbols[atom_idx] == 'C':
            if atom_idx in aromatic_carbons:
                shift_figure[atom_idx] = 2
            elif atom_idx in double_bonded_carbons:
                shift_figure[atom_idx] = 3
            else:
                shift_figure[atom_idx] = 1
        else:
            shift_figure[atom_idx] = 0
    
    return shift_figure



def find_atoms_in_same_ring(atom_index :int, smiles : str):
    """
    Finds all atoms in the same aromatic ring(s) as a specified atom in a molecule.

    This function takes an atom index and a SMILES string representing a molecule,
    and returns a sorted list of indices of all atoms that are in the same aromatic ring(s)
    as the specified atom.

    Parameters:
    atom_index (int): The index of the atom for which to find ring neighbors.
    smiles (str): A valid SMILES string representing the molecule.

    Returns:
    list: A sorted list of indices of atoms in the same aromatic ring(s) as the specified atom.
          Returns "Invalid SMILES" if the input SMILES string is invalid.

    Raises:
    ValueError: If the atom index is out of range of the molecule's atoms.

    Example:
    >>> find_atoms_in_same_ring(1, 'c1ccccc1CCC')
    [0, 1, 2, 3, 4, 5]
    """
    # Load the molecule from SMILES
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        return "Invalid SMILES"
    
    # Sanitize the molecule to ensure ring info is accurate
    Chem.SanitizeMol(mol)
    
    # Get ring info
    ring_info = mol.GetRingInfo()
    
    # List to store indices of atoms in the same ring(s) as the specified atom
    atoms_in_same_rings = set()

    # Check each ring to see if it contains the specified atom
    for ring in ring_info.AtomRings():
        if atom_index in ring and mol.GetAtomWithIdx(atom_index).GetIsAromatic():
            # Add all atoms in this ring to the set
            atoms_in_same_rings.update(ring)

    # Return the result as a sorted list, but you might want to remove the original atom from the list
    return sorted(atoms_in_same_rings)



def atom_branches(atom_index, smiles :str):
    """
    Finds all branches (substructures) starting from a specified atom in a molecule.

    This function takes an atom index and a SMILES string representing a molecule,
    and returns a list of dictionaries. Each dictionary represents a branch originating
    from the specified atom, with keys as atom indices and values as the depth (number
    of bonds) from the central atom.

    Parameters:
    atom_index (int): The index of the atom from which to find branches.
    smiles (str): A valid SMILES string representing the molecule.
    Returns:
    list: A list of dictionaries, where each dictionary represents a branch. The keys
          are atom indices, and the values are the depth (number of bonds) from the 
          central atom. Each branch dictionary includes the central atom with a depth of 0.

    Raises:
    ValueError: If the input SMILES string is invalid or the atom index is out of range.

    Example:
    >>> atom_branches(0, 'CCO')
    [{0: 0, 1: 1}, {0: 0, 2: 1}]
    """
    # Load the molecule from the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError("Invalid SMILES string")

    if atom_index < 0 or atom_index >= mol.GetNumAtoms():
        raise ValueError("Atom index is out of range")
    # Load the molecule from the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    
    # Initialize a list to store the branches as dictionaries
    branches = []

    # Get the atom object from the molecule
    central_atom = mol.GetAtomWithIdx(atom_index)

    # Function to find branches using depth-first search
    def dfs(atom, current_branch, visited_atoms, depth):
        if atom.GetIdx() in visited_atoms:
            return

        visited_atoms.add(atom.GetIdx())

        
        # Record the depth (number of bonds) from the central atom
        current_branch[atom.GetIdx()] = depth

        if atom.GetIsAromatic and atom.GetSymbol()=='C':
            aromatic_start = atom.GetIdx()
            if atom.GetIdx() in find_atoms_in_same_ring(aromatic_start, smiles):
                return
            
        for neighbor in atom.GetNeighbors():
            if neighbor.GetIdx() not in visited_atoms:
                dfs(neighbor, current_branch, visited_atoms, depth + 1)
    
    # Explore each neighboring atom
    for neighbor in central_atom.GetNeighbors():
        branch = {}
        visited_atoms = set([atom_index])  # Initialize with the central atom
        dfs(neighbor, branch, visited_atoms, 1)
        branches.append(branch)
    for branch in branches:
        branch[central_atom.GetIdx()]= 0

    return branches





def map(idx : int, smiles : str):
    """
    Maps the atoms and their depths (number of bonds) from a specified atom in a molecule.

    This function takes an atom index and a SMILES string representing a molecule,
    and returns a dictionary mapping each atom index to its depth (number of bonds)
    from the specified atom.

    Parameters:
    idx (int): The index of the central atom from which to map depths.
    smiles (str): A valid SMILES string representing the molecule.

    Returns:
    dict: A dictionary where keys are atom indices, and values are the depths (number
          of bonds) from the specified central atom.

    Raises:
    ValueError: If the input SMILES string is invalid or the atom index is out of range.

    Example:
    >>> map(0, 'CC(C)CC(=O)CC')
    {0: 1, 1: 0, 2: 1, 3: 1, 4: 2, 5: 3, 6: 3, 7: 4}
    """
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError("Invalid SMILES string")

    if idx < 0 or idx >= mol.GetNumAtoms():
        raise ValueError("Atom index is out of range")
    d = {}
    a_b = atom_branches(idx, smiles)
    for b in a_b:
        for a in list(b.keys()):
            d[a] = b[a]
    return d




def branches_to_smiles(branches, smiles : str):
    """
    Generates SMILES representations of substructures (branches) in a molecule.

    This function takes a list of dictionaries representing branches in a molecule,
    where each dictionary maps atom indices to their depths (number of bonds) from a central atom.
    It then generates SMILES representations for each branch and returns them as a list.

    Parameters:
    branches (list): A list of dictionaries where each dictionary represents a branch.
                     The keys are atom indices, and the values are the depths (number of bonds)
                     from a central atom.
    smiles (str): A valid SMILES string representing the molecule.

    Returns:
    list: A list of SMILES representations of substructures (branches) in the molecule with aromatic
    rings cut off. 

    Raises:
    ValueError: If the input SMILES string is invalid.

    Example:
    >>> smiles = 'CC(C)CC(=O)CCc1ccccc1'
    >>> branches = atom_branches(1, smiles)
    >>> branch_smiles = branches_to_smiles(branches, smiles)
    >>> print(branch_smiles)
    ['CC', 'CC', 'cCCC(=O)CC']
    """
    mol = Chem.MolFromSmiles(smiles)
    def get_connected_subgraph_smiles(mol, atom_indices):
        # Create a new editable molecule
        emol = Chem.EditableMol(Chem.Mol())
        idx_map = {}
        
        # Add atoms
        for idx in atom_indices:
            #If the atom hasn't any branches, i.e monoheavy atom molecules.
            if not atom_indices:
                return []
            new_idx = emol.AddAtom(mol.GetAtomWithIdx(idx))
            idx_map[idx] = new_idx
        
        # Add bonds
        for idx in atom_indices:
            for nbr in mol.GetAtomWithIdx(idx).GetNeighbors():
                if nbr.GetIdx() in atom_indices and idx < nbr.GetIdx():
                    bond = mol.GetBondBetweenAtoms(idx, nbr.GetIdx())
                    emol.AddBond(idx_map[idx], idx_map[nbr.GetIdx()], order=bond.GetBondType())

        new_mol = emol.GetMol()
        return Chem.MolToSmiles(new_mol)

    branch_smiles = []
    for branch in branches:
        atom_indices = list(branch.keys())
        smiles_branch = get_connected_subgraph_smiles(mol, atom_indices)
        branch_smiles.append(smiles_branch)

    return branch_smiles




def clean_aromatics(central_atom : int, smiles :str):
    mol = Chem.MolFromSmiles(smiles)
    updated_branch_smiles = []
    #Someone could build upon this to work on aromatic rings and use this for the identification of functional 
    #groups for aromatic carbons and thus do NMR for molecules with two aromatic rings.
    if not mol.GetAtomWithIdx(central_atom).GetIsAromatic:
        return
    else:
        branch_smiles = branches_to_smiles(atom_branches(central_atom, smiles), smiles)
        for smile in branch_smiles:
            smiles_updated = ''
            i = 1
            for element in smile:
                if element == 'c':
                    smiles_updated += f'c{i}ccccc{i}'
                    i +=1
                else:  
                    smiles_updated += element 
            smiles_updated = canonicalize_smiles(smiles_updated)
            updated_branch_smiles.append(smiles_updated)
    return updated_branch_smiles




dict_shifts_identify = {1 : {
    "C": 0.54,        # Assuming 'alkyl' as a generic single carbon alkyl group. Also, it is the difference between methane and ethane.
    "C=C": 1.33,      # Generic double bond
    "C#C": 1.52,      # Generic triple bond
    "C=C[C]": 1.52,   # Alkyl on a double bond
    "C#Cc1ccccc1": 1.77,  # Phenyl on a triple bond
    "c1ccccc1": 1.85, # Phenyl
    "F": 3.15,
    "Cl": 2.48,
    "Br": 2.29,
    "O": 2.46,        # Hydroxy group
    #"OC": 2.27,       # Alkoxy group
    "Oc1ccccc1": 2.89, # Phenoxy group
    #"OC(=O)C": 2.99,    # OCoalkyl, assumed
    "OC(=O)c1ccccc1": 3.23, # OCophenyl
    "N": 1.69,        # Amino group
    "NC": 1.60,       # Alkylamino
    "N(C)C" : 1.41,   # Dialkylamine (assuming simple dialkyl)
    "[NH]c1ccccc1" : 2.15,
    "[NH2+]C": 2.31,   # Protonated amine  
    "N(C)c1ccccc1": 2.39, # N(alkyl)phenyl
    "[NH3+]": 2.31,
    "[NH2+]C": 2.31,    # Assuming amine with substitution
    "[NH+](C)C": 2.46,
    "[N+](C)(C)C": 2.56, # Tertiary amine
    "[NH]C(=O)C": 2.23, # Assuming another tertiary amine configuration
    "N(C)C(=O)C": 2.23,
    "NC(=O)c1ccccc1": 2.33, # Amide with phenyl
    "S": 1.63,        # Thiol group
    "SC": 1.66,       # Alkylthio
    "Sc1ccccc1": 1.92, # Phenylthio
    "C(=O)": 1.58,   # Coalkyl, assumed ketone
    "C(=O)c1ccccc1": 2.08, # COPhenyl, assumed ketone with phenyl
    "C(=O)O": 1.44,   # Carboxylic acid
    "C(=O)OC": 1.49,  # Ester with alkyl
    "C(=O)Oc1ccccc1": 1.74, # Ester with phenyl
    "C(=O)N": 1.39,   # Amide
    "C(=O)NC": 1.39,  # Amide with alkyl substitution
    "C(=O)N(C)C": 1.99, # Amide with dialkyl substitution
    "C(=O)Nc1ccccc1": 1.59, # Amide with phenyl
    "C#N": 1.73       # Nitrile
}, 
 2 : {'[H]': (0, 0, 0),
 'C': (-0.2, -0.12, -0.22),
 'CC': (-0.14, -0.06, -0.17),
 'C(C)C': (-0.13, -0.08, -0.18),
 'C(C)(C)C': (0.02, -0.08, -0.21),
 'CCl': (0.0, 0.0, 0.0),
 'C(F)(F)F': (0.32, 0.14, 0.2),
 'C(Cl)(Cl)Cl': (0.64, 0.13, 0.1),
 'CO': (-0.07, -0.07, -0.07),
 'C=C': (0.06, -0.03, -0.1),
 'C=CC1=CC=CC=C1': (0.15, -0.01, -0.16),
 'C#C': (0.15, -0.02, -0.01),
 'C#CC1=CC=CC=C1': (0.19, 0.02, 0.0),
 'c1ccccc1': (0.37, 0.2, 0.1),
 'F': (-0.26, 0.0, -0.2),
 'Cl': (0.03, -0.02, -0.09),
 'Br': (0.18, -0.08, -0.04),
 'I': (0.39, -0.21, 0.0),
 'O': (-0.56, -0.12, -0.45),
 'OC': (-0.48, -0.09, -0.44),
 'OCC': (-0.46, -0.1, -0.43),
 'OC1=CC=CC=C1': (-0.29, -0.05, -0.23),
 'O=C(OC)': (-0.25, 0.03, -0.13),
 'O=C(Oc1ccccc1)': (-0.09, 0.09, -0.08),
 'OS(=O)(=O)C': (-0.05, 0.07, -0.01),
 'N': (-0.75, -0.25, -0.65),
 'NC': (-0.8, -0.22, -0.68),
 'N(C)C': (-0.66, -0.18, -0.67),
 'NC(=O)C': (0.12, -0.07, -0.28),
 'N(C)C(=O)C': (-0.16, 0.05, -0.02),
 'NN': (-0.6, -0.08, -0.55),
 'N=Nc1ccccc1': (0.67, 0.2, 0.2),
 '[N+](=O)[O-]': (0.58, 0.31, 0.37),
 '[N+](=O)(O)[O-]': (0.95, 0.26, 0.38),
 'S': (-0.08, -0.16, -0.22),
 'SC': (-0.08, -0.1, -0.24),
 'Sc1ccccc1': (0.06, -0.09, -0.15),
 'S(=O)(=O)OC': (0.6, 0.26, 0.33),
 'S(=O)(=O)Cl': (0.76, 0.35, 0.45),
 'C=O': (0.56, 0.22, 0.29),
 'C(=O)C': (0.62, 0.14, 0.21),
 'C(=O)CC': (0.63, 0.13, 0.2),
 'COC(C)(C)C': (0.44, 0.05, 0.05),
 'C(=O)c1ccccc1': (0.47, 0.13, 0.22),
 'C(=O)O': (0.85, 0.18, 0.27),
 'C(=O)OC': (0.71, 0.11, 0.21),
 'C(=O)OCC(C)C': (0.7, 0.09, 0.19),
 'C(=O)Oc1ccccc1': (0.9, 0.17, 0.27),
 'C(=O)N': (0.61, 0.1, 0.17),
 'C(=O)Cl': (0.84, 0.22, 0.36),
 'C(=O)Br': (0.8, 0.21, 0.37),
 'C=Nc1ccccc1': (0.6, 0.2, 0.2),
 'C#N': (0.36, 0.18, 0.28),
 'Si(C)(C)C': (0.22, -0.02, -0.02),
 'P(O)(OC)OC': (0.48, 0.16, 0.24)}, 
3: {
 '[H]': (0, 0, 0),
 'C': (0.45, -0.22, -0.28),
 'C1CCCCC1': (0.69, -0.25, -0.28),
 'Cc1ccccc1': (1.05, -0.29, -0.32),
 'CF': (0.7, 0.11, -0.04), 'CCl': (0.7, 0.11, -0.04),'CBr': (0.7, 0.11, -0.04),
 'C(F)F': (0.7, 0.11, 0.34),
 'C(F)(F)F': (0.66, 0.61, 0.32),
 'CO': (0.64, -0.01, -0.02),
 'CN': (0.58, -0.1, -0.08),
 'CS': (0.71, -0.13, -0.22),
 'CC(=O)': (0.69, -0.08, -0.06),
 'C=C': (0.47, 0.38, 0.12),
 'C=CC=CC': (1.24, 0.02, 0.05),
 'c1ccccc1': (1.6, -0.05, 0.05),
 'C#C' : (0.47, 0.38, 0.12),
 'c1ccccc1O': (1.65, 0.19, 0.09),
 'F': (1.54, -0.4, -1.02),
 'Cl': (1.08, 0.18, 0.13),
 'Br': (1.07, 0.45, 0.55),
 'I': (1.14, 0.81, 0.88),
 'OC': (1.22, -1.07, -1.21),
 'OC=C': (1.21, -0.6, -1.0),
 'OC(=O)C': (2.11, -0.35, -0.64),
 'N': (0.8, -0.26, -0.21),
 'NC': (0.8, -1.26, -1.21),
 'N(C)C': (0.8, -1.26, -1.21),
 'N(C=C)C': (1.17, -0.53, -0.99),
 'NC(=O)C': (2.08, -0.57, -0.72),
 'N=Nc1ccccc1': (2.39, 1.11, 0.67),
 'N(=O)(=O)': (1.87, 1.3, 0.62),
 'SC': (1.11, -0.29, -0.13),
 'S(=O)(OC)': (1.27, 0.67, 0.41),
 'S(=O)(=O)C': (1.55, 1.16, 0.93),
 'SC(=O)C': (1.41, 0.06, 0.02),
 'SC#N': (0.94, 0.45, 0.41),
 'S(F)(F)(F)(F)F': (1.68, 0.61, 0.49),
 'C=O': (1.06, 0.91, 0.74),
 'C(=O)O': (0.8, 0.98, 0.32),
 'C(=O)OC': (0.78, 1.01, 0.46),
 'C(=O)N(C)C': (1.37, 0.98, 0.46),
 'C(=O)Cl': (1.11, 1.46, 1.01),
 'C#N': (0.27, 0.75, 0.55),
 'P(OCC)OCC': (0.66, 0.88, 0.67),
 'OP(OCC)OCC': (1.33, -0.34, -0.66),
 'C(Cl)Cl': (0.7, 0.11, 0.34),
 'C(Br)Br': (0.7, 0.11, 0.34)}}





dict = builtins.dict
sorted_data = {k: dict(sorted(dict_shifts_identify[k].items(), key=lambda item: len(item[0]), reverse=True)) for k in range(1, 4)}



def shift_0(idx, smiles: str):
    """
    Determines the chemical shift of a specific atom in a molecule, but only if the atom has a hydrogen bonded to it.

    This function takes an atom index and a SMILES string representing a molecule,
    and returns the chemical shift of the specified atom based on predefined criteria.

    Parameters:
    idx (int): The index of the atom for which to determine the chemical shift.
    smiles (str): A valid SMILES string representing the molecule.

    Returns:
    float or str: The chemical shift of the specified atom if it has a hydrogen bonded to it.
                  If the criteria for determining the chemical shift is not met, or if the atom
                  does not have a hydrogen bonded to it, returns a string indicating the situation.

    Example:
    >>> shift_0(2, 'CCO')
    2.46
    >>> shift_0(2, 'CCP')
    An error occurred, no values of shift_0 for atom with symbol: P
    >>> shift_0(2, 'CCF')
    An error occurred, atom with symbol: F does not have a bonded hydrogen
    """
    mol = Chem.MolFromSmiles(smiles)
    atom = mol.GetAtomWithIdx(idx)
    symbol = atom.GetSymbol()

    # Check if the atom has at least one hydrogen bonded to it
    dict_Hs = multiplicity(smiles)[1]

    if idx not in dict_Hs.keys():
        return f'An error occurred, atom with symbol: {symbol} does not have a bonded hydrogen'

    if symbol == 'O':
        return 2.46  # Based on images of the chemical shift of the hydroxyl group in ethanol.
    elif symbol == 'N':
        return 2.4  # Based on images of the chemical shift of methanamine. Even though typically it is between 0.5 and 5 ppm.
    elif symbol == 'S':
        return 0.72  # It is the value of chemical shift of hydrogen sulfide.
    return f'An error occurred, no values of shift_0 for atom with symbol: {symbol}'





def search_algo(idx : int, smiles : str):
    """
    Searches for substructures in a molecule based on a specified atom index.

    This function takes an atom index and a SMILES string representing a molecule,
    and searches for substructures in the molecule based on predefined criteria.
    It returns a list of dictionaries, where each dictionary contains functional groups
    found within a certain distance from the specified atom.

    Parameters:
    idx (int): The index of the central atom for searching substructures.
    smiles (str): A valid SMILES string representing the molecule.

    Returns:
    list: A list of dictionaries, where each dictionary represents functional groups found
          within a certain distance from the specified central atom. The keys are distances
          from the central atom, and the values are lists of SMILES strings representing
          the found functional groups.

    Example:
    >>> smiles = 'CCCC(=O)O'
    >>> search_algo(0, smiles)
    [{1: ['C'], 2: ['C'], 3: ['C(=O)O']}]
    """
    mol = Chem.MolFromSmiles(smiles)
    List = []
    if len(mol.GetAtomWithIdx(idx).GetNeighbors()) >= 2: # or idx == 0:
        branches = clean_aromatics(idx, smiles) #Gets all the branches of the central atom.
        
        for branch in branches:                     #Looks through all of them individually.
            #in the case of a aldehyde:
            if branch == 'C=O':
                continue
            new_mol = Chem.MolFromSmiles(branch)
            dict_smiles_found = {}                  
            map_dist = map(0, branch)    #Hold the distance from the central atom (which is zero because of the branches) to all other atoms.
            for atom in range(0,max(list(map_dist.values()))+1):
                dict_smiles_found[atom] = []        #It will hold every functional group and the distance from central atom. 
            #Every branch starts at index 0. 
            atoms_seen = set([0])
            #Looks through all functional groups and finds the shortest distance from the central atom.
            for smiles_data in list(sorted_data[1].keys()):  
                mol_data = Chem.MolFromSmiles(smiles_data)
                gsm = new_mol.GetSubstructMatches(mol_data)
                if gsm!=():
                    for resemblance in gsm:
                        dict = {}
                        for a in resemblance:
                            if a in list(map_dist.keys()):
                                dict[a] = map_dist[a]
                                y_min = min(list(dict.values()))
                        if not any(resemblance[a] in atoms_seen for a in range(0,len(resemblance))):
                            atoms_seen.update(a for a in resemblance)
                            dict_smiles_found[y_min].append(smiles_data)  


                          
            del dict_smiles_found[0]
            Final_dict = {key: value for key, value in dict_smiles_found.items() if value != []}
            List.append(Final_dict)
    else:
        #Deals with the case of methane.
        if len(smiles) > 1:
            branch = clean_aromatics(idx, smiles)[0] #Gets only one branch from the central atom.
            #in the case of formaldehyde:
            if branch == 'C=O':
                return 0
            new_mol = Chem.MolFromSmiles(branch)
            dict_smiles_found = {}                  
            map_dist = map(0, branch)    #Hold the distance from the central atom (which is zero because of the branches) to all other atoms.
            del map_dist[0]
            for atom in range(0,max(list(map_dist.values()))+1):
                dict_smiles_found[atom] = []        #It will hold every functional group and the distance from central atom. 
            #Every branch starts at index 0. 
            atoms_seen = set([0])
            #Looks through all functional groups and finds the shortest distance from the central atom.
            for smiles_data in list(sorted_data[from_mol_to_shift_figure(smiles)[idx]].keys()):  
                mol_data = Chem.MolFromSmiles(smiles_data)
                gsm = new_mol.GetSubstructMatches(mol_data)
                if gsm != ():
                    for resemblance in gsm:
                        dict = {}
                        for a in resemblance:
                            if a in list(map_dist.keys()):
                                dict[a] = map_dist[a]
                                y_min = min(list(dict.values()))
                            else:
                                pass

                        if not any(resemblance[a] in atoms_seen for a in range(0,len(resemblance))):
                            atoms_seen.update(a for a in resemblance)
                            dict_smiles_found[y_min].append(smiles_data)                             
        else:
            return []
        del dict_smiles_found[0]
        Final_dict = {key: value for key, value in dict_smiles_found.items() if value != []}
        List.append(Final_dict)
        
    return List



def shift_1(idx : int, smiles :str):
    """
    Calculates the chemical shift of the hydrogens of a specified alkyl carbon in a molecule 
    based on nearby functional groups.

    This function takes an atom index and a SMILES string representing a molecule, and calculates
    the chemical shift of the specified atom based on nearby functional groups within a certain
    radius. It considers the effect of each functional group on the chemical shift, dividing it
    by the square of the distance from the functional group to the specified atom.

    Parameters:
    idx (int): The index of the atom for which to calculate the chemical shift.
    smiles (str): A valid SMILES string representing the molecule.

    Returns:
    float: The calculated chemical shift of the specified atom.

    Example:
    >>> shift_1(1, 'CCCC=O')
    1.475
    """
    shift_dict = {}
    shift = 0
    Idx_groups = search_algo(idx, smiles)
    table = from_mol_to_shift_figure(smiles)[idx]
    #Considers all functional groups within a radius of 4 and divides the effect with a ratio of distance squared.
    for dict in Idx_groups:
        for dist, groups in dict.items():
            if dist <= 4:
                for group in groups:
                    shift += sorted_data[table][group]*(1/(dist)**2)

    return shift



def process_smiles(smiles):
    """
    Process a SMILES string by converting it to its canonical form and removing slashes and hyphens.

    Parameters:
        smiles (str): A SMILES string of a molecule.

    Returns:
        str: The processed canonical SMILES string with slashes and hyphens removed.
    """
    if not isinstance(smiles, str) or not smiles:
        raise ValueError("Input must be a non-empty string.")
    
    try:
        # Convert the SMILES string to a molecule object
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            raise ValueError("Invalid SMILES string.")
        
        # Generate the canonical SMILES string from the molecule object
        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
        
        # Remove slashes and hyphens
        processed_smiles = canonical_smiles.replace('/', '').replace('\\', '').replace('-', '')
        
        return processed_smiles
    except Exception as e:
        raise ValueError(f"An error occurred while processing the SMILES string: {str(e)}")



def find_aromatic_carbon_indexes(smiles):
    """
    Identify and isolate the indexes of aromatic carbon atoms in a SMILES string.

    Parameters:
        smiles (str): A SMILES string of a molecule.

    Returns:
        list of lists: Sorted lists of indexes for aromatic carbons in each ring, based on their positions in the original SMILES string.
    """
    if not isinstance(smiles, str) or not smiles:
        raise ValueError("Input must be a non-empty string.")

    # Parse the SMILES string to create a molecule
    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        raise ValueError("Invalid SMILES string")

    # Create a mapping of the atom indexes to their positions in the original SMILES string
    atom_positions = []
    atom_counter = 0
    for i, char in enumerate(smiles):
        if char == 'c':
            atom_positions.append(i)
            atom_counter += 1
        else:
            atom_counter += 1
            pass  # adjust for ring closure digits

    return atom_positions



def sort_sublists(main_list):
    # Sort each sublist individually
    sorted_main_list = [sorted(sublist) for sublist in main_list]
    return sorted_main_list



def find_aromatic_substituents_single(og_string, indexes):
    """
    Identify the substituents for a single aromatic cycle in a given string.

    Parameters:
        og_string (str): The original string representing the molecule.
        indexes (list of int): List of 6 indexes indicating positions of aromatic carbons in the cycle.

    Returns:
        list of str: List of 6 substituents corresponding to the positions around the aromatic cycle.
    """
    substrings = ['','','','','','']
    index_to_delete = indexes[5]+1
    s = og_string[:index_to_delete] + og_string[index_to_delete + 1:]
    for i in range(0, indexes[0]):
        substrings[0] += s[i]
    if substrings[0] == '':
        substrings[0] = 'H'
    for y in range(indexes[1]+1, indexes[2]):
        substrings[1] += s[y]
    if substrings[1] == '':
        substrings[1] = 'H'
    for m in range(indexes[2]+1, indexes[3]):
        substrings[2] += s[m]
    if substrings[2] == '':
        substrings[2] = 'H'
    for n in range(indexes[3]+1, indexes[4]):
        substrings[3] += s[n]
    if substrings[3] == '':
        substrings[3] = 'H'
    for l in range(indexes[4]+1, indexes[5]):
        substrings[4] += s[l]
    if substrings[4] == '':
        substrings[4] = 'H'
    for w in range(indexes[5]+1, len(s)):
        substrings[5] += s[w]
    if substrings[5] == '':
        substrings[5] = 'H'
    return substrings


def find_substituents_aromatic_multiple(s, indexes):
    """
    Find substituents for multiple aromatic cycles in a given string.

    Parameters:
        s (str): The original string representing the molecule.
        indexes (list of list of int): List of lists where each sublist contains 6 indexes 
                                       indicating positions of aromatic carbons in each cycle.

    Returns:
        list of list of str: List of lists where each sublist contains 6 substituents 
                             corresponding to the positions around each aromatic cycle.
    """
    if not isinstance(s, str) or not s:
        raise ValueError("Input must be a non-empty string.")
    
    result = []
    subs = []
    for i in range(len(indexes)):
        aromatic_sub = find_aromatic_substituents_single(s, indexes[i])
        subs.append(aromatic_sub)
    return subs



def process_smiles_lists(smiles_lists):
    """
    Remove unnecessary parentheses from SMILES strings and canonicalize the SMILES strings,
    ensuring that 'H' is maintained as is.

    Parameters:
        smiles_lists (list of list of str): A list of lists containing SMILES strings.

    Returns:
        list of list of str: A list of lists with SMILES strings processed to remove unnecessary
                             parentheses and canonicalized, with 'H' left unchanged.
    """
    processed_lists = []
    for sublist in smiles_lists:
        processed_sublist = []
        for smile in sublist:
            # Ensure the input is a string
            if not isinstance(smile, str):
                raise TypeError("All elements of the input must be strings.")
            # Remove parentheses if the string starts with '(' and ends with ')'
            if smile.startswith('(') and smile.endswith(')'):
                smile = smile[1:-1]
            # Canonicalize the SMILES string
            if smile == 'H':
                processed_sublist.append('H')
            else:
                try:
                    mol = Chem.MolFromSmiles(smile)
                    if mol is not None:  # Ensure mol is valid
                        canonical_smiles = Chem.MolToSmiles(mol, canonical=True)
                        processed_sublist.append(canonical_smiles)
                    else:
                        processed_sublist.append(smile)  # Return original if RDKit can't parse it
                except:
                    processed_sublist.append(smile)  # Handle any other exceptions
        processed_lists.append(processed_sublist)
    return processed_lists



data_dict = {
 'C':5,
 'CC':4,
 'C(C)C':3,
 'C(C)(C)C':2,
 'CCl':2,
 'C(F)(F)F':1,
 'C(Cl)(Cl)Cl':1,
 'CO':3,
 'C=C':3,
 'C=CC1=CC=CC=C1':1,
 'C#C':1,
 'C#CC1=CC=CC=C1':1,
 'c1ccccc1':2,
 'C1CCCCC1':2,
 'F':4,
 'Cl':4,
 'Br':4,
 'I':4,
 'O':5,
 'OC':3,
 'OCC':2,
 'OC1=CC=CC=C1':1,
 'O=C(OC)':2,
 'O=C(Oc1ccccc1)':1,
 'OS(=O)(=O)C':1,
 'N':5,
 'NC':4,
 'N(C)C':3,
 'NC(=O)C':2,
 'N(C)C(=O)C':1,
 'NN':1,
 'N=Nc1ccccc1':1,
 '[N+](=O)[O-]': 2,'[O][N+]=O':2,
 '[N+](=O)(O)[O-]':1,
 'S':3,
 'SC':2,
 'Sc1ccccc1':1,
 'S(=O)(=O)OC':1,
 'S(=O)(=O)Cl':1,
 'C=O':4,
 'C(=O)C':3,
 'C(=O)CC':2,
 'COC(C)(C)C':1,
 'C(=O)c1ccccc1':1,
 'C(=O)O':3,
 'C(=O)OC':2,
 'C(=O)OCC(C)C':1,
 'C(=O)Oc1ccccc1':1,
 'C(=O)N':1,
 'C(=O)Cl':1,
 'C(=O)Br':1,
 'C=Nc1ccccc1':1,
 'C#N':1,
 '[Si](C)(C)C':1,
 'P(O)(OC)OC':1}

#9) Value dictionary for shift calculation
aromatic_dict = {'H': (0, 0, 0),
 'C': (-0.2, -0.12, -0.22),
 'CC': (-0.14, -0.06, -0.17),
 'C(C)C': (-0.13, -0.08, -0.18),
 'C(C)(C)C': (0.02, -0.08, -0.21),
 'CCl': (0.0, 0.0, 0.0),
 'C(F)(F)F': (0.32, 0.14, 0.2),
 'C(Cl)(Cl)Cl': (0.64, 0.13, 0.1),
 'CO': (-0.07, -0.07, -0.07),
 'C=C': (0.06, -0.03, -0.1),
 'C=CC1=CC=CC=C1': (0.15, -0.01, -0.16),
 'C#C': (0.15, -0.02, -0.01),
 'C#CC1=CC=CC=C1': (0.19, 0.02, 0.0),
 'c1ccccc1': (0.37, 0.2, 0.1),
 'F': (-0.26, 0.0, -0.2),
 'Cl': (0.03, -0.02, -0.09),
 'Br': (0.18, -0.08, -0.04),
 'I': (0.39, -0.21, 0.0),
 'O': (-0.56, -0.12, -0.45),
 'OC': (-0.48, -0.09, -0.44),
 'OCC': (-0.46, -0.1, -0.43),
 'OC1=CC=CC=C1': (-0.29, -0.05, -0.23),
 'O=C(OC)': (-0.25, 0.03, -0.13),
 'O=C(Oc1ccccc1)': (-0.09, 0.09, -0.08),
 'OS(=O)(=O)C': (-0.05, 0.07, -0.01),
 'N': (-0.75, -0.25, -0.65),
 'NC': (-0.8, -0.22, -0.68),
 'N(C)C': (-0.66, -0.18, -0.67),
 'NC(=O)C': (0.12, -0.07, -0.28),
 'N(C)C(=O)C': (-0.16, 0.05, -0.02),
 'NN': (-0.6, -0.08, -0.55),
 'N=Nc1ccccc1': (0.67, 0.2, 0.2),
 '[N+](=O)[O-]': (0.58, 0.31, 0.37), '[O][N+]=O':(0.58, 0.31, 0.37),
 '[N+](=O)(O)[O-]': (0.95, 0.26, 0.38),
 'S': (-0.08, -0.16, -0.22),
 'SC': (-0.08, -0.1, -0.24),
 'Sc1ccccc1': (0.06, -0.09, -0.15),
 'S(=O)(=O)OC': (0.6, 0.26, 0.33),
 'S(=O)(=O)Cl': (0.76, 0.35, 0.45),
 'C=O': (0.56, 0.22, 0.29),
 'C(=O)C': (0.62, 0.14, 0.21),
 'C(=O)CC': (0.63, 0.13, 0.2),
 'COC(C)(C)C': (0.44, 0.05, 0.05),
 'C(=O)c1ccccc1': (0.47, 0.13, 0.22),
 'C(=O)O': (0.85, 0.18, 0.27),
 'C(=O)OC': (0.71, 0.11, 0.21),
 'C(=O)OCC(C)C': (0.7, 0.09, 0.19),
 'C(=O)Oc1ccccc1': (0.9, 0.17, 0.27),
 'C(=O)N': (0.61, 0.1, 0.17),
 'C(=O)Cl': (0.84, 0.22, 0.36),
 'C(=O)Br': (0.8, 0.21, 0.37),
 'C=Nc1ccccc1': (0.6, 0.2, 0.2),
 'C#N': (0.36, 0.18, 0.28),
 '[Si](C)(C)C': (0.22, -0.02, -0.02),
 'P(O)(OC)OC': (0.48, 0.16, 0.24)}



def delete_atoms(mol, atom_indices):
    """
    Delete specified atoms from an RDKit molecule.

    Parameters:
    mol (rdkit.Chem.Mol): The molecule from which atoms should be deleted.
    atom_indices (list of int): List of atom indices to delete.

    Returns:
    rdkit.Chem.Mol: The molecule after deletion, or None if the resulting molecule is invalid.
    """
    if mol is None or atom_indices is None or not atom_indices:
        return mol  # Return the original molecule if no indices to delete or mol/atom_indices is None

    editable_mol = Chem.EditableMol(mol)
    try:
        for index in sorted(atom_indices, reverse=True):
            editable_mol.RemoveAtom(index)
        new_mol = editable_mol.GetMol()
        Chem.SanitizeMol(new_mol)
        
        # Check if the molecule is fragmented and handle it
        fragments = Chem.GetMolFrags(new_mol, asMols=True)
        if len(fragments) > 1:
            new_mol = Chem.CombineMols(*fragments)  # Combine fragments if there are multiple
        
        return new_mol
    except Exception as e:
        print(f"Error during atom deletion: {e}")
        return None
    


def search_fctgrps_substituents(target_smiles, smiles_dict):
    """
    Search the mol object for any SMILES encodings present in the dictionary
    in the order of priority defined, repeatedly until no more matches are found.

    Parameters:
    target_smiles (str): A SMILES string representing the target molecule.
    smiles_dict (dict): A dictionary where keys are SMILES strings and values are priority (lower value means higher priority).

    Returns:
    list of dict: A list of dictionaries where each dictionary has the matched SMILES string and the index at which the match was found,
                  or an empty list if no match is found.
    """
    target_mol = Chem.MolFromSmiles(target_smiles)
    if target_mol is None:
        print("Invalid target SMILES string.")
        return []

    matches = []
    atom_map = list(range(target_mol.GetNumAtoms()))  # Keep track of original indices

    # Sort the dictionary by priority (values)
    sorted_smiles = sorted(smiles_dict.items(), key=lambda item: item[1])
    
    while True:
        found_match = False
        for smiles, _ in sorted_smiles:
            query_mol = Chem.MolFromSmiles(smiles)
            if query_mol:
                match_indices = target_mol.GetSubstructMatch(query_mol)
                if match_indices:
                    # Record the match with the original index
                    original_index = atom_map[match_indices[0]]
                    matches.append({smiles: original_index})
                    
                    # Create a copy of the target molecule and delete the matched atoms
                    new_mol = delete_atoms(target_mol, match_indices)
                    
                    if new_mol:
                        target_mol = new_mol  # Continue searching in the modified molecule
                        # Update the atom map by removing the matched indices
                        for index in sorted(match_indices, reverse=True):
                            del atom_map[index]
                        found_match = True
                        break  # Restart the search for the highest priority match in the modified molecule
        
        if not found_match:
            break
    
    return matches



def search_multiple_aromatics(substituents, smiles_dict):
    """
    Searches for the functional groups in smiles_dict in the substituents list of list, works for the case of multiple aromatic rings.

    Parameters:
    substituents (list of list of str): List of lists containing SMILES strings or 'H'.
    smiles_dict (dict): Dictionary where keys are SMILES strings and values are priorities.

    Returns:
    list of list of dict: A list of lists containing dictionaries with matched SMILES strings and their indices.
    """
    results = []  # Start with an empty list of lists

    if not substituents or not smiles_dict:
        return results  # Handle empty inputs

    for i in range(len(substituents)):
        # Ensure results has enough sublists
        while len(results) <= i:
            results.append([])
        for y in range(len(substituents[i])):
            # Ensure each sublist in results has enough elements
            while len(results[i]) <= y:
                results[i].append([{'H': 0}])
            
            if substituents[i][y] == 'H':
                pass
            else:
                result_list_dict = search_fctgrps_substituents(substituents[i][y], smiles_dict)
                # If there are results, append them; otherwise, keep the placeholder
                if result_list_dict:
                    results[i][y] = result_list_dict
                else:
                    results[i][y] = [{'No match': None}]

    return results





def replace_keys(dict_list, value_dict, tuple_index):
    """
    Replace functional group smiles string in the list of matched dictionaries with the value of shift for that specific group at the correct position relative to the hydrogen.

    Parameters:
    dict_list (list of dict): A list of dictionaries where each dictionary has a single key-value pair.
    value_dict (dict): A dictionary where keys are the original keys to be replaced and values are tuples containing new keys.
    tuple_index (int): The index of the tuple in value_dict to use for replacement.

    Returns:
    list of dict: A new list of dictionaries with keys replaced based on value_dict and tuple_index.
    """
    new_dict_list = copy.deepcopy(dict_list)
    
    for dict_item in new_dict_list:
        key, value = next(iter(dict_item.items()))
        if key in value_dict:
            new_key = value_dict[key][tuple_index]
            del dict_item[key]
            dict_item[new_key] = value
            
    return new_dict_list

def sum_distance(dict_list):
    """
    Calculate the sum of key values from a list of dictionaries. If the index is 0, the key value is added directly.
    Otherwise, the key value is divided by the square of the index before being added to the total sum.

    Parameters:
    dict_list (list of dict): A list of dictionaries where keys are numeric values and values are indices.

    Returns:
    float: The total sum calculated based on the provided rules.
    """
    if dict_list is None:
        return 0

    total_sum = 0
    for dict_item in dict_list:
        for key, index in dict_item.items():
            if index == 0:  # If the index is 0, add the key value directly
                total_sum += key
            else:  # Otherwise, add the key value divided by the square of the index
                total_sum += key / (index ** 2)
    return total_sum

def calculation_shift_aromatic(values_dict, matches_one_ring):
    """
    Calculate the aromatic shift values for a ring structure based on substituents.

    Parameters:
    values_dict (dict): A dictionary where keys are substituent names and values are tuples representing their values at different positions.
    matches_one_ring (list of list of dict): A list of lists where each sublist represents a position on the ring, and each dictionary contains a single key-value pair with the substituent and its index.

    Returns:
    list: A list containing the calculated shift values or 'sub' for each position in the ring.
    """
    final_shift = ['','','','','','']
    n = len(matches_one_ring)
    for i in range(len(matches_one_ring)):
        sublist = matches_one_ring[i]
        first_dict = sublist[0]
        first_key = list(first_dict.keys())[0]
        if first_key == 'H':
            ValueDistance_iplus1 = replace_keys(matches_one_ring[(i+1)%n], values_dict, 0)
            ValueDistance_iplus2 = replace_keys(matches_one_ring[(i+2)%n], values_dict, 1)
            ValueDistance_iplus3 = replace_keys(matches_one_ring[(i+3)%n], values_dict, 2)
            ValueDistance_imins1 = replace_keys(matches_one_ring[(i-1)%n], values_dict, 0)
            ValueDistance_imins2 = replace_keys(matches_one_ring[(i-2)%n], values_dict, 1)
            Value_i_plus_1 = sum_distance(ValueDistance_iplus1)
            Value_i_plus_2 = sum_distance(ValueDistance_iplus2)
            Value_i_plus_3 = sum_distance(ValueDistance_iplus3)
            Value_i_mins_1 = sum_distance(ValueDistance_imins1)
            Value_i_mins_2 = sum_distance(ValueDistance_imins2)
            Total_shift_H = round(7.26 + Value_i_plus_1 + Value_i_plus_2 + Value_i_plus_3 + Value_i_mins_1 + Value_i_mins_2, 5)
            final_shift[i] = Total_shift_H
        else:
            final_shift[i] = 'sub'
    return final_shift



def multiple_shift_rings(values_dict, matches):
    """
    Calculate the aromatic shift values for multiple rings based on substituents.

    Parameters:
    values_dict (dict): A dictionary where keys are substituent names and values are tuples representing their values at different positions.
    matches (list of list of list of dict): A list of lists, where each sublist represents a ring and contains lists of dictionaries with substituent matches.

    Returns:
    list of list: A list of lists containing the calculated shift values or 'sub' for each ring.
    """
    final_shift_list = []
    for t in range(len(matches)):
        final_shift_list.append('')

    for y in range(len(matches)):
        ring_y_list = calculation_shift_aromatic(values_dict, matches[y])
        final_shift_list[y] = ring_y_list

    return final_shift_list



def match_smiles_and_mol(smiles: str, aromatic_carbon_indices: list):
    """
    Match the aromatic carbon indices from a SMILES string to the corresponding Mol object indices.

    Parameters:
    smiles (str): A SMILES string representing the molecule.
    aromatic_carbon_indices (list of int): A list of indices of aromatic carbons in the SMILES string.

    Returns:
    list of tuple: A list of tuples where each tuple contains a Mol object index and the corresponding SMILES index.

    Raises:
    ValueError: If the SMILES string is invalid or if the number of provided aromatic carbon indices
                does not match the number of aromatic carbons in the molecule.
    """
    if not smiles:
        return []

    # Convert the SMILES string to a Mol object
    mol = Chem.MolFromSmiles(smiles)
    
    if mol is None:
        raise ValueError("Invalid SMILES string")
    
    # Create a mapping from SMILES index to Mol object index
    mol_indices = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.GetIsAromatic()]
    
    # Verify we have the same number of aromatic carbons
    if len(mol_indices) != len(aromatic_carbon_indices):
        raise ValueError("The number of provided aromatic carbon indices does not match the number of aromatic carbons in the molecule.")
    
    aromatic_carbons = []
    
    # Iterate through the provided aromatic carbon indices
    for smiles_index in aromatic_carbon_indices:
        mol_index = mol_indices[aromatic_carbon_indices.index(smiles_index)]
        aromatic_carbons.append((mol_index, smiles_index))
    
    return aromatic_carbons



def create_shift_dict(shift_list, index_list):
    """
    Create a dictionary with as key the index of the aromatic carbon, that has a hydrogen, in the mol object, and as 
    associated value has the shift value for the hydrogen attached to it.

    Parameters:
    index_list (list of tuple of int): each tuple contains the aromatic carbon indexes in the mol object and smiles respectively.
    shift (list of float): A list of shifts of aromatic hydrogens in the molecule.

    Returns:
    result_dict (int,float): a dictionary as said in the function explanation.

    Raises:
    ValueError: If the SMILES string is invalid or if the number of provided aromatic carbon indices
                does not match the number of aromatic carbons in the molecule.
    """
    result_dict = {}
    for idx, (key, _) in enumerate(index_list):
        if idx < len(shift_list) and isinstance(shift_list[idx], float):
            result_dict[key] = shift_list[idx]
    return result_dict



def main_aromatic(smiles):
    """
    Input: 
        smiles (str): The SMILES string of the molecule.
    
    Output: 
        dict: A dictionary with aromatic carbon indexes as keys and corresponding shift values as values.
    """
    if not smiles:
        return {}

    try:
        # Clean the SMILES string
        clean_smiles = process_smiles(smiles)
    except ValueError as e:
        print(f"Error processing SMILES: {e}")
        return {}

    try:
        # Find indexes of aromatic carbons
        list_c = find_aromatic_carbon_indexes(clean_smiles)
        list_format = [list_c]
        
        # Find substituents on aromatic carbons
        list_subs = find_substituents_aromatic_multiple(clean_smiles, list_format)
        
        # Process the list of substituents
        clean_list_lists = process_smiles_lists(list_subs)
        
        # Identify functional groups
        fct_grp = search_multiple_aromatics(clean_list_lists, data_dict)
        
        # Calculate shift values for aromatic hydrogens
        shift_list = multiple_shift_rings(aromatic_dict, fct_grp)
        
        # Match SMILES and molecule object
        equivalence = match_smiles_and_mol(clean_smiles, list_c)
        
        # Create the shift dictionary
        result_dict = create_shift_dict(shift_list[0], equivalence)
        
        return result_dict
    except Exception as e:
        print(f"An error occurred: {e}")
        return {}
    


def shift(idx : int, smiles : str):
    """
    Calculates the chemical shift of a specified atom in a molecule based on its environment.

    This function calculates the chemical shift of the specified atom in a molecule represented
    by a SMILES string. It considers the atom's local environment, such as nearby functional groups
    or molecular structure, to determine the chemical shift.

    Parameters:
    idx (int): The index of the atom for which to calculate the chemical shift.
    smiles (str): A valid SMILES string representing the molecule.

    Returns:
    float: The calculated chemical shift of the specified atom.

    Example:
    >>> shift(1, 'CC=O')
    9.54
    >>> shift(0, 'CC=O')
    1.81
    """
    table = from_mol_to_shift_figure(smiles)[idx]
    mol = Chem.MolFromSmiles(smiles)
    atom = mol.GetAtomWithIdx(idx)
    
    if table == 1:
        #In the case of an aldehyde.
        if atom.GetSymbol() == 'C':  
            # Iterate over the bonds of the carbon atom
            for bond in atom.GetBonds():
                # Check if the bond is a double bond
                if bond.GetBondType() == Chem.rdchem.BondType.DOUBLE:
                    # Get the other atom in the bond
                    other_atom = bond.GetOtherAtom(atom)
                    # Check if the other atom is oxygen
                    if other_atom.GetSymbol() == 'O':  
                        #+9 is arbitrary but hydrogen aldehydes are approximately in this spectrum of the shift.
                        return round(shift_1(idx, smiles) + 9 ,4) 
        return round(shift_1(idx, smiles) + 0.23,4)
        """
    Strategy changed during the developpment of the code. The shift of aromatic hydrogens are all calculated at the same time with main_aromatic().
    if table == 2:
        return round(shift_2(idx, smiles) + 7.26, 4)
        
    This part could not be done due to lack of time.
    if table == 3:
        #return round(shift_3(idx, smiles) + 5.25, 4)
        """
    if table == 0:
        return round(shift_0(idx, smiles),4)



def info_from_multiplicity(multiplicity):
    """
    Provides information on the expected shape and intensity of peaks in NMR spectra based on multiplicity.

    This function takes the multiplicity of a peak in an NMR spectrum as input and returns two lists:
    - The first list contains the positions where the peaks should be located.
    - The second list contains the relative intensities of the peaks.

    Parameters:
    multiplicity (int): The multiplicity of the peak in the NMR spectrum.

    Returns:
    tuple: A tuple containing two lists:
           - The first list contains the positions of the peaks.
           - The second list contains the relative intensities of the peaks.

    Example:
    >>> info_from_multiplicity(1)
    ([-0.5, 0.5], [1, 1])
    """
    if multiplicity == 0:
        return [0], [1]
    if multiplicity == 1:
        return [-0.5, 0.5], [1,1]
    if multiplicity == 2:
        return [-1, 0, 1], [0.5, 1, 0.5]
    if multiplicity == 3:
        return list(np.arange(-1.5, 1.6)), [0.33, 1, 1, 0.33]
    if multiplicity == 4:
        return list(np.arange(-2,2.1)), [0.25, 0.5, 1, 0.5, 0.25]
    if multiplicity == 5:
        return list(np.arange(-2.5,2.6)), [0.33, 0.625, 1, 1, 0.625, 0.33]
    if multiplicity == 6:
        return list(np.arange(-3,3.1)), [0.09, 0.42, 0.91, 1, 0.91, 0.42, 0.09]
    if multiplicity == 7:
        return list(np.arange(-3.5, 3.6)), [0.06, 0.35, 0.67, 1, 1, 0.67, 0.35, 0.06]
    if multiplicity == 8:
        return list(np.arange(-4, 4.1)), [0.03, 0.23, 0.45, 0.70, 1, 0.7, 0.45, 0.23, 0.03]
    if multiplicity == 9:
        return list(np.arange(-4.5, 4.6)), [0.05, 0.15, 0.58, 0.79, 1, 1, 0.79, 0.58, 0.15, 0.05]






def norm(x, a , x0, sigma): 
    """
    Calculate the value of a normalized Gaussian function at a given point.

    Parameters:
    x (float): The point at which to evaluate the function.
    a (float): The amplitude of the Gaussian function.
    x0 (float): The mean or center of the Gaussian function.
    sigma (float): The standard deviation, controlling the width of the Gaussian curve.

    Returns:
    float: The value of the normalized Gaussian function at the given point `x`.

    Formula:
    The normalized Gaussian function is defined as:
        f(x) = a * exp(-(x - x0)^2 / (2 * sigma^2))

    Example:
    >>> norm(0, 1, 0, 1)
    0.3989422804014327
    """
    return a*np.exp(-(x-x0)**2/(2*sigma**2)) 




def get_name_from_smiles(smiles: str):
    """
    Get the name of a compound from its SMILES string using PubChemPy.

    Parameters:
    smiles (str): The SMILES string of the compound.

    Returns:
    str: The name of the compound obtained from PubChem.

    Raises:
    ValueError: If no compound is found for the given SMILES.
    Exception: If an error occurs during the search process.
    """
    try:
        # Search PubChem for the compound by its SMILES string
        compounds = pcp.get_compounds(smiles, 'smiles')
        if compounds:
            # Return the name of the first matching compound
            return compounds[0].iupac_name
        else:
            raise ValueError("No compound found for the given SMILES.")
    except Exception as e:
        raise Exception(f"An error occurred: {e}")
    


def is_valid_smiles(smiles: str):
    """
    Check if a string is a valid SMILES.

    Parameters:
    smiles (str): The input string to check.

    Returns:
    bool: True if the string is a valid SMILES, False otherwise.
    """
    mol = Chem.MolFromSmiles(smiles)
    return mol is not None




def get_smiles_from_input(input_string: str):
    """
    Convert an input string to a valid RDKit molecule SMILES.

    The input can be a SMILES string or a chemical name.

    Parameters:
    input_string (str): The input string, which can be either a SMILES or a chemical name.

    Returns:
    str: The SMILES string representing the input molecule.

    Raises:
    ValueError: If the input string cannot be converted to a molecule.
    """
    if is_valid_smiles(input_string):
        return input_string
    else:
        try:
            # Treat as chemical name and convert to SMILES then to mol
            smiles = get_smiles(input_string)
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"Cannot convert chemical name to molecule: {input_string}")
            return smiles
        except Exception as e:
            raise ValueError(f"An error occurred: {e}")

def has_aromatic_ring(smiles: str) -> bool:
    """
    Checks if the given SMILES string contains a 6-carbon aromatic ring.
    
    Args:
        smiles (str): The SMILES string of the molecule.
        
    Returns:
        bool: True if there is a 6-carbon aromatic ring, False otherwise.
    """
    # Parse the SMILES string to create a molecule object
    molecule = Chem.MolFromSmiles(smiles)
    
    # Ensure the molecule is valid
    if molecule is None:
        raise ValueError("Invalid SMILES string provided.")
    
    # Find all aromatic rings in the molecule
    aromatic_rings = [ring for ring in Chem.GetSymmSSSR(molecule) if all(molecule.GetAtomWithIdx(idx).GetIsAromatic() for idx in ring)]
    
    # Check if any aromatic ring has exactly 6 carbon atoms
    for ring in aromatic_rings:
        if len(ring) == 6 and all(molecule.GetAtomWithIdx(idx).GetSymbol() == 'C' for idx in ring):
            return True
            
    return False


def NMR(name: str):
    """
    Plot the NMR spectrum of a molecule based on its name or SMILES representation.

    Parameters:
    name (str): The name or SMILES representation of the molecule.

    Returns:
    tuple: A tuple containing the matplotlib.pyplot object and RDKit Mol object.

    Raises:
    Exception: If an error occurs during processing.
    """
    try:        
        # The SMILES of the molecule is acquired.
        crude_smiles = get_smiles_from_input(name)

        # Rearrange the molecule SMILES.
        smiles = canonicalize_smiles(crude_smiles)
        # The IUPAC name is acquired to show it on the graph.
        if is_valid_smiles(name):
            chemical_name = get_name_from_smiles(smiles)
        else:
            chemical_name = name

        # The SMILES is turned into a Mol object.
        mol = Chem.MolFromSmiles(smiles)

        # The multiplicity of each hydrogen is calculated.
        mult, dict_Hs = multiplicity(smiles)
        dict_shift = {}
        # The shift of each hydrogen is calculated.
        for atom in mol.GetAtoms():
            if atom.GetIdx() in list(dict_Hs.keys()) and not atom.GetIsAromatic():
                idx = atom.GetIdx()
                idx_shift = shift(idx, smiles)
                dict_shift[idx] = idx_shift
        if has_aromatic_ring(smiles):
            dict_shift_aromatics = main_aromatic(smiles)
            dict_shift.update(dict_shift_aromatics)
        # The graph is plotted.
        plt.figure(figsize=(8, 5), dpi=150)
        plt.xlabel('δ (ppm)', fontsize=12)
        plt.ylabel('Number of hydrogens', fontsize=12)
        plt.title(f'NMR spectrum of {chemical_name}', fontsize=14, weight='bold')
        plt.grid(False)

        # Add a subtle background color
        plt.gca().set_facecolor('#f9f9f9')

        # Remove spines for cleaner appearance
        plt.gca().spines['top'].set_visible(False)
        plt.gca().spines['right'].set_visible(False)

        # Spacing between the peaks.
        spacing = 0.06

        # How wide should the peaks be calculated for (it is wise for it to be small as bigger values will
        # yield very similar results).
        spread = 0.5

        # How wide the line should be.
        lw = 0.25

        # How wide the peaks should be
        sigma = 0.005

        x_lists = []
        y_lists = []

        for atom in list(mult.keys()):
            r, h = info_from_multiplicity(mult[atom])
            # Considers each peak one after the other.
            for peak in range(0, len(r)):
                # Where each individual peak should be.
                x_basis = np.arange(dict_shift[atom] + r[peak] * spacing - spread, dict_shift[atom] + r[peak] * spacing + spread, 0.0005)
                x_lists.append(x_basis)
                # Calculates each peak individually. The intensity of the biggest peak directly proportional to
                # the number of hydrogens on a heavy atom.
                y = norm(x_basis, 1, dict_shift[atom] + r[peak] * spacing, sigma) * h[peak] * dict_Hs[atom]
                y_lists.append(y)

        summed_y = defaultdict(float)

        for x_list, y_list in zip(x_lists, y_lists):
            for x, y in zip(x_list, y_list):
                summed_y[x] += y

        # Sort the summed results by x values
        summed_x = sorted(summed_y.keys())
        summed_y_values = [summed_y[x] for x in summed_x]
        if summed_y_values:
            ymax = max(summed_y_values)
            xmax = max(summed_x)
        #If no hydrogens
        else:
            ymax = 1
            xmax = 12

        # Plot the summed data
        plt.plot(summed_x, summed_y_values, linewidth=lw, color='black')

        # Set x-axis and y-axis limits. The x axis is constant except for abnormally big shifts.
        if xmax > 11:
            plt.xlim(xmax + 0.1 * xmax, 0)
        else:
            plt.xlim(12, 0)
        plt.ylim(0, ymax + 0.1 * ymax)

        # Draws the molecule considered so it is easy to know if the getsmilesfromsame function messed up or not.
        # It is also handy so as to use the function Show(int) where it will turn red the peak which corresponds
        # to the atom chosen.
        IPythonConsole.drawOptions.addAtomIndices = True
        return plt, mol
    except Exception as e:
        raise Exception(f"An error occurred: {e}")


def Show(name : str, z : int):
    """
    Display the NMR spectrum of a molecule with a highlighted peak for the specified hydrogen index.

    Parameters:
    name (str): The name or SMILES representation of the molecule.
    z (int): The index of the hydrogen atom to highlight in the spectrum.

    Returns:
    tuple: A tuple containing the matplotlib.pyplot object and RDKit Mol object with the highlighted atom.

    Raises:
    ValueError: If the specified index `z` does not correspond to a hydrogen atom.
    """
    # Generate the general NMR plot using the Plot function.
    plt, mol = NMR(name)

   # The SMILES of the molecule is acquired.
    crude_smiles = get_smiles_from_input(name)
    # Rearrange the molecule SMILES.
    smiles = canonicalize_smiles(crude_smiles)
    # The SMILES is turned into a Mol object.
    mol = Chem.MolFromSmiles(smiles)
    mult, dict_Hs = multiplicity(smiles)
    dict_shift = {}

    for atom in mol.GetAtoms():
        if atom.GetIdx() in list(dict_Hs.keys()) and not atom.GetIsAromatic():
                idx = atom.GetIdx()
                idx_shift = shift(idx, smiles)
                dict_shift[idx] = idx_shift
    if has_aromatic_ring(smiles):
        dict_shift_aromatics = main_aromatic(smiles)
        dict_shift.update(dict_shift_aromatics)

    if z not in dict_shift.keys():
        raise ValueError(f'The index {z} does not have a hydrogen')

    # Highlight the specified peak.
    spacing = 0.06
    spread = 0.5
    lw = 0.25
    sigma = 0.005

    ymax = 0
    for i in list(dict_shift.keys()):
        if i == z:
            r, h = info_from_multiplicity(mult[i])
            for j in range(0, len(r)):
                x_basis = np.arange(dict_shift[i] + r[j] * spacing - spread, dict_shift[i] + r[j] * spacing + spread, 0.0005)
                y = norm(x_basis, 1, dict_shift[i] + r[j] * spacing, sigma) * h[j] * dict_Hs[i]
                plt.plot(x_basis, y, color='r', linewidth=lw)

    # Redraw the molecule with highlighted atom indices.
    Zs = get_keys_from_value(dict_shift, dict_shift[z])
    mol.__sssAtoms = [z for z in Zs]

    # Show the plot
    plt.show()
    return plt, mol




def get_keys_from_value(d, value):
    """Return a list of keys from the dictionary `d` that have the specified `value`."""
    return [key for key, val in d.items() if val == value]




