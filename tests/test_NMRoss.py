import pytest
from nmross.NMRoss import *

def test_get_smiles_butanol():
    expected_smiles = 'CCCCO'  # Known SMILES for butanol
    result = get_smiles('butanol')
    assert result == expected_smiles, f"Expected {expected_smiles}, but got {result}"

def test_get_smiles_non_existent():
    expected_message = 'No results found for the given molecule name.'
    result = get_smiles('non_existent_molecule')
    assert result == expected_message, f"Expected '{expected_message}', but got '{result}'"

def test_get_smiles_error_handling():
    # Assuming we want to check if an invalid input is handled properly
    result = get_smiles(None)
    assert result.startswith("An error occurred:"), f"Expected error message, but got '{result}'"

if __name__ == "__main__":
    test_get_smiles_butanol()
    test_get_smiles_non_existent()
    test_get_smiles_error_handling()
    print("All tests passed.")




def test_multiplicity_butanol():
    smiles = 'CCCCO'
    expected_result = (
        {0: 2, 1: 5, 2: 4, 3: 2, 4: 0},
        {0: 3, 1: 2, 2: 2, 3: 2, 4: 1}
    )
    result = multiplicity(smiles)
    assert result == expected_result, f"Expected {expected_result}, but got {result}"

if __name__ == "__main__":
    test_multiplicity_butanol()
    print("All tests passed.")





def test_from_mol_to_shift_figure():
    smiles = 'c1ccc(C=CCO)cc1CCOCCCO'
    expected_output = {
        0: 2, 1: 2, 2: 2, 4: 3, 5: 3, 6: 1, 7: 0,
        8: 2, 10: 1, 11: 1, 13: 1, 14: 1, 15: 1, 16: 0
    }
    result = from_mol_to_shift_figure(smiles)
    assert result == expected_output, f"Expected {expected_output}, but got {result}"

if __name__ == "__main__":
    test_from_mol_to_shift_figure()
    print("All tests passed.")





def test_find_atoms_in_same_ring():
    atom_index = 7
    smiles = 'OCCCOCCc1cccc(C=CCc2ccccc2)c1'
    expected_output = [7, 8, 9, 10, 11, 21]
    result = find_atoms_in_same_ring(atom_index, smiles)
    assert result == expected_output, f"Expected {expected_output}, but got {result}"

if __name__ == "__main__":
    test_find_atoms_in_same_ring()
    print("All tests passed.")





def test_atom_branches():
    atom_index = 1
    smiles = 'CC(C)CC(=O)CC'
    expected_output = [
        {0: 1, 1: 0},  # Branch starting from atom 1 going to atom 0
        {2: 1, 1: 0},  # Branch starting from atom 1 going to atom 2
        {3: 1, 4: 2, 5: 3, 6: 3, 7: 4, 1: 0}  # Branch starting from atom 1 going to atom 7
    ]
    result = atom_branches(atom_index, smiles)
    assert result == expected_output, f"Expected {expected_output}, but got {result}"

if __name__ == "__main__":
    test_atom_branches()
    print("All tests passed.")



def test_map():
    idx = 1
    smiles = 'CC(C)CC(=O)CC'
    expected_output = {0: 1, 1: 0, 2: 1, 3: 1, 4: 2, 5: 3, 6: 3, 7: 4}
    result = map(idx, smiles)
    assert result == expected_output, f"Expected {expected_output}, but got {result}"

if __name__ == "__main__":
    test_map()
    print("All tests passed.")





def test_branches_to_smiles():
    smiles = 'CC(C)CC(=O)CCc1ccccc1'
    branches = atom_branches(1, smiles)
    expected_output = ['CC', 'CC', 'cCCC(=O)CC']
    result = branches_to_smiles(branches, smiles)
    assert sorted(result) == sorted(expected_output), f"Expected {expected_output}, but got {result}"

if __name__ == "__main__":
    test_branches_to_smiles()
    print("All tests passed.")




def test_clean_aromatics():
    smiles = 'CCC(=O)N(c1c(CCC)cccc1)C1C(Cl)CN(CCc2cc(I)cc(CCO)c2)CC1'
    expected_output = ['CC', 'CC(=O)N(c1ccccc1)C1CCN(CCc2ccccc2)CC1Cl']
    result = clean_aromatics(1, smiles)
    assert result == expected_output, f"Expected {expected_output}, but got {result}"

if __name__ == "__main__":
    test_clean_aromatics()
    print("All tests passed.")




def test_shift_oxygen_with_hydrogen():
    assert shift_0(2, 'CCO') == 2.46

def test_shift_nitrogen_with_hydrogen():
    assert shift_0(0, 'N') == 2.4

def test_shift_sulfur_with_hydrogen():
    assert shift_0(0, 'S') == 0.72

def test_no_hydrogen_bonded():
    expected_message = 'An error occurred, atom with symbol: F does not have a bonded hydrogen'
    assert shift_0(2, 'CCF') == expected_message

def test_no_value_for_symbol():
    expected_message = 'An error occurred, no values of shift_0 for atom with symbol: P'
    assert shift_0(2, 'CCP') == expected_message

if __name__ == "__main__":
    test_shift_oxygen_with_hydrogen()
    test_shift_nitrogen_with_hydrogen()
    test_shift_sulfur_with_hydrogen()
    test_no_hydrogen_bonded()
    test_no_value_for_symbol()
    print("All tests passed.")





def test_search_algo_with_central_carbon():
    smiles = 'CCCC(=O)O'
    result = search_algo(0, smiles)
    expected = [{1: ['C'], 2: ['C'], 3: ['C(=O)O']}]
    assert result == expected, f"Expected {expected}, but got {result}"

def test_search_algo_no_neighbors():
    smiles = 'C'
    result = search_algo(0, smiles)
    expected = []
    assert result == expected, f"Expected {expected}, but got {result}"

def test_search_algo_with_formaldehyde():
    smiles = 'C=O'
    result = search_algo(0, smiles)
    expected = 0
    assert result == expected, f"Expected {expected}, but got {result}"

if __name__ == "__main__":
    test_search_algo_with_central_carbon()
    test_search_algo_no_neighbors()
    test_search_algo_with_formaldehyde()
    print("All tests passed.")




def test_shift_1_butyric_acid():
    smiles = 'CCCC(=O)O'
    result = shift_1(0, smiles)
    expected = 0.835
    assert abs(result - expected) < 1e-3, f"Expected {expected}, but got {result}"

def test_shift_1_propanol():
    smiles = 'CCCO'
    result = shift_1(1, smiles)
    expected = 1.695  # Example expected value, adjust as needed based on actual data
    assert abs(result - expected) < 1e-3, f"Expected {expected}, but got {result}"

def test_shift_1_ethanol():
    smiles = 'CCO'
    result = shift_1(1, smiles)
    expected = 3.0  # Example expected value, adjust as needed based on actual data
    assert abs(result - expected) < 1e-3, f"Expected {expected}, but got {result}"

def test_shift_1_methanol():
    smiles = 'CO'
    result = shift_1(0, smiles)
    expected = 2.46  # Example expected value, adjust as needed based on actual data
    assert abs(result - expected) < 1e-3, f"Expected {expected}, but got {result}"

def test_shift_1_no_neighbors():
    smiles = 'C'
    result = shift_1(0, smiles)
    expected = 0.0  # Example expected value, adjust as needed based on actual data
    assert abs(result - expected) < 1e-3, f"Expected {expected}, but got {result}"

if __name__ == "__main__":
    test_shift_1_butyric_acid()
    test_shift_1_propanol()
    test_shift_1_ethanol()
    test_shift_1_methanol()
    test_shift_1_no_neighbors()
    print("All tests passed.")




def test_valid_smiles():
    assert process_smiles("C/C=C/C") == "CC=CC"

def test_valid_smiles_with_hyphen():
    assert process_smiles("C/C=C/C-C") == "CC=CCC"

def test_invalid_smiles():
    try:
        process_smiles("InvalidSMILES")
    except ValueError:
        pass
    else:
        assert False, "Expected ValueError"

def test_non_string_input():
    try:
        process_smiles(12345)
    except ValueError:
        pass
    else:
        assert False, "Expected ValueError"

def test_empty_string():
    try:
        process_smiles("")
    except ValueError:
        pass
    else:
        assert False, "Expected ValueError"

if __name__ == "__main__":
    test_valid_smiles()
    test_valid_smiles_with_hyphen()
    test_invalid_smiles()
    test_non_string_input()
    test_empty_string()
    print("All tests passed.")




def test_find_aromatic_carbon_indexes_benzene():
    assert find_aromatic_carbon_indexes("c1ccccc1") == [0, 2, 3, 4, 5, 6]

def test_find_aromatic_carbon_indexes_toluene():
    assert find_aromatic_carbon_indexes("c1ccccc1C") == [0, 2, 3, 4, 5, 6]

def test_invalid_smiles():
    try:
        find_aromatic_carbon_indexes("InvalidSMILES")
    except ValueError:
        pass
    else:
        assert False, "Expected ValueError"

def test_non_string_input():
    try:
        find_aromatic_carbon_indexes(12345)
    except ValueError:
        pass
    else:
        assert False, "Expected ValueError"

def test_empty_string():
    try:
        find_aromatic_carbon_indexes("")
    except ValueError:
        pass
    else:
        assert False, "Expected ValueError"

if __name__ == "__main__":
    test_find_aromatic_carbon_indexes_benzene()
    test_find_aromatic_carbon_indexes_toluene()
    test_invalid_smiles()
    test_non_string_input()
    test_empty_string()
    print("All tests passed.")




def test_single_aromatic_cycle():
    assert find_substituents_aromatic_multiple("c1ccccc1", [[0, 2, 3, 4, 5, 6]]) == [['H', 'H', 'H', 'H', 'H', 'H']]

def test_single_aromatic_cycle_with_substituents():
    assert find_substituents_aromatic_multiple("c1ccccc1C", [[0, 2, 3, 4, 5, 6]]) == [['H', 'H', 'H', 'H', 'H', 'C']]

def test_multiple_aromatic_cycles():
    assert find_substituents_aromatic_multiple("c1ccccc1c2ccccc2", [[0, 2, 3, 4, 5, 6], [8, 10, 11, 12, 13, 14]]) == [['H', 'H', 'H', 'H', 'H', 'c2ccccc2'], ['c1ccccc1', 'H', 'H', 'H', 'H', 'H']]

def test_empty_string():
    try:
        find_substituents_aromatic_multiple("", [[0, 2, 3, 4, 5, 6]])
    except ValueError:
        pass
    else:
        assert False, "Expected ValueError"

def test_non_string_input():
    try:
        find_substituents_aromatic_multiple(12345, [[0, 2, 3, 4, 5, 6]])
    except ValueError:
        pass
    else:
        assert False, "Expected ValueError"

if __name__ == "__main__":
    test_single_aromatic_cycle()
    test_single_aromatic_cycle_with_substituents()
    test_multiple_aromatic_cycles()
    test_empty_string()
    test_non_string_input()
    print("All tests passed.")





def test_canonicalize_smiles():
    assert process_smiles_lists([['C1CCC1', 'H', 'C2CCC2']]) == [['C1CCC1', 'H', 'C1CCC1']]

def test_invalid_smiles():
    assert process_smiles_lists([['(Invalid)', 'H', 'C2CCC2']]) == [['Invalid', 'H', 'C1CCC1']]

def test_empty_string():
    assert process_smiles_lists([['', 'H', 'C2CCC2']]) == [['', 'H', 'C1CCC1']]

def test_empty_list():
    assert process_smiles_lists([[]]) == [[]]

def test_non_string_input():
    try:
        process_smiles_lists([['C1CCC1', 'H', 12345]])
    except TypeError:
        pass
    else:
        assert False, "Expected TypeError"

if __name__ == "__main__":
    test_canonicalize_smiles()
    test_invalid_smiles()
    test_empty_string()
    test_empty_list()
    test_non_string_input()
    print("All tests passed.")




def test_invalid_deletion():
    mol = Chem.MolFromSmiles("CCO")
    new_mol = delete_atoms(mol, [10])  # Invalid index
    assert new_mol is None

def test_empty_atom_indices():
    mol = Chem.MolFromSmiles("CCO")
    new_mol = delete_atoms(mol, [])  # No atoms to delete
    assert new_mol is not None
    assert Chem.MolToSmiles(new_mol) == "CCO"

def test_none_molecule():
    new_mol = delete_atoms(None, [1])  # None molecule
    assert new_mol is None

if __name__ == "__main__":
    test_invalid_deletion()
    test_empty_atom_indices()
    test_none_molecule()
    print("All tests passed.")





def test_valid_target_smiles():
    target_smiles = "CCOCCO"
    smiles_dict = {"CCO": 1, "O": 2}
    matches = search_fctgrps_substituents(target_smiles, smiles_dict)
    expected_matches = [{"CCO": 0}, {"CCO": 3}]
    assert matches == expected_matches, f"Expected {expected_matches}, but got {matches}"

def test_invalid_target_smiles():
    target_smiles = "InvalidSMILES"
    smiles_dict = {"CCO": 1, "O": 2}
    matches = search_fctgrps_substituents(target_smiles, smiles_dict)
    assert matches == [], f"Expected [], but got {matches}"

def test_empty_smiles_dict():
    target_smiles = "CCOCCO"
    smiles_dict = {}
    matches = search_fctgrps_substituents(target_smiles, smiles_dict)
    assert matches == [], f"Expected [], but got {matches}"

def test_no_matches():
    target_smiles = "CCOCCO"
    smiles_dict = {"N": 1, "Cl": 2}
    matches = search_fctgrps_substituents(target_smiles, smiles_dict)
    assert matches == [], f"Expected [], but got {matches}"

def test_partial_matches():
    target_smiles = "CCOCCO"
    smiles_dict = {"CC": 1, "OCC": 2}
    matches = search_fctgrps_substituents(target_smiles, smiles_dict)
    expected_matches = [{"CC": 0}, {"CC": 3}]
    assert matches == expected_matches, f"Expected {expected_matches}, but got {matches}"

def test_full_overlap():
    target_smiles = "CCOCCO"
    smiles_dict = {"CCOCCO": 1}
    matches = search_fctgrps_substituents(target_smiles, smiles_dict)
    expected_matches = [{"CCOCCO": 0}]
    assert matches == expected_matches, f"Expected {expected_matches}, but got {matches}"

if __name__ == "__main__":
    test_valid_target_smiles()
    test_invalid_target_smiles()
    test_empty_smiles_dict()
    test_no_matches()
    test_partial_matches()
    test_full_overlap()
    print("All tests passed.")




def test_valid_substituents():
    substituents = [['c1ccccc1', 'H', 'c1ccccc1C']]
    smiles_dict = {'c1ccccc1': 1, 'C': 2}
    results = search_multiple_aromatics(substituents, smiles_dict)
    expected_results = [[[{'c1ccccc1': 0}], [{'H': 0}], [{'c1ccccc1': 0}, {'C': 6}]]]
    assert results == expected_results, f"Expected {expected_results}, but got {results}"

def test_empty_substituents():
    substituents = []
    smiles_dict = {'c1ccccc1': 1, 'C': 2}
    results = search_multiple_aromatics(substituents, smiles_dict)
    assert results == [], f"Expected [], but got {results}"

def test_no_matches():
    substituents = [['c1ccccc1', 'H', 'c1ccccc1C']]
    smiles_dict = {'N': 1, 'Cl': 2}
    results = search_multiple_aromatics(substituents, smiles_dict)
    expected_results = [[[{'No match': None}], [{'H': 0}], [{'No match': None}]]]
    assert results == expected_results, f"Expected {expected_results}, but got {results}"

def test_mixed_substituents():
    substituents = [['c1ccccc1', 'H'], ['C1=CC=CC=C1']]
    smiles_dict = {'c1ccccc1': 1, 'C1=CC=CC=C1': 2}
    results = search_multiple_aromatics(substituents, smiles_dict)
    expected_results = [[[{'c1ccccc1': 0}], [{'H': 0}]], [[{'c1ccccc1': 0}]]]
    assert results == expected_results, f"Expected {expected_results}, but got {results}"

if __name__ == "__main__":
    test_valid_substituents()
    test_empty_substituents()
    test_no_matches()
    test_mixed_substituents()
    print("All tests passed.")





def test_basic_shift():
    values_dict = {
        'H': (0, 0, 0),
        'Cl': (1, 0.5, 0.2),
        'NO2': (2, 1, 0.5)
    }
    matches_one_ring = [
        [{'H': 0}],
        [{'Cl': 1}],
        [{'NO2': 2}],
        [{'H': 0}],
        [{'Cl': 1}],
        [{'NO2': 2}]
    ]
    result = calculation_shift_aromatic(values_dict, matches_one_ring)
    expected = [9.51, 'sub', 'sub', 9.51, 'sub', 'sub']
    assert result == expected, f"Expected {expected}, but got {result}"

def test_no_substituents():
    values_dict = {
        'H': (0, 0, 0)
    }
    matches_one_ring = [
        [{'H': 0}],
        [{'H': 0}],
        [{'H': 0}],
        [{'H': 0}],
        [{'H': 0}],
        [{'H': 0}]
    ]
    result = calculation_shift_aromatic(values_dict, matches_one_ring)
    expected = [7.26, 7.26, 7.26, 7.26, 7.26, 7.26]
    assert result == expected, f"Expected {expected}, but got {result}"

def test_mixed_substituents():
    values_dict = {
        'H': (0, 0, 0),
        'Cl': (1, 0.5, 0.2),
        'NO2': (2, 1, 0.5)
    }
    matches_one_ring = [
        [{'H': 0}],
        [{'H': 0}],
        [{'Cl': 1}],
        [{'H': 0}],
        [{'NO2': 2}],
        [{'H': 0}]
    ]
    result = calculation_shift_aromatic(values_dict, matches_one_ring)
    expected = [8.01, 8.385, 'sub', 8.76, 'sub', 7.96]
    assert result == expected, f"Expected {expected}, but got {result}"

def test_empty_ring():
    values_dict = {
        'H': (0, 0, 0),
        'Cl': (1, 0.5, 0.2)
    }
    matches_one_ring = []
    result = calculation_shift_aromatic(values_dict, matches_one_ring)
    expected = ['', '', '', '', '', '']
    assert result == expected, f"Expected {expected}, but got {result}"

if __name__ == "__main__":
    test_basic_shift()
    test_no_substituents()
    test_mixed_substituents()
    test_empty_ring()
    print("All tests passed.")





def test_no_substituents():
    values_dict = {
        'H': (0, 0, 0)
    }
    matches = [
        [
            [{'H': 0}],
            [{'H': 0}],
            [{'H': 0}],
            [{'H': 0}],
            [{'H': 0}],
            [{'H': 0}]
        ]
    ]
    result = multiple_shift_rings(values_dict, matches)
    expected = [
        [7.26, 7.26, 7.26, 7.26, 7.26, 7.26]
    ]
    assert result == expected

def test_empty_matches():
    values_dict = {
        'H': (0, 0, 0),
        'Cl': (1, 0.5, 0.2)
    }
    matches = []
    result = multiple_shift_rings(values_dict, matches)
    expected = []
    assert result == expected

if __name__ == "__main__":
    test_no_substituents()
    test_empty_matches()
    print("All tests passed.")





def test_valid_smiles_and_indices():
    smiles = "c1ccccc1"
    aromatic_carbon_indices = [0, 1, 2, 3, 4, 5]
    result = match_smiles_and_mol(smiles, aromatic_carbon_indices)
    expected = [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5)]
    assert result == expected, f"Expected {expected}, but got {result}"

def test_invalid_smiles():
    smiles = "InvalidSMILES"
    aromatic_carbon_indices = [0, 1, 2, 3, 4, 5]
    try:
        match_smiles_and_mol(smiles, aromatic_carbon_indices)
    except ValueError:
        pass
    else:
        assert False, "Expected ValueError"

def test_mismatched_indices():
    smiles = "c1ccccc1"
    aromatic_carbon_indices = [0, 1, 2, 3]  # Not enough indices
    try:
        match_smiles_and_mol(smiles, aromatic_carbon_indices)
    except ValueError:
        pass
    else:
        assert False, "Expected ValueError"

def test_empty_smiles():
    smiles = ""
    aromatic_carbon_indices = []
    result = match_smiles_and_mol(smiles, aromatic_carbon_indices)
    expected = []
    assert result == expected, f"Expected {expected}, but got {result}"

def test_no_aromatic_carbons():
    smiles = "CCO"
    aromatic_carbon_indices = []
    result = match_smiles_and_mol(smiles, aromatic_carbon_indices)
    expected = []
    assert result == expected, f"Expected {expected}, but got {result}"

if __name__ == "__main__":
    test_valid_smiles_and_indices()
    test_invalid_smiles()
    test_mismatched_indices()
    test_empty_smiles()
    test_no_aromatic_carbons()
    print("All tests passed.")





def test_toluene():
    expected_result = {2: 7.06, 3: 7.14, 4: 7.04, 5: 7.14, 6: 7.06}
    result = main_aromatic("Cc1ccccc1")
    assert result == expected_result, f"Expected {expected_result}, but got {result}"

def test_benzene():
    expected_result = {0: 7.26, 1: 7.26, 2: 7.26, 3: 7.26, 4: 7.26, 5: 7.26}
    result = main_aromatic("c1ccccc1")
    assert result == expected_result, f"Expected {expected_result}, but got {result}"

def test_aniline():
    expected_result = {2: 6.51, 3: 7.01, 4: 6.61, 5: 7.01, 6: 6.51}
    result = main_aromatic("Nc1ccccc1")
    assert result == expected_result, f"Expected {expected_result}, but got {result}"

def test_invalid_name():
    expected_result = {}  # Expected result for an invalid molecule name
    result = main_aromatic("invalid_molecule")
    assert result == expected_result, f"Expected {expected_result}, but got {result}"

def test_empty_name():
    expected_result = {}  # Expected result for an empty molecule name
    result = main_aromatic("")
    assert result == expected_result, f"Expected {expected_result}, but got {result}"

if __name__ == "__main__":
    test_toluene()
    test_benzene()
    test_aniline()
    test_invalid_name()
    test_empty_name()
    print("All tests passed.")



def test_shift_aldehyde():
    smiles = 'CC=O'
    result = shift(1, smiles)
    expected = 9.54
    assertAlmostEqual(result, expected, places=2)

def test_shift_methyl_group():
    smiles = 'CC=O'
    result = shift(0, smiles)
    expected = 1.81
    assertAlmostEqual(result, expected, places=2)

def test_shift_hydroxyl_group():
    smiles = 'CCO'
    result = shift(2, smiles)
    expected = 2.46  # Example expected value, adjust as needed based on actual data
    assertAlmostEqual(result, expected, places=2)

def test_shift_ammonia():
    smiles = 'CCN'
    result = shift(2, smiles)
    expected = 2.4  # Example expected value, adjust as needed based on actual data
    assertAlmostEqual(result, expected, places=2)

def test_shift_hydrogen_sulfide():
    smiles = 'CCS'
    result = shift(2, smiles)
    expected = 0.72  # Example expected value, adjust as needed based on actual data
    assertAlmostEqual(result, expected, places=2)

def assertAlmostEqual(actual, expected, places):
    assert round(actual - expected, places) == 0, f"Expected {expected}, but got {actual}"

if __name__ == "__main__":
    test_shift_aldehyde()
    test_shift_methyl_group()
    test_shift_hydroxyl_group()
    test_shift_ammonia()
    test_shift_hydrogen_sulfide()
    print("All tests passed.")





def test_multiplicity_0():
    positions, intensities = info_from_multiplicity(0)
    assert positions == [0]
    assert intensities == [1]

def test_multiplicity_1():
    positions, intensities = info_from_multiplicity(1)
    assert positions == [-0.5, 0.5]
    assert intensities == [1, 1]

def test_multiplicity_2():
    positions, intensities = info_from_multiplicity(2)
    assert positions == [-1, 0, 1]
    assert intensities == [0.5, 1, 0.5]

def test_multiplicity_3():
    positions, intensities = info_from_multiplicity(3)
    assert positions == list(np.arange(-1.5, 1.6))
    assert intensities == [0.33, 1, 1, 0.33]

def test_multiplicity_4():
    positions, intensities = info_from_multiplicity(4)
    assert positions == list(np.arange(-2, 2.1))
    assert intensities == [0.25, 0.5, 1, 0.5, 0.25]

def test_multiplicity_5():
    positions, intensities = info_from_multiplicity(5)
    assert positions == list(np.arange(-2.5, 2.6))
    assert intensities == [0.33, 0.625, 1, 1, 0.625, 0.33]

def test_multiplicity_8():
    positions, intensities = info_from_multiplicity(8)
    assert positions == list(np.arange(-4, 4.1))
    assert intensities == [0.03, 0.23, 0.45, 0.70, 1, 0.7, 0.45, 0.23, 0.03]

if __name__ == "__main__":
    test_multiplicity_0()
    test_multiplicity_1()
    test_multiplicity_2()
    test_multiplicity_3()
    test_multiplicity_4()
    test_multiplicity_5()
    test_multiplicity_8()
    print("All tests passed.")





def test_norm():
    # Test parameters
    a = 1
    x0 = 0
    sigma = 1
    step = 0.5
    x_values = np.arange(-1, 1.1, step)
    expected_values = [norm(x, a, x0, sigma) for x in x_values]

    # Expected results calculated using the Gaussian function formula
    expected_results = [
        0.6065306597126334,  # norm(-1, 1, 0, 1)
        0.8824969025845955,  # norm(-0.5, 1, 0, 1)
        1.0,                 # norm(0, 1, 0, 1)
        0.8824969025845955,  # norm(0.5, 1, 0, 1)
        0.6065306597126334   # norm(1, 1, 0, 1)
    ]

    # Assert each calculated value is approximately equal to the expected result
    for calculated, expected in zip(expected_values, expected_results):
        assertAlmostEqual(calculated, expected, places=7)


if __name__ == "__main__":
    test_norm()
    print("All tests passed.")





def test_get_name_from_smiles():
    # Test with a valid SMILES
    smiles_valid = 'CCO'
    assert get_name_from_smiles(smiles_valid) == 'ethanol'

def test_is_invalid_smiles():
    # Test with an invalid SMILES
    smiles_invalid = 'invalid_smiles'
    assert not is_valid_smiles(smiles_invalid)
    
def test_is_valid_smiles():
    # Test with a valid SMILES
    smiles_valid = 'CCO'
    assert is_valid_smiles(smiles_valid)

    # Test with an invalid SMILES
    smiles_invalid = 'invalid_smiles'
    assert not is_valid_smiles(smiles_invalid)

def test_get_smiles_from_input():
    # Test with a valid SMILES
    input_valid_smiles = 'CCO'
    assert get_smiles_from_input(input_valid_smiles) == 'CCO'

    # Test with a valid chemical name
    input_valid_name = 'Ethanol'
    assert get_smiles_from_input(input_valid_name) == 'CCO'

    # Test with an invalid input
    input_invalid = 'invalid_input'
    try:
        get_smiles_from_input(input_invalid)
        assert False, "Expected ValueError"
    except ValueError:
        pass

if __name__ == "__main__":
    test_get_name_from_smiles()
    test_is_invalid_smiles()
    test_is_valid_smiles()
    test_get_smiles_from_input()
    print("All tests passed.")




def test_nmr_function():
    # Test with a valid SMILES representation
    smiles = 'CCO'
    plt, mol = NMR(smiles)
    assert plt is not None
    assert mol is not None

if __name__ == "__main__":
    test_nmr_function()
    print("All tests passed.")





def test_get_keys_from_value():
    # Test case with a dictionary having multiple keys with the same value
    d = {'a': 1, 'b': 2, 'c': 1, 'd': 3}
    value = 1
    expected_result = ['a', 'c']
    assert get_keys_from_value(d, value) == expected_result

    # Test case with an empty dictionary
    d = {}
    value = 5
    expected_result = []
    assert get_keys_from_value(d, value) == expected_result

    # Test case with a dictionary having keys but no matching value
    d = {'x': 1, 'y': 2, 'z': 3}
    value = 5
    expected_result = []
    assert get_keys_from_value(d, value) == expected_result

if __name__ == "__main__":
    test_get_keys_from_value()
    print("All tests passed.")



def test_show_function():
    # Test case with a valid molecule name and valid hydrogen index
    name = 'CCO'
    z = 1
    plt, mol = Show(name, z)
    assert plt is not None
    assert mol is not None

    # Test case with a valid molecule name but invalid hydrogen index
    name = 'CCO'
    z = 5  # Index 5 doesn't have a hydrogen in this molecule
    try:
        Show(name, z)
        assert False, "Expected ValueError"
    except ValueError:
        pass

    # Test case with an invalid molecule name
    name = 'Invalid_Molecule_Name'
    z = 1
    try:
        Show(name, z)
        assert False, "Expected Exception"
    except Exception:
        pass

if __name__ == "__main__":
    test_show_function()
    print("All tests passed.")
