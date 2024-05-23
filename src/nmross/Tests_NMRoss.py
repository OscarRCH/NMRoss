import unittest

class TestSmilesFetcher(unittest.TestCase):

    def test_get_smiles_butanol(self):
        expected_smiles = 'CCCCO'  # Known SMILES for butanol
        result = get_smiles('butanol')
        self.assertEqual(result, expected_smiles)

    def test_get_smiles_non_existent(self):
        expected_message = 'No results found for the given molecule name.'
        result = get_smiles('non_existent_molecule')
        self.assertEqual(result, expected_message)

    def test_get_smiles_error_handling(self):
        # Assuming we want to check if an invalid input is handled properly
        result = get_smiles(None)
        self.assertTrue(result.startswith("An error occurred:"))




class TestMultiplicity(unittest.TestCase):

    def test_multiplicity_butanol(self):
        smiles = 'CCCCO'
        expected_result = (
            {0: 2, 1: 5, 2: 4, 3: 2, 4: 0},
            {0: 3, 1: 2, 2: 2, 3: 2, 4: 1}
        )
        result = multiplicity(smiles)
        self.assertEqual(result, expected_result)




class TestFromMolToShiftFigure(unittest.TestCase):

    def test_from_mol_to_shift_figure(self):
        smiles = 'c1ccc(C=CCO)cc1CCOCCCO'
        expected_output = {
            0: 2, 1: 2, 2: 2, 4: 3, 5: 3, 6: 1, 7: 0,
            8: 2, 10: 1, 11: 1, 13: 1, 14: 1, 15: 1, 16: 0
        }
        result = from_mol_to_shift_figure(smiles)
        self.assertEqual(result, expected_output)




class TestFindAtomsInSameRing(unittest.TestCase):

    def test_find_atoms_in_same_ring(self):
        atom_index = 7
        smiles = 'OCCCOCCc1cccc(C=CCc2ccccc2)c1'
        expected_output = [7, 8, 9, 10, 11, 21]
        result = find_atoms_in_same_ring(atom_index, smiles)
        self.assertEqual(result, expected_output)




class TestAtomBranches(unittest.TestCase):

    def test_atom_branches(self):
        atom_index = 1
        smiles = 'CC(C)CC(=O)CC'
        expected_output = [
            {0: 1, 1: 0},  # Branch starting from atom 1 going to atom 0
            {2: 1, 1: 0},  # Branch starting from atom 1 going to atom 2
            {3: 1, 4: 2, 5: 3, 6: 3, 7: 4, 1: 0}  # Branch starting from atom 1 going to atom 7
        ]
        result = atom_branches(atom_index, smiles)
        self.assertEqual(result, expected_output)




class TestMapFunction(unittest.TestCase):

    def test_map(self):
        idx = 1
        smiles = 'CC(C)CC(=O)CC'
        expected_output = {0: 1, 1: 0, 2: 1, 3: 1, 4: 2, 5: 3, 6: 3, 7: 4}
        result = map(idx, smiles)
        self.assertEqual(result, expected_output)




class TestBranchesToSmiles(unittest.TestCase):

    def test_branches_to_smiles(self):
        smiles = 'CC(C)CC(=O)CCc1ccccc1'
        branches = atom_branches(1, smiles)
        expected_output = ['CC', 'CC', 'cCCC(=O)CC']
        result = branches_to_smiles(branches, smiles)
        self.assertEqual(sorted(result), sorted(expected_output))




class TestCleanAromatics(unittest.TestCase):
    def test_clean_aromatics(self):
        smiles = 'CCC(=O)N(c1c(CCC)cccc1)C1C(Cl)CN(CCc2cc(I)cc(CCO)c2)CC1'
        expected_output = ['CC', 'CC(=O)N(c1ccccc1)C1CCN(CCc2ccccc2)CC1Cl']
        result = clean_aromatics(1, smiles)
        self.assertEqual(result, expected_output)




class TestShift0(unittest.TestCase):

    def test_shift_oxygen_with_hydrogen(self):
        self.assertEqual(shift_0(2, 'CCO'), 2.46)

    def test_shift_nitrogen_with_hydrogen(self):
        self.assertEqual(shift_0(0, 'N'), 7.8)

    def test_shift_sulfur_with_hydrogen(self):
        self.assertEqual(shift_0(0, 'S'), 0.72)

    def test_no_hydrogen_bonded(self):
        self.assertEqual(shift_0(2, 'CCF'), 'An error occurred, atom with symbol: F does not have a bonded hydrogen')

    def test_no_value_for_symbol(self):
        self.assertEqual(shift_0(2, 'CCP'), 'An error occurred, no values of shift_0 for atom with symbol: P')




class TestSearchAlgo(unittest.TestCase):

    def test_search_algo_with_central_carbon(self):
        smiles = 'CCCC(=O)O'
        result = search_algo(0, smiles)
        expected = [{1: ['C'], 2: ['C'], 3: ['C(=O)O']}]
        self.assertEqual(result, expected)

    def test_search_algo_no_neighbors(self):
        smiles = 'C'
        result = search_algo(0, smiles)
        expected = []
        self.assertEqual(result, expected)

    def test_search_algo_with_formaldehyde(self):
        smiles = 'C=O'
        result = search_algo(0, smiles)
        expected = 0
        self.assertEqual(result, expected)



class TestShift1(unittest.TestCase):

    def test_shift_1_butyric_acid(self):
        smiles = 'CCCC=O'
        result = shift_1(1, smiles)
        expected = 1.475
        self.assertAlmostEqual(result, expected, places=3)

    def test_shift_1_propanol(self):
        smiles = 'CCCO'
        result = shift_1(1, smiles)
        expected = 1.695  # Example expected value, adjust as needed based on actual data
        self.assertAlmostEqual(result, expected, places=3)

    def test_shift_1_ethanol(self):
        smiles = 'CCO'
        result = shift_1(1, smiles)
        expected = 3.0  # Example expected value, adjust as needed based on actual data
        self.assertAlmostEqual(result, expected, places=3)

    def test_shift_1_methanol(self):
        smiles = 'CO'
        result = shift_1(0, smiles)
        expected = 2.46  # Example expected value, adjust as needed based on actual data
        self.assertAlmostEqual(result, expected, places=3)

    def test_shift_1_no_neighbors(self):
        smiles = 'C'
        result = shift_1(0, smiles)
        expected = 0.0  # Example expected value, adjust as needed based on actual data
        self.assertAlmostEqual(result, expected, places=3)




class TestProcessSmiles(unittest.TestCase):
    def test_valid_smiles(self):
        self.assertEqual(process_smiles("C/C=C/C"), "CC=CC")
    
    def test_valid_smiles_with_hyphen(self):
        self.assertEqual(process_smiles("C/C=C/C-C"), "CC=CCC")
    
    def test_invalid_smiles(self):
        with self.assertRaises(ValueError):
            process_smiles("InvalidSMILES")
    
    def test_non_string_input(self):
        with self.assertRaises(ValueError):
            process_smiles(12345)
    
    def test_empty_string(self):
        with self.assertRaises(ValueError):
            process_smiles("")




class TestFindAromaticCarbonIndexes(unittest.TestCase):
    def test_valid_smiles(self):
        self.assertEqual(find_aromatic_carbon_indexes("c1ccccc1"), [[0, 2, 3, 4, 5, 6]])
    
    def test_valid_smiles_with_substituent(self):
        self.assertEqual(find_aromatic_carbon_indexes("c1ccccc1C"), [[0, 2, 3, 4, 5, 6]])
    
    def test_invalid_smiles(self):
        with self.assertRaises(ValueError):
            find_aromatic_carbon_indexes("InvalidSMILES")
    
    def test_non_string_input(self):
        with self.assertRaises(ValueError):
            find_aromatic_carbon_indexes(12345)
    
    def test_empty_string(self):
        with self.assertRaises(ValueError):
            find_aromatic_carbon_indexes("")




class TestFindSubstituentsAromaticMultiple(unittest.TestCase):
    def test_single_aromatic_cycle(self):
        self.assertEqual(
            find_substituents_aromatic_multiple("c1ccccc1", [[0, 2, 3, 4, 5, 6]]),
            [['H', 'H', 'H', 'H', 'H', 'H']]
        )
    
    def test_single_aromatic_cycle_with_substituents(self):
        self.assertEqual(
            find_substituents_aromatic_multiple("c1ccccc1C", [[0, 2, 3, 4, 5, 6]]),
            [['H', 'H', 'H', 'H', 'H', 'C']]
        )
    
    def test_multiple_aromatic_cycles(self):
        self.assertEqual(
            find_substituents_aromatic_multiple("c1ccccc1c2ccccc2", [[0, 2, 3, 4, 5, 6], [8, 10, 11, 12, 13, 14]]),
            [['H', 'H', 'H', 'H', 'H', 'c2ccccc2'], ['c1ccccc1', 'H', 'H', 'H', 'H', 'H']]
        )
    
    def test_empty_string(self):
        with self.assertRaises(ValueError):
            find_substituents_aromatic_multiple("", [[0, 2, 3, 4, 5, 6]])
    
    def test_non_string_input(self):
        with self.assertRaises(ValueError):
            find_substituents_aromatic_multiple(12345, [[0, 2, 3, 4, 5, 6]])




class TestProcessSmilesLists(unittest.TestCase):
    def test_canonicalize_smiles(self):
        self.assertEqual(
            process_smiles_lists([['C1CCC1', 'H', 'C2CCC2']]),
            [['C1CCC1', 'H', 'C1CCC1']]
        )
    
    def test_invalid_smiles(self):
        self.assertEqual(
            process_smiles_lists([['(Invalid)', 'H', 'C2CCC2']]),
            [['Invalid', 'H', 'C1CCC1']]
        )
    
    def test_empty_string(self):
        self.assertEqual(
            process_smiles_lists([['', 'H', 'C2CCC2']]),
            [['', 'H', 'C1CCC1']]
        )
    
    def test_empty_list(self):
        self.assertEqual(
            process_smiles_lists([[]]),
            [[]]
        )
    
    def test_non_string_input(self):
        with self.assertRaises(TypeError):
            process_smiles_lists([['C1CCC1', 'H', 12345]])




class TestDeleteAtoms(unittest.TestCase):  
    def test_invalid_deletion(self):
        mol = Chem.MolFromSmiles("CCO")
        new_mol = delete_atoms(mol, [10])  # Invalid index
        self.assertIsNone(new_mol)
    
    def test_empty_atom_indices(self):
        mol = Chem.MolFromSmiles("CCO")
        new_mol = delete_atoms(mol, [])  # No atoms to delete
        self.assertIsNotNone(new_mol)
        self.assertEqual(Chem.MolToSmiles(new_mol), "CCO")
    
    def test_none_molecule(self):
        new_mol = delete_atoms(None, [1])  # None molecule
        self.assertIsNone(new_mol)




class TestSearchFctgrpsSubstituents(unittest.TestCase):
    def test_valid_target_smiles(self):
        target_smiles = "CCOCCO"
        smiles_dict = {"CCO": 1, "O": 2}
        matches = search_fctgrps_substituents(target_smiles, smiles_dict)
        expected_matches = [{"CCO": 0}, {"CCO": 3}]
        self.assertEqual(matches, expected_matches)
    
    def test_invalid_target_smiles(self):
        target_smiles = "InvalidSMILES"
        smiles_dict = {"CCO": 1, "O": 2}
        matches = search_fctgrps_substituents(target_smiles, smiles_dict)
        self.assertEqual(matches, [])
    
    def test_empty_smiles_dict(self):
        target_smiles = "CCOCCO"
        smiles_dict = {}
        matches = search_fctgrps_substituents(target_smiles, smiles_dict)
        self.assertEqual(matches, [])
    
    def test_no_matches(self):
        target_smiles = "CCOCCO"
        smiles_dict = {"N": 1, "Cl": 2}
        matches = search_fctgrps_substituents(target_smiles, smiles_dict)
        self.assertEqual(matches, [])
    
    def test_partial_matches(self):
        target_smiles = "CCOCCO"
        smiles_dict = {"CC": 1, "OCC": 2}
        matches = search_fctgrps_substituents(target_smiles, smiles_dict)
        expected_matches = [{"CC": 0}, {"CC": 3}]
        self.assertEqual(matches, expected_matches)
    
    def test_full_overlap(self):
        target_smiles = "CCOCCO"
        smiles_dict = {"CCOCCO": 1}
        matches = search_fctgrps_substituents(target_smiles, smiles_dict)
        expected_matches = [{"CCOCCO": 0}]
        self.assertEqual(matches, expected_matches)



class TestSearchMultipleAromatics(unittest.TestCase):
    def test_valid_substituents(self):
        substituents = [['c1ccccc1', 'H', 'c1ccccc1C']]
        smiles_dict = {'c1ccccc1': 1, 'C': 2}
        results = search_multiple_aromatics(substituents, smiles_dict)
        expected_results = [[[{'c1ccccc1': 0}], [{'H': 0}], [{'c1ccccc1': 0}, {'C': 6}]]]
        self.assertEqual(results, expected_results)
    
    def test_empty_substituents(self):
        substituents = []
        smiles_dict = {'c1ccccc1': 1, 'C': 2}
        results = search_multiple_aromatics(substituents, smiles_dict)
        self.assertEqual(results, [])
    
    def test_no_matches(self):
        substituents = [['c1ccccc1', 'H', 'c1ccccc1C']]
        smiles_dict = {'N': 1, 'Cl': 2}
        results = search_multiple_aromatics(substituents, smiles_dict)
        expected_results = [[[{'No match': None}], [{'H': 0}], [{'No match': None}]]]
        self.assertEqual(results, expected_results)
    
    def test_mixed_substituents(self):
        substituents = [['c1ccccc1', 'H'], ['C1=CC=CC=C1']]
        smiles_dict = {'c1ccccc1': 1, 'C1=CC=CC=C1': 2}
        results = search_multiple_aromatics(substituents, smiles_dict)
        expected_results = [[[{'c1ccccc1': 0}], [{'H': 0}]], [[{'c1ccccc1': 0}]]]
        self.assertEqual(results, expected_results)




class TestCalculationShiftAromatic(unittest.TestCase):
    def test_basic_shift(self):
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
        self.assertEqual(result, expected)
    
    def test_no_substituents(self):
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
        self.assertEqual(result, expected)
    
    def test_mixed_substituents(self):
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
        self.assertEqual(result, expected)
    
    def test_empty_ring(self):
        values_dict = {
            'H': (0, 0, 0),
            'Cl': (1, 0.5, 0.2)
        }
        matches_one_ring = []
        result = calculation_shift_aromatic(values_dict, matches_one_ring)
        expected = ['', '', '', '', '', '']
        self.assertEqual(result, expected)




class TestMultipleShiftRings(unittest.TestCase):
    def test_no_substituents(self):
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
        self.assertEqual(result, expected)

    def test_empty_matches(self):
        values_dict = {
            'H': (0, 0, 0),
            'Cl': (1, 0.5, 0.2)
        }
        matches = []
        result = multiple_shift_rings(values_dict, matches)
        expected = []
        self.assertEqual(result, expected)




class TestMatchSmilesAndMol(unittest.TestCase):
    def test_valid_smiles_and_indices(self):
        smiles = "c1ccccc1"
        aromatic_carbon_indices = [0, 1, 2, 3, 4, 5]
        result = match_smiles_and_mol(smiles, aromatic_carbon_indices)
        expected = [(0, 0), (1, 1), (2, 2), (3, 3), (4, 4), (5, 5)]
        self.assertEqual(result, expected)

    def test_invalid_smiles(self):
        smiles = "InvalidSMILES"
        aromatic_carbon_indices = [0, 1, 2, 3, 4, 5]
        with self.assertRaises(ValueError):
            match_smiles_and_mol(smiles, aromatic_carbon_indices)
    
    def test_mismatched_indices(self):
        smiles = "c1ccccc1"
        aromatic_carbon_indices = [0, 1, 2, 3]  # Not enough indices
        with self.assertRaises(ValueError):
            match_smiles_and_mol(smiles, aromatic_carbon_indices)
    
    def test_empty_smiles(self):
        smiles = ""
        aromatic_carbon_indices = []
        result = match_smiles_and_mol(smiles, aromatic_carbon_indices)
        expected = []
        self.assertEqual(result, expected)
    
    def test_no_aromatic_carbons(self):
        smiles = "CCO"
        aromatic_carbon_indices = []
        result = match_smiles_and_mol(smiles, aromatic_carbon_indices)
        expected = []
        self.assertEqual(result, expected)




class TestMainAromatic(unittest.TestCase):
    def test_toluene(self):
        expected_result = {2: 7.06, 3: 7.14, 4: 7.04, 5: 7.14, 6: 7.06}
        result = main_aromatic("Cc1ccccc1")
        self.assertEqual(result, expected_result)

    def test_benzene(self):
        expected_result = {0: 7.26, 1: 7.26, 2: 7.26, 3: 7.26, 4: 7.26, 5: 7.26}  
        result = main_aromatic("c1ccccc1")
        self.assertEqual(result, expected_result)

    def test_aniline(self):
        expected_result = {2: 6.51, 3: 7.01, 4: 6.61, 5: 7.01, 6: 6.51}  
        result = main_aromatic("Nc1ccccc1")
        self.assertEqual(result, expected_result)

    def test_invalid_name(self):
        expected_result = {}  # Expected result for an invalid molecule name
        result = main_aromatic("invalid_molecule")
        self.assertEqual(result, expected_result)

    def test_empty_name(self):
        expected_result = {}  # Expected result for an empty molecule name
        result = main_aromatic("")
        self.assertEqual(result, expected_result)




class TestShiftFunction(unittest.TestCase):

    def test_shift_aldehyde(self):
        smiles = 'CC=O'
        result = shift(1, smiles)
        expected = 9.54
        self.assertAlmostEqual(result, expected, places=2)

    def test_shift_methyl_group(self):
        smiles = 'CC=O'
        result = shift(0, smiles)
        expected = 1.81
        self.assertAlmostEqual(result, expected, places=2)

    def test_shift_hydroxyl_group(self):
        smiles = 'CCO'
        result = shift(2, smiles)
        expected = 2.46  # Example expected value, adjust as needed based on actual data
        self.assertAlmostEqual(result, expected, places=2)

    def test_shift_ammonia(self):
        smiles = 'CCN'
        result = shift(2, smiles)
        expected = 7.8  # Example expected value, adjust as needed based on actual data
        self.assertAlmostEqual(result, expected, places=2)

    def test_shift_hydrogen_sulfide(self):
        smiles = 'CCS'
        result = shift(2, smiles)
        expected = 0.72  # Example expected value, adjust as needed based on actual data
        self.assertAlmostEqual(result, expected, places=2)




class TestInfoFromMultiplicity(unittest.TestCase):

    def test_multiplicity_0(self):
        positions, intensities = info_from_multiplicity(0)
        self.assertEqual(positions, [0])
        self.assertEqual(intensities, [1])

    def test_multiplicity_1(self):
        positions, intensities = info_from_multiplicity(1)
        self.assertEqual(positions, [-0.5, 0.5])
        self.assertEqual(intensities, [1, 1])

    def test_multiplicity_2(self):
        positions, intensities = info_from_multiplicity(2)
        self.assertEqual(positions, [-1, 0, 1])
        self.assertEqual(intensities, [0.5, 1, 0.5])

    def test_multiplicity_3(self):
        positions, intensities = info_from_multiplicity(3)
        self.assertEqual(positions, list(np.arange(-1.5, 1.6)))
        self.assertEqual(intensities, [0.33, 1, 1, 0.33])

    def test_multiplicity_4(self):
        positions, intensities = info_from_multiplicity(4)
        self.assertEqual(positions, list(np.arange(-2, 2.1)))
        self.assertEqual(intensities, [0.25, 0.5, 1, 0.5, 0.25])

    def test_multiplicity_5(self):
        positions, intensities = info_from_multiplicity(5)
        self.assertEqual(positions, list(np.arange(-2.5, 2.6)))
        self.assertEqual(intensities, [0.33, 0.625, 1, 1, 0.625, 0.33])
        
    def test_multiplicity_8(self):
        positions, intensities = info_from_multiplicity(8)
        self.assertEqual(positions, list(np.arange(-4, 4.1)))
        self.assertEqual(intensities, [0.03, 0.23, 0.45, 0.70, 1, 0.7, 0.45, 0.23, 0.03])




class TestNormFunction(unittest.TestCase):

    def test_norm(self):
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
            self.assertAlmostEqual(calculated, expected, places=7)




class TestCompoundFunctions(unittest.TestCase):
    def test_get_name_from_smiles(self):
        # Test with a valid SMILES
        smiles_valid = 'CCO'
        self.assertEqual(get_name_from_smiles(smiles_valid), 'ethanol')

    def test_is_invalid_smiles(self):
        # Test with an invalid SMILES
        smiles_invalid = 'invalid_smiles'
        self.assertFalse(is_valid_smiles(smiles_invalid))
        
    def test_is_valid_smiles(self):
        # Test with a valid SMILES
        smiles_valid = 'CCO'
        self.assertTrue(is_valid_smiles(smiles_valid))

        # Test with an invalid SMILES
        smiles_invalid = 'invalid_smiles'
        self.assertFalse(is_valid_smiles(smiles_invalid))

    def test_get_smiles_from_input(self):
        # Test with a valid SMILES
        input_valid_smiles = 'CCO'
        self.assertEqual(get_smiles_from_input(input_valid_smiles), 'CCO')

        # Test with a valid chemical name
        input_valid_name = 'Ethanol'
        self.assertEqual(get_smiles_from_input(input_valid_name), 'CCO')

        # Test with an invalid input
        input_invalid = 'invalid_input'
        with self.assertRaises(ValueError):
            get_smiles_from_input(input_invalid)




class TestNMRFunction(unittest.TestCase):
    
    def test_nmr_function(self):
        # Test with a valid SMILES representation
        smiles = 'CCO'
        plt, mol = NMR(smiles)
        self.assertIsNotNone(plt)
        self.assertIsNotNone(mol)




class TestGetKeysFromValue(unittest.TestCase):

    def test_get_keys_from_value(self):
        # Test case with a dictionary having multiple keys with the same value
        d = {'a': 1, 'b': 2, 'c': 1, 'd': 3}
        value = 1
        expected_result = ['a', 'c']
        self.assertEqual(get_keys_from_value(d, value), expected_result)

        # Test case with an empty dictionary
        d = {}
        value = 5
        expected_result = []
        self.assertEqual(get_keys_from_value(d, value), expected_result)

        # Test case with a dictionary having keys but no matching value
        d = {'x': 1, 'y': 2, 'z': 3}
        value = 5
        expected_result = []
        self.assertEqual(get_keys_from_value(d, value), expected_result)


class TestShowFunction(unittest.TestCase):

    def test_show_function(self):
        # Test case with a valid molecule name and valid hydrogen index
        name = 'CCO'
        z = 1
        plt, mol = Show(name, z)
        self.assertIsNotNone(plt)
        self.assertIsNotNone(mol)

        # Test case with a valid molecule name but invalid hydrogen index
        name = 'CCO'
        z = 5  # Index 5 doesn't have a hydrogen in this molecule
        with self.assertRaises(ValueError):
            Show(name, z)

        # Test case with an invalid molecule name
        name = 'Invalid_Molecule_Name'
        z = 1
        with self.assertRaises(Exception):
            Show(name, z)