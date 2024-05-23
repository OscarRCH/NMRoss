"""
Gives an approximate 1H NMR when given a SMILES or an IUPAC name of a molecule with maximum one aromatic ring and no double bonds.
"""

# Ensure compatibility with future annotations (useful for type hints)
from __future__ import annotations

# Import necessary libraries and modules
import pubchempy as pcp
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import IPythonConsole
import numpy as np
import builtins
import copy
import matplotlib.pyplot as plt
from collections import defaultdict
from rdkit import RDLogger

# Enable SVG drawing in IPython Console
IPythonConsole.ipython_useSVG = True

# Suppress RDKit warnings
lg = RDLogger.logger()
lg.setLevel(RDLogger.CRITICAL)

# Define the version of the package
__version__ = "0.0.1"
