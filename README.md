![Project Logo](assets/banner.png)

![Coverage Status](assets/coverage-badge.svg)

<h1 align="center">
NMRoss
</h1>

<br>


Gives an approximate 1H NMR when given a smiles or a IUPAC name of a molecule with maximum one aromatic ring and no double bonds.

## `👷‍♂️:` Installation

## 1. Fork the Repository

First, fork the repository to your own GitHub account. This creates a copy of the repository under your GitHub account.

## 2. Clone the Repository

Clone the forked repository to your local machine: replace 'your-username' with your username.
```
git clone https://github.com/your-username/nmross.git
cd nmross
```

## 3. Create a Virtual Environment

Create a new virtual environment to keep your dependencies isolated: and activate the environment.
```
conda create -n nmross python=3.10
conda activate nmross
```
## 4. Install dependencies via conda

To avoid any build issues when installing some of the necessary packages for installation, it is easier to install them using conda.
```
conda install -c conda-forge numpy matplotlib rdkit
```




## 5. Install the Package Locally

With the virtual environment activated, install the package locally using 'pip install .': this installs the package that are in the current folder so make sure you are in the environment nmross.
```
pip install .
```
## 6. Install JupyterLab 

Install jupyter lab in the nmross environment. You can launch jupyter lab by copying the second line.

```
pip install jupyter lab
jupyter lab
```


## `🧠:` Usage


The two main functions of this module are NMR(name : str) and Show(name : str, index : int).

The function NMR takes as an input either the IUPAC name of the compound or the smiles representation of it. This function outputs a plot of the NMR of the molecule and the Mol object of the molecule. The plot appears automatically when using NMR(name), but if needed there is a Mol object that allows to see which atom on the molecule has which index.

The function Show takes as inputs either the IUPAC name of the compound or the smiles representation of it and an index in the molecule. The function outputs a plot of the 1H NMR of the molecule with the peak of the atom chosen hihlighted and a Mol object with all equivalent hydrogens highlighted. The function only ever plots a NMR spectrum if the atom with the index chosen has at least one hydrogen. The plot appears automatically when using the Show() method but the Mol object needs to be stocked in a variable to be shown.

One who only ever wants to know the 1H NMR of ethanol would use the package like this:

```python
from nmross.NMRoss import NMR
plt, mol = NMR('CCO')
mol
```

One who would like to know which peak is the one of the terminal carbon would use the package like this. First:

```python
from nmross.NMRoss import NMR
from nmross.NMRoss import Show
plt, mol = NMR('CCO')
mol
```

The second line serves to identify how the program gives the indices to the atoms in the molecule. One of these indices can then be used during the following step. In this case, the terminal carbon is the one that is meant to be identified. Thus, the next steps are:
```python
plt, mol = Show('CCO', 0)
mol
```
This shows that the primary carbon has hydrogens which are non equivalent to other hydrogens in the molecule.


Click to add a cell.




## `💥:` Development installation

Initialize Git (only for the first time). 

Note: You should have create an empty repository on `https://github.com:CedricRossboth/nmross`.

```
git init
git add * 
git add .*
git commit -m "Initial commit" 
git branch -M main
git remote add origin git@github.com:CedricRossboth/nmross.git 
git push -u origin main
```

Then add and commit changes as usual. 

To install the package, run

```
(nmross) $ pip install -e ".[test,doc]"
```

### Run tests and coverage

```
(conda_env) $ pip install tox
(conda_env) $ tox
```



