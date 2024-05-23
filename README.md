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

```python
from nmross.NMRoss import NMR
from nmross.NMRoss import Show

NMR("molecule name or smiles")
Show("molecule name or smiles",atom index)
```

You can quickly copy these lines to get a complete view of the package's functionalities. The `NMR` function generates an image of the NMR spectrum along with an image of the molecule to help visualise it. The `Show` does the exact same, however, it adds the option to select a carbon index in the molecule, that will be highlighted in red on the molecule image, aswell as highlight the peaks of the hydrogens on that carbon in the NMR spectrum.




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



