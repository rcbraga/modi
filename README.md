# MODelability Index (MODI)

MODelability Index (MODI) that estimates the feasibility of obtaining predictive QSAR models (correct classification rate above 0.7) for a binary data set of bioactive compounds. MODI is defined as an activity class-weighted ratio of the number of nearest-neighbor pairs of compounds with the same activity class versus the total number of pairs.

# System requirements
This program will run on all operating systems.

# How to install
To run MODI, you need to:
* Clone this repository
* setup Python and third-party dependencies.
* (optionally) install this package.

## Cloning this repository
See "Clone" button on this page for further information

## Setting up Python and third-party dependencies
Our package has been developed and tested using Python 3 and the following
versions of the third-party packages:
* rdkit
* numpy
* pandas
* copy
* tqdm

## Chemical descriptors and fingerprints implementations
Our package has been developed and tested using :
* maacs 
* Morgan (radius 2 and nBits 2048)
* morgan4 (radius 4 and nBits 2048)
* morgan_chiral2 (radius 2 and nBits 2048)
* morgan_chiral2 (radius 4 and nBits 2048)
* SiRMS (Simplex representation of molecular structure - a chemoinformatic tool for calculation of simplex (fragment) descriptors)

# Example runs
Note: You need to be running the following commands after you
```cd {THIS_REPOSITORY}``` if you have not installed this package.

To see the help on how to run the program go to examples folder:
``` 
Tutorial_classification.ipynb
```


# Publication
[Data Set Modelability by QSAR](https://pubs.acs.org/doi/10.1021/ci400572x)

# Acknowledgements


# Contact
For any suggestions, comments, issues, please contact us using "rodolpho-at-insilicall-dot-com"