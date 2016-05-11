# DWT-toolbox source codes
This directory contains all the source codes, both in C++ and Python, that are used in the DWT-toolbox. For compiling the
C++ codes, consult with the `../manual.txt`.

## Compiling the C++ codes
A `Makefile` was made that is helps for compiling the C++ codes and producing the executables. For compiling, simply
execute the following command in the terminal:

```
make all
```

Note that first the Boost library will be uncompressed, and consequently the source codes will be compiled. After running
this command, the following executables will be found in the `Source` directory.

* DWT_TFBS_prediction: For scanning and calculating TF binding score for a set of DNA sequences under the DWT model.
* positional_dependency_posterior: For calculating the posterior of positional dependency between all pairs of positions within
TF binding sites (TFBSs).
* motevo: A program for TFBS identification under the classical position specific weight matrix (PSWM) model. For more information
see the original [paper](http://www.ncbi.nlm.nih.gov/pubmed/22334039).


## Python scripts
We have implemented several functionalities in Python (ver 2.70). This includes:

* diLogo.py: For generating diLogo, a graphical representation of binding specificity under the positional dependency
model.
* DWT_TFBS_prediction.py: A wrapper script for running the C++ code for TFBS identification. The output of the program will
be a set of TFBSs with information such as their position within the sequence and the posterior probability of TF binding.
* fitting_DWT_model.py: It is a useful program for fitting a DWT model given a set of DNA sequences and an initial PSWM model.
* generate_DWT_model.py: For generating a DWT flat file, given a set of k-mers. This program basically calculate the di-nucleotide
  frequencies between all pairs of positions


For running this scripts, it is required to have  `pip` amd  `argparse` installed. To that end, user can simply execute
the following commands:

```
pip install pyx==0.12.1
pip install argparse
```

For more information see `manual.txt`.


## Contacts and credits
The DWT model is developed in C++ by Saeed Omidi and Erik van Nimwegen.
For citing our work, please check out our website:
http://nimwegenlab.org/

Email:
saeed.omidi@gmail.com ,
erik.vannimwegen@unibas.ch

On Twitter:
@NimwegenLab ,
@SaeedOmidi

Address:
Core program Computational and Systems Biology,
Biozentrum, University of Basel,
Klingelbergstrasse 50-70,
4056 Basel, Switzerland


## License
This code is distributed under the terms of the MIT license. Please see `LICENSE.md` for more information.

