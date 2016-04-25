# DWT-toolbox: A collection of programs for analyzing positional dependency within TF binding sites
Here you find a number of useful programs for analyzing TF binding specificity under the DWT model. We have implemented
several functionalities of the DWT model in a form of a website which can be found at [dwt.unibas.ch](https://dwt.unibas.ch).

Here we are providing a set of tools and scripts that can be used offline. For each tool, there is a directory that in
addition to the source code and examples, there is a `README` and a `manual.txt` file. The content of each directory is
independent from the other directories. The `README` and `manual.txt` files explain the content of the directory and give
guidelines on running the programs. Moreover, the content of `Examples` directory  under each directory demonstrates the
usage of the program.

More specifically, following tools are provided in this repository: 
* DWT_model: This directory contains a C++ implementation of the DWT model, which users can scan DNA sequences for identifying binding sites under the DWT model. 
* Positional_Dependency_Posterior: In this directory, a C++ implementation can be found that is useful for calculating the posterior probabilities of positional dependencies within TF binding sites. 
* diLogo: A Python script for generating graphical representation of the positional dependency model, provided a DWT flat file.
* Fitting_DWT_model: It contains a program, implemented in Python, for fitting a DWT model given a set of DNA sequences and an initial position-speicific weight matrix (PSWM) model.
This program uses other programs in this pipeline, therefore, it's recommended first to compile the other programs.

For compiling the C++ codes, we have created `Makefile` under `Source` directories. In Linux and Mac, by simply running  `make` in the `Source` directory, the executable of the program will be created. 

For implementing the C++ codes, we have used several functions from the [Boost library](http://www.boost.org/). We have put a compatible version of the Boost directory under the `lib` directory, which needs to be first un-tarred. This is done by running the following command in the command line under the `lib` directory:

```
tar -xvf boost_1_41_0.tar.gz
```

## Dependencies for the diLogo.py 
In order to run the diLogo script, it is required to have a version of Python 2.7. In addition to that, the script uses `pyx` and `argparse` libraries. For installing these libraries, we use `pip` as it's demonstrated bellow: 

```
pip install pyx==0.12.1
pip install argparse
```
Note that the compatible version of the `pyx` library for Python 2.7 is 0.12.1. Hence, in above command we specifically installed the compatible version of `pyx`.


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
