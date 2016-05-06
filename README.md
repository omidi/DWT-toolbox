# DWT-toolbox: A collection of programs for analyzing positional dependency within TF binding sites
Here you find a number of useful programs for analyzing TF binding specificity under the DWT model. We have implemented
several functionalities of the DWT model in a form of a website which can be found at [dwt.unibas.ch](https://dwt.unibas.ch).

Here we are providing a set of tools and scripts that can be used offline. All these tools are described in detail in the
`manual.txt` file. Moreover, the content of `Examples` directory includes several sub-directory that demonstrate the usage of
each tools separately.

More specifically, following tools are provided in this repository:

* DWT_model: A C++ implementation of the DWT model, which users can scan DNA sequences for identifying binding sites under the DWT model.
* DWT_model: A Python wrapper for predicting TFBS within a set of DNA sequences. It reports the position of binding sites
 within the DNA sequences as well as their DWT score and posterior of binding.
* Positional_Dependency_Posterior: A C++ implementation can be found that is useful for calculating the posterior probabilities
of positional dependencies within TF binding sites.
* diLogo: A Python script for generating diLogo, a graphical representation of the positional dependency model, provided a DWT flat file.
* Fitting_DWT_model: A Python implementation for a pipeline that can be used for fitting the DWT models, given a set of
DWT sequences and an initial position specific weight matrix (PSWM) model. This program uses other scripts and executables, therefore,
it's first required to make sure all the other programs in the pipeline are functional.

For compiling the C++ codes, we have created `Makefile` under `Source` directory. In Linux and Mac, by simply running  `make`
in the `Source` directory, the executable of the program will be created.

For implementing the C++ codes, we have used several functions from the [Boost library](http://www.boost.org/). We have
put a compatible version of the Boost directory in the `Source` directory, which needs to be first un-tarred. This is easily
done by executing the following command in `Source` directory:

```
tar -xvf boost_1_41_0.tar.gz
```

## Dependencies for the diLogo.py 
In order to run the diLogo script, it is required to have a version of Python 2.7. In addition to that, the script uses
`pyx` and `argparse` libraries. For installing these libraries, we use `pip` as it's demonstrated bellow:

```
pip install pyx==0.12.1
pip install argparse
```

For more information see `manual.txt`.

Note that the compatible version of the `pyx` library for Python 2.7 is 0.12.1. Hence, in above command we explicitly specifically
the required version of `pyx`.


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
