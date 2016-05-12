# DWT-toolbox: A collection of programs for analyzing positional dependency within TF binding sites
Here you find a number of useful programs for analyzing TF binding specificity under the dinucleotide weight tensor (DWT) model. Under the DWT model, unlike the classical position specific weight matrix (PSWM) model, we allow positional dependency between any pairs of positions and use this information for describing the binding specificity of a TF. Within this toolbox, we implemented several functionalities that are useful for analyzing TF binding under the DWT model. Moreover, a website has been implemented that offers similar functionalities. Here is the address of our website: [http://dwt.unibas.ch/fcgi/dwt](http://dwt.unibas.ch/fcgi/dwt).

In this bundle, we tried to describe the usage of these tools in detail. Therefore, users are recommended to consult with the `manual.txt` file. Moreover, the content of `Examples` directory includes several sub-directory that each of them demonstrate the usage of each tool separately.

More specifically, users can find the following functionalities: 

* DWT model: Using this feature of the toolbox, users can identify TF binding sites (TFBSs). It is a C++ implementation of the DWT model, which scans DNA sequences for identifying TFBSs under the DWT model. This code comes with a useful Python wrapper script that can be used for running the C++ code and calculating the posterior probability of the TF binding.
* Posterior of positional dependency: For calculating the probability of positional dependency between any pairs of positions, we have implemented a C++ program. Given a set of binding sites, this program can give an indea about which pairs of positions are highly dependent. 
* diLogo: A Python script for generating 'diLogo' which is a a graphical representation of the positional dependency model. DiLogo is the cousin of the well-known sequence logo, which includes information on positional dependency and the di-nucleotide tendency at the dependent positions. This script receives, as input, a DWT flat file and posteriors of positional dependency between all pairs of positions. 
* Fitting a DWT model: Sometimes we would like to fit a DWT model to a set of observed TF binding DNA sequences. For example, peaks that are extracted from a ChIP-seq assay. This Python implementation can be used for this purpose. To run the script, users need a set of DNA sequences (FASTA format) and an initial PSWM model. The latter can be extracted from known PSWM databases, such as [JASPAR](http://jaspar.genereg.net/) or [SwissRegulon](http://swissregulon.unibas.ch/fcgi/sr/swissregulon). 

For compiling the C++ codes, we created `Makefile` under `Source` directory. In Linux and Mac, it's fairly easy: just run  `make all` command while in the `Source` directory and the binaries for the programs will be created. Here, we have explicitly assumed that the `g++` is already installed in the machine. For more information see [here](https://gcc.gnu.org/). 

For implementing the C++ codes, we made use of several functions from the [Boost library](http://www.boost.org/). We
put a compatible version of the Boost directory in the `Source` directory, which needs to be first uncompressed. This is easily
done by executing the following command in `Source` directory:

```
tar -xvf boost_1_41_0.tar.gz
```

## Dependencies for the diLogo.py 
In order to run the diLogo script, it is required to have Python 2.7 already installed. In addition to that, the script uses
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
[@NimwegenLab](https://twitter.com/NimwegenLab) , 
[@SaeedOmidi](https://twitter.com/saeedomidi)

Address:
Core program Computational and Systems Biology, 
Biozentrum, University of Basel, 
Klingelbergstrasse 50-70, 
4056 Basel, Switzerland


## License 
This code is distributed under the terms of the MIT license. Please see `LICENSE.md` for more information. 
