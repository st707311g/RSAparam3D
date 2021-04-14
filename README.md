# RSAparam3D: A module of RSAtrace3D, an RSA vectorization software, for 3D RSA phenotyping

![top image](figures/top_image.jpg) 

XXX

## Installation

### A Backbone Software RSAtrace3D

Clone RSAtrace3D repository if it is not already cloned:

    $ git clone https://github.com/st707311g/RSAtrace3D.git

For detailed installation instructions, please refer to [here](https://github.com/st707311g/RSAtrace3D) and make sure you have installed it correctly.

### RSAparam3D

Move into the `mod` directory in the RSAtrace3D root directory.

    $ cd RSAtrace3D/mod

Clone RSAparam3D repository.

    $ git clone https://github.com/st707311g/RSAparam3D.git

Move back to the the RSAtrace3D root directory.

    $ cd .. 

Make sure you have installed RSAparam3D correctly.

    $ python list_modules.py | grep RSAparam3D
    RootTraitBackbone, <class 'mod.RSAparam3D.RSAparam3D.Root_RSAparam3D'>
    RSATraitBackbone, <class 'mod.RSAparam3D.RSAparam3D.RSA_RSAparam3D'>
    ExtensionBackbone, <class 'mod.RSAparam3D.RSAparam3D_simulatior.ExtensionRSAparam3D'>

You will find three modules, namely `Root_RSAparam3D`, `RSA_RSAparam3D`, and `ExtensionRSAparam3D`.


## Citation

Please cite the following article:

XXX

## Copyright

National Agriculture and Food Research Organization (2021)

## Project homepage
https://rootomics.dna.affrc.go.jp/en/
