# Multi-slice beam propagation scattering model for 3D refractive-index reconstruction of multiple-scattering objects

Please cite the following paper when using the code:

[High-resolution 3D refractive index microscopy of multiple-scattering samples from intensity images](https://www.osapublishing.org/optica/abstract.cfm?uri=optica-6-9-1211)

Shwetadwip Chowdhury, Michael Chen, Regina Eckert, David Ren, Fan Wu, Nicole Repina, and Laura Waller

## Contents
1. [Data](#data)
2. [Usage](#usage)
3. [Resources](#Resources)
4. [Updates](#updates)

## Data
Data can be found at this link: [3D Intensity-based ODT data](https://drive.google.com/drive/folders/19eQCMjTtiK8N1f1nGtXlfXkEa8qL6kDl?usp=sharing)

README.txt contains information about the organization of the data in the C. elegans data folders.

## Usage 
0. Clone this repo: ```git clone https://github.berkeley.edu/Waller-Lab/multi-slice.git```
1. Download the UNLocBoX regularization toolbox, an extensive open-source compilation of several optimizers utilizing proximity operators. Make sure to put this toolbox into your download path. This toolbox can be found at this link: [UNLocBoX toolbox](https://epfl-lts2.github.io/unlocbox-html/).
2. Downlod the data (from above) into your download path.
3. Run main_phantom.m to run MSBP forward-model and inverse recontruction on simulated phantom
4. Run main_data.m to run MSBP inverse recontruction on C. elegan data

## Resources
Python implementation of Multi-slice as well [Multi-layer Born](https://www.osapublishing.org/optica/abstract.cfm?uri=optica-7-5-394) scattering models can be found [here](https://github.com/Waller-Lab/multi-layer-born).

## Updates
08/27/2020:
1. Added first version of code to repo. Slightly updated version compared to what is used in the paper.
