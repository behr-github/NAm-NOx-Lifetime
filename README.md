# NAm-NOx-Lifetime
Code used in Laughner and Cohen 2019 to study changing NOx lifetime in North America

## Fair use policy
If you use part of this code in your research, please cite it with the appropriate DOI in any publications resulting from its use.
If this code is critical to the conclusions of a paper using it, please reach out to Ron Cohen (rccohen@berkeley.edu) or 
Josh Laughner (jlaughner@berkeley.edu) and inquire if offering coauthorship would be appropriate.

## Fair warning
This repository represents ~2 years worth of analysis. Not all of the methods contained within have been updated to use
the final version of the data. If you have questions, please reach out to one of the above contacts.

## Installation
Some parts of the code require Python functionality that in turn relies on Cython components that must be compiled and
installed. This can be handled by running `BEHR_initial_setup.m` in `BEHR-core-utils`. If it fails, it will provide instructions
on how to try it manually. We find that configuring Matlab to use a Python 3 installation works much more reliably than Python 2.
See [Matlab's pyversion function](https://www.mathworks.com/help/matlab/ref/pyversion.html?s_tid=srchtitle) to switch your Python 
interpreter.

Included in this repository is a Matlab package for parsing WRF KPP files into a chemical mechanism for analysis. This part is 
not pertinent to the analysis in Laughner and Cohen 2019, but if you try to use some of the methods in `misc_wrf_lifetime_analysis`,
you may need to set it up. It will require compiling the TUV model ([provided by NCAR](https://www2.acom.ucar.edu/modeling/tropospheric-ultraviolet-and-visible-tuv-radiation-model)
found in `WRF-Matlab-Parsing/tuv5.2_source/` (type `make` in that directory, required a Fortran compiler).
