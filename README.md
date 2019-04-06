# NAm-NOx-Lifetime
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.26295879.svg)](https://doi.org/10.5281/zenodo.2629579)

Code used in Laughner and Cohen 2019 to study changing NOx lifetime in North America

## Fair use policy
If you use part of this code in your research, please cite it with the appropriate DOI (see the above badge) in any publications resulting from its use.
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
you may need to set it up. It will require compiling the TUV model ([provided by NCAR](https://www2.acom.ucar.edu/modeling/tropospheric-ultraviolet-and-visible-tuv-radiation-model))
found in `WRF-Matlab-Parsing/tuv5.2_source/` (type `make` in that directory, required a Fortran compiler).

## Getting started

I encourage you to begin by reviewing the plots made in `lifetime_paper_figs.m`. These are the main figures in the paper
and will point you to the most useful parts of the analysis code. The EMG fit files are provided at https://doi.org/10.6078/D1RQ4V.
The .mat files should be placed in `BEHR-emissions-analysis/Workspaces/EMGFits`. (The netCDF files provided in that repository
are not used by this code, but are simply provided for non-Matlab users.)

If you are interested in generating your own line densities, the functions relevant to that are in 
`BEHR-emissions-analysis/misc_emissions_analysis.m`. These require that you have run `BEHR_initial_setup.m`
in `BEHR-core-utils` and correctly set up the paths `behr_mat_dir` and `wrf_profiles` in the `behr_paths.m`
file it creates (in `BEHR-core-utils/Utils/Constants`).

* `make_location_winds_file` - creates the initial `locs` structure with day-by-day winds from WRF
files.
* `make_rotated_line_densities` - makes NO2 line densities from the BEHR product, assuming that
you have BEHR data in .mat files. These can be created using `behrhdf2struct` in `BEHR-core-utils/Utils`
and saving the resulting `Data` struct (it must be a variable named `Data`) in a .mat file with the name
given by `behr_filename` (in `BEHR-core-utils/Utils/Constants`) for that date, profile temporal resolution
("daily") and region ("US"). 
* `make_emg_fits` - fits the line densities with EMG functions.
