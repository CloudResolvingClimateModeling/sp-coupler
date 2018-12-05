# Superparameterization of OpenIFS with DALES
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1968305.svg)](https://doi.org/10.5281/zenodo.1968305)


This repository contains a script for running the global atmospheric model [OpenIFS](https://confluence.ecmwf.int/display/OIFS/OpenIFS+Home)
coupled to local cloud-resolving LES simulations. The LES used is [DALES](https://github.com/dalesteam/dales),
the Dutch Atmospheric Large Eddy Simulation.


## Authors

Fredrik Jansson (CWI, Amsterdam),
Gijs van den Oord (Netherlands e-Science center, Amsterdam),
Inti Pelupessy (Netherlands e-Science center, Amsterdam),
Pier Siebesma (TU Delft and KNMI),
Daan Crommelin (CWI, Amsterdam),


## Singularity image

For easy setup of the superparameterized simulation, we provide a
Singularity recipe. This recipe can be used to build a Singularity
container including everything required to run the simulation.
See the [Singularity documentation](https://www.sylabs.io/docs/).


* Install singularity on a computer where you have root access.

* Build the image. This step requires root access.
```
sudo singularity build sp.img Singularity
```
The build procedure will ask for a user name and password for the OpenIFS git repository at ECMWF,
to download the modified OpenIFS.


* To run the simulation, launch a shell inside the container. This step does not require root access,
and can be done on a different machine.
```
singularity shell sp.img
```
By default Singularity mounts the user's home directory inside the image. If you have the sp-coupler directory somewhere in your home directory,
the singularity shell will be opened there.

Run the example simulation with
```
./run_T21_sockets.sh
```


### Limitations of the Singularity setup

It's unclear whether the Singularity image supports running on multiple nodes. AMUSE launches the workers using MPI_COMM_SPAWN,
and this may not work over multiple nodes in this setup. For large runs, we recommend a manual installation for now.


## Example case

This repository contains a small example case which can be run on a single workstation, with OpenIFS on a T21 grid coupled to two DALES models. 

`oifs-input/` contains the files required to run OpenIFS for the small T21 grid. This is the standard OpenIFS test case bundled with OpenIFS itself.

`dales-input/` contains files required for DALES. This is a case with 64 x 64 x 160 grid points. The horizontal resolution can easily be changed by editing the file `namoptions.001`.


`run_T21_mpi.sh` run example simulation using MPI. For simulations using one or more computer nodes.

`run_T21_sockets.sh` run example simulation using the AMUSE sockets channel. For simulations that fit within one node.

`run_T21_nospawn.sh` run example simulation with work-around for MPI that does not support spawn. Experimental, provided as-is.


In the Singularity image, the sockets variant works immediately. The MPI variant requires the following command to load the openMPI module:
```
eval `/usr/bin/modulecmd sh load mpi/openmpi-x86_64`
```

The full run of 100 time steps took about 13h on a quad-core workstation (i7-4790).


## Model settings

Model settings and input data are provided in three places:
* an OpenIFS input directory, containing an initial state, and model settings in fort.4
* a DALES input directory, containing model settings in namoptions.001
* options for the model coupling, provided on the command line of the coupling script. For a list of them, run `./spmaster.py --help`

For a small example, see `run_T21_sockets.sh`.


## Results

All model output is organized in an output directory:
```
dales-work-nnn/   DALES work directory, one for each LES instance.
                    surf_xy*.nc  surface fields: liquid water path, rain water path, total water path, accumulated surface rain
                    cross*.nc    cross section fields of the LES volume
les-input         copy of the DALES input files.
oifs-work         OpenIFS work directory, contains output from the global model, mainly in GRIB format.
spifs.nc          netCDF file containing vertical profiles and tendencies for the superparameterized columns.
timing.txt        CPU time statistics per time step for all models.
```

OpenIFS and DALES can be configured as usual with their respective input files, in particular the type and frequency of the output they provide.
See the model documentation for details.


## Requirements and manual installation procedure

For initial tests, we recommend trying the Singularity image, since it simplifies the installation.
The singularity recipe in the file `Singularity` can also be used as instructions for a manual setup.

For a manual setup, the following tools and libraries are required:

* C and Fortran compilers, e.g. gcc and gfortran
* make
* cmake
* netCDF4
* eccodes or gribapi
* MPI
* mpi4py
* the following Python modules:
```
pip install --upgrade mercurial moviepy f90nml numpy scipy matplotlib nose h5py docutils netCDF4 shapely psutil
```

Next, install the following programs, in this order:

* AMUSE http://amusecode.org/
* OMUSE https://bitbucket.org/omuse/omuse
  The OMUSE Makefiles downloads and builds the two models. 
    * OpenIFS  (note: requires username/password from ECMWF)
    * DALES    https://github.com/CloudResolvingClimateModeling/dales


Note that OpenIFS might require several environment variables to be set both at compilation and at runtime.
See [the OpenIFS manual](https://confluence.ecmwf.int/display/OIFS/OpenIFS+User+Guides).

Once the above are installed, you will need to add the python modules to your PYTHONPATH:
```
export PYTHONPATH=<AMUSE clone path>/src:<SP-coupler clone path>/splib
```
to run the main driver script in this repo, spmaster.py. To view all the superparametrization options and configurations (e.g. the choice of the superparametrized region), type
```
./spmaster.py --help
```

## License

The code in this repository is available under the Apache 2.0 license.

DALES and OpenIFS have their own licenses.

