# Superparameterization of OpenIFS with DALES

This repository contains a script for running the global atmospheric model OpenIFS
coupled to local cloud-resolving LES simulations. The LES used is DALES, the Dutch Atmospheric Large Eddy Simulation.

## Authors

Fredrik Jansson (CWI, Amsterdam)
Gijs van den Oord (Netherlands e-Science center, Amsterdam)
Inti Pelupessy (Netherlands e-Science center, Amsterdam)
Pier Siebesma (TU Delft and KNMI)
Daan Crommelin (CWI, Amsterdam)


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


`oifs-input/` contains the files required to run OpenIFS for the small T21 grid. This is the standard OpenIFS test case bundled with OpenIFS itself.

`dales-input/` contains files required for DALES. This is a case with 64 x 64 x 160 grid poins. The horizontal resolution can easily be changed by editing the file `namoptions.001`.


`run_T21_mpi.sh` run example simulation using MPI. For simulations using one or more computer nodes.

`run_T21_sockets.sh` run example simulation using the AMUSE sockets channel. For simulations that fit within one node.

`run_T21_nospawn.sh` run example simulation with work-around for MPI that does not support spawn. Experimental, provided as-is.


In the Singularity image, the sockets variant works immediately. The mpi variant requires the following command to load the openMPI module:
```
eval `/usr/bin/modulecmd sh load mpi/openmpi-x86_64`
```

## Model settings

Model settings and input data are provided in three places:
* an OpenIFS input directory, containing an initial state, and model settings in fort.4
* a Dales input directory, containing model settings in namoptions.001
* options for the model coupling, provided on the command line of the coupling script. For a list of them, run `./spmaster.py --help`

For a small example, see `run_T21_sockets.sh`.

## Results

All model output is organized in an output directory:
```
dales-work-nnn/   Dales work directory, one for each LES instance.
les-input         Dales input files.
oifs-work         OpenIFS work directory, contains output from the global model.
spifs.nc          netCDF file containing vertical profiles and tendencies for the superparameterized columns.
timing.txt        CPU time statistics per time step for all models.
```


## Requirements and manual installation procedure

For initial tests, we recommend trying the Singularity image, since it simplifies the installation.
The Singularity recipe in the file Singularity can also be used as instructions for a manual setup.

For a manual setup, the following are required:

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

* AMUSE http://amusecode.org/
* OMUSE https://bitbucket.org/omuse/omuse
The OMUSE Makefiles downloads and builds the two models:
* OpenIFS  (note: requires username/password from ECMWF)
* DALES    https://github.com/CloudResolvingClimateModeling/dales 

Once the above are installed, the model coupling scripts in this repository can run.
* sp-coupler (this repository)  
