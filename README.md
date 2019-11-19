# Superparameterization of OpenIFS with DALES
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1968305.svg)](https://doi.org/10.5281/zenodo.1968305)


This repository contains a script for running the global atmospheric model [OpenIFS](https://confluence.ecmwf.int/display/OIFS/OpenIFS+Home)
coupled to local cloud-resolving LES simulations. The LES used is [DALES](https://github.com/dalesteam/dales),
the Dutch Atmospheric Large Eddy Simulation.

A description of the coupling procedure and simulation results are given in 
[Jansson, F., van den Oord, G., Pelupessy, I., Grönqvist, J. H., Siebesma, A. P., & Crommelin, D. (2019). Regional superparameterization in a global circulation model using large eddy simulations. Journal of Advances in Modeling Earth Systems, 11](https://doi.org/10.1029/2018MS001600)

Interfaces to the models are built with [OMUSE](https://bitbucket.org/omuse/omuse/src/default/)
The interfaces are documented in the [OMUSE documentation](https://omuse.readthedocs.io/en/latest/)

## Authors

Fredrik Jansson (CWI, Amsterdam),
Gijs van den Oord (Netherlands e-Science center, Amsterdam),
Inti Pelupessy (Netherlands e-Science center, Amsterdam),
Maria Chertova (Netherlands e-Science center, Amsterdam),
Pier Siebesma (TU Delft and KNMI),
Daan Crommelin (CWI, Amsterdam),

## License

The code in this repository is available under the Apache 2.0 license.

DALES and OpenIFS have their own licenses.


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


### Format of spifs.nc

The output file spifs.nc contains vertical profiles of model variables and superparameterization tendencies
for every superparameterized grid point and global model time step.
The data is organized in groups according to the grid point where the model is located,
for example all data for the DALES at grid point 888 is located in the group 888/ in the netCDF file.
In general, variables in upper case relate to the global model, and variables in lower case relate to the local model.
Forcings *on* the global model are denoted e.g. f_T, and on the local model f_thl.

#### Vertical coordinates

Profiles in the local model use `zf`, in the root group of the file, as vertical coordinate. These are constant in time and the same for all the local models.
For the global model, the vertical coordinate is `Zf`, which depends on both the grid point and time (because the global model's
levels are not on a fixed height but defined by pressure, they vary in time and space).

#### Variables

The most important variables are summarized below.

================  ======= ===========================================================================
OpenIFS Variable  Unit    Description
================  ======= ===========================================================================
lat, lon          degrees grid point coordinates
U,V               m/s     velocity components in x, y directions
T                 K       temperature
SH                kg/kg   specific humidity (i.e. water vapor, not cloud condensate)
QL                kg/kg   specific cloud condensate, liquid
QI                kg/kg   specific cloud condensate in the form of ice
QT                kg/kg   total specific humidity, SH+QL+QI
Pf                Pa      pressure
A                 -       cloud fraction
f_U, f_V          m/s^2   forcings on global model
f_T               K/s
f_SH, f_QL, f_QI  kg/kg/s
================  ======= ===========================================================================

================  ======= ===========================================================================
DALES Variable    Unit    Description
================  ======= ===========================================================================
u, v              m/s     velocity components in x, y directions
thl               K       liquid water potential temperature
qt                kg/kg   total specific humidity
ql                kg/kg   condensed water specific humidity
wthl              K m/s   surface heat flux
wqt               m/s     surface moisture flux
f_u, f_v          m/s^2   forcings on local model
f_thl             K/s
f_qt              kg/kg/s
================  ======= ===========================================================================



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

# Installation notes for specific systems

## Installation on Arch Linux

When configuring OMUSE, one must explicitly specify python2, since the default is python3.

```
cd amuse
PYTHON=python2 ./configure --with-netcdf=/usr/
make framework

export DOWNLOAD_CODES=all

cd src/omuse/community/dales
make

cd ../oifs
make
```

## Installation on Fedora

Fedora's netcdf require some extra settings, becuse the module files
and .inc files are in different places. We specify the module path
with FCFLAGS: Another issue seen on Fedora is that make in the dales
directory fails with `build.py: error: No module named
dalesreader`. One solution is to add . to PYTHONPATH. This seems to confuse mercurial though.


```
FCFLAGS=-I/usr/lib64/gfortran/modules ./configure --with-netcdf=/usr
make framework

export DOWNLOAD_CODES=all

export PYTHONPATH=$PYTHONPATH:.  # for dalesreader to be found when creating the interface code
cd src/omuse/community/dales
make

cd ../oifs
make
```


## Installation on ECMWF Cray system


### Initial setup

#### Load modules
```
prgenvswitchto intel

module load python/2.7.12-01
module load netcdf4/4.4.1
module load cmake/3.12.0
module load git
module load eccodes
```

#### Other settings
```
# https proxy
export https_proxy=proxy:2222

export AMUSE_DIR=$PERM/2019/amuse/
export PYTHONPATH=$PYTHONPATH:$AMUSE_DIR/src/

source $PERM/meteo/bin/activate

# OpenIFS compilation options
export OIFS_COMP=intel
export OIFS_BUILD=craynomp

# Cray setup: all compilers are invoked with these names:
export OIFS_FC=ftn
export OIFS_CC=cc

export OIFS_GRIB_API_DIR=$ECCODES_DIR
export OIFS_GRIB_API_LIB="-L $ECCODES_LIB_DIR -leccodes_f90 -leccodes"
export OIFS_GRIB_API_INCLUDE="-I $ECCODES_INCLUDE_DIR"

export FCFLAGS="-convert big_endian"

# On the Cray, we don't want any linking flags for Lapack
# they are included when using the Cray compiler wrappers
export OIFS_LAPACK_LIB=" "

# DALES compilation options
export SYST=ECMWF-intel
export DALES_FCFLAGS="-g -traceback -O3 -r8 -xHost -fpp"
#these flags apply to the interface only

```


#### virtual Python environment

```
pip install --user virtualenv
PATH=$PATH:~/.local/bin/

cd $PERM
virtualenv meteo
source $PERM/meteo/bin/activate
pip install --upgrade mercurial moviepy f90nml numpy scipy matplotlib nose h5py docutils netCDF4 shapely psutil
```                     

#### mpi4py on ECMWF

Since mid-2018 the mpi4py installed with the python modules at ECMWF no longer works. It can be installed manually from source.
This should be done with the same set of compilers and modules loaded as used for everything else.

* activate the virtual python environment, and with the intel compiler and our modules loaded.
```
cd $PERM
wget https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-3.0.0.tar.gz -O mpi4py-3.0.0.tar.gz
tar zxf mpi4py-3.0.0.tar.gz
cd mpi4py-3.0.0

# add an enry for the Cray system in mpi.cfg
cat >> mpi.cfg <<EOF
[cray]
mpicc = cc
mpicxx = CC
extra_link_args = -shared
EOF

python setup.py build --mpi=cray
python setup.py install 

cd $PERM
```

##### Notes
FJ tried to compile mpi4py with the gnu compiler (`prgenvswitchto gnu`). Compilation seemed OK, but python segfaulted when testing the coupled system. Compiling mpi4py with the intel compiler seems to work - no module changes needed.

[Source for mpi4py instructions](http://jaist-hpc.blogspot.com/2015/02/mpi4py.html)

The following instructions install omuse and amuse sibe by side in the directory $PERM/2019/.
Then a symlink in amuse/src is created, to omuse/src/omuse, so that the path amuse/src/omuse/community still works.

### OMUSE

```
cd $PERM/2019
hg clone --insecure https://bitbucket.org/omuse/omuse
```


### Amuse

```
git clone https://github.com/fjansson/amuse
cd amuse
git checkout spawnless

cd src
ln -s $PERM/2019/omuse/src/omuse omuse
# so that the old path amuse/src/omuse/community still works

cd ..
```

This version is our own no-spawn fork for use at ECMWF. Elsewhere, the official amuse can be used:
<https://github.com/amusecode/amuse/> .

  
```
#make AMUSE find the right python:
export PYTHON=python

./configure FC=ftn CC=cc --with-netcdf=`nc-config --prefix`
# some libraries will not be found, e.g. gsl. This is OK 

make framework
```


### OpenIFS and DALES
OpenIFS and DALES can be cloned using the OMUSE make file.

```
export DOWNLOAD_CODES=all
# DOWNLOAD_CODES=all will checkout entire repo with ssh, intended for developers of the components.
# DOWNLOAD_CODES=latest will (shallow) checkout latest revision only
# DOWNLOAD_CODES=<anything else> will (shallow) checkout release tag spifs_v1.0.0

export AMUSE_DIR=$PERM/2019/amuse/
export PYTHONPATH=$PYTHONPATH:$AMUSE_DIR/src/
export PATH=$PATH:$AMUSE_DIR/bin/
```

```
cd community/dales
make
cd ../..

cd community/oifs
make
# note: this downloads OpenIFS, which requires ECMWF credentials

```

