
## Requirements

cmake, mpi4py, netCDF4, and the following python modules:

```
pip install --upgrade mercurial moviepy f90nml numpy scipy matplotlib nose h5py docutils netCDF4 shapely psutil
```

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


## Example case


`oifs-input/` contains the files required to run OpenIFS for the small T21 grid. This is the standard OpenIFS test case bundled with OpenIFS itself.

`dales-input/` contains files required for DALES. This is a case with 64 x 64 x 160 grid poins. The horizontal resolution can easily be changed by editing the file `namoptions.001`.


`run_T21_mpi.sh` run example simulation using MPI. For simulations using one or more computer nodes.

`run_T21_sockets.sh` run example simulation using the AMUSE sockets channel. For simulations that fit within one node.

`run_T21_nospawn.sh` run example simulation with work-around for MPI that does not support spawn. Experimental, provided as-is.

