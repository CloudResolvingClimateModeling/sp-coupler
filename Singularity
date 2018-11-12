Bootstrap: docker
From: centos:7

%setup
mkdir ${SINGULARITY_ROOTFS}/opt/splib

%files
splib/*.py      /opt/splib/
spmaster.py     /opt/
oifs-input/     /opt/
dales-input/    /opt/



%post
yum -y update
yum install -y epel-release
yum groupinstall -y "Development Tools"
yum install -y git mercurial gcc-gfortran cmake python-devel python-pip python-wheel wget openmpi-devel mpi4py-openmpi netcdf-devel netcdf-fortran-devel fftw-devel gmp-devel mpfr-devel gsl-devel atls atlas-devel blas-devel lapack-devel perl-Digest-MD5 perl-Time-Piece perl-IO-Compress
pip install --upgrade --ignore-installed pip setuptools
pip install moviepy f90nml numpy scipy matplotlib nose h5py docutils netCDF4 shapely psutil cython

cd /opt

wget https://software.ecmwf.int/wiki/download/attachments/3473437/grib_api-1.21.0-Source.tar.gz?api=v2 -O grib_api-1.21.0-Source.tar.gz
tar -xzf grib_api-1.21.0-Source.tar.gz
cd grib_api-1.21.0-Source
mkdir build
cd build
cmake ..
make
make install
export GRIB_API_DIR=$PWD
cd /opt
git clone  -b release-11.1 --single-branch --depth=1 https://github.com/amusecode/amuse.git 
cd amuse
export PYTHON=python
export MODULEPATH=/etc/modulefiles
eval `/usr/bin/modulecmd sh load mpi/openmpi-x86_64`
./configure FC=gfortran FCFLAGS="-I/usr/include -I/usr/lib64/gfortran/modules"
make framework
cd src

hg clone -r v1.1 https://bitbucket.org/omuse/omuse
cd omuse


cd community/oifs
export OIFS_GRIB_API_DIR=$GRIB_API_DIR
export DOWNLOAD_CODES=0
git clone https://git.ecmwf.int/scm/~g.vandenoord_esciencecenter.nl/oifs40r1-lib.git oifs-repo
make

cd ../dales
export DOWNLOAD_CODES=1
make


%environment
PYTHONPATH=/opt/amuse/src/:/opt/splib/
LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/local/lib/
#  LD_LIBRARY_PATH=/.singularity.d/libs:/usr/local/lib/
 GRIB_SAMPLES_PATH=/usr/local/share/grib_api/ifs_samples/grib1_mlgrib2
MODULEPATH=/etc/modulefiles

export PYTHONPATH GRIB_SAMPLES_PATH LD_LIBRARY_PATH MODULEPATH



