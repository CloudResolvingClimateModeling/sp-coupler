#!/bin/bash

# This is a simple spifs run which fits on a single PC.
# OpenIFS with the T21 grid and 2 LES instances close to Barbados.
# This version uses the Amuse MPI channel.
# Demands an MPI implementation that supports comm_spawn.
# (Normally they do, Cray doesn't)

# number of GCM(OpenIFS) processes
N_GCM=1

# number of LES instances
N_LES=2

# number of LES processes per instance
LES_PROCS=1

OIFSDIR=oifs-input
DALESDIR=dales-input
OUT=output


mpiexec -n 1 python ./spmaster.py --steps 100  \
--poly        20 -50  10 -50   10 -40   20 -40	\
--gcmprocs $N_GCM --numles $N_LES --lesprocs $LES_PROCS  \
--gcmdir=$OIFSDIR --gcmexp=TEST \
--lesdir=$DALESDIR --odir=$OUT --cplsurf \
--channel=mpi 

# optional LES spinup: total spin-up time(s), number of steps, strength of the forcing towards the GCM profile
# --spinup 14400 --spinup_steps 16 --spinup_forcing 0.25





