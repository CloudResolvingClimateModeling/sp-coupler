#!/bin/bash

# This is a simple spifs run which fits on a single PC.
# OpenIFS with the T21 grid and 2 LES instances close to Barbados.
# This version uses the Amuse sockets channel - which at the moment works on
# a single node but not over networked nodes.

# number of GCM(OpenIFS) processes
N_GCM=1

# number of LES instances
N_LES=2

# number of LES processes per instance
LES_PROCS=1


OIFSDIR=oifs-input
DALESDIR=dales-input
OUT=output


python2 ./spmaster.py --steps 100  \
--poly        20 -50  10 -50   10 -40   20 -40	\
--gcmprocs $N_GCM --numles $N_LES --lesprocs $LES_PROCS  \
--gcmdir=$OIFSDIR --gcmexp=TEST \
--lesdir=$DALESDIR --odir=$OUT --cplsurf \
--channel=sockets 

# optional LES spinup: total spin-up time(s), number of steps, strength of the forcing towards the GCM profile
# --spinup 14400 --spinup_steps 16 --spinup_forcing 0.25







