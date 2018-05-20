#!/bin/bash

# This is a simple spifs run which fits on a single PC.
# OpenIFS with the T21 grid and 2 LES instances close to Barbados.
# This version uses the Amuse sockets channel - which at the moment works on
# a single node but not over networked nodes.

N_GCM=2
N_LES=2
LES_PROCS=2

OIFSDIR=oifs-input
DALESDIR=dales-input
OUT=output


python2 ./spmaster.py --steps 5  \
--poly        20 -50  10 -50   10 -40   20 -40	\
--gcmprocs $N_GCM --numles $N_LES --lesprocs $LES_PROCS  \
--channel=sockets \
--gcmdir=$OIFSDIR --gcmexp=TEST \
--lesdir=$DALESDIR --odir=$OUT --cplsurf


# --spinup 14400 --spinup_steps 16





