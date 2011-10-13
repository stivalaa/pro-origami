#!/bin/bash
#
# File:    build_omega_db_pbs_script_minlen4.sh
# Author:  Alex Stivala
# Created: November 2009
#
# PBS script for building omega+distance matrix database min SSE length 4
# (assumes distance matrix already build by build_tableaux_db_pbs_script.sh)
#
# $Id$
 	
#PBS -N build_tableaux_db

#PBS -l walltime=23:0:0

#PBS -l nodes=1

module load python/2.6.2-gcc

cd $PBS_O_WORKDIR
set CONV_RSH = ssh

INPUT_PDB_ROOT=/usr/local/ASTRAL/pdbstyle-sel-gs-bib-95-1.75
OUTPUT_TABLEAUX_DIR=/home/alexs/tableauxdb/minlen4-ASTRAL-sel-gs-bib-95-1.75
TABLEAUX_PICKLE=$OUTPUT_TABLEAUX_DIR/omegadb.pickle
DISTMATRIX_PICKLE=$OUTPUT_TABLEAUX_DIR/distmatrixdb.pickle
TABLEAUXDB_ASCII=$OUTPUT_TABLEAUX_DIR/omegadistmatrixdb.ascii
OPTIONS="-m 4 -t dssp -3 -5 -p none"

buildtableauxdb.py -n $OPTIONS $INPUT_PDB_ROOT $TABLEAUX_PICKLE      
#buildtableauxdb.py -d $OPTIONS $INPUT_PDB_ROOT $DISTMATRIX_PICKLE

#wait
if [ $? -ne 0 ]; then
	echo "build failed"
	exit 1
fi

convdbomega2.py $TABLEAUX_PICKLE $DISTMATRIX_PICKLE > $TABLEAUXDB_ASCII

times

