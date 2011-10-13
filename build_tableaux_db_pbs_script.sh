#!/bin/bash
#
# File:    build_tableaux_db_pbs_script.sh
# Author:  Alex Stivala
# Created: July 2009
#
# PBS script for building tableaux+distance matrix database
#
 	
#PBS -N build_tableaux_db

#PBS -l walltime=23:0:0

#PBS -l nodes=2

module load python/2.6.2-gcc

cd $PBS_O_WORKDIR
set CONV_RSH = ssh

INPUT_PDB_ROOT=/usr/local/ASTRAL/pdbstyle-sel-gs-bib-95-1.75
OUTPUT_TABLEAUX_DIR=/home/alexs/tableauxdb/ASTRAL-sel-gs-bib-95-1.75
TABLEAUX_PICKLE=$OUTPUT_TABLEAUX_DIR/tableauxdb.pickle
DISTMATRIX_PICKLE=$OUTPUT_TABLEAUX_DIR/distmatrixdb.pickle
TABLEAUXDB_ASCII=$OUTPUT_TABLEAUX_DIR/tableauxdistmatrixdb.ascii
OPTIONS="-t dssp -3 -5 -p none"

buildtableauxdb.py $OPTIONS $INPUT_PDB_ROOT $TABLEAUX_PICKLE      &
buildtableauxdb.py -d $OPTIONS $INPUT_PDB_ROOT $DISTMATRIX_PICKLE

wait
if [ $? -ne 0 ]; then
	echo "build failed"
	exit 1
fi

convdb2.py $TABLEAUX_PICKLE $DISTMATRIX_PICKLE > $TABLEAUXDB_ASCII

times

