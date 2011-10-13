#!/bin/bash
#
# File:    sort_tableaux_db_pbs_script.sh
# Author:  Alex Stivala
# Created: February 2010
#
# PBS script for sorting tableaux+distance matrix database
#
 	
#PBS -N sort_tableaux_db

#PBS -l walltime=23:0:0

#PBS -l nodes=1

module load python/2.6.2-gcc

cd $PBS_O_WORKDIR
set CONV_RSH = ssh


OUTPUT_TABLEAUX_DIR=/home/alexs/tableauxdb/ASTRAL-sel-gs-bib-95-1.75
TABLEAUX_PICKLE=$OUTPUT_TABLEAUX_DIR/tableauxdb.pickle
DISTMATRIX_PICKLE=$OUTPUT_TABLEAUX_DIR/distmatrixdb.pickle
TABLEAUXDB_ASCII=$OUTPUT_TABLEAUX_DIR/tableauxdistmatrixdb.sorted.ascii

convdb2.py -s $TABLEAUX_PICKLE $DISTMATRIX_PICKLE > $TABLEAUXDB_ASCII

times

