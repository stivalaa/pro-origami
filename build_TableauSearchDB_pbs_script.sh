#!/bin/bash
#
# File:    build_TableauSearchDB_pbs_script.sh
# Author:  Alex Stivala
# Created: November 2009
#
# PBS script for building database of tableaux (numeric) and distance
# matrices for (modified) version of Arun's TableauSearch
# (TableauComparer). Uses 2 CPUs - simultaneously build tableaux and
# distance matrix
#
# requires PATH and PYTHONPATH already set up in environment 
#
# $Id$
 	
#PBS -N build_tableausearch_db
#PBS -l walltime=8:0:0
#PBS -l nodes=2

# input ASTRAL pdbstyle hierarchy
ASTRAL=/usr/local/ASTRAL/pdbstyle-sel-gs-bib-95-1.75

# output directory
TABSEARCHDB=${HOME}/TableauSearchDB

TABOPTS="-e -3 -5 -t dssp -p none"


cd $PBS_O_WORKDIR
set CONV_RSH = ssh

for i in `find $ASTRAL -name \*.ent`
do
  pytableaucreate.py $TABOPTS $i > $TABSEARCHDB/`basename $i`.angles &
  pytableaucreate.py -d $TABOPTS  $i > $TABSEARCHDB/`basename $i`.distmatrix
  wait
done

find $TABSEARCHDB -name \*.angles > ${TABSEARCHDB}/dbfnames.txt

times

