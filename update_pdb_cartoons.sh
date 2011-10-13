#!/bin/sh
#
# update_pdb_cartoons.sh - generate or regenerate cartoons for new
#                          or updated PDB entries after the local PDB
#                          in the 'divided' directory hierarchy
#                          has been updated by update_pdb.sh 
#                          (files pdb1qlp.ent.gz under directory ql/ etc.)
#
# File:    update_pdb_cartoons.sh
# Author:  Alex Stivala
# Created: September 2010
#
# Usage: update_pdb_cartoons.sh
#
# uses the make_cartoons.sh script for each PDB entry to handle
# ptgraph2.py, Dunnart, svg2png and re-running to minimize overlaps etc.
#
# $Id$
#

export PATH=/local/munk/proorigami-prod/dunnart/trunk:/local/munk/proorigami-prod/ptgraph:${PATH}

# location of PDB hierarchy (input)
PDBROOT=/local/munk/data/pdb/pdb

# output directory hierarchy root (.svg, .ps. output files created here)
# WARNING: existing files for updated pdb entries overwritten if they exist
OUTROOT=/local/munk/data/proorigami_output_pdb-rendered

# timestamp file used to find new/modified files since last run
TIMESTAMP_FILE=${OUTROOT}/last_update_pdb_cartoons

# lockfile used to indicate this sript is currently running
LOCK_FILE=${OUTROOT}/update_pdb_cartoons_running

# working directory
TEMPDIR=/var/tmp/updatepdb$$

if [ $# -ge 1 ]; then
    echo "Usage: $0" >&2
    exit 1
fi

echo "$0: started at " `date`

if [ -f $LOCK_FILE ]; then
  echo "lock file $LOCK_FILE exists, script exiting" >&2
  exit 1
fi 

mkdir $TEMPDIR


date >> ${LOCK_FILE}

trap "rm -f ${LOCK_FILE}" 0

version=`ptgraph2.py -z`

oldir=`pwd`
cd $TEMPDIR
find ${PDBROOT} -type f -name pdb\* -newer ${TIMESTAMP_FILE}.prev | while read pdbfile
do
  echo "$0: processing new/updated file $pdbfile"
  pdbname=`basename $pdbfile`
  div=`expr substr $pdbname 5 2`
  pdbid=`expr substr $pdbname 4 4 | tr a-z A-Z`
  make_cartoon.sh $pdbfile
  mv ${pdbid}*.svg ${pdbid}*.png ${OUTROOT}/${div}
done

cd $oldir
mv ${TIMESTAMP_FILE} ${TIMESTAMP_FILE}.prev
date > ${TIMESTAMP_FILE}
rm -f ${LOCK_FILE}
rm -r ${TEMPDIR}


echo "$0: finished at " `date`
times

