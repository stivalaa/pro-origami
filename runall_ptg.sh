#!/bin/sh
#
# runall_ptg.sh - run ptgraph2.py on all locally downloaded PDB files
#
# Usage: runall_ptg.sh [ ptgraph2 options... ]
#
# $Id$
#
# Now that there is ptgrunall_divpdb.sh for 'production' run over entire
# local 'divided' PDB structure, this script is used for quick testing
# by using the small (about 100) set of PDB files under 
# /local/charikar/astivala/pdb
# Running on this subset takes about 1/4 hour, as opposed to a couple of 
# weeks for the whole PDB.
#

myname=`basename $0`
mydir=$PWD

# location of PDB files (input)
##PDBDIR=/local/charikar/astivala/pdb
PDBDIR=/home/astivala/pdb

# output directory (.svg, .ps. output files created here)
# WARNING: contents overwritten if they exist
OUTDIR=output

# ptgraph executable
PTGRAPH=$mydir/ptgraph2.py

#default ptgraph2 options
DEFAULT_PTGOPTS="-r35 -c -t stride -k purple -uw -l crossing:black,red,green,navy,blue -b sequential -j -e auto -f auto -o gradient -q"

if [ $# -ge 1 ]; then
    ptopts=$*
else
    ptopts=${DEFAULT_PTGOPTS}
fi


if [ ! -d ${OUTDIR} ]; then
    mkdir ${OUTDIR}
fi


version=`${PTGRAPH} -z`
echo "${myname}: ptgraph2 is ${version}"
echo "${myname}: ptgraph2 options are: ${ptopts}"
echo "${myname}: reading PDB files from ${PDBDIR}"
echo "${myname}: writing output files to ${OUTDIR}"
echo
cd ${OUTDIR}
for pdbfile in ${PDBDIR}/*.pdb
do
  echo "${myname}: processing `basename ${pdbfile}`"
  if [ ! -z "${ptopts}" ]; then
      ${PTGRAPH} ${ptopts} ${pdbfile}
  else
      ${PTGRAPH} ${pdbfile}
  fi
  echo
done
