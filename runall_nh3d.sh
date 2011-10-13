#!/bin/sh
#
# runall_nh3d.sh - run ptgraph2.py on all Nh3D PDB files
#
# Usage: runall_ptg.sh [ ptgraph2 options... ]
#
# $Id$
#
#
# Run ptgraph2.py for all the PDB files in the
# Nh3D data set (Thiruv et al 2005 BMC Struct Biol 5:12)
#

myname=`basename $0`
mydir=$PWD

# location of PDB files (input)
PDBDIR=/local/charikar/Nh3D/v3.0

# output directory (.svg, .ps. output files created here)
# WARNING: contents overwritten if they exist
OUTDIR=nh3d_output

# ptgraph executable
PTGRAPH=$mydir/ptgraph2.py

ptopts=$*

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
  pdbid=`grep 'REMARK   The PDBID of the protein: ' ${pdbfile} | awk '{print $7}'`
  echo "${myname}: processing `basename ${pdbfile}` for ${pdbid}"
  cathid=`basename ${pdbfile} .pdb`
  tmppdbfile=/tmp/${cathid}.pdb
  echo "HEADER   ${pdbid} ${cathid}" > ${tmppdbfile}
  grep -v '^HEADER' ${pdbfile} >>${tmppdbfile}
  if [ ! -z "${ptopts}" ]; then
      ${PTGRAPH} ${ptopts} ${tmppdbfile}
  else
      ${PTGRAPH} ${tmppdbfile}
  fi
  echo
  rm ${tmppdbfile}
done

