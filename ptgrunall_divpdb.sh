#!/bin/sh
#
# runall_divpdb.sh - run ptgraph2.py on all locally downloaded PDB files
#                    in the 'divided' directory hierarchy
#                    (files pdb1qlp.ent.gz under directory ql/ etc.)
#
# Note that are over 47000 files in this hierarchy so don't expect this
# to finish quickly, and make sure there is enough output space available.
#
# Usage: runall_divpdb.sh [ -o ] [ ptgraph2 options... ]
#
# If -o is specified, overwrite output files if they exist.
# The default behaviour is to leave existing files and only write output
# files that don't exist; useful for restarting the job if it is stopped
# partway through.
#
# $Id$
#

myname=`basename $0`
mydir=$PWD

# location of PDB hierarchy (input)
##PDBROOT=/local/munk/data/ASTRAL/pdbstyle-sel-gs-bib-95-1.75
PDBROOT=/local/munk/data/pdb/pdb/

# output directory hierarchy root (.svg, .ps. output files created here)
# WARNING: contents overwritten if they exist
OUTROOT=/local/munk/proorigami-test/output_pdb

# ptgraph executable
PTGRAPH=$mydir/ptgraph2.py

#default ptgraph2 options (must match tableau+distmatrix db so autoselect
# query SSE works correctly on QP and SA tableau searching on the webserver)
#DEFAULT_PTGOPTS="-r35 -t dssp -k purple -uw -l crossing:black,red,green,navy,blue -b sequential -j -e auto -f auto -o gradient"
DEFAULT_PTGOPTS="-r35 -t dssp -k purple -l crossing:black,red,green,navy,blue -b sequential -j -e auto -f auto -o gradient -p cath:/local/munk/proorigami-test/cath/CathDomall.v3.3.0"

overwrite=0
if [ $# -gt 0 ]; then
  if [ $1 = "-o" ]; then
      overwrite=1
      shift 1
  fi
fi


if [ $# -ge 1 ]; then
    ptopts=$*
else
    ptopts=${DEFAULT_PTGOPTS}
fi

echo "${myname}: started at " `date`

version=`${PTGRAPH} -z`

if [ ! -d ${OUTROOT} ]; then
    mkdir ${OUTROOT}
fi

echo "${myname}: reading PDB files from under ${PDBROOT}"
echo "${myname}: writing output files under ${OUTROOT}"
if [ $overwrite -eq 1 ]; then
    echo "${myname}: overwriting existing output files"
else 
    echo "${myname}: leaving existing output files"
fi
echo "${myname}: ptgraph2 version is: ${version}"
echo "${myname}: ptgraph2 options are: ${ptopts}"
echo
cd ${OUTROOT}
for divdir in ${PDBROOT}/??
do
  if [ ! -d ${divdir} ]; then
      continue
  fi
  ddname=`basename ${divdir}`
  if [ ! -d ${OUTROOT}/$ddname ]; then
      mkdir ${OUTROOT}/$ddname
  fi
  cd ${OUTROOT}/$ddname
  # ptgraph2 handles both pdb????.ent and pdb????.ent.gz and d??????.ent files
  for pdbfile in ${divdir}/*.ent*
  do
    if [ $overwrite -eq 1 ]; then
        doit=1
    else
        bname=`basename $pdbfile`
        if [ `expr substr $bname 1 3` == "pdb" ]; then
          pdbid=`expr substr $bname 4 4 | tr '[a-z]' '[A-Z]'`
        else
          pdbid=`expr substr $bname 1 7`
        fi
        if [ -f ${pdbid}.svg -o -f ${pdbid}-1.svg ]; then
            doit=0
        else
            doit=1
        fi
    fi

    if [ $doit -eq 1 ]; then
        echo "${myname}: processing `basename ${pdbfile}`"
        if [ ! -z "${ptopts}" ]; then
            ${PTGRAPH} ${ptopts} ${pdbfile}
        else
            ${PTGRAPH} ${pdbfile}
        fi
    else
        echo "${myname} output for pdbid $pdbid exists, skipping"
    fi
    echo
  done
done
echo "${myname}: finished at " `date`
