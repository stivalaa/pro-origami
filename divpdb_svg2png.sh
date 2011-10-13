#!/bin/sh
#
# divpdb_svg2png.sh - run svg2png.sh on all ptgraph2+dunnart .svg output files
#                     created by dunnart_runall_divpdb.sh
#
# Run svg2png.sh mode over all SVG files created from the dividied PDB
# hierarchy by ptgrunall_divpdb.sh plus dunnart_runall_divpdb.sh,
# Note that are over 100 000 files in this hierarchy so don't expect this
# to finish quickly. 
#
# Usage: dvipdb_svg2png.sh
#
# $Id$
#
myname=`basename $0`
mydir=$PWD

# input/output directory hierarchy root (output root of ptgrunall_divpdb.sh)
# .svg files are left as is, and each gets a .png counterpart created.
PTGROOT=/local/munk/data/proorigami_output_pdb-rendered

SVG2PNG=svg2png.sh

if [ $# -ne 0 ]; then
    echo "Usage: ${myname}" >&2
    exit 1
fi

echo "${myname}: started at " `date`

echo "${myname}: processing SVG files under ${PTGROOT}"
echo
for divdir in ${PTGROOT}/??
do
  if [ ! -d ${divdir} ]; then
      continue
  fi
  cd ${divdir}
  for svgfile in *.svg
  do
    echo "${myname}: processing `basename ${svgfile}`"
    ${SVG2PNG} ${svgfile} > /dev/null 2>&1
    if [ $? -ne 0 ]; then
        echo "${myname}: error processing `basename ${svgfile}`"
    fi
  done
done
echo "${myname}: finished at " `date`

