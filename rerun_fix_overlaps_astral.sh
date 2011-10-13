#!/bin/sh
#
# rerun_fix_overlaps_astral.sh - use make_cartoons.sh to rerun ptgraph2
#                    and Dunnart on the output directory hierarchy from
#                    ptgrunall_divpdb.sh followed by dunnart_runall_divpdb.sh
#                    for those SVG files that have a nonzero number of
#                    overlaps as reported by Dunnart in the SVG.
#                    
#
# Usage: rerun_fix_overlaps_astral.sh
#
# $Id$
#

myname=`basename $0`

# location of ASTRAL SCOP hierarchy (input)
ASTRALROOT=/local/munk/data/ASTRAL/pdbstyle-sel-gs-bib-95-1.75

# output directory hierarchy root (.svg, .ps. output files created here)
# WARNING: overwrites .svg, .png for cartoons that have overlaps
#OUTROOT=/local/munk/proorigami-test/output-rendered
OUTROOT=/local/munk/var/tmp/output-rendered

if [ $# -ne 0 ]; then
    echo usage: $myname >&2
    exit 1
fi

echo "${myname}: PDB files read from under ${ASTRALROOT}"
echo "${myname}: writing output files under ${OUTROOT}"
echo "${myname}: started at " `date`

for svgfile in `find ${OUTROOT} -name \*.svg -exec grep -H count \{} \; | fgrep -v '"0"'|awk '{print $1}' | sed 's/:$//'`
do
    odir=`dirname $svgfile`
    scopsid=`basename $svgfile .svg`
    div=`echo $scopsid | cut -c3-4`
    pdbfile=${ASTRALROOT}/${div}/${scopsid}.ent
    ( cd $odir ;  make_cartoon.sh $pdbfile )
    
done

echo "${myname}: finished at " `date`

