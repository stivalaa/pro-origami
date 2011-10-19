#!/bin/sh
#
# rerun_fix_bad_svg_pdb.sh -
#
#  rerun make_cartoon.sh on DB hierarchy cartoon output to fix bad 
#  SVG (Dunnart not run/crashed) caused by bug in old version of 
#  make_cartoons.sh that caused good domains to be overwritten by
# bad SVG when other domains need larger gapsize (in these cases,
# which ammounted to about 19% of all domains inthe database,
# the SVG file is wrong (not processed by Dunnart) but the PNG
# file still OK when getting cartoon from database).
#                    
#
# WARNING: this will overwrite SVG and PNG files with new versions
# under the OUTROOT hierarchy
#
# Usage: rerun_fix_bad_svg_pdb.sh
#
# $Id$
#

myname=`basename $0`

# location of PDB hierarchy (input)
PDBROOT=/local/munk/data/pdb/pdb

# output directory hierarchy root (.svg, .ps. output files created here)
# WARNING: overwrites .svg, .png for cartoons with unprocessed svg
OUTROOT=/local/munk/data/proorigami_output_pdb-rendered
##OUTROOT=/var/tmp/proorigami_output_pdb-rendered

if [ $# -ne 0 ]; then
    echo usage: $myname >&2
    exit 1
fi

echo "${myname}: PDB files read from under ${PDBROOT}"
echo "${myname}: writing output files under ${OUTROOT}"
echo "${myname}: started at " `date`

for divdir in ${OUTROOT}/??
do
  for svgfile in ${divdir}/*.svg
  do
    grep 'for use with Dunnart' $svgfile >/dev/null 2>&1; 
    if [ $? -eq 0 ]; then
      odir=`dirname $svgfile`
      svgpdbid=`basename $svgfile .svg`
      pdbid=`expr substr ${svgpdbid} 1 4 | tr A-Z a-z`
      div=`expr substr ${pdbid} 2 2`
      pdbfile=${PDBROOT}/${div}/pdb${pdbid}.ent.gz
      echo "* Reprocessing: ${pdbfile}"
      ( cd $odir ;  make_cartoon.sh $pdbfile )
    fi
  done  
done

echo "${myname}: finished at " `date`
times

