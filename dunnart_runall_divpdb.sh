#!/bin/sh
#
# dunnart_runall_divpdb.sh - run dunnart on all ptgraph2 .svg output files
#                            created by ptgrunall_divpdb.sh
#
# Run Dunnart in batch mode over all SVG files created from the dividied PDB
# hierarchy by ptgrunall_divpdb.sh, using the constraints etc. in those
# files to replace them with full SVG files with constraints satsified,
# objects in place, connectors routed, etc.
# Note that are over 100 000 files in this hierarchy so don't expect this
# to finish quickly. 
# NOTE: This modifies the .svg files - it does NOT make a new copy - 
# maybe you should back the whole hierarchy up before starting this script.
#
# Usage: dunnart_runall_divpdb.sh
#
#
# Now the hacky thing. Since Dunnart requires a display (haven't
# made the -b 'batch mode' really a batch mode yet), we need to have
# the environment variable SDL_VIDEODRIVER set to 'dummy' to use the SDL
# null video device driver so it doens't try to open X11 display etc.
#
# $Id$
#
myname=`basename $0`
mydir=$PWD

# input/output directory hierarchy root (output root of ptgrunall_divpdb.sh)
# WARNING: this overwrites the .svg files in place.
##PTGROOT=/local/cluster/users/astivala/ptgoutput
PTGROOT=/local/munk/proorigami-test/output_pdb-rendered

# Dunnart executable 
DUNNART=/local/munk/proorigami-test/dunnart/trunk/dunnart

# time in seconds to limit Dunnart CPU time to before assuming failure
DUNNART_CPUTIME_LIMIT=180

# connector nudging distance
NUDGEDIST=4

# See header comments - this uses the null video device driver in SDL
export SDL_VIDEODRIVER=dummy

if [ $# -ne 0 ]; then
    echo "Usage: ${myname}" >&2
    exit 1
fi

if [ ! -x ${DUNNART} ]; then
    echo "${myname}: cannot find dunnart executable" >&2
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
  for svgfile in ${divdir}/*.svg
  do
    echo "${myname}: processing `basename ${svgfile}`"
    (ulimit -t ${DUNNART_CPUTIME_LIMIT} ; ${DUNNART} -b -y -w ${NUDGEDIST} -z 100 ${svgfile} > /dev/null)
    if [ $? -ne 0 ]; then
        echo "${myname}: error processing `basename ${svgfile}`"
    fi
  done
done
echo "${myname}: finished at " `date`
