#!/bin/sh
#
# renderall.sh - Run Dunnart on local test PDB files output (created by
#                runall_ptg.sh) and render to PNG
#
# Run Dunnart in batch mode over all SVG files created from the local PDB
# test files by runall_ptg.sh, using the constraints etc. in those
# files to replace them with full SVG files with constraints satsified,
# objects in place, connectors routed, etc.
# Then convert to PNG with svg2png.sh
#
# Usage: renderall.sh [outdir]
#
# Ruin this script from where the output directory is under.
# (ie same as runall_ptg.sh) If outdir is not specified, ./output
# is used.
#
# Becuase of the way Dunnart currently works, Dunnart
#  must be run from the Dunnart
# directory. (E.g. /local/charikar/proorigami/dunnart).
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

# input/output directory (output of ptgrunall.sh)
# WARNING: this overwrites the .svg files in place.
PTGROOT=${mydir}/output

# Dunnart build directory
#DUNNART_DIR=${HOME}/dunnart/trunk
DUNNART_DIR=${HOME}//home/astivala/dunnart.local.r2794-headtaillabels

# Dunnart executable (must be in DUNNART_DIR)
DUNNART=dunnart

# time in seconds to limit Dunnart CPU time to before assuming failure
DUNNART_CPUTIME_LIMIT=150


# See header comments - this uses the null video device driver in SDL
export SDL_VIDEODRIVER=dummy

# SVG to PNG convertor
SVG2PNG=svg2png.sh

if [ $# -gt 1 ]; then
    echo "Usage: ${myname} [outdir]" >&2
    exit 1
fi

if [ $# -eq 1 ]; then
    outdir=${mydir}/$1
else
    outdir=${PTGROOT}
fi

if [ ! -x ${DUNNART_DIR}/${DUNNART} ]; then
    echo "${myname}: cannot find dunnart executable" >&2
    exit 1
fi

cd ${DUNNART_DIR}

echo "${myname}: started at " `date`

echo "${myname}: processing SVG files under ${outdir}"
echo
for svgfile in ${outdir}/*.svg
do
  echo "${myname}: processing `basename ${svgfile}`"
  #(ulimit -t ${DUNNART_CPUTIME_LIMIT} ; ${DUNNART} -b -y -w5 -z100 ${svgfile} > /dev/null)
  (ulimit -t ${DUNNART_CPUTIME_LIMIT} ; ${DUNNART} -b -y  ${svgfile} > /dev/null)
  if [ $? -ne 0 ]; then
      echo "${myname}: error processing `basename ${svgfile}`"
  fi
done
echo "${myname}: finished Dunnart at " `date`

echo "${myname}: converting SVG to PNG"
cd ${outdir}
for svgfile in *.svg
do
  echo "${myname}: converting `basename ${svgfile}`"
  ${SVG2PNG} ${svgfile} > /dev/null 2>&1
  if [ $? -ne 0 ]; then
      echo "${myname}: error converting `basename ${svgfile}`"
  fi
done
echo "${myname}: finished at " `date`

