#!/bin/sh
#
# editsvglayer.sh - Edit the layer attribute in Dunnart SVG in order
#                   to hide connector overlap highlight
#
# Usage: editsvglayer.sh svgfile
#
#
# Alex Stivala
# April 2010
#
# $Id$
#


if [ $# -ne 1 ]; then
    echo "Usage: $0 svgfile" >&2
    exit 1
fi

svgfile=$1

infile=${svgfile}.tmp.$$
outfile=${svgfile}
mv ${svgfile} ${infile}
sed 's/\(.*id="layer1".*\)display:inline\(.*\)/\1display:none\2/' < ${infile} > ${outfile}
rm ${infile}

