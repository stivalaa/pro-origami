#!/bin/sh
#
# list_overlap_diagrams.sh - list SCOP sids of diagrams with overlapping
#                            connector count not 0
#
#
# Usage: list_overlaps_diagrams.sh dirname
#      
#         dirname is root of directory hierarhcy containing all the diagrams
#         e.g. as created by dunnart_runall_divpdb.sh
#  
# $Id$

if [ $# -ne 1 ]; then
	echo "Usage: $0 dirname" >&2
	exit 1
fi
outdir=$1

find $outdir -name \*.svg -exec grep -H count \{} \; | \
              fgrep -v '"0"'|awk '{print $1}' | awk -F/ '{print $NF}' |cut -c1-7

