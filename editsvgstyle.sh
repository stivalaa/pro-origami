#!/bin/sh
#
# editsvgstyle.sh - Edit the CSS style attribute in Dunnart SVG in order
#                   to change the connector stroke width.
#
# Usage: editsvgstyle.sh linewidth svgfile
#
#         linewdith is the connector stroke width in pixels
#         svgfile is the SVG file to edit NB this file is changed
#
# Alex Stivala
# April 2008
#
# $Id$
#

fontsize=11

if [ $# -ne 2 ]; then
    echo "Usage: $0 linewdith svgfile" >&2
    exit 1
fi

linewidth=$1
svgfile=$2

# linewidth2 is larger linewdith for highlighting when cursor over it
linewidth2=`expr $linewidth \* 2 + 1`

styleattr=`printf \
"  <style type=\"text/css\"> "\
".shape   { "\
"stroke:black; "\
"stroke-width:1px; "\
"fill:#f0f0d2; "\
"} "\
".connector { "\
"fill:none; "\
"stroke:black; "\
"stroke-width:%dpx; "\
"stroke-linecap:butt; "\
"stroke-linejoin:miter; "\
"stroke-opacity:1; "\
"} "\
".connector:hover { "\
"stroke-width:%dpx; "\
"} "\
".overlaphighlight { "\
"fill:none; "\
"stroke:#9100ff; "\
"stroke-width:15px; "\
"stroke-linecap:round; "\
"stroke-linejoin:round; "\
"stroke-opacity:0.61254614; "\
"} "\
".cluster { "\
"stroke:#60cdf3; "\
"stroke-linejoin:round; "\
"stroke-width:10; "\
"stroke-miterlimit:4; "\
"fill:#60cdf3; "\
"opacity:0.33333333; "\
"} "\
".tLabel { "\
"text-anchor:middle; "\
"text-align:center; "\
"} "\
".tsLabel { "\
"font-size:%dpx; "\
"font-style:normal; "\
"font-weight:normal; "\
"fill:black; "\
"fill-opacity:1; "\
"stroke:none; "\
"font-family:DejaVu Sans; "\
"} "\
".tshtLabel { "\
"font-size:%dpx; "\
"font-style:italic; "\
"font-weight:normal; "\
"fill:black; "\
"fill-opacity:1; "\
"stroke:none; "\
"font-family:DejaVu Sans; "\
"} "\
".frLabel { "\
"font-size:%dpx; "\
"font-style:normal; "\
"font-weight:normal; "\
"fill:black; "\
"fill-opacity:1; "\
"stroke:none; "\
"stroke-width:1px; "\
"stroke-linecap:butt; "\
"stroke-linejoin:miter; "\
"stroke-opacity:1; "\
"font-family:DejaVu Sans; "\
"font-stretch:normal; "\
"font-variant:normal; "\
"text-anchor:middle; "\
"text-align:center; "\
"progression-align:center; "\
"writing-mode:lr; "\
"line-height:125%%; "\
"} "\
'</style>\n' \
$linewidth $linewidth2 $fontsize $fontsize`

infile=${svgfile}.tmp.$$
outfile=${svgfile}
mv ${svgfile} ${infile}
sed "/<style /c\\${styleattr}"  < ${infile} > ${outfile}
rm ${infile}

 
