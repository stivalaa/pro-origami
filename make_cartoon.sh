#!/bin/sh
#
# make_cartoon.sh - make SVG and PNG cartoon for a PDB file
#
# Usage: make_cartoon.sh pdbfile
#
# Runs ptgraph2.py, Dunnart and svg2png.sh to make the .svg and .png
# versions of the cartoon for given PDB file. The PDB file can be 
# in any directory, but output is in cwd.
#
#
# $Id$
#

#PTGRAPH2_OPTIONS="-r35 -t dssp -k purple -l crossing:black,red,green,navy,blue -b sequential -j -e auto -f auto -o gradient -p ddomain"
PTGRAPH2_OPTIONS="-r35 -t dssp -k purple -l crossing:black,red,green,navy,blue -b sequential -j -e auto -f auto -o gradient -p cath:/local/munk/proorigami-test/cath/CathDomall.v3.3.0"

#DUNNART=$HOME/dunnart/trunk/dunnart
DUNNART=/local/munk/proorigami-prod/dunnart/trunk/dunnart

# time in seconds to limit Dunnart CPU time to before assuming failure
DUNNART_CPUTIME_LIMIT=150

DUNNART_OPTS="-b -y -w4 -z100"

#  this uses the null video device driver in SDL
export SDL_VIDEODRIVER=dummy

if [ $# -ne 1 ]; then
    echo "usage $0 pdbfile" >&2
    exit 1
fi
pdbfile=$1

tmpfile=/var/tmp/mkptg$$

ptgraph2.py $PTGRAPH2_OPTIONS $pdbfile > $tmpfile

if [ $? -ne 0 ]; then
    echo "ptgraph2.py failed"
    rm $tmpfile
    exit 1
fi

filelist=`grep 'writing file' $tmpfile |  cut -d' ' -f3`

for svgfile in $filelist ; do
    (ulimit -t ${DUNNART_CPUTIME_LIMIT} ; $DUNNART $DUNNART_OPTS ${svgfile} ) > /dev/null 2>&1 
    if [ $? -ne 0 ]; then
        echo dunnart failed on $svgfile >&2
     else
        cp $svgfile ${svgfile}.orig
    fi
done
restore_list=""
for svgfile in $filelist ; do
        overlapcount=`grep count $svgfile | sed 's/.*count="\([0-9]*\).*/\1/'`
        orig_overlapcount=$overlapcount
        echo "overlap count for $svgfile (default) is $overlapcount"
        # rerun with progressively larger gap sizes until we get no overlaps
        # or gap size is 'too big'
        gapsize=56  # starts at default 55
        while [ $overlapcount -gt 0 -a $gapsize -le 80 ];
        do
            ptgraph2.py $PTGRAPH2_OPTIONS -g $gapsize $pdbfile > $tmpfile
            if [ $? -ne 0 ]; then
                echo "ptgraph2.py failed"
                rm $tmpfile
                exit 9
            fi
            (ulimit -t ${DUNNART_CPUTIME_LIMIT} ; $DUNNART $DUNNART_OPTS ${svgfile} ) > /dev/null 2>&1 
            if [ $? -ne 0 ]; then
                echo dunnart failed on $svgfile >&2
            fi
            overlapcount=`grep count $svgfile | sed 's/.*count="\([0-9]*\).*/\1/'`
            echo "overlap count for $svgfile (-g $gapsize) is $overlapcount"
            gapsize=`expr $gapsize + 1`
        done
        # if there was no improvement, just use default gapsize version
        if [ $overlapcount -ge $orig_overlapcount ]; then
            old_overlapcount=`grep count ${svgfile}.orig | sed 's/.*count="\([0-9]*\).*/\1/'`
            echo "going back to to default gapsize (overlap count $old_overlapcount) for $svgfile"
            restore_list="${restore_list} ${svgfile}"
        fi
done

for restore_file in $restore_list ; do
    cp ${restore_file}.orig ${restore_file}
done

for svgfile in $filelist ; do
    svg2png.sh $svgfile
    rm ${svgfile}.orig
done

rm $tmpfile

