#!/usr/bin/env python
###############################################################################
#
# convdbomega2.py - Convert numeric tableaux db plus Numeric distance matrix
#              to single ASCII format file.
#
# File:    convdbomega2.py
# Author:  Alex Stivala
# Created: November 2009
#
#
# Usage:
#    convdbomega2.py inputnumericdb inputdistmatrixdb > outputfile
#
# (output is to stdout)
#
# Requires the Numeric library.
#
# $Id$
#
###############################################################################

import sys
import pickle
import numpy.oldnumeric as Numeric
from ptutils import isNaN


"""
This script converts Omega matrix database in pickled Numeric.array format
and distance matrix database in pickled Numeric.array format
built by buildtableauxdb.py to a simple fixed field width ASCII
format useful for parsing by other programs (especially FORTRAN).
There must be both a tableau and distance matrix for each entry (i.e.
the two input databases contain data for same identifiers).

The format of the tableau input is pickled Numeric.array objects built by
buildtableauxdb.py, and the distance matrix input is pickled Numeric.array
objects built by buildtableauxdb.py (with -d option).


The format of the 'database' is a text file with an entry for each
structure.
The first line of an entry is the identifier and
order of tableau (i.e. dimension of square array), then
each subsequent row is a row of the tableau, lower triangle
only (since it is symmetric).
The diagonal entries are meaningless (self-angle) in tableaux,
and are included instead to specify the SSE type, with
the following codes:

0.000 beta strand
1.000 alpha helix
2.000 pi helix
3.000 3_10 helix

Width of identifier is 8 chars, blank padded on right,
width of order is 4 digits, blank padded on left.
There is a single space between identifier and order.
Each entry in Omega matrix is in radians in [-pi, pi] format 
F6.3 with a space between each on a line, and one line
per row of matrix.

Following the tableau is the distance matrix.
Each row is a row of the distance matrix, lower triangle
only (since it is symmetric).
The diagonal entries are meaningless (self-distance)
and are included instead to specify the SSE type, with
the same codes as the Omega matrix above.

Each entry in matrix is in Angstroms format
F6.3 with a space between each on a line, and one line
per row of matrix.
NB any NaN values are converted to 0.000 in the output.


E.g.:

 /local/charikar/astivala/tableauxdb/astral/tableauxdb.numeric.dmat.ascii

 D1UBIA_    6
  0.000 
  2.650  0.000
 -1.170  2.150  1.000
  2.040 -1.140  2.080  0.000
 -1.260  1.560 -1.110  2.990  0.000
 -0.590  2.100 -1.230  2.570 -0.720  0.000
  0.000 
  4.501  0.000 
 11.662 10.386  1.000 
 10.588 13.738 11.815  0.000 
 15.025 18.692 17.143  6.466  0.000 
  7.549 11.072 12.248  4.583  9.903  0.000 

There is a blank line between each entry.

Note any NaN values are converted to 0.0 in the output.
"""

def isNaN(x):
    """
    Test if supplied float is an IEEE not-a-number (NaN).
    For some reason Python does not hav a function to do this,
    and nor does Numeric (although numpy and scipy have support for it).
    
    Parameters:
        x - float to test for NaN

    Return value:
        True if x is NaN, else False.
    """
    # NaN is the only float value that is not equal to itself (IEEE
    # standard)
    if x != x:
        return True
    else:
        return False
    

def usage(prog):
    """
    print usage message and exit
    """
    sys.stderr.write("Usage: " + prog + " inputtableauxdb inputdistmatrixdb > outputfile\n")
    sys.exit(1)
    

def main():
    """
    main for convdbomega2.py - load Omega Numeric.array pickle and
    Numeric.array distance matrix pickle and output as ascii
    """
    if len(sys.argv) != 3:
        usage(sys.argv[0])

    dbfile = sys.argv[1]
    distmatrixdbfile = sys.argv[2]

    db = pickle.load(open(dbfile))
    distmatrixdb = pickle.load(open(distmatrixdbfile))
    first = True
    for pdbid,dbtablist in db.iteritems():
        tabnum = 0
        while tabnum < len(dbtablist):
            omega = dbtablist[tabnum]
            n = Numeric.shape(omega)[0]
            name = pdbid
            if len(dbtablist) > 1:
                name += str(tabnum)
            try:
                distmatrix = distmatrixdb[pdbid][tabnum]
            except KeyError:
                sys.stderr.write('ERROR: no distance matrix for id ' +
                    pdbid + ' - skipped\n')
                tabnum += 1
                continue
            if len(distmatrix) != n:
                sys.stderr.write('ERROR: dist matrix order ' + 
                    str(len(distmatrix)) + ' but tableau order ' +
                    str(n) + ' for id ' + pdbid + 
                    ' - skipped\n')
                tabnum += 1
                continue
            if not first:
                sys.stdout.write('\n')
            else:
                first = False
            sys.stdout.write('%6s %4d\n' % (name, n))
            for i in xrange(n):
                for j in xrange(i+1):
                    if isNaN(omega[i,j]):
                        angle = 0.0
                    else:
                        angle = omega[i,j]
                    sys.stdout.write('%6.3f ' % angle)
                sys.stdout.write('\n')
            for i in xrange(n):
                for j in xrange(i+1):
                    if isNaN(distmatrix[i,j]):
                        dist = 0.0
                    else:
                        dist = distmatrix[i,j]
                    sys.stdout.write('%6.3f ' % dist)
                sys.stdout.write('\n')
            tabnum += 1

            
if __name__ == "__main__":
    main()
