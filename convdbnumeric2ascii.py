#!/usr/bin/env python
###############################################################################
#
# convdbnumeric2ascii.py - Convert numeric tableaux db to ASCII format
#
# File:    convdbnumeric2ascii.py
# Author:  Alex Stivala
# Created: July 2008
#
#
# Usage:
#    convdbnumeric2ascii.py inputdb > outputfile
#
# (output is to stdout)
#
# This script is standalone, ie does not import any of the ptgraph
# etc. modules. It requires only the Numeric python library.
#
# $Id$
#
###############################################################################

import sys
import pickle
import numpy.oldnumeric as Numeric


"""
This script converts a Numeric (Omega matrix) tableaux database as
built by buildtableauxdb.py -n to a simple fixed field width ASCII
format useful for parsing by other programs (especially FORTRAN).

The format of the input is pickled Numeric.array values saved by
buildtableauxdb.py.

The format of the output 'database' is a text file with an entry for each
structure. The first line of an entry is the identifier and
order of tableau (i.e. dimension of square Omega matrix), then
each subsequent row is a row of the Omega matrix, lower triangle
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

E.g.:

D1UBIA_    6
 0.000 
 2.650  0.000
-1.170  2.150  1.000
 2.040 -1.140  2.080  0.000
-1.260  1.560 -1.110  2.990  0.000
-0.590  2.100 -1.230  2.570 -0.720  0.000

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
    sys.stderr.write("Usage: " + prog + " inputdb > outputfile\n")
    sys.exit(1)
    

def main():
    """
    main for convdbnumeric2ascii.py - load Numeric pickle and output as ascii
    """
    if len(sys.argv) != 2:
        usage(sys.argv[0])

    dbfile = sys.argv[1]
    db = pickle.load(open(dbfile))
    for pdbid,dbtablist in db.iteritems():
        tabnum = 0
        while tabnum < len(dbtablist):
            omega = dbtablist[tabnum]
            n = Numeric.shape(omega)[0]
            name = pdbid
            if len(dbtablist) > 1:
                name += str(tabnum)
            sys.stdout.write('%6s %4d\n' % (name, n))
            for i in xrange(n):
                for j in xrange(i+1):
                    if isNaN(omega[i,j]):
                        angle = 0.0
                    else:
                        angle = omega[i,j]
                    sys.stdout.write('%6.3f ' % angle)
                sys.stdout.write('\n')
            sys.stdout.write('\n')
            tabnum += 1

            
if __name__ == "__main__":
    main()
