#!/usr/bin/env python
###############################################################################
#
# convdbpacked2ascii.py - Convert PTTableauPacked tableaux db to ASCII format
#
# File:    convdbpacked2ascii.py
# Author:  Alex Stivala
# Created: July 2008
#
#
# Usage:
#    convdbpacked2ascii.py inputdb > outputfile
#
# (output is to stdout)
#
# $Id$
#
###############################################################################

import sys
import pickle
from pttableau import PTTableauPacked


"""
This script converts tableaux database in pickled PTTableauPacked format
built by buildtableauxdb.py to a simple fixed field width ASCII
format useful for parsing by other programs (especially FORTRAN).

The format of the input is pickled PTTableauPacked objects built by
buildtableauxdb.py.


The format of the 'database' is a text file with an entry for each
structure.
The first line of an entry is the identifier and
order of tableau (i.e. dimension of square array), then
each subsequent row is a row of the tableau, lower triangle
only (since it is symmetric).
The diagonal entries are meaningless (self-angle) in tableaux,
and are included instead to specify the SSE type, with
the following codes:

e     beta strand
xa    alpha helix
xi    pi helix
xg    3_10 helix

Width of identifier is 8 chars, blank padded on right,
width of order is 4 digits, blank padded on left.
There is a single space between identifier and order.
Each entry in tableau is two characters, with a space betwen
each on a line, and one line
per row of matrix.

E.g.:

/local/charikar/astivala/tableauxdb/astral/tableauxdb.ascii
 T F
D1UBIA_    6
e  
OT e  
LE RT xa 
RT LE RT e  
LE RD LE OT e  
PE RT LE OT PE e  

There is a blank line between each entry.

"""

def usage(prog):
    """
    print usage message and exit
    """
    sys.stderr.write("Usage: " + prog + " inputdb > outputfile\n")
    sys.exit(1)
    

def main():
    """
    main for convdbpacked2ascii.py - load PTTableauPacked pickle and 
    output as ascii
    """
    if len(sys.argv) != 2:
        usage(sys.argv[0])

    dbfile = sys.argv[1]
    db = pickle.load(open(dbfile))
    for pdbid,dbtablist in db.iteritems():
        tabnum = 0
        while tabnum < len(dbtablist):
            tableau = dbtablist[tabnum]
            n = len(tableau)
            name = pdbid
            if len(dbtablist) > 1:
                name += str(tabnum)
            sys.stdout.write('%6s %4d\n' % (name, n))
            for i in xrange(n):
                for j in xrange(i+1):
                    sys.stdout.write(tableau[(i,j)] + ' ')
                sys.stdout.write('\n')
            sys.stdout.write('\n')
            tabnum += 1

            
if __name__ == "__main__":
    main()
