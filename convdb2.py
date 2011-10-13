#!/usr/bin/env python
###############################################################################
#
# convdb2.py - Convert PTTableauPacked tableaux db plus Numeri distance matrix
#              to single ASCII format file.
#
# File:    convdb2.py
# Author:  Alex Stivala
# Created: August 2008
#
#
# Usage:
#    convdb2.py [-l] [-s] inputtableauxdb inputdistmatrixdb [< inputlist] > outputfile
#
# (output is to stdout)
# Input on stdin is list (one per line) of identifiers to include in output,
# if -l i specified.
# Sort the db by tableau size if -s is specified.
#
# Requires the Numeric library, as well as the pttableau module.
#
# $Id$
#
###############################################################################

import sys
import getopt
import pickle
import numpy.oldnumeric as Numeric
from pttableau import PTTableauPacked
from ptutils import isNaN


"""
This script converts tableaux database in pickled PTTableauPacked format
and distance matrix database in pickled Numeric.array format
built by buildtableauxdb.py to a simple fixed field width ASCII
format useful for parsing by other programs (especially FORTRAN).
There must be both a tableau and distance matrix for each entry (i.e.
the two input databases contain data for same identifiers).

The format of the tableau input is pickled PTTableauPacked objects built by
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

Following the tableau is the distance matrix.
Each row is a row of the distance matrix, lower triangle
only (since it is symmetric).
The diagonal entries are meaningless (self-distance)
and are included instead to specify the SSE type, with
the following codes:

0.000 beta strand
1.000 alpha helix
2.000 pi helix
3.000 3_10 helix

Each entry in matrix is in Angstroms format
F6.3 with a space between each on a line, and one line
per row of matrix.
NB any NaN values are converted to 0.000 in the output.


E.g.:

/local/charikar/astivala/tableauxdb/astral/tableauxdb.ascii
 T F
D1UBIA_    8
e  
OT e  
LE RT xa 
PD OS RD xg 
RT LE RT LS e  
LE RD LE LS OT e  
RT LS LS RD PE OS xg 
PE RT LE RD OT PE RT e  
 0.000 
 4.501  0.000 
 1.662 10.386  1.000 
16.932 17.644  9.779  3.000 
10.588 13.738 11.815 10.527  0.000 
15.025 18.692 17.143 15.341  6.466  0.000 
15.298 17.276 16.276 20.075 13.264 11.610  3.000 
 7.549 11.072 12.248 12.446  4.583  9.903 15.689  0.000 

There is a blank line between each entry.

"""



def sizecmp(dbent1, dbent2):
    """
    Comparison function for (pdbid,dbtablist)
    tuples used to sort database by size
    """
    tab1 = dbent1[1][0]
    tab2 = dbent2[1][0]
    if len(tab1) < len(tab2):
        return -1
    elif len(tab1) > len(tab2):
        return 1
    else:
        return 0
    

def usage(prog):
    """
    print usage message and exit
    """
    sys.stderr.write("Usage: " + prog + " [-ls] inputtableauxdb inputdistmatrixdb[< inputlist]  > outputfile\n")
    sys.stderr.write("   -l : read list of SIDs on stdin\n")
    sys.stderr.write("   -s : sort database by size\n")
    sys.exit(1)
    

def main():
    """
    main for convdb2.py - load PTTableauPacked pickle and Numeric.array
    distance matrix pickle and 
    output as ascii
    """
    use_sidlist = False
    do_sort = False

    try:
        opts, args = getopt.getopt(sys.argv[1:], "ls")
    except getopt.GetoptError:
        usage(sys.argv[0])
        
    for opt,arg in opts:
        if opt == "-l":  # read list of identifeirs on stdin
            use_sidlist = True
        elif opt == "-s": # sort db  by size
            do_sort = True
        else:
            usage(sys.argv[0])
    
    if len(args) != 2:
        usage(sys.argv[0])

    dbfile = args[0]
    distmatrixdbfile = args[1]

    if use_sidlist:
        # make dictionary of identifiers from list on stdin
        sid_dict = dict([(x.strip(), None) for x in sys.stdin.readlines()])

    sys.stderr.write('loading tableaux...\n')
    db = pickle.load(open(dbfile))
    sys.stderr.write('loading distance matrices...\n')
    distmatrixdb = pickle.load(open(distmatrixdbfile))
    sys.stderr.write('writing ASCII format tableaux+distmatrices...\n')
    first = True
    count = 0
    total_count = 0
    dblist = [ (pdbid,dbtablist) for (pdbid,dbtablist) in db.iteritems()
               if len(dbtablist) > 0] # remove those with no tableaux
    if do_sort:
        sys.stderr.write("sorting database...\n")
        dblist.sort(cmp=sizecmp)
        
    for pdbid,dbtablist in dblist:
        if use_sidlist and not sid_dict.has_key(pdbid):
            continue
        tabnum = 0
        while tabnum < len(dbtablist):
            tableau = dbtablist[tabnum]
            n = len(tableau)
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
            if len(distmatrix) != len(tableau):
                sys.stderr.write('ERROR: dist matrix order ' + 
                    str(len(distmatrix)) + ' but tableau order ' +
                    str(len(tableau)) + ' for id ' + pdbid + 
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
                    sys.stdout.write(tableau[(i,j)] + ' ')
                sys.stdout.write('\n')
            for i in xrange(n):
                for j in xrange(i+1):
                    if isNaN(distmatrix[i,j]):
                        dist = 0.0
                    else:
                        dist = distmatrix[i,j]
                    sys.stdout.write('%6.3f ' % dist)
                sys.stdout.write('\n')
            total_count += 1
            tabnum += 1
        count += 1
    sys.stderr.write('wrote %d tableaux+distmatrices for %d entries\n' 
                       %  (total_count, count))

            
if __name__ == "__main__":
    main()
