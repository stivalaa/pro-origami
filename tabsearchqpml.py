#!/usr/bin/env python
###############################################################################
#
# tabsearchqpml - tableau search as relaxed QP solved using MATLAB spsolqp
#
# File:    tabsearchqpml.py
# Author:  Alex Stivala
# Created: May 2008
#
# $Id$
#
#
###############################################################################

"""
Search a tableaux database and allow using the relaxed quadratic
programming formulation as per MNAligner (see MNAligner directory and
MATLAB code therein).

This module uses Prof. Y. Ye's SPSOLQP MATLAB implementation of interior
point method for the solution (see spsolqp.m for references) so we
need to interface to MATLAB. This is done with the mlabwrap module
(http://mlabwrap.sourceforge.net/), so we require the python libraries:

. mlablwrap
  http://mlabwrap.sourceforge.net/
  v1.0 was used.

. numpy
  http://numpy.scipy.org/
  v1.01 was used
  
MATLAB Version 7.2.0.283 (R2006a) was used.

TODO: reimplement the nonconvenex qudratic programming solver in
FORTRAN or C and use that instead for greater efficiency and removing
necessity of MATLAB and interface to it.
"""
import os,sys
import getopt
import pickle
from time import strftime,localtime

from numpy import *
from mlabwrap import mlab

from pttableau import PTTableauPacked
from tableaubuild import get_tableaux
from tabsearchlib import *

#-----------------------------------------------------------------------------
#
# Constants
#
#-----------------------------------------------------------------------------

# Location of mfiles e.g. fspsolqp.m etc.
# FIXME should do this in some better way
MFILES_PATH = "/home/charikar/pgrad/astivala/phd/MNAligner"

#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------
 
def tabmatch_qp_matlab(tab_a, tab_b, use_numeric, disallow_type_mismatch,
                       use_ordering = False):
    """
    Match two tableaux using relaxed quadratic programming solved with
    fast interior point method via MATLAB SPSOLQP implementation.

    Parameters:
       tab_a - PTTableauPacked (or Numeric.array if use_numeric) for tableau A
       tab_b - PTTableauPacked (or Numeric.array if use_numeric) for tableau B
       use_numeric - if True use Numeric.array Omega matrices not tableaux
       disallow_type_mistmatch - if True give effectively infinite penalty
                      to matches between SSEs of different type (strand/helix)
       use_ordering - If True give effectively infinte penalty to matches
                     between SSEs not mainiting sequence order between
                     the tableuax i.e. if i < k and j >= l for i,k
                     indices in tab_a and j,l indices in tab_b.

    Return value:
     tuple (matchmat, ObjVal) where
       ObjVal   - final objective value
    """

    # This is python version of TableauMatch.m - hopefully setting up
    # arrays In Python then passing to MATLAB is not as slow as setting
    # it all up in MATLAB (but who knows unti I try it?)
    DimA = len(tab_a)
    DimB = len(tab_b)

    # Q: Sparse symmetric objective matrix.
    # A: Sparse constraint left-hand matrix
    # b: constraint right-hand column vector
    # c: objective column vector
    
    b=zeros((DimA+DimB,1))
    c=zeros((DimA*DimB+DimA+DimB,1))
    Q=zeros((DimA*DimB+DimA+DimB, DimA*DimB+DimA+DimB))
    A=zeros((DimA+DimB, DimA*DimB+DimA+DimB))

    for i in xrange(DimA):
        for j in xrange(DimB):
            A[i, i*DimB+j] = 1

    for i in xrange(DimB):
        for j in xrange(DimA):
            A[DimA+i, j*DimB+i] = 1

    for i in xrange(DimA+DimB):
        A[i, DimA*DimB+i] = 1
        b[i, 0] = 1

    for i in xrange(DimA):
        for j in xrange(DimB):
            for k in xrange(DimA):
                for l in xrange(DimB):
                    if i != k and j != l:  # diagonals are SSE type not angle
                        if (disallow_type_mismatch and
                            tab_a[(i,i)] != tab_b[(j,j)]):
                            # severely penalize matches between SSEs
                            # that are not of the same type (strand/helix)
                            Q[i*DimB+j, k*DimB+l] = -999
                        else:
                            if (use_ordering and  
                                ((i < k and j > l) or
                                 (i > k and j < l))): # SSEs out of order
                                Q[i*DimB+j, k*DimB+l] = -999
                            elif use_numeric:
                                Q[i*DimB+j, k*DimB+l] = \
                                      score_numeric(tab_a[(i,k)], tab_b[(j,l)])
                            else:
                                Q[i*DimB+j, k*DimB+l] = \
                                     score_discrete(tab_a[(i,k)], tab_b[(j,l)])

    Q = -Q
    c = -c

    # Removed matchmat and y output parameter, maybe faster not having them
    #matchmat - matching matrix where (i,j) is 1.0 for SSE i matching SSE j
    #     x, y, obhis = mlab.fspsolqp(Q, A, b, c, nout=3)
    #     matchmat = reshape(x[:DimA*DimB],(DimB,DimA))
    #     return (matchmat, obhis[0,shape(obhis)[1]-1])
    
    try:
        obhis  = mlab.f1spsolqp(Q, A, b, c) # this interface returns only obhis
    except:
        # sometimes we get errors like:
        #
        # File "/usr/local/apps/python-2.5.1//lib/python/mlabwrap.py",
        # line 508, in _do handle_out(mlabraw.eval(self._session,
        # '[%s]=%s;' % (", ".join(resSL), cmd))) mlabraw.error:
        # Undefined function or variable "y".
        #
        # but only occasionally, and not consistently.
        # So we'll just give up when it happen and hope it doesn't happen
        # too much.
        sys.stderr.write('Exception caught, return NaN:\n' + str(sys.exc_info()[0]) +'\n')
        return float('NaN')
        

    return obhis[0,shape(obhis)[1]-1]


    
def tabsearch_qp_matlab(tableaux_db, qtab, use_numeric, disallow_type_mismatch,
                        use_ordering):
    """
    Search tableaux database (hash table of PTTableauPacked) for matches
    to tableau tab using relaxed quadratic programming solved with
    fast interior point method via MATLAB SPSOLQP implementation.
    This is currently doing nothing but computing match score for
    each tableau in db and building list of them.

    Currently the results are just lines on stdout, each line is

    SCOPsid score

    e.g.:
    
    d1qlpa_ -80.9998932482
    d1q9ha_ -73.999871453
    ...

    Parameters:
       tableaux_db -  dict of { pdbid : [PTTableauPacked list] }
                           or { pdbid : [Numeric.array list] } if use_numeric
       qtab - PTTableauPacked representation of tableau to search for
       use_numeric - if True use numeric Omega matrix rather than tableau
       disallow_type_mistmatch - if True give effectively infinite penalty
                      to matches between SSEs of different type (strand/helix)
       use_ordering - If True give effectively infinte penalty to matches
                     between SSEs not mainiting sequence order between
                     the tableuax i.e. if i < k and j >= l for i,k
                     indices in tab_a and j,l indices in tab_b.
                      

    Return value:
       list of (score, pdbid) tuples.

    """
    # TODO: more intelligent procesing of list so instead of just putting
    # all in list and sorting later, throwing away those
    # that cannot be in top n or something.
    num_tableaux = 0
    if verbose:
        sys.stderr.write('searching db...\n')
    scorelist = [] # list of (score, pdbid) tuples
    for pdbid,dbtablist in tableaux_db.iteritems():
        for dbtab in dbtablist:
            objval = tabmatch_qp_matlab(qtab, dbtab, use_numeric,
                                        disallow_type_mismatch,
                                        use_ordering)
            sys.stdout.write(pdbid + ' ' +  str(objval) + '\n')
            num_tableaux += 1
            scorelist.append((objval, pdbid))
    if verbose:
        sys.stderr.write('searched %d tableaux.\n' % num_tableaux)
    return scorelist



############################################################################
#
#  Testing functions
#
#

# neat way to get a small subset of db (scop domains named 'd?a????'):
#   smalldb = dict(itertools.ifilter(lambda kv : kv[0][2]=="a" ,db.iteritems()))


def test_qpml():
    """
    Run on two small tableaux to test.
    """
    
    ta = get_tableaux('/local/charikar/astivala/pdb/1HH1.pdb', 'dssp',
                      'none',False,False,use_numeric=True)[1][0]
    tb = get_tableaux('/local/charikar/astivala/pdb/1CD8.pdb', 'dssp',
                      'none',False,False, use_numeric=True)[1][0]
    
    objval = tabmatch_qp_matlab(ta, tb, use_numeric=True,
                                disallow_type_mismatch = True,
                                use_ordering = True)
    print 'objval = ',objval


#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------

def usage(progname):
    """
    Print usage message and exit
    """
    sys.stderr.write("Usage: " + progname + " [-ontv] "
                     "<dbname> <querytableau>\n")
    sys.exit(1)


def main():
    """
    main for tabsearchqpml

    Usage: tabsearchqpml [-ntv] <dbname> <querytableau>

    
    -n use Numeric.array Omega matrix rather than discrete tableau
       (both db and query must therefore be numeric not tableaux).

    -o disallow matches between SSEs that are not in the same sequence
       order in the two tableaux. (Constrain matches so that order
       is preserved).

    -t disallow matches between SSEs that are not the same type
       (i.e. between strands and helices, alpha and pi helices, etc.)

    -v turns on debug output to stderr

    <querytableau> is a packed format [sub-]tableau as created by the -o option
                   of pytableaucreate or Omega matrix for -n

    <dbname> is a tableaux database as creaetd by buildtableaudb

    Note that the same options should be used for building the tableaux
    database and the query tableau (e.g. if one has -35 for pi/310 helices
    then so should the other, otherwise search is likely to be less
    successful).
    
    """
    global verbose
    verbose = False
    use_numeric = False
    disallow_type_mismatch = False
    use_ordering = False

    mlab.addpath(MFILES_PATH) #FIXME shoudl set this up some other way

    try:
        opts,args = getopt.getopt(sys.argv[1:], "notv?")
    except:
        usage(os.path.basename(sys.argv[0]))
    for opt,arg in opts:
        if opt == "-n": # numeric
            use_numeric = True
        elif opt == "-o": # ordering constraints
            use_ordering = True
        elif opt == "-t": # disallow matches between differnt type SSEs
            disallow_type_mismatch = True
        elif opt == "-v": # verbose
            verbose = True # this module only
        else:
            usage(os.path.basename(sys.argv[0]))

    if len(args) != 2:
        usage(os.path.basename(sys.argv[0]))

    dbfilename = args[0]
    query_tabfile = args[1]

    qtabp = pickle.load(open(query_tabfile))
    if not use_numeric:
        if not isinstance(qtabp, PTTableauPacked):
            sys.stderr.write('invalid tableau file ' + query_tabfile + '\n')
            sys.exit(1)
    db = load_db(dbfilename, use_numeric, verbose)

    timestamp = strftime("%d%b%Y %H:%M:%S", localtime())
    sys.stdout.write('# ' + ' '.join(sys.argv) + '\n')
    sys.stdout.write('# $Id$\n'.replace('$',' '))
    sys.stdout.write('# '+ timestamp + '\n')
   

    #writes results to stdout
    tabsearch_qp_matlab(db, qtabp, use_numeric, disallow_type_mismatch,
                        use_ordering)
    
if __name__ == "__main__":
    main()
