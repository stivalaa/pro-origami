#!/usr/bin/env python
###############################################################################
#
# tabsearchdp - TableauSearch using dynamic programming 'alignment-like' method
#
# File:    tabsearchdp.py
# Author:  Alex Stivala
# Created: June 2008
#
# $Id$
#
#
###############################################################################

"""
Search a tableaux database with the alignment-like dynamic programming
approach as per TableauSearch in Konagurthu et al (2008) (section 6).

"""
import os,sys
import getopt
import pickle
from time import strftime,localtime
import resource # for getrusage()

import numpy.oldnumeric as Numeric

from pttableau import PTTableauPacked
from tableaubuild import get_tableaux
from tabsearchlib import *

#-----------------------------------------------------------------------------
#
# Constants
#
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
#
# Class definitions
#
#-----------------------------------------------------------------------------


class SeqAlignNW:
    """
    Sequence alignment using the Needleman-Wunsch (1970) algorithm
    for global alignment. This is a generic implementation so that
    sequences are lists (of any type) and the substitution 'matrix'
    is actually a dictionary hashed by tuples of the type used in the
    sequences, or alternatively a function taking such a tuple as argument.
    """

    #
    # member functions
    #
    
    def __init__(self, seqA, seqB,
                 subst_matrix,
                 subst_function,
                 gap_penalty):
        """
        Construct a SeqAlignNW object for the supplied sequences and
        substitution matrix and gap penalty.

        Parameters:   seqA         - Sequence A to align.
                      seqB         - Sequence B to align.
                      subst_matrix - The substitution matrix. Actually a
                                     dictionary keyed by a tuple of types
                                     in seqA and seqB.
                      subst_function - May be used instead of subst_matrix.
                                      A function taking 2 arguments a,b
                                      where a and b are types
                                      in the sequences, and returning
                                      substitution score for (a,b).
                      gap_penalty  - The cost (<0) of an insertion/deletion.

       Raises excptions:
            ValueError if both subst_matrix and subst_function are provided.
        """

        self.seqA = seqA
        self.seqB = seqB
        if subst_matrix and subst_function:
            raise ValueError('only one of substition matrix or function allowed')
        self.subst_matrix = subst_matrix
        self.subst_function = subst_function # one of matrix or function is None
        self.gap_penalty = gap_penalty
        self.n1 = len(seqA) # length of first sequence
        self.n2 = len(seqB) # length of 2nd sequence

            
        # The H matrix is implemented as a dictionary keyed by (i,j)
        # tuples where i is the end of a subsequence in the
        # first sequence and j is the  end of a subsequence
        # in the second sequence.
        # I.e. keyed by (i,j) for sequence A 0..i and sequenc B 0..j.
        
        self.H = {}

#         # for bottom-up version:
        
#         # we only need to initialize H0_k0 = H0_0l = 0 for 0 <= k <= n
#         # and 0 <= l <= m but easier to init all to 0 using numpy arrays
#         self.H0  = Numeric.zeros((len(seqA)+1, len(seqB)+1))


    def dynprogH(self, i, j):
        """
        The dynamic programming (top-down memoization implementation) of
        the Needleman-Wunsch H matrix for subsequence 0..i in sequence A
        and 0..j in sequence B .

        The H matrix is implemented as a dictionary keyed by (i,j).

        Parameters:
            i  - end position in seqA
            j  - end position in seqB

        Uses members:
          read/write:
            H   - The H dynamic programming matrix (as a dict), 2 dimensional
                   as described above.
           readonly:
            seqA         - first sequence
            seqB         - second sequence
            subst_matrix - subsitution score matrix (as a dictionary)
            subst_function - or subsitutino score function
            gap_penalty  - gap cost

        Return value:
            The value of the H matrix at (i,j).
            
        """
        assert(i >= 0)
        assert(i <= len(self.seqA))
        assert(j >= 0)
        assert(j <= len(self.seqB))


        # if the value here has already been computed then return it

        if self.H.has_key((i,j)):
            return self.H[(i,j)]


        #
        # The initialization case: zero for zero-length sequences
        #

        if i == 0 or j == 0:
            score = 0
            self.H[(i,j)] = score
            return score
            
        #
        # For the recursive case we
        # find the maximum over the 3 cases:
        #   1. seqA[i] aligned with seqB[j]
        #   2. gap in sequence A
        #   3. gap in sequence B
        #
        if self.subst_function:
            subst_score = self.subst_function(self.seqA[i-1], self.seqB[j-1])
        else:
            subst_score = self.subst_matrix[(self.seqA[i-1], self.seqB[j-1])]
            
        aligned = self.dynprogH(i-1, j-1) + subst_score
        gapA = self.dynprogH(i, j-1) + self.gap_penalty
        score = max(aligned, gapA)
        gapB = self.dynprogH(i-1, j) + self.gap_penalty
        score = max(aligned, gapA, gapB)
        self.H[(i,j)] = score
        return score



    def backtraceH(self, i, j):
        """
        Get the list of aligned positions by tracing back through the
        dynamic programming hash table H

        Parameters:
             i  - the i co-ord in H at which to start the backtrace
             j  - the j co-ord in H at which to start the backtrace

        Uses members:
          readonly:
             H - the dynamic programming hashtable as computed
                  by dynprogH()
             seqA - sequence A
             seqB - sequence B
             subst_matrix - subsitution score matrix (as a dictionary)
             subst_function - or substitution score function
             gap_penalty  - gap cost
             

        Return value:
          list of (i,j) tuples where each (i,j) means position i in seqA
          is aligned with position j in seqB, sorted by i,j ascending
        """

        assert(i > 0)
        assert(i <= len(self.seqA))
        assert(j > 0)
        assert(j <= len(self.seqB))

               
        epsilon = self.epsilon
        
        btlist = []

        while i > 0 and j > 0:
            assert(self.H.has_key((i,j)))
            if self.subst_function:
                subst_score = self.subst_function(self.seqA[i-1],self.seqB[j-1])
            else:
                subst_score = self.subst_matrix[(self.seqA[i-1],self.seqB[j-1])]
            if self.H.has_key((i-1, j-1)) and \
                abs(self.H[(i,j)] -
                   (self.H[(i-1, j-1)] +
                    subst_score)) < epsilon:
                btlist.append((i-1, j-1)) # case 1: aligned at i-1 and j-1
                i -= 1
                j -= 1
            elif self.H.has_key((i, j-1)) and \
                 abs(self.H[(i, j)] -
                     (self.H[(i, j-1)] + self.gap_penalty)) < epsilon:
                j -= 1 # case 2: gap in seqA
            elif self.H.has_key((i-1, j)) and \
                 abs(self.H[(i, j)] -
                     (self.H[(i-1, j)] + self.gap_penalty)) < epsilon:
                i -= 1 # case 3: gap in sewqB
            else:
                print 'backtrace failure',i,j
                assert(False)

        btlist.reverse()
        return btlist



    def build_alignment(self, btlist):
        """
        Build the alignment of the two sequences by building two strings
        from the two sequences with '-' inserted for the gaps.
        Since this is really meant for local alignment, a different character
        '~' is used for padding at start where best local match is found
        these don't count as gap penalties.


        Parameters:
           btlist - the list of (i,j) tuples meaning seqA[i] is aligned
                    with seqB[i]. Must be sorted by i (and therefore j)
                    ascending (as returned from backtrace()).

        Uses members:
           seqA - sequence A
           seqB - sequence B

        Return value:
           tuple (alnA, alnB) where alnA is the sequence seqA with gap
           characters '-' inserted as per the supplied aligned pair list
           and similarly for alnB for seqB.
        """

        alnA = ""
        alnB = ""
        previ = 0
        prevj = 0
        for (i,j) in btlist:
            i += 1
            j += 1
            gaplenA = j - prevj
            gaplenB = i - previ

            if previ == 0 and prevj == 0:
                gapchar = '~'
            else:
                gapchar = '-'
                
            if gaplenA > 0:
                alnA += (gaplenA-1) * gapchar
                alnB += self.seqB[prevj : prevj+gaplenA-1]

            if gaplenB > 0:
                alnB += (gaplenB-1) * gapchar
                alnA += self.seqA[previ : previ+gaplenB-1]

            alnA += self.seqA[i-1]
            alnB += self.seqB[j-1]

            previ = i
            prevj = j

        # handle between the last pair in the list and the end of the sequences
        gaplenB = len(self.seqB) - prevj
        gaplenA = len(self.seqA) - previ
        alnA += self.seqA[previ : previ+gaplenA]
        alnB += self.seqB[prevj : prevj+gaplenB]
        if gaplenA < gaplenB:
            alnA += (gaplenB-gaplenA) * '-'
        else:
            alnB += (gaplenA-gaplenB) * '-'
    

        return (alnA, alnB)


#     def bu_dynprogH(self):
#         """
#         Traditional dynamic programming (bottom-up iterative implementation) of
#         the 2-dimensional Needleman-Wunsch H matrix.


#         Uses members:
#           read/write:
#             H0    - The H dynamic programming matrix (bottom-up version)

#            readonly:
#             seqA         - first sequence
#             seqB         - second sequence
#             subst_matrix - The substitution matrix. Actually a
#                            dictionary keyed by a tuple of types
#                            in seqA and seqB.
#             subst_function - May be used instead of subst_matrix.
#                             A function taking 2 arguments a,b
#                             where a and b are types
#                             in the sequences, and returning
#                             substitution score for (a,b).
#             gap_penalty  - gap cost

#         Return value:
#             None. H0 contains the Needleman-Wunsch d.p. matrix.

#         Precondition:
#             H0 is initialized to H0_k0 = H0_0l = 0
#                                  for 0 <= k <= n and 0 <= l <= m
            
#         """


#         #
#         # find the maximum over the 3 cases:
#         #   1. seqA[i-1] aligned with seqB[j-1]
#         #   2. gap in sequence A
#         #   3. gap in sequence B
#         #

#         for i in range(1, len(self.seqA)+1):
#             for j in range(1, len(self.seqB)+1):
#                 if self.subst_function:
#                     subst_score = self.subst_function(self.seqA[i-1],
#                                                       self.seqB[j-1])
#                 else:
#                     subst_score = self.subst_matrix[(self.seqA[i-1],
#                                                      self.seqB[j-1])]
#                 aligned = self.H0[i-1, j-1] + subst_score
#                 gapA = self.H0[i, j-1] + self.gap_penalty
#                 gapB = self.H0[i-1, j] + self.gap_penalty

#                 score = max(aligned, gapA, gapB)
#                 self.H0[i,j] = score
                

# NB bottom-up version actually SLOWER than recursive version!

    def computedp(self):
        """
        Compute the final score from the dynamic programming by calling
        dynprogH() to compute the value at (n1,n2).

        Parameters: None.

        Return value: Value of H(n1, n2), the final dp result.

        Uses members: n1, n2 (readonly)

        """
        return self.dynprogH(self.n1, self.n2)
    

    ##########################################################################

    #
    # member data
    #

    # epsilon is the tolerance for comparing float point numbers for equality
    # needed for the traceback procedure
    epsilon = 1.0e-10 # TODO: find out if numpy has some builtin for this


#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------


def is_tableau_code(tabcode):
    """
    Return True iff two character string is a valie tableau code.
    """
    if ( (tabcode[0] not in ['L','R','P','O'] or
          tabcode[1] not in ['E','D','S','T']) and
         tabcode != 'HH' and tabcode != 'KK' ):
         return False
    else:
        return True

def score_discrete_ssetype(tabcode_a, tabcode_b):
    """
    Return the tableau matching score between the two tableaux entries
    tabcode_a amd tabcode_b, as per Kamat et al (2008) (s5.1).,
    with effecitvely negative infinity score for SSE type mismatch
    """
    if is_tableau_code(tabcode_a) and is_tableau_code(tabcode_b):
        return score_discrete(tabcode_a, tabcode_b)
    else:
        if tabcode_a != tabcode_b:
            return  -99999
        else:
            return 0

def score_numeric_deg_ssetype(omega_a, omega_b):
    """
    Return the tableau matching score between two Omega matrix entries
    omega_a and omega_b, as per Kamat et al (2008),
    with effiectvely negative infinty score for SSE type mismatch

    Parameters:
        omega_a - angle in (-pi, pi]
        omega_b - angle in (-pi, pi]
    Return value:
        score betweem omega_a and omega_b
    """
    if (omega_a not in [0,1,2,3] and omega_b not in [0,1,2,3]):
        return score_numeric_deg(omega_a, omega_b)
    else:
        if omega_a != omega_b:
            return  -99999
        else:
            return 0


def tabmatch_dp(tab_a, tab_b, use_numeric, disallow_type_mismatch=False):
    """
    Match two tableaux using the TableauSearch dp (Konagurthu et al 2008).

    Parameters:
       tab_a - PTTableauPacked (or Numeric.array if use_numeric) for tableau A
       tab_b - PTTableauPacked (or Numeric.array if use_numeric) for tableau B
       use_numeric - if True use Numeric.array Omega matrices not tableaux
       disallow_type_mistmatch - if True give effectively infinite penalty
                      to matches between SSEs of different type (strand/helix)

    Return value:
       score of matching the two tableaux.
    """
    
    DimA = len(tab_a)
    DimB = len(tab_b)

    if use_numeric:
        if disallow_type_mismatch:
            score_func = score_numeric_deg_ssetype
        else:
            score_func = score_numeric_deg
        gap_cost = -20
    else:
        if disallow_type_mismatch:
            score_func = score_discrete_ssetype
        else:
            score_func = score_discrete
        gap_cost = -1

    
    # first we build a scoring matrix by comparing every row in tab_a
    # with every row in tab_b, treating each row as a sequence for
    # alignment.
    score_matrix = Numeric.empty((DimA, DimB), 'd')
    for i in xrange(DimA):
        if use_numeric:
            seq_ai = list(tab_a[i])
        else:
            seq_ai = tab_a.getrow(i)
        for j in xrange(DimB):
            if use_numeric:
                seq_bi = list(tab_b[j])
            else:
                seq_bi = tab_b.getrow(j)
            nw = SeqAlignNW(seq_ai, seq_bi, None, score_func, gap_cost)
            score_matrix[i,j] = nw.computedp()


    #print score_matrix
    
    # now do a sequence alignment of sequences of SSEs using score_matrix
    # just built, which has a score for each pair of SSEs,
    # as the substitution matrix.
    # We can do this just by using sequences of integers reprsenting
    # index of SSEs, which are used pairwise as key to score_matrix
    # ie. aligning SSE i in sequence A with SSE j in sequence B has
    # substitution score score_matrix[i,j]
    seq_a = range(DimA)
    seq_b = range(DimB)
    nw = SeqAlignNW(seq_a, seq_b, score_matrix, None, gap_cost)
    score = nw.computedp()

    #print nw.backtraceH(DimA, DimB)

    return score

            
def tableausearch_dp(tableaux_db, qtab, use_numeric, disallow_type_mismatch):
    """
    Search tableaux database (hash table of PTTableauPacked) for matches
    to tableau tab using the Konagurthu et al (2008) dynamic programming
    TableauSearch.
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
            objval = tabmatch_dp(qtab, dbtab, use_numeric,
                                        disallow_type_mismatch)
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

def xmatch_score(a, b):
    """
    Simple scoring function: 1 for same value, else 0
    """
    if a == b:
        return 1
    else:
        return 0
    

def test_seqalign():
    """
    Print an alignment of two sequences to show how it works.
    """
    s1 = 'GGGGGACCCTTTACAAC'
    s2 = 'GGGCACCCCAAC'
    gapcost = -1
    nw = SeqAlignNW(s1,s2,None,xmatch_score,gapcost)
    score = nw.computedp()
    print 'score = ', score
    (a,b) = nw.build_alignment(nw.backtraceH(len(s1), len(s2)))
    print a + '\n' + b

    # recalculate score to verify
    
    rscore = 0
    assert(len(a) == len(b))
    for i in range(len(a)):
        if a[i] == '-' or b[i] == '-':
            rscore += gapcost
        elif a[i] == b[i]:
            rscore += 1
    print 'recalculated score = ', rscore
    

def test_dp():
    """
    Run on two small tableaux to test.
    """
    
    ta = get_tableaux('/local/charikar/astivala/pdb/d1ubia_.ent', 'dssp',
                      'none',True,True,use_numeric=True)[1][0]
    tb = get_tableaux('/local/charikar/ASTRAL/pdbstyle-1.73/xd/d1xd3b_.ent',
                      'dssp',
                      'none',True,True, use_numeric=True)[1][0]
    
    score = tabmatch_dp(ta, tb, use_numeric=True,
                         disallow_type_mismatch = False)
    print 'score = ',score


#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------

def usage(progname):
    """
    Print usage message and exit
    """
    sys.stderr.write("Usage: " + progname + " [-ntv] "
                     "<dbname> <querytableau>\n")
    sys.exit(1)


def main():
    """
    main for tabsearchdp

    Usage: tabsearchdp [-ntv] <dbname> <querytableau>

    
    -n use Numeric.array Omega matrix rather than discrete tableau
       (both db and query must therefore be numeric not tableaux).

    -t disallow matches between SSEs that are not the same type
       (i.e. between strands and helices, alpha and pi helices, etc.)

    -v turns on debug output to stderr

    <querytableau> is a packed format [sub-]tableau as created by the -o option
                   of pytableaucreate (or Omega matrix for -n)

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

    try:
        opts,args = getopt.getopt(sys.argv[1:], "ntv?")
    except:
        usage(os.path.basename(sys.argv[0]))
    for opt,arg in opts:
        if opt == "-n": # numeric
            use_numeric = True
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
    sys.stdout.write('# $Id$\n'.replace('$', ' '))
    sys.stdout.write('# '+ timestamp + '\n')


    #writes results to stdout
    tableausearch_dp(db, qtabp, use_numeric, disallow_type_mismatch)

    # time.clock() is too problematic (CLOCKS_PER_SEC etc.) so we will
    # use getrusage() instead
    # (Although note this time includes time spent parsing files, etc.)
    user_cpu_time = (resource.getrusage(resource.RUSAGE_SELF).ru_utime + 
                     resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime)
    sys.stdout.write('# User cpu seconds: %.1f\n' % user_cpu_time)
    

if __name__ == "__main__":
    main()
