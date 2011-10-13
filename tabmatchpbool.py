#!/usr/bin/env python
###############################################################################
#
# tabmatchpbool - formulate tableau matching as pseudo-Boolean problem
#
# File:    tabmatchpbool.py
# Author:  Alex Stivala
# Created: May 2008
#
# $Id$
#
#
###############################################################################

"""
Given two input PDB files, run MiniSat+ or CPLEX on a psuedo-Boolean problem
formulation of the tableau matching problem between the two tableaux
for those proteins, based on the formulation of the problem as QP
and ILP in:

Konagurthu, Stuckey and Lesk 2008 'Structural search and retrieval using
a tableau representation of protein folding patterns' 
Bioinformatics 24(5):645-651.

The problem forumation is in format for input to MinSat+
(http://minisat.se/MiniSat+.html). The reference for MinSat+ is
Een and Sorensson (2006) 'Translating Pseuso-Boolean Constraints into SAT'
Journal on Satisfiability, Boolean Modelling and Computation 2:1-25

The program 'minisat+' must be in the PATH, it is run using popen().
The output is parsed and put to stdout as a list of tuples (i,j) where
SSE i in tableau a matches SSE j in tableau b (zero-based).

Example usage:

 tabmatchpbool.py 1QLP.pdb 7API.pdb

Filenames may be either in the format above or the pdbq1lp.pdb format.
Compressed pdb files are supported (gzip) (e.g. pdb1qlp.ent.gz).
See usage description in docstring for main().

It is written in Python and depends on some Python libraries:

. BioPython (including Bio.PDB)
  http://www.biopython.org

  Reference for Bio.PDB is:
  Hamelryck and Manderick 2003 "PDB parser and structure class implemented
  in Python" Bioinformatics 19:2308-2310

  which in turn depends on Numeric
  http://sourceforge.net/projects/numpy


Developed on Linux 2.6.9 (x86_64) with Python 2.5.1
and BioPython 1.43 with Numeric 24.2
"""

import warnings # so we can suppress the annoying tempnam 'security' warning
import sys,os
import getopt
import re
import numpy.oldnumeric as Numeric
from Bio.PDB import *

import ptsecstruct
from ptnode import ptnode_set_verbose
from ptdomain import *
from ptutils import cleanup_tmpdir
import getdomains
from tableaubuild import TableauBuild

#-----------------------------------------------------------------------------
#
# Class definitions 
#
#-----------------------------------------------------------------------------


#
# Empty classes for exceptions
#

class NoSSE_Exception(Exception): # raised when no helices or strands found
    pass
            
#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------



def get_tableau(pdb_filename,
                pdb_struct,
                secstruct_program,
                domain_program,
                include_310_helices = False,
                include_pi_helices = False,
                sse_id_list = None):
    """
    For the supplied filemame, read PDB format data from that file
    and return tableau for that structre. Note for multidomains,
    returns a list of tableaux (so for single domain, return value is a list
    of one PTTableau).

    Paramteters:
       pdb_filename - filename of PDB file to read
       pdb_struct - Bio.PDB parsed PDB structure
       secstruct_program - secondary structure definition program
                       ('stride' or 'dssp' or 'pdb') to use.
       domain_progam - domain decompositino method ('ddomain','cath', etc.)
       include_310_helices - if True, include 3_10 helices in the graph
       include_pi_helices - if True, include pi helices in the graph
       sse_id_list - list of ints representing SSE sequential id numbers
                     to include in tableau. Default None.
                     When None, all SSEs are included.

    Return value: List of PTTableau objects, one per domain.
    """
    (pdbid,suffix) = os.path.splitext(os.path.basename(pdb_filename))
    pdbid = pdbid.upper()
    if len(pdbid) >= 6 and pdbid[:3] == "PDB":
        pdbid = pdbid[3:7]

    if secstruct_program == "pdb":
        secstruct = ptsecstruct.read_secstruct_from_pdb_file(pdb_filename)
        if secstruct != None:
            secstruct.pdb_header = pdb_struct.header['head']
        else:
            secstruct_program = "dssp"
            sys.stderr.write('WARNING: error with HELIX or SHEET cards in PDB'
                             ': ' + secstruct_program +
                             ' will be used instead\n')
    else:
        secstruct = None

    if secstruct == None:
        # read secondary structure information from STRIDE or DSSP
        if secstruct_program == "stride":
            secstruct = ptsecstruct.read_secstruct_from_stride(pdb_filename)
        elif secstruct_program == "dssp":
            secstruct = ptsecstruct.read_secstruct_from_dssp(pdb_filename)
        else:
            assert(False)


    if domain_program != None:
        domain_list = getdomains.get_domains(domain_program,
                                             pdbid, pdb_filename, pdb_struct)
    else:
        domain_list = [PTDomain(None, None)] # one-domain protein, no further info

    tableau_list = []
    for domain in domain_list:
        ptg = TableauBuild(pdb_struct, pdbid,
                           include_310_helices, include_pi_helices)
        # build tableaubuild object from secondary structure
        try:
            ptg.build_graph_from_secstruct(secstruct, domain)
        except NoSSE_Exception:
            sys.stderr.write('WARNING: No helices or strands found in ' +
                             pdbid +
                             ': skipping domain\n')
            continue


        if verbose:
            for nodelist in ptg.iter_chains():
                for node in nodelist:
                    sys.stderr.write(str(node) + '\n')

        # if list of int SSE sequential ids supplied, convert to list of
        # PTNode objects
        if sse_id_list:
            try:
                ptnode_list = [ptg.seqnum2node[sse_id] for sse_id in sse_id_list]
            except KeyError,k:
                sys.stderr.write("SSE sequential id " + str(k)
                                 + " does not exist\n")
                sys.exit(1)
        else:
            ptnode_list = None

        # NB must ensure we do not use the HH and KK codes, which would
        # break the tableau scoring function score()
        ptg.build_tableau(pdbid, domain, ptnode_list, use_hk=False)
        
        tableau_list.append(ptg.tableau)
    return tableau_list


def score(tabcode_a, tabcode_b):
    """
    Return the tableau matching score between the two tableaux entries
    tabcode_a amd tabcode_b, as per Kamat et al (2008) (s5.1). This
    score is -2 if the tableau entries are equal, -1 if they are equal
    in only one position, else 2.
    NB these are negated since MinSat+ minimizes objective functoin.
    """
    if tabcode_a[0] == tabcode_b[0]:
        if tabcode_a[1] == tabcode_b[1]:
            return -2
        else:
            return -1
    elif tabcode_a[1] == tabcode_b[1]:
        return -1
    else:
        return 2

def linear_objective_function_coefficients(tableau_a, tableau_b):
    """
    Build the integer line objective function, a summation over i,k in
    tab_a and j,l in tab_b of m^2*n^2 terms where m is order of tableau a
    and n is order of tableau b. Each term is the score function at (i,j,k,l)
    multiplied by the indicator variable x_{ijkl}, where x_{ijkl} (boolean)
    is 1 if SSE i in A is matched with SSE j in B AND SSE k in A is matched
    with SSE l in B, else 0.

    WARNING: Note the large (quartic) number of terms: eg this comes to 810 000
    terms if m = n = 30.

    This is a generator function that yields one coefficient at a time.
    
    Parameters:
        tableau_a - PTTableau object for protein A
        tableau_b - PTTableau object for protein B

    Return value:
         (generator function; yields int at a time from:)
         list of integers, where each is coefficient for indicator
         variable x_r where 1 <= r <= m^2*n^2 
    """
    m = len(tableau_a)
    n = len(tableau_b)
    r = 0
    for i in xrange(m):
        for j in xrange(n):
            for k in xrange(m):
                for l in xrange(n):
                    yield score(tableau_a[(i, k)], tableau_b[(j, l)])
                    assert (r == get_x_indicator_variable_index(i,j,k,l,m,n))
                    r += 1


def get_x_indicator_variable_index(i, j, k, l, m, n):
    """
    Map the i,j,k,l indices to the sequential indicator variable index
    as generated by linear_objective_function_coefficients().
    This is basically the (4-dimensional) 'array equation' (as per
    row-major arrays in C for example).

    Note that for MiniSat+, the variables are juist indexed sequentially
    (from 1), and we are mapping the x_{ijkl} to x_r for
    0 <= r < m^2*n^2 variables,
    This function gets the sequential index for an x_{ijkl} variable.

    Parameters:
        i, j, k, l - indices for indicator variables as per loop in
                     linear_objective_function_coefficients()
        m - order of tableau a (0 <= i,k < m)
        n - order of tableau b (0 <= j,l < n)

    Return value:
        index r of indicator variable x_{r} as generated by
        linear_objective_function_coefficients()
        
    """
    return  ((i*n + j)*m + k)*n + l 


def get_y_indicator_variable_index(i, j, m, n):
    """
    Map the i,j indices to the sequential indicator variable index
    for the y_{ij} variable.

    This is basically the (2-dimensional) 'array equation' (as per
    row-major arrays in C for example).

    Note that for MiniSat+, the variables are juist indexed sequentially
    and we are mapping the y_{ij} to y_r for 0 <= r < m*n variables.
    This function gets the sequential index for a y_{ij} variable.

    Parameters:
        i, j   - indices for y indicator variable
        m - order of tableau a (0 <= i,k < m)
        n - order of tableau b (0 <= j,l < n)

    Return value:
        index r of indicator variable y_{r} corresponding to y_{ij}
        
    """
    return i*n + j


def get_ij_from_index(r, m, n):
    """
    The inverse of get_y_indicator_variable_index(): given the indiator
    variable index, return the (i ,j) pair for y_{ij} to which it
    corresponds. So when we get a solution from MiniSat+, the variable
    index given to this function gets converted to (i,j) our y_{ij}, meaning
    that if that indiator variable is 1, SSE i in tableau a is matched with
    SSE j in tableau b.

    Parameters:
        index r of indicator variable x_{r} 
        m - order of tableau a (0 <= i,k < m)
        n - order of tableau b (0 <= j,l < n)

    Return value:
        i, j   - indices for y indicator variable corresponding to
                 MiniSat+ variable y_r
    """
    i = r / n
    j = r % n
    return (i,j)
     
    
def write_linear_objective_function(tableau_a, tableau_b, fh,
                                    format='minisat+'):
    """
    Build the integer line objective function, a summation over i,k in
    tab_a and j,l in tab_b of m^2*n^2 terms where m is order of tableau a
    and n is order of tableau b. Each term is the score function at (i,j,k,l)
    multiplied by the indicator variable x_{ijkl}, where x_{ijkl} (boolean)
    is 1 if SSE i in A is matched with SSE j in B AND SSE k in A is matched
    with SSE l in B, else 0.
    For MiniSat+, these indicator variables x are numbered 1..r where r=m^2*n^2.

    WARNING: Note the large (quartic) number of terms: eg this comes to 810 000
    terms if m = n = 30.

    Parameters:
        tableau_a - PTTableau object for protein A
        tableau_b - PTTableau object for protein B
        fh - open for write filehandle to write objective
             function to.
        format - 'minisat+' or 'cplex'. Default 'minisat+'.

    Return value:
        Index of last coefficient (starts at 1) = m^2*n^2
        which is index of last indicator variable used.

    Raises exceptions:
         ValueError on uknown format.
    """
    r = 0 
    for c in linear_objective_function_coefficients(tableau_a, tableau_b):
        if format == 'minisat+':
            fh.write("%+d" % c + "*" + "x" + str(r) + " ")
        elif format == 'cplex':
            fh.write("%+d" % c + "x" + str(r) + " ")
        else:
            raise ValueError('unknown format ' + format + '\n')
        r += 1
    if format == 'minisat+':
        fh.write(";\n")
    else:
        fh.write("\n")
    return r-1


def write_atmostone_constraints(m, n, fh, format='minisat+'):
    """
    Write the constraints to ensure that each SSE is matched with at most
    one SSE in the other tableau i.e.
    \sum {i=1}_m y_{ij} <= 1
    here the y_{ij} are for MiniSat+ new variables so we start them at y_r

    Parameters:
        m - order of tableau a (0 <= i,k < m)
        n - order of tableau b (0 <= j,l < n)
       fh - open for write filehandle to write MiniSat+ format constraints to
       format - 'minisat+' (deafult) or 'cplex'.
       
    Return value:
       None. fh is written to.

    Raises exceptions:
         ValueError on uknown format.
       
    """
#     for i in xrange(m):
#         for j in xrange(n):
#             fh.write("* +1*y" + str(i) + ',' + str(j) + ' ')
#         fh.write("<= +1;\n")
#     for j in xrange(n):
#         for i in xrange(m):
#             fh.write("* +1*y" + str(i) + ',' + str(j) + ' ')
#         fh.write("<= +1;\n")
                     
    for i in xrange(m):
        for j in xrange(n):
            if format == 'minisat+':
                fh.write("+1*y" + str(get_y_indicator_variable_index(i,j,m,n))+" ")
            elif format == 'cplex':
                fh.write("+1y" + str(get_y_indicator_variable_index(i,j,m,n))+" ")
            else:
                raise ValueError('unknown format ' + format + '\n')
        if format == 'minisat+':
            fh.write("<= +1;\n")
        else:
            fh.write("<= 1\n")
    for j in xrange(n):
        for i in xrange(m):
            if format == 'minisat+':
                fh.write("+1*y" + str(get_y_indicator_variable_index(i,j,m,n))+" ")
            elif format == 'cplex':
                fh.write("+1y" + str(get_y_indicator_variable_index(i,j,m,n))+" ")
            else:
                raise ValueError('unknown format ' + format + '\n')
        if format == 'minisat+':
            fh.write("<= +1;\n")
        else:
            fh.write("<= +1\n")


def write_conjunction_constraints(m, n, fh, format='minisat+'):
    """
    Write the constraints to ensure that the introduced indicator variables
    x_{ijkl} values do not exceed the values of the 'original' indicator
    variables y_{ij} from whch they were derived as the conjunction
    x_{ijkl} = y_{ij} /\ y_{kl}.
    These constraints are:
       x_{ijkl} <= y_{ij},    0 <= i,k < m, 0 <= j,l < n
       x_{ijkl} <= y_{kl},    0 <= i,k < m, 0 <= j,l < n

    re-written for pseudo-boolean MiniSat+ constraint format as:

       x_{ijkl} - y_{ij} <= 0,    0 <= i,k < m, 0 <= j,l < n
       x_{ijkl} - y_{kl} <= 0,    0 <= i,k < m, 0 <= j,l < n
    

    WARNING: Note the large (quartic) number of constraints.

    Parameters:
        m - order of tableau a (0 <= i,k < m)
        n - order of tableau b (0 <= j,l < n)
       fh - open for write filehandle to write constraints to
       format - 'minisat+' (default) or 'cplex'.
       
    Return value:
       None. fh is written to.

    Raises exceptions:
         ValueError on uknown format.
       
    """
#     for i in xrange(m):
#         for j in xrange(n):
#             for k in xrange(m):
#                 for l in xrange(n):
#                     fh.write("* +1*x" +
#                              str(i)+','+str(j)+','+str(k)+','+str(l) +
#                              " -1*y" +
#                              str(i)+','+str(j) +
#                              " <= 0;\n")
#                     fh.write("* +1*x" +
#                              str(i)+','+str(j)+','+str(k)+','+str(l) +
#                              " -1*y" +
#                              str(k)+','+str(l) +
#                              " <= 0;\n")
    
    for i in xrange(m):
        for j in xrange(n):
            for k in xrange(m):
                for l in xrange(n):
                    if format == 'minisat+':
                        fh.write("+1*x" +
                                 str(get_x_indicator_variable_index(i,j,k,l,m,n)) +
                                 " -1*y" +
                                 str(get_y_indicator_variable_index(i,j, m, n)) +
                                 " <= 0;\n")
                        fh.write("+1*x" +
                                 str(get_x_indicator_variable_index(i,j,k,l,m,n)) +
                                 " -1*y" +
                                 str(get_y_indicator_variable_index(k,l, m, n)) +
                                 " <= 0;\n")
                    elif format == 'cplex':
                        fh.write("+1x" +
                                 str(get_x_indicator_variable_index(i,j,k,l,m,n)) +
                                 " -1y" +
                                 str(get_y_indicator_variable_index(i,j, m, n)) +
                                 " <= 0\n")
                        fh.write("+1x" +
                                 str(get_x_indicator_variable_index(i,j,k,l,m,n)) +
                                 " -1y" +
                                 str(get_y_indicator_variable_index(k,l, m, n)) +
                                 " <= 0\n")

                    else:
                        raise ValueError('unknown format ' + format + '\n')


def write_explicit_xy_constraints(m, n, fh, format='minisat+'):
    """
    Write constraints to push x_ijkl constraint to 1 when y_jk and
    y_kl 1 (constraint (12)). These constraints are:

    y_{ij} + y_{kl} <= x_{ijkl} + 1,    0 <= i,k < m, 0 <= j,l < n

    re-written for pseudo-boolean MiniSat+ constraint format as:

    y_{ij} + y_{kl} - x_{ijkl} <= 1,    0 <= i,k < m, 0 <= j,l < n
    

    Parameters:
        m - order of tableau a (0 <= i,k < m)
        n - order of tableau b (0 <= j,l < n)
       fh - open for write filehandle to write MiniSat+ format constraints to
       format - 'minisat+' (default) or 'cplex'.
       
    Return value:
       None. fh is written to.

    Raises exceptions:
         ValueError on uknown format.
       
       
    """
#     for i in xrange(m):
#         for j in xrange(n):
#             for k in xrange(m):
#                 for l in xrange(n):
#                     fh.write("* +1*y" +
#                              str(i) + ',' + str(j) +
#                              " +1*y" +
#                              str(k) + ',' + str(l) +
#                              " -1*x" +
#                              str(i)+','+str(j)+','+str(k)+','+str(l) +
#                              " <= +1;\n")
    
    for i in xrange(m):
        for j in xrange(n):
            for k in xrange(m):
                for l in xrange(n):
                    if format == 'minisat+':
                        fh.write("+1*y" +
                                 str(get_y_indicator_variable_index(i, j, m, n)) +
                                 " +1*y" +
                                 str(get_y_indicator_variable_index(k, l, m, n)) +
                                 " -1*x" +
                                 str(get_x_indicator_variable_index(i,j,k,l,m,n)) +
                                 " <= +1;\n")
                    elif format == 'cplex':
                        fh.write("+1y" +
                                 str(get_y_indicator_variable_index(i, j, m, n)) +
                                 " +1y" +
                                 str(get_y_indicator_variable_index(k, l, m, n)) +
                                 " -1x" +
                                 str(get_x_indicator_variable_index(i,j,k,l,m,n)) +
                                 " <= +1\n")
                    else:
                        raise ValueError('unknown format ' + format + '\n')




def write_ordering_constraints(m, n, fh, format='minisat+'):
    """
    Write constraints to preserve order of SSEs (constraint (9)).
    These constraints are:

    y_{ij} + y_{kl} <= 1,     0 <= i < k < m, 0 <= l < j < n

    Parameters:
        m - order of tableau a (0 <= i,k < m)
        n - order of tableau b (0 <= j,l < n)
       fh - open for write filehandle to write MiniSat+ format constraints to
       format - 'minisat+' (default) or 'cplex'.
       
    Return value:
       None. fh is written to.

    Raises exceptions:
         ValueError on uknown format.
    """
    for i in xrange(m):
        for k in xrange(i+1, m):
            for l in xrange(n):
                for j in xrange(l+1, n):
                    if format == 'minisat+':
                        fh.write("+1*y" +
                                 str(get_y_indicator_variable_index(i, j, m, n)) +
                                 " +1*y" +
                                 str(get_y_indicator_variable_index(k, l, m, n)) +
                                 " <= 1;\n")
                    elif format == 'cplex':
                        fh.write("+1y" +
                                 str(get_y_indicator_variable_index(i, j, m, n)) +
                                 " +1y" +
                                 str(get_y_indicator_variable_index(k, l, m, n)) +
                                 " <= 1\n")
                    else:
                        raise ValueError('unknown format ' + format + '\n')

                            



def write_minisatplus_spec(tableau_a, tableau_b, fh,
                           use_explicit_xy_constraint,
                           use_ordering_constraint):
    """
    Write the MiniSat+ problem specification for pseudo-Boolean satisfaction
    formulation of tableau matching.

    Parameters:
        tableau_a - PTTableau object for protein A
        tableau_b - PTTableau object for protein B
        fh - open for write filehandle to write MiniSat+ problem spec to.
        use_explicity_xy_constraint - if True, include constraint to
                       push x_ijkl constraint to 1 when y_jk and y_kl 1
                       (constraint (12))
        use_ordering_constraint - if True, include constraint to preserve
                       order of SSEs.
    
    """
    m = len(tableau_a)
    n = len(tableau_b)
    fh.write("min: ")
    r = write_linear_objective_function(tableau_a, tableau_b, fh)
    assert(r == m**2 * n**2 - 1)
    write_atmostone_constraints(m, n, fh)
    write_conjunction_constraints(m, n, fh)
    if use_explicit_xy_constraint:
        write_explicit_xy_constraints(m, n, fh)
    if use_ordering_constraint:
        write_ordering_constraints(m, n, fh)


def write_cplex_spec(tableau_a, tableau_b, fh,
                     use_explicit_xy_constraint,
                     use_ordering_constraint):
    """
    Write the CPLEX  problem specification for integer linear program
    formulation of tableau matching.

    Parameters:
        tableau_a - PTTableau object for protein A
        tableau_b - PTTableau object for protein B
        fh - open for write filehandle to write MiniSat+ problem spec to.
        use_explicity_xy_constraint - if True, include constraint to
                       push x_ijkl constraint to 1 when y_jk and y_kl 1
                       (constraint (12))
        use_ordering_constraint - if True, include constraint to preserve
                       order of SSEs.
    
    """
    # TODO:
    # Currently we write a file of commands for the interactive optimizer,
    # should probably use the .lp file format instead (or do it properly
    # with one of the CPLEX APIs)
    m = len(tableau_a)
    n = len(tableau_b)
    fh.write("enter tabmatch\n")
    fh.write("minimize\n")
    r = write_linear_objective_function(tableau_a, tableau_b, fh,
                                        format='cplex')
    assert(r == m**2 * n**2 - 1)
    fh.write("st\n")
    write_atmostone_constraints(m, n, fh, format='cplex')
    write_conjunction_constraints(m, n, fh, format='cplex')
    if use_explicit_xy_constraint:
        write_explicit_xy_constraints(m, n, fh, format='cplex')
    if use_ordering_constraint:
        write_ordering_constraints(m, n, fh, format='cplex')
    fh.write('end\n')
    fh.write('change problem milp\n')
    for i in xrange(r):
        fh.write('change type x' + str(i) + ' b\n')
    for i in xrange(m*n):
        fh.write('change type y' + str(i) + ' b\n')
    fh.write('optimize\n')
    fh.write('display solution objective\n')
    fh.write('display solution variables -\n')



def parse_minisat_result(fh):
    """
    Parse the MiniSat+ output, or actually just the (extremely long!) line
    starting with 'v ' that lists all variables like x1 for true -x1 for false
    e.g.
    
    v x1 -x2 -x3 -x4 -x5 -x6 

    The result returned from here is a list of integers where the MiniSat+
    x variable is true, e.g. just [1] in the example above.
    
    Parameters:
       fh - open for reading filehandle of MiniSat+ output

    Return value:
       list of integers where corresponding variable is True by MiniSat+
    """
    truelist = []
    for line in fh:
        if line[:2] == "v ":
            values = line.split()
            for v in values[1:]:
                if v[0] == "y":  # The y indicator variables are what we want
                    index = int(v[1:])
                    truelist.append(index)
                elif v[0] == "x":
                    pass # not doing anything with the x indicator variables
                elif v[0] != "-":
                    sys.stderr.write("cannot parse value " + v + "\n")
        elif verbose:
            sys.stderr.write(line)
    return truelist

def convert_truelist_to_ij(truelist, m, n):
    """
    Given a list of MiniSat+ indicator variables that are assigned True in
    the solution, as retured by parse_minisat_result(), return a list
    of (i,j) tuples where each (i,j) is the y_{ij} index correspdonding
    to the MiniSat+ indicator variable index in the truelist.

    Parameters:
       truelist - list of MiniSat+ y indicator variable indices assigned True
        m - order of tableau a (0 <= i,k < m)
        n - order of tableau b (0 <= j,l < n)
       
    Return value:
       list of (i,j) tuples sa described above.
    """
    return [get_ij_from_index(r,m,n) for r in truelist]


def tabmatch_minisatplus(tableau_a, tableau_b, use_explicit_xy_constraint,
                         use_ordering_constraint):
    """
    Build the tableau matching problem as pseudo-boolean SAT in
    MinSat+ format, run MiniSat+ on it, and parse result.  The result
    is a list of (i,j) tuples where each (i,j) is the y_{ij} index
    correspdonding to the MiniSat+ indicator variable index in the
    truelist.

    Warning: the tempfile used here can be very large, due to very large
    number of constraints and variables introduced.
    
    Parameters:
        tableau_a - PTTableau object for protein A
        tableau_b - PTTableau object for protein B
        use_explicit_xy_constraint - if True, include constraint to
                  explicitly push x_{ijkl} to 1 when y_{ij} and y_{kl}
                  are 1.
        use_ordering_constraint - if True, include constraint to preserve
                  order of SSEs.
        
    Return value:
       list of (i,j) tuples sa described above.
    
    """
    m = len(tableau_a)
    n = len(tableau_b)
    # minisat+ cannot read from stdin, need a temp file
    tmpfile = os.tempnam(None, "msin")
    tmpfd = open(tmpfile, 'w')
    try:
        write_minisatplus_spec(tableau_a, tableau_b, tmpfd,
                               use_explicit_xy_constraint,
                               use_ordering_constraint)
        tmpfd.close()
        if verbose:
            sys.stderr.write("running minisat+...")
        fd = os.popen("minisat+ " + tmpfile)
        truelist = parse_minisat_result(fd)
        ijlist = convert_truelist_to_ij(truelist, m, n)
    finally:
        os.unlink(tmpfile)
    return ijlist
    
    
#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------


def usage(progname):
    """
    Print usage message and exit
    """
    sys.stderr.write("Usage: " + progname +
            " [-35cov] [-t struct_prog] "
            "[-p domain_prog] [-x sse_num_list] [-y sse_num_list] <PDBfile1> <PDBfile2>\n")
    sys.stderr.write("  -3 include 3_10 helices\n")
    sys.stderr.write("  -5 include pi helices\n")
    sys.stderr.write("  -p domain decomposition method/db\n"
                     "     valid values are none (default), "
                     "ddomain, cath:cdffile, pdomains:pdomainsfile\n")
    sys.stderr.write("  -t struct_prog : use struct_prog define " \
                     "secondary structure\n")
    sys.stderr.write("       supported is 'pdb' (default) or 'stride' or 'dssp'\n")
    sys.stderr.write("  -l use CPLEX rather than MiniSat+\n")
    sys.stderr.write("  -c include explicit constraint on introduced conjunction"
                      " indicator variables\n")
    sys.stderr.write("  -o include ordinering constraints\n")
    sys.stderr.write("  -s do not run MiniSat+, just write spec to stdout\n")
    sys.stderr.write("  -x sse_num_list : specifies comma-separated list of "
                     "SSE sequential numbers to include in the PDBFile1 tableau\n")
    sys.stderr.write("  -y sse_num_list : specifies comma-separated list of "
                     "SSE sequential numbers to include in the PDBFile2 tableau\n")
    sys.stderr.write("  -s sse_num_list : specifies comma-separated list of "
                     "SSE sequential numbers to include in the tableau\n")
    sys.stderr.write("  -v print verbose debugging messages to stderr\n")
    sys.exit(1)


def main():
    """
    main for tabmatchpbool.py

    Usage: tabmatchpbool [-35lcosv][-t structprog] [-p domainprog]
                         [-x sse_num_list] [-y sse_num_list]
                         <PDBfile1> <PDBfile2>


    -3 specifies to include 3_10 helices in the diagram. Default is only
       alpha helices.

    -5 specifies to include pi helices in the diagram. Defaul is only
       alpha helices.

    -p specify the domain decomposition method.
       Valid values are 'none' (default), 'ddomain', 'cath:cdf_filename'.

    -t specifies the secondary structure assignment program to use.
       Currently suppoed is 'pdb' and 'dfh,ssp' and 'stride'. Default 'pdb'.

    -l use ILP formulation for CPLEX (MILP) rather than pseudo-Boolean with
       MiniSat+
       
    -c include constraint explicitly push x_{ijkl} to 1 when y_{ij} and y_{kl}
       are 1.

    -o include constraints so that order of SSEs is preserved.

    -s do not run MiniSat+ or CPLEX, just write problem specification to stdout.

    -x sse_num_list specifies a comman-separated
       list of SSE sequential ids to build the
       tableau for PDBfile1. SSE sequential id's start at 1 and go from N to C
       terminus. E.g. -s1,5,8 includes only the 1st, 5th and 8ths SSEs.
       Numbers do not restart at chains (but do restart in each domain).
       These nubmers are those assigned by 'ptgraph2 -b sequential' option.

       TODO: this currently does not make sense when multiple domains
       are being procssed, this option applies to each domain.

    -y sse_num_list specifies a comman-separated
       list of SSE sequential ids to build the
       tableau for PDBfile2. SSE sequential id's start at 1 and go from N to C
       terminus. E.g. -s1,5,8 includes only the 1st, 5th and 8ths SSEs.
       Numbers do not restart at chains (but do restart in each domain).
       These nubmers are those assigned by 'ptgraph2 -b sequential' option.

       TODO: this currently does not make sense when multiple domains
       are being procssed, this option applies to each domain.

    -v specifies verbose mode: debugging output is written to stderr.
    """
    global verbose
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "35lcosp:t:s:x:y:v?")
    except getopt.GetoptError:
        usage(os.path.basename(sys.argv[0]))

    valid_secstruct_programs = ["dssp", "stride", "pdb"]
    valid_domain_programs = getdomains.valid_domain_programs + [r"none"]
    valid_domain_programs_re = [ re.compile(re_str) for re_str in
                                 valid_domain_programs ]

    verbose = False # global (python globals are only 'global' to module though)
    secstruct_program = "pdb"
    include_310_helices = False
    include_pi_helices = False
    domain_program = "none"
    sse_id_lists = [None, None]
    use_explicit_xy_constraint = False
    use_ordering_constraint = False
    run_solver = True
    use_cplex = False

    for opt,arg in opts:
        if opt == "-3":   # include 3_10 helices
            include_310_helices = True
        elif opt == "-5": # include pi helices
            include_pi_helices = True
        elif opt == "-p": # domain parsing program
            domain_program = None
            for valid_domarg_re in valid_domain_programs_re:
                if valid_domarg_re.match(arg):
                    domain_program = arg
                    break
            if domain_program == None:
                sys.stderr.write("valid values for -p are: " +
                                 str(valid_domain_programs) + "\n")
                usage(sys.argv[0])
        elif opt == "-t":
            if arg not in valid_secstruct_programs:
                sys.stderr.write("valid values for -t are: " +
                                 str(valid_secstruct_programs) + "\n")
                usage(sys.argv[0])
            secstruct_program = arg
        elif opt == "-l":   # use CPLEX not MiniSat+
            use_cplex = True
        elif opt == "-c":   # explicit conjunction constraint is to be use
            use_explicit_xy_constraint = True
        elif opt == "-o":   # use ordering constraint
            use_ordering_constraint = True
        elif opt == "-s":   # write problem formulation to stdout, don't run sat
            run_solver = False
        elif opt == "-x" or opt == "-y":
            sse_id_list_str = arg.split(',')
            sse_id_list = []
            sse_id_uniq_dict = {} # { id : True } just for checking all unique
            for sse_id_str in sse_id_list_str:
                if sse_id_str.isdigit():
                    if sse_id_uniq_dict.has_key(int(sse_id_str)):
                        sys.stderr.write("duplicate SSE sequential number "  +
                                         sse_id_str + "\n")
                        usage(sys.argv[0])
                    sse_id_uniq_dict[int(sse_id_str)] = True
                    sse_id_list.append(int(sse_id_str))
                else:
                    sys.stderr.write("not a valid SSE sequential number '" +
                                     sse_id_str + "'\n")
                    usage(sys.argv[0])
            sse_id_list.sort() # ensure SSEs are in order
            if opt == "-x":
                sse_id_lists[0] = list(sse_id_list)
            else:
                sse_id_lists[1] = list(sse_id_list)
        elif opt == "-v": # verbose
            verbose = True # this module only
            ptnode_set_verbose(True) # ptnode module
            ptsecstruct.ptsecstruct_set_verbose(True) # ptsecstruct module
            ptdomain_set_verbose(True) # ptdomain module
        else:
            usage(sys.argv[0])

    if len(args) != 2:
        usage(os.path.basename(sys.argv[0]))
    pdb_filenames = args

    tableaux = []
    for (pdb_filename, sse_id_list) in zip(pdb_filenames, sse_id_lists):
        # check for compressed files. We only support gzip (.gz)
        # Note we are not using the zlib or GzipFile python modules
        # since we are calling to external programs which require the
        # file uncompressed themsevles anyway so we'll just run gzip
        # to uncompress the file to a temporary directory.
        pdb_file_basename = os.path.basename(pdb_filename)
        (name,extension) = os.path.splitext(pdb_file_basename)
        if extension == '.gz':
            TMPDIR = os.tempnam(None, "ptgz")
            os.mkdir(TMPDIR)
            tmp_pdbfilename = os.path.join(TMPDIR, name)
            os.system("gzip " + pdb_filename + " -d -c > " + tmp_pdbfilename)
            our_pdb_filename = tmp_pdbfilename
            used_tmp_file = True
        else:
            our_pdb_filename = pdb_filename
            used_tmp_file = False

        try:
            pdbid = name.upper()
            if len(pdbid) >= 6 and pdbid[:3] == "PDB":
                pdbid = pdbid[3:7]
            # parse PDB file
            pdb_parser = PDBParser()
            pdb_struct = pdb_parser.get_structure(pdbid, our_pdb_filename)
            # create the Tableau
            tableau_list = get_tableau(our_pdb_filename,
                                       pdb_struct,
                                       secstruct_program,
                                       domain_program,
                                       include_310_helices,
                                       include_pi_helices,
                                       sse_id_list)
            if len(tableau_list) > 1:
                # TODO: handle multiple domains in at least one of the
                # PDB files. But for now OK to just use SCOP/ASTRAL inputs
                # for single domain or do no domain decomp for whole PDB
                # entry.
                sys.stderr.write("ERROR: Only single domains supported for now:"
                                 " retry without domain decomposition\n")
                sys.exit(1)
            else:
                tableaux.append(tableau_list[0])
        finally:
            if used_tmp_file:
                cleanup_tmpdir(TMPDIR)

    if verbose:
        sys.stderr.write('A:\n')
        sys.stderr.write(str(tableaux[0]))
        sys.stderr.write('\nB:\n')
        sys.stderr.write(str(tableaux[1]))


    if run_solver:
        if use_cplex:
            # TODO run CPLEX
            sys.stderr.write('running CPLEX not yet implemented\n')
        else:
            ijlist = tabmatch_minisatplus(tableaux[0], tableaux[1],
                                          use_explicit_xy_constraint,
                                          use_ordering_constraint)
            sys.stdout.write(str(ijlist))
            sys.stdout.write('\n')
    else:
        if use_cplex:
            write_cplex_spec(tableaux[0], tableaux[1], sys.stdout,
                               use_explicit_xy_constraint,
                               use_ordering_constraint)
        else:
            write_minisatplus_spec(tableaux[0], tableaux[1], sys.stdout,
                                   use_explicit_xy_constraint,
                                   use_ordering_constraint)
    
    
if __name__ == "__main__":
    warnings.filterwarnings('ignore', 'tempnam', RuntimeWarning) 
    main()
