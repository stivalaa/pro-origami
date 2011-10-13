###############################################################################
#
# pttableau.py - Object to represent protein tableaux and functions to
#                parse output of TableauCreator program into tableau object.
#
#
# File:    pttableau.py
# Author:  Alex Stivala
# Created: October 2007
#
#
# $Id$
#
###############################################################################
"""
Ths module contains routines to generate protein tableaux using the
axes fitted to SSEs by functions in the ptnode.py module, and the relative
angle calculate in that module.

This module also contains functions to parse the output of
Arun Konargurthu's TableauCreator program, which was used before
tableaux were re-implemented internally to this code.
NOTE: not yet published or available (as of October 2007).
IMPORTANT: this requires my modified versin of TableauCreator,
and also a patched version of Bio.PDB file PDBIO.py :
context diff for patching it is Bio.PDB.PDBIO.py.diff
This patch was made relative to BioPython release 1.43.

Tableaux are described by Kamat and Lesk 2007
'Contact Patterns Between Helices and Strands of Sheet Define Protein
 Folding Patterns' Proteins 66:869-876
and Lesk 2003 'From Electrons to Proteins and Back Again'
Int. J. Quant. Chem. 95:678-682
and Lesk 1995 'Systematic representation of folding patterns'
J. Mol. Graph. 13:159-164.

Two classes are provided, PTTableau and PTTableauPacked. The latter is
a more compact format based on LAPACK style symmetric matrix packed
array storage, useful for holding whole databse of tableaux in memory
(see buildtableauxdb.py) and dumping/loading it.  They can be used
interchangeably when it comes to getting and setting with []
(__getitem__ and __setitem__). (They should probably both inherit from
some tableau base class to make this explicity, but it doesn't really
matter with Python 'duck typing' of which this is a use (abuse?).
"""

import os,sys
from math import pi
import numpy.oldnumeric as Numeric
from Bio.PDB import *

from ptnode import *
from ptdomain import PTDomain
from ptsecstruct import pdb_chainid_to_stride_chainid
from ptutils import cleanup_tmpdir

#-----------------------------------------------------------------------------
#
# Module globals
#
#-----------------------------------------------------------------------------


# constants

TABLEAU_MIN_HELIX_LEN  = 4 # min length of a helix in TableauCreator is 4
TABLEAU_MIN_STRAND_LEN = 2 # min length of a strand in TableauCreator is 2

# global variables

verbose = False

#-----------------------------------------------------------------------------
#
# Class definitions 
#
#-----------------------------------------------------------------------------

class PTTableau:
    """
    The PTTableau class is a protein tableau, as per Kamat and Lesk 2007
    'Contact Patterns Between Helices and Strands of Sheet Define Protein
    Folding Patterns' Proteins 66:869-876.

    The tableau is a 2 dimensional symmetric matrix indexed by SSEs
    in the protein where each entry
    is a two character code representing the angle between those SSEs.
    (See paper(s) for details).

    We implement it as a mapping container, i.e. using __getitem__ and
    __setitem__ so that elements can bet get/set with dictionary/array
    type syntax e.g. tableau[(helix1,strand2)]. (index is a tuple of SSEs
    represented by PTNode objects)- NB must have in parens to ensure tuple.
    This is implemented with a standard dictionary object,
    and since it is symmetric, only one copy is stored, the one where
    i < j in (i,j) index; however either can be get/set, they are swapped
    internally if necessary. Accessing (i,i) returns
    'xa', 'xi', 'xg' for respectively alpha,pi,310 helices and
    'e ' for strand.


    """

    def __init__(self, nodelist):
        """
        Intialize a PTTableau with no tableau entries set yet.

        Parameters:
           nodelist - list of nodes that will be in the tableau, ins
                      residue sequence number order.
        """
        self.tabdict = {} # { (res1, res2) : code }; see class documentation\
        self.nodelist = nodelist # ptnodes with those not in tableau removed
        

    def __str__(self):
        """
        Return string representation of the tableau; we will write a full matrix
        just like TableauCreator does.
        """
        s = ""
        for sse1 in self.nodelist:
            for sse2 in self.nodelist:
                try:
                    s += self[(sse1, sse2)] + ' '
                except KeyError:
                    s += "?? "
            s += '\n'
        return s
    
        

    #
    # methods defined to implement container type
    #

    def __len__(self):
        """
        Return number of SSEs in the tableau

        Parameters: None
        Return value: Number of SSEs in nodelist for building tableau
        """
        return len(self.nodelist)


    def __getitem__(self, ssepair):
        """
        Return the entry in the tableau for the pair of SSEs
        (sse1, sse2) where sse1 and sse2 are PTNode objects;
        or if ssepair is (i,j) where i,j are integers, the corresponding
        tableau entry for (nodelist[i], nodelist[j]).

        Parameters:
           ssepair - tuple (sse1,sse2) (PTNode objects) or
                     tuple (i,j) (integers) to look up tableau entry for
        Return value:
           two character tableau string e.g. 'RD' or 'HH', or '  ' (2 spaces).
           On the main diagonal (self-orientation) since this has no menaing
           we return a (two-char) encoding of the SSE type instead:
           'xa', 'xi', 'xg' for respectively alpha,pi,310 helices and
           'e ' for strand.

        Raises Exceptions:
            TypeError if ssepair is not PTNode pair or int pair.
        """
        ssespec1 = ssepair[0]
        ssespec2 = ssepair[1]
        if isinstance(ssespec1, PTNode) and isinstance(ssespec2, PTNode):
            sse1 = ssespec1
            sse2 = ssespec2
        elif isinstance(ssespec1, int) and isinstance(ssespec2, int):
            sse1 = self.nodelist[ssespec1]
            sse2 = self.nodelist[ssespec2]
        else:
            raise TypeError("bad tuple type in PTTableau getitem")
        if sse1 == sse2:
            if isinstance(sse1, PTNodeHelix):
                if sse1.get_type() == "ALPHA":
                    return "xa"
                elif sse1.get_type() == "PI":
                    return "xi"
                elif sse1.get_type() == "310":
                    return "xg"
                else:
                    return "??" # should not happen
            elif isinstance(sse1, PTNodeStrand):
                return "e "
            else:
                return "??" # should not happen
        elif sse1 < sse2:
            ssepair = (sse1,sse2)
        else:
            ssepair = (sse2,sse1)
        return self.tabdict[ssepair]


    def __setitem__(self, ssepair, tabcode):
        """
        Set the entry in the tableau for the pair of SSEs (sse1,sse2)
        specified as the key (ssepair) parameter to the tabcode value.

        Parameters:
           ssepair - tuple (sse1,sse2) to set.
           tabccode - two character tableau string e.g. 'RD' or 'HH', or '  '.

        Return value: None

        Raises exceptions:
           TypeError  if tabcode is not a valid 2 char uppercase string or '  '
        """
        if len(tabcode) != 2 or not tabcode.isupper() and not tabcode.isspace():
            raise TypeError("bad tableau code '" + tabcode + "'\n")

        if (tabcode[0] not in ['L','R','P','O'] or   \
            tabcode[1] not in ['E','D','S','T']) and \
           tabcode != 'HH' and tabcode != 'KK':
            raise TypeError("bad tableau code '" + tabcode + "'\n")

        sse1 = ssepair[0]
        sse2 = ssepair[1]
        if sse1 == sse2:
            return
        elif sse1 < sse2:
            ssepair = (sse1, sse2)
        else:
            ssepair = (sse2, sse1)
        self.tabdict[ssepair] = tabcode
        

    # have not implemented: __delitem__, __iter__, __contains__


    # TODO: work out how to implement things like tab[2:] to get row 2,
    # ust like Numeric.array etc.
    def getrow(self, i):
        """
        Return a row of the tableau as a list of tableau codes.

        Parameters:
           i - row to get 0 <= i < len(self)

        Return value:
           list of two-character tableau codes for row i.
        """
        return [self[(i,j)] for j in xrange(len(self))]
            

class PTTableauPacked:
    """
    The PTTableauPacked class is a compact representation
    of a  protein tableau, as per Kamat and Lesk 2007
    'Contact Patterns Between Helices and Strands of Sheet Define Protein
    Folding Patterns' Proteins 66:869-876.

    The tableau is a 2 dimensional symmetric matrix indexed by SSEs
    in the protein where each entry
    is a two character code representing the angle between those SSEs.
    (See paper(s) for details).

    We implement it as a mapping container, i.e. using __getitem__ and
    __setitem__ so that elements can bet get/set with dictionary/array
    type syntax e.g. tableau[(1,2)]. (index is a pair of sequential
    SSE numbers, from 0 to n-1 where n is the order of tableau ie number
    of SSEs) - NB must have in parens to ensure tuple.

    This is the compact respresentation, storing tableau as simply
    a linear string of two-character tableau codes, in the same as as the
    LAPACK 'packed' format for triangular/symmetric arrays. i.e.
    each column of the matrix is stored in sequence.
    We could save even more space by using only 4 bits for each tableau
    code (since there are only 16 possible codes), but in Python
    it doesn't really make sense to try to be so efficient - but
    we are trying to save space to some degree so that the entire ASTRAL
    PDB non-redundant set or similar can be loaded as tableaux in
    memory.
    As it happens, strings in python don't even support item assignemnt,
    so we have to store it as a list anyway  i.e ['xa','OT',...]
    instead of 'xaOT...'

    Unlike PTTableau, this format contains just the tableau codes
    and diagonal SSE type entries, i.e. just character data. there
    are no PTNode object references or anything, so it is simple and
    quick to dump/load with Python pickle module (or similar) with
    no need to build all sorts of other objects (PTNode, Bio.PDB.Structure,
    etc.).

    Accessing (i,i) returns 'xa', 'xi', 'xg' for respectively
    alpha,pi,310 helices and 'e ' for strand.
    """

    def __init__(self, tableau):
        """
        Intialize a PTTableauPacked given an already built tableaux in
        the full PTTableau format.

        Parameters:
            tableau - an already built PTTableau object
        """
        #self.n = 3 
        #self.uplist = ['xa','OS','e ','OT','PE','xa']
        self.n = len(tableau)  # order of tableau (number of SSEs)
        self.uplist = []       # packed format of matrix upper triangle
                               # NB COLUMN-MAJOR (LAPACK style)
        for j in range(self.n):
            for i in range(j+1):
                try:
                    tabcode = tableau[(i,j)]
                except:
                    tabcode = '??'
                self.uplist.append(tabcode)
        assert(len(self.uplist) ==  self.n * (self.n + 1) / 2)

    def __str__(self):
        """
        Return string representation of the tableau; we will write a full matrix
        just like TableauCreator does.
        """
        s = ""
        for i in range(self.n):
            for j in range(self.n):
                s += self[(i,j)] + ' '
            s += '\n'
        return s
        

    #
    # methods defined to implement container type
    #

    def __len__(self):
        """
        Return number of SSEs respresented the tableau

        Parameters: None
        Return value: order of tableau
        """
        return self.n
    

    def __getitem__(self, ssepair):
        """
        Return the entry in the tableau for the pair of SSEs
        (i,j) where i,j are integers, 0 <= i,j < n.

        Parameters:
           ssepair - tuple (i,j) (integers) to look up tableau entry for
        Return value:
           two character tableau string e.g. 'RD' or 'HH', or '  ' (2 spaces).
           On the main diagonal (self-orientation) since this has no menaing
           we return a (two-char) encoding of the SSE type instead:
           'xa', 'xi', 'xg' for respectively alpha,pi,310 helices and
           'e ' for strand.

        Raises Exceptions:
            TypeError if ssepair is not PTNode pair or int pair.
        """
        i = ssepair[0] + 1
        j = ssepair[1] + 1 # more convenient to have 1 < i,j <= n internally
        if j < i:
            tmp = i
            i = j
            j = tmp
        r = i + j*(j-1)/2
        r -= 1             # back to zero-based for list indexing
        return self.uplist[r]


    def __setitem__(self, ssepair, tabcode):
        """
        Set the entry in the tableau for the pair of SSEs (sse1,sse2)
        specified as the key (ssepair) parameter to the tabcode value.

        Parameters:
           ssepair - tuple (i,j) to set, 0 <= i,j < n.
           tabccode - two character tableau string e.g. 'RD' or 'HH',
                      or 'xa','e ' etc. if i==j for SSE type on diagonal.

        Return value: None

        Raises exceptions:
           TypeError  if tabcode is not a valid 2 char uppercase string
                       or lowercase type code for diagonal.
        """
        i = ssepair[0]
        j = ssepair[1]
        if i == j:
            if tabcode not in ['xa','xi','xg','e ']:
                raise TypeError('bad tableau sse type code ' + tabcode + '\n')
        else:
            if len(tabcode) != 2 or not tabcode.isupper():
                raise TypeError("bad tableau code '" + tabcode + "'\n")
            
            if ( (tabcode[0] not in ['L','R','P','O'] or
                  tabcode[1] not in ['E','D','S','T']) and
                 tabcode != 'HH' and tabcode != 'KK' ):
                raise TypeError("bad tableau code '" + tabcode + "'\n")

            if j < i:
                tmp = i
                i = j
                j = tmp
        i += 1
        j += 1  # more convenient to have 1 < i,j <= n internally
        r = i + j*(j-1)/2  # location in packed rep assuming each entry len 1
        r -= 1             # back to zero-based for list indexing
        self.uplist[r] = tabcode
        
    # have not implemented: __delitem__, __iter__, __contains__
    
    # TODO: work out how to implement things like tab[2:] to get row 2,
    # ust like Numeric.array etc.
    def getrow(self, i):
        """
        Return a row of the tableau as a list of tableau codes.

        Parameters:
           i - row to get 0 <= i < len(self)

        Return value:
           list of two-character tableau codes for row i.
        """
        # TODO: we could do this more efficiently for packed tableau
        return [self[(i,j)] for j in xrange(len(self))]



#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

def angle_to_tabcode(omega):
    """
    Convert an angle (radians in (-pi, pi]) to a two-character tableau
    code (double quadrant encoding) as described in the papers cited
    at top of module.

    Parmaters:
       omega - relative angle to encode
    Return value:
       two-character tableua code (OS, PD, etc.)
    Raises exceptions:
       ValueError if angle is out of range 
    """
    if omega > -pi/4 and omega <= pi/4:
        tabcode = 'P'                           # parallel
    elif omega > pi/4 and omega <= 3*pi/4:
        tabcode = 'R'                           # crossing-right
    elif (omega > 3*pi/4 and omega <= pi) or omega > -pi and omega <= -3*pi/4:
        tabcode = 'O'                           # antiparallel (opposite)
    elif omega > -3*pi/4 and omega <= -pi/4:
        tabcode = 'L'                           # crossing-left
    else:
        raise ValueError('bad omega value ' + str(omega) + '\n')

    if omega > 0 and omega <= pi/2:
        tabcode += 'D'                          # dinner
    elif omega > pi/2 and omega <= pi:
        tabcode += 'T'                          # tea
    elif omega > -pi and omega <= -pi/2:
        tabcode += 'S'                          # supper
    elif omega > -pi/2 and omega <= 0:
        tabcode += 'E'                          # elevenses
    else:
        raise ValueError('bad omega value ' + str(omega) + '\n')

    return tabcode



def compute_tableau(ptnode_list, pdb_structure, use_hk=True):
    """
    Build a PTTableau object for the tableau by computing relative angles
    between all SSEs in the ptnode_list.

    Parameters:
        ptnode_list - list of PTNode objects (ie iterable of PTNode)
                         representing the SSEs (helices,strands) the
                         tabelau is for.
        pdb_structure - parsed Bio.PDB structure
        use_hk - If True, use the HH and KK codes for respectively
                 antiparallel and parallel strands of the same sheet.
                 Default True.

    Return value:
       PTTableau object with entry for each pair of SSEs.
    """
    tableau = PTTableau(ptnode_list)
    for i in range(len(ptnode_list)):
        for j in range(i+1, len(ptnode_list)):
            omega = ptnode_list[i].relative_angle(ptnode_list[j], pdb_structure)
            if omega != None:
                try:
                    tabcode = angle_to_tabcode(omega)
                except ValueError:
                    sys.stderr.write('WARNING: catch bad tableau angle, seting Parallel (%d,%d)\n' %(i,j))
                    tabcode = "PE" # NaN -> 0.0 -> parallel: should not happen but does e.g. d7pcka_ -35
                # set tabcode to HH for antiparallel strands and
                # KK for parallel strands
                if (use_hk and
                    isinstance(ptnode_list[i], PTNodeStrand) and
                    isinstance(ptnode_list[j], PTNodeStrand) and
                    (ptnode_list[i].get_sheet_id() != None and
                     ptnode_list[i].get_sheet_id() ==
                     ptnode_list[j].get_sheet_id())):
                    if tabcode[0] == 'O':
                        tableau[(ptnode_list[i], ptnode_list[j])] = 'HH'
                    elif tabcode[0] == 'P':
                        tableau[(ptnode_list[i], ptnode_list[j])] = 'KK'
                    else:
                        tableau[(ptnode_list[i], ptnode_list[j])] = tabcode
                else:
                    tableau[(ptnode_list[i], ptnode_list[j])] = tabcode

    if verbose:
        sys.stderr.write(str(tableau))

    return tableau


def compute_omega_matrix(ptnode_list, pdb_structure):
    """
    Return the omega (relative angle, in radians) matrix as a 2D Numeric.array
    by computing relative angles between all SSEs in the ptnode_list

    Parameters:
        ptnode_list - list of PTNode objects (ie iterable of PTNode)
                         representing the SSEs (helices,strands) the
                         tabelau is for.
        pdb_structure - parsed Bio.PDB structure

    Return value:
       Numeric.array square symmetric (order length of ptnode_list) where
       each entry is relative angle between SSEs in radians.
       Main diagonal entries set to 0.
    
    """
    n = len(ptnode_list)
    omega_array = Numeric.zeros((n, n), Numeric.Float)
    for i in range(n):
        for j in range(i+1, n):
            omega = ptnode_list[i].relative_angle(ptnode_list[j], pdb_structure)
            if omega == None:
                omega_array[i, j] = float('NaN')
            else:
                omega_array[i, j] = omega
            omega_array[j, i] = omega_array[i, j]

    # set the diagonal as follows:
    # 0.00 for strand
    # 1.00 for alpha helix
    # 2.00 for pi helix
    # 3.00 for 3_10 helix
    for i in range(n):
        if isinstance(ptnode_list[i], PTNodeHelix):
            if ptnode_list[i].get_type() == "ALPHA":
                v = 1.00
            elif ptnode_list[i].get_type() == "PI":
                v = 2.00
            elif ptnode_list[i].get_type() == "310":
                v = 3.00
            else:
                pass # should not happen
        elif isinstance(ptnode_list[i], PTNodeStrand):
            v = 0.00
        omega_array[i,i] = v

    return omega_array

    
#-----------------------------------------------------------------------------
#
# Classes and functions for running external TableauCreator
#
#-----------------------------------------------------------------------------


# Inherit from the PDBIO.Select class for writing only parts of PDB to file
# See the Bio.PDB documentation: biopython-1.43/Doc/biopdb_faq.pdf
class DomainSelect(Select):
    """
    The DomainSelect class inherits from the PDBIO.Select class
    and overrides function to select only certain residues for writing
    ATOM records in the domain we are interested in to the
    simplified PDB file for TableauCreator.

    See the Bio.PDB documentation by Thomas Hamelryck:
      biopython-1.43/Doc/biopdb_faq.pdf
    """
    def __init__(self, domain):
        """
        Constructor for the DomainSelect class, sets the domain member
        used to accept only residues in that domain.
        Parameters:
           domain - ptdomain object of domain to select residues from
        """
        self.domain = domain

    def __repr__(self):
        """
        Overrides the base __repr__ to write out the domain we have
        """
        return "<DomainSelect: " + str(self.domain) + ">"
    
    def accept_residue(self, residue):
        """
        overrides the base accept_residue() function to accept only
        residues in our domain. Also reject HETATMS.
        Paramteters:
           residue - Bio.PDB Residue object of residue to test
        Return value:
           1 to accept residue, 0 to reject.
        """
        chain = residue.get_parent()
        chainid = pdb_chainid_to_stride_chainid(chain.get_id())
        # id of a residue in Bio.PDB is tuple (hetatm, resseqnum, icode)
        # so we choose those where chain and residue PDB number
        # is in the domain.
        resnum = residue.get_id()[1]
        if (self.domain.is_in_domain(chainid, resnum) and
            residue.get_id()[0] == ' '):
            return 1
        else:
            return 0



def get_tableau_from_pdbstruct(pdbid, domain,
                               pdb_structure, ptnode_list):
    """
    Build a PTTableau object for the tableau by first creating a
    simple PDB file with only the ATOM records for residues in the
    domain we are processing, and also a .SSEsInfo file containing the
    secnodary structure assignments we already have, then running
    TableauCreator on it (using our simple PDB file and SSEsInfo) and
    parsing the output.

    Parameters:
        pdbid - PDB identifier of the strucutre
        domain - The PTDomain object for our current domain
        pdb_structure - parsed Bio.PDB structure
        ptnode_list - list of PTNode objects (ie iterable of PTNode)
                         representing the SSEs (helices,strands) the
                         tabelau is for.
    Return value:
       PTTableau object built from TableauCreator output

    """
    TMPDIR = os.tempnam(None, "pttabin")
    os.mkdir(TMPDIR)
    try:
        filename = pdbid
        if domain.domainid != None:
            filename += '-' + domain.domainid
        filename += '.pdb'
        domain_pdb_filename = os.path.join(TMPDIR, filename)
        io = PDBIO()
        io.set_structure(pdb_structure)
        io.save(domain_pdb_filename, DomainSelect(domain))

        ssesinfo_filename = os.path.join(TMPDIR, filename + ".input-SSEsInfo")
        write_ssesinfo(ssesinfo_filename, ptnode_list)

        tableau =  read_tableau_from_tableaucreator(domain_pdb_filename,
                                                    ptnode_list,
                                                    ssesinfo_filename)
        os.unlink(domain_pdb_filename)
        os.unlink(ssesinfo_filename)
    finally:
        cleanup_tmpdir(TMPDIR)
    return tableau

    
def read_tableau_from_tableaucreator(pdb_filename, ptnode_list,
                                     ssesinfo_filename):
    """
    Run Arun's TableauCreator program on the supplied pdb_filename
    using SSEsInfo file.

    Parameters:
       pdb_filename - PDB file to run TableauCreator on
       ptnode_list - list of PTNode objects (ie iterable of PTNode)
                         representing the SSEs (helices,strands) the
                         tabelau is for.
       ssesinfo_filename - filename of the .SSEsInfo file that was written
                           to define SSEs for TableauCreator.
    Return value:
       PTTableau object built from TableauCreator output
       
    NB: TableauCreator is not yet published or available (October 2007)
    and I am using a private version which Arun sent me, which I modified
    to add the -s option to use STRIDE rather than DSSP
    and to have the -i option to parse .SSEsInfo files.
    """

    # TableauCreator needs an output directory where it writes all its
    # intermediate/output files, only puts progress information/errors
    # to stdout/stderr.

    tmpdir = os.tempnam(None, "pttab")
    os.mkdir(tmpdir)
    command = "TableauCreator "
    command += "-i " + ssesinfo_filename + " "
    command += pdb_filename + " " + tmpdir
    command += " >/dev/null"
    if verbose:
        sys.stderr.write("running '" + command + "'...")
    os.system(command)
    if verbose:
        sys.stderr.write("done\n")
    # output files are:
    #   <pdbfilename>.angles
    #   <pdbfilename>.SSEsInfo
    #   <pdbfilename>.stride or <pdbfilename>.dssp
    #   <pdbfilename>.tableau
    outfile_prefix = os.path.join(tmpdir, os.path.basename(pdb_filename))
    if not os.path.isfile(os.path.join(tmpdir, "TABCREATE_OK")):
        sys.stderr.write("ERROR: TableauCreator failed\n")
        cleanup_tmpdir(tmpdir)
        return None
    # Now the tricky thing is TableauCreator indexes its matrix just with
    # purely sequential numbers from 0 (as conventional)
    # assuming all SSEs in one domain and in fact one chain
    # (so we handle this by creating our own simple PDB file with only
    # ATOM records for our current domain, and only one TER record on
    # end so chains concatenated effectively).
    # And also (as in comments above functions) we have the dodginess of
    # doing the same thing in different ways in multiple places (DSSP/STRIDE
    # parsing, PDB parsing, etc.).
    # So let's check that the TableauCreator SSE info lines up with ours
    # (otherwise we can't use the tableau data).

    # parse the SSEsInfo file and check lines up with ptnodes,
    # returns list of ptnodes corresponding to Tableau entries (may be shorter
    # than our input node list; some removed as no equivalent in tableua).
    nodelist = parse_tableaucreator_ssesinfo(outfile_prefix + '.SSEsInfo',
                                             ptnode_list)
    if nodelist != None:
        tableau_filename = outfile_prefix + ".tableau"
        tableau = parse_tableaucreator_output(tableau_filename, nodelist)
        if tableau != None:
            if verbose:
                sys.stderr.write(str(tableau))
        else:
            sys.stderr.write('WARNING: problem parsing TableauCreator output;\n'
                         '         tableau information will not be used\n')
    else:
        sys.stderr.write('WARNING: problem with TableauCreator output;\n'
                         '         tableau information will not be used\n')
        tableau = None

    cleanup_tmpdir(tmpdir)
    return tableau
 

def parse_tableaucreator_ssesinfo(filename, nodelist):
    """
    Parse the .SSEsInfo file created by TableauCreator and check that
    it lines up with our SSE info in the form of the list of helix/strand
    PTNodes.

    Parameters:
      filename - filename of the .SSEsInfo file
      nodelist - list of PTNodes, in order of residue sequence number

    Return value:
      nodelist where nodes with no tableau entry removed (too short) if
      OK (they line up) else None (different number of nodes,
      residue sequence numbers/types don't match, etc).
    """

    # first remove all nodes with len < TABLEAU_MIN_SSE_LEN, since Tableau
    # Creator won't have entries for them; we will have to set them
#     ptnodelist = [ node for node in nodelist \
#                    if ( (isinstance(node, PTNodeHelix) and
#                          node.get_span() >= TABLEAU_MIN_HELIX_LEN) or
#                         (isinstance(node, PTNodeStrand) and
#                          node.get_span() >= TABLEAU_MIN_STRAND_LEN) )
#                  ]

    ptnodelist = nodelist 
    # FIXME: no longer need this filtering of too short SSEs,
    # now that .SSEsInfo input is being used
    if len(ptnodelist) != len(nodelist):
        sys.stderr.write('WARNING: no tableau entry for '
                         + str(len(nodelist)-len(ptnodelist)) +
                         ' nodes due to length to small' 
                         + '\n')
    
    fh = open(filename)
    # first line is number of SSEs, subsequent lines are 
    # type start_resnum end_resnum
    # where type is E or H (DSSP code) and resnums are PDB residue numbers
    # 
    numlines = int(fh.readline())
    if numlines < 2:
        sys.stderr.write('ERROR: bad SSEsInfo data\n')
        fh.close()
        return None
    linenum = 1
    sseinfo = [] # tuple (type, start, end)
    line = fh.readline()
    while line != "":
        fields = line.split()
        sseinfo.append((fields[0], int(fields[1]), int(fields[2])))
        linenum += 1
        line = fh.readline()
    fh.close()
    if len(sseinfo) != numlines:
        sys.stderr.write('ERROR: TableauCreator SSEsInfo file specified ' \
                         + str(numlines) + ' entries but ' \
                         + str(len(sseinfo)) + ' read\n')
        return None
    if len(sseinfo) != len(ptnodelist):
        sys.stderr.write('ERROR: TableauCreator SSEsInfo file has ' \
                         + str(numlines) + ' entries but ' \
                         + 'we have ' + str(len(ptnodelist)) + ' SSEs\n')
        return None
    for i in range(len(sseinfo)):
        sse = sseinfo[i]
        ptnode = ptnodelist[i]
        if sse[0] == 'H' and not isinstance(ptnode, PTNodeHelix) or \
           sse[0] == 'E' and not isinstance(ptnode, PTNodeStrand) or \
           sse[1] != ptnode.get_start_res_seq() or \
           sse[2] != ptnode.get_end_res_seq():
            sys.stderr.write('ERROR: TableauCreator SSEInfo entry ' + \
                             str(sse) + ' does not match node ' +
                             str(ptnode) + '\n')
            return None
    return ptnodelist


def parse_tableaucreator_output(filename, nodelist):
    """
    Parse the .tableau file created by TableauCreator

    Parameters:
       filename - filename of the .tableau file
       nodelist - list of PTNodes, in order of residue sequence number

    Return value:
      PTTableau for the tableau parsed
    """
    tableau = PTTableau(nodelist)
    # First line of file is number of SSEs (order of square matrix)
    # The whole matrix is in the file (it is symmetric), diagonal elements
    # and other unset (ie non-contact) elements are set to '--', or '  '
    # (two spaces) when cannot be calculated.
    fh = open(filename)
    numlines = int(fh.readline())
    if numlines < 2:
        sys.stderr.write('ERROR: bad tableau data\n')
        fh.close()
        return None
    linenum = 1
    line = fh.readline()
    i = 0
    while line != "":
        # we will use the fact that fields are fixed length rather than
        # splitting on space separator as fields may be set to '  ' (2 spaces)
        # when error calculating (helices too short etc.) in TableauCreator.
        # fields are two chars with two spaces between each
        node_i = nodelist[i]
        if len(line) < len(nodelist) * 4:
            sys.stderr.write('ERROR: bad line in tableau; line too short:\n')
            sys.stderr.write(line)
            fh.close()
            return None
        for j in range(i+1, len(nodelist)): # no need to store both i,j and j,i
            node_j = nodelist[j]
            col = 4*j  # each field as 2 char tabcode then 2 spaces
            tabcode = line[col:col+2]
            if tabcode !=  '  ' and tabcode != '--':
                tableau[(node_i, node_j)] = line[col:col+2]
        i += 1
        linenum += 1
        line = fh.readline()
    fh.close()
    return tableau
        


def write_ssesinfo(filename, nodelist):
    """
    Write a TableauCreator .SSEsInfo file describing the SSE asignments
    we have to the specified filename. This is so that we avoid having
    TableauCreator re-run DSSP or STRIDE for the assignments, which is
    inefficient and leads to inconsistencies. TableauCreator has been
    modified to be able to read this .SSEsInfo file instead, allowing
    the same assignments we have to be re-used by TableauCreator.

    WARNING: file is overwritten if it exists
    
    Parameters:
       filenanme - filename to write SSEsInfo to
       nodelist - list of PTNodes defining the SSEs

    Return value: None
    """
    # The format of the .SSEsInfo file is (see writeTableauAnglesSSEinfo())
    # that the first line has number of records and each subsequent
    # line (record)
    # is whitespace-separated:
    # dssp-code start end chainid
    # e.g.
    # H   10  21  A
    # Only H and E codes are used.
    # blank chainid is not allowed, '-' used instead.
    fh = open(filename, 'w')
    fh.write(str(len(nodelist)) + "\n")
    for node in nodelist:
        if isinstance(node, PTNodeHelix):
            typecode = 'H'
        elif isinstance(node, PTNodeStrand):
            typecode = 'E'
        else:
            assert(False)
        fh.write(typecode + " " + str(node.get_start_res_seq()) + " " +
                 str(node.get_end_res_seq()) + " " + node.get_chainid() + "\n")
    fh.close()

    
def pttableau_set_verbose(verb):
    """
    set the module global verbose flag in this module to supplied value
    Parameters: verb - True (for verbose output) or False
    Return value: None
    Uses globals: verbose (in this module)
    """
    global verbose
    verbose = verb
    
