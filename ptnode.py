###############################################################################
#
# ptnode.py - node (strand, helix, terminus) classes for ptgraph.
#
# File:    ptnode.py
# Author:  Alex Stivala
# Created: July 2007
#
# $Id$
#
# A PTnode is a node in the protein topology graph. It represents a
# secondary structure which is an alpha helix or beta sheet.
#
# Since these are simple (linear) secondary structures, the
# PTnode will have two
# edges, one from the previous (towards N-terminal) PTnode and one to
# the next (towards C-terminal) PTnode. So we just keep these nodes
# in a list, the ordering of which is sufficient for this purpose.
#
# Specific PTNode types (Helix, Strand, Terminus) are derived from PTNode.
#
#
###############################################################################

import sys
from math import pi,atan2,acos

from numpy.oldnumeric import array
from numpy.oldnumeric.linear_algebra import singular_value_decomposition  # part of Numeric
from Bio.PDB import *

from ptsecstruct import stride_chainid_to_pdb_chainid
from geometry import *
from ptmfile import mfile_write_strand, mfile_write_helix
from ptutils import get_int_icode,biopdbresid_to_pdbresseq,char_if_not_blank

#-----------------------------------------------------------------------------
#
# Module constants
#
#-----------------------------------------------------------------------------

ALPHA = 100    # multiplier of dircos to get second point on line
EPSILON = 1e-4 # epsilon for testing closeness to 0 or pi in fit_axis

#-----------------------------------------------------------------------------
#
# Module globals
#
#-----------------------------------------------------------------------------

verbose = False
global issued_warning
issued_warning = {} # warning for (chainid,res_seq_num) not found

#-----------------------------------------------------------------------------
#
# Class definitions
#
#-----------------------------------------------------------------------------

class PTNode:
    """
    A PTnode is a node in the protein topology graph. It represents a
    secondary structure which is an alpha helix or beta strand.

    Since these are simple (linear) secnodary structures, the
    PTnode will have two
    edges, one from the previous (towards N-terminal) PTnode and one to
    the next (towards C-terminal) PTnode. So we just keep these nodes
    in a list, the ordering of which is sufficient for this purpose.

    Specific PTNode types (Helix, Strand, Terminus) are derived from this
    class.

    The 'rich comparison' operators (<, <=, <, etc.) are overridden
    so that they are integer comparison on the start residue sequence number.
    Note equality and inequality are not overriden.

    """
    def __init__(self, nodeid, seqnum, start_res_seq, end_res_seq, chainid,
                 domainid,
                 pdb_residue_list, pdb_resid_dict,
                 fake_resids = False):

        """
        Construct a PTNode with supplied nodeid and type, and empty
        hydrogen bond list.

        Parameters:
           nodeid   -  unique identifier for this node (string)
           seqnum   -  sequence number of this node, note need not be unique
           start_res_seq - smallest PDB residue sequence number in the structure
           end_res_seq - largest PDB residue sequence number in the structure
             Note these residue sequence 'numbers' are strings, may have
             insertion codes.
          chainid = chain identifier for this SSE
          domainid = domain identifier for domain this SSE is in
          pdb_residue_list - list of all residues (for all chains) in the protein
          pdb_resid_dict -  dict of { {chainid,pdb_resseq) : seqindx }
                               where chainid and pdb_resseq make up
                               the PDB residue identifier, the pdb_resseq
                               being string resnum+icode if any e.g.
                               '60' or '60A', seqindx is the indiex
                               into sequential list of all residues
                                pdb_residue_list.
        fake_resids - Default False. IF True, do not look up resids in
                      pdb_resid_dict as they do not really exist; used for
                      terminus nodes.

          
        Exceptions:
           raises ValueEror if end_res_seq < start_res_seq

        """
        
        if ( not fake_resids ):
            try:
                start_res_indx = pdb_resid_dict[(chainid, start_res_seq)]
            except KeyError:
                if not issued_warning.has_key((chainid, start_res_seq)):
                    sys.stderr.write('WARNING: residue ' + start_res_seq +
                                     ' (chain ' +
                                     chainid +
                                     ') not found. May be HETATM.\n')
                    issued_warning[(chainid,start_res_seq)] = True
                while not pdb_resid_dict.has_key((chainid, start_res_seq)):
                    start_res_seq = str(get_int_icode(start_res_seq)[0] + 1)
                start_res_indx = pdb_resid_dict[(chainid, start_res_seq) ]
            try:
                end_res_indx = pdb_resid_dict[(chainid, end_res_seq)]
            except KeyError:
                if not issued_warning.has_key((chainid, end_res_seq)):
                    sys.stderr.write('WARNING: residue ' + end_res_seq +
                                     ' (chain ' +
                                     chainid +
                                     ') not found. May be HETATM.\n')
                    issued_warning[(chainid,end_res_seq)] = True
                while not pdb_resid_dict.has_key((chainid, end_res_seq)):
                    end_res_seq = str(get_int_icode(end_res_seq)[0] - 1)
                end_res_indx = pdb_resid_dict[(chainid, end_res_seq)]
            
            if ( end_res_indx < start_res_indx ):
                raise ValueError("end residue seqnum " + str(end_res_indx) + " < " +
                                 "start residue seqnum " + str(start_res_indx))
        
        self.nodeid = nodeid
        self.seqnum = seqnum # assigned sequentially, maybe not unique,
                             # used so strands can have a sequence number
                             # and helices also (both start at 1)
        self.start_res_seq = start_res_seq # start PDB residue seq num
        self.end_res_seq = end_res_seq # end PDB residue seq num
        self.chainid = chainid  # id of the chain this node is in
        self.pdb_residue_list = pdb_residue_list
        self.pdb_resid_dict = pdb_resid_dict
        self.hydrogen_bond_list = []   # list of (PTnode, r1, r2, dist) tuples
        self.sideways = False  # True if element to be
                               # drawn left-right not up-down
                               # use get/set_sideways()

        self.reversed = False # True if strand is to be drawn in reverse
                              # direction. This is set based on the
                              # parallel/antiparaell relationshiup to
                              # neighbour strand (so once the first
                              # in a sheet is set to False, all the others
                              # are set based on this and following the
                              # bridge relationships and their 'N' or 'P' flag).
                              # Also used for helices, where the tableau
                              # information is used i.e. orienation as
                              # parallel or antiparallel.
                              # use get/set_reversed()

        self.residue_list = [] # list of Bio.PDB Residue objects in this SSE.
                               # Built and returned by get_residue_list()
        self.resname_list = [] # list of 3 letter residue names in this SSE
                               # built by build_resname_sequence()
        self.resid_list = []   # list of string PDB residue sequence numbers
                               # built by build_resname_sequence()

        self.residue_ordinal_dict = None
                    #  dict { pdb_resseq : ordinal } mapping
                    #   the PDB sequence number string to integer orinal position
                    #   in sequence. Used to memoize this funtino so fast time it
                    #   is called the dict is built, subsequent calls look up in
                    #   dictionary.
                    # built get by get_residue_ordinal()

        self.fake_resids = fake_resids

        self.seq_component_id = None # intger sequence connect componented
                                     # id used for the 'fold' color scheme:
                                     # all strands in a sheet that are
                                     # connected in sequence with no
                                     # other SSEs in between have
                                     # the same id
                                     # similar for helices in clusters if
                                     # used.
                                     # use get/set_seq_component_id()

        self.domain_id = domainid  # domain identifier for domain this SSE is in
            

        # The following are only used by domainfc.py not ptgraph2.py:

        self.closenode_list=[]# list of (node, dist) tuples for
                              # PTNodes within some distance threshold
                              # of this node
        self.endchain = False # True if first or last in chain.
                              # use get/set_endchain()


        
        

    def __str__(self):
        """
        Return String representation of the node as 'TYPE id'
        """
        return "PTNode" + " " + self.nodeid

    def __lt__(self, other):
        """
        Compare based on chainid and start residue sequence number
        """
        assert(isinstance(other, PTNode))
        return (self.pdb_resid_dict[(self.chainid, self.start_res_seq)] <
                self.pdb_resid_dict[(other.chainid,other.start_res_seq)])

    def __le__(self, other):
        """
        Compare based on chainid and start residue sequence number
        """
        assert(isinstance(other, PTNode))
        return (self.pdb_resid_dict[(self.chainid, self.start_res_seq)] <=
                self.pdb_resid_dict[(other.chainid,other.start_res_seq)])


    def __gt__(self, other):
        """
        Compare based on chainid and start residue sequence number
        """
        assert(isinstance(other, PTNode))
        return (self.pdb_resid_dict[(self.chainid, self.start_res_seq)] >
                self.pdb_resid_dict[(other.chainid,other.start_res_seq)])

    def __ge__(self, other):
        """
        Compare based on chainid and start residue sequence number
        """
        assert(isinstance(other, PTNode))
        return (self.pdb_resid_dict[(self.chainid, self.start_res_seq)] >=
                self.pdb_resid_dict[(other.chainid,other.start_res_seq)])

    def get_chainid(self):
        """
        Return the chain identifier in this structure node
        Parameters: None
        Uses membe data (readonly):
           chainid
        Return value:
           The chain identifier of this structure node
        """
        return self.chainid

    def get_start_res_seq(self):
        """
        Return the lowest residue sequence number in this structure node

        Parameters: None

        Uses member data: (readonly)
           start_res_seq

        Return value:
           The residue sequence number at the start of this structure
           NB this is a pdb residue sequence nubmer so it is a string
           and may have an insertion code
        """
        return self.start_res_seq


    def get_end_res_seq(self):
        """
        Return the highest residue sequence number in this structure node

        Parameters: None

        Uses member data: (readonly)
           end_res_seq

        Return value:
           The residue sequence number at the end of this structure
           NB this is a pdb residue sequence number so it is a astring and
           may have an insertion code.
        """
        return self.end_res_seq


    def get_span(self):
        """
        Return the length in residues spanned by this structure node

        Parameters: None

        Uses member data: (readonly)
           end_res_seq, start_res_seq

        Return value: length in resisudes spanned by this node, ie the
             number of residues in this node.
        """
        return len(self.get_residue_list())

    def is_in_interval(self, res_seq_num):
        """
        Return True iff the supplied residue sequence number is contained
        in the interval spanned by this PTNode.

        Parameters:
            res_seq_num - PDB resisude sequence number to test (assumed in
                           this chain)

        Uses member data (readonly):
            chainid,start_res_seq, end_res_seq - start/end res seq num and
                              chainid of this node
            pdb_resid_dict - dict mapping chainid,resnum+icode to sequential
                            index

        Return value:
            True if res_seq_num is >= start_res_seq and <= end seq num
            else False
        """
        if self.fake_resids:
            return False
        try:
            res_indx = self.pdb_resid_dict[(self.chainid, res_seq_num)]
        except KeyError:
            if not issued_warning.has_key((self.chainid, res_seq_num)):
                sys.stderr.write('WARNING: residue ' + res_seq_num + ' (chain ' +
                                 self.chainid + ') not found. May be HETATM.\n')
                issued_warning[(self.chainid,res_seq_num)] = True
            return False
        return \
           ( res_indx >= self.pdb_resid_dict[(self.chainid,self.start_res_seq)]
             and
             res_indx <= self.pdb_resid_dict[(self.chainid,self.end_res_seq)] )


    def add_hbond(self, other_node, resnum1, resnum2, dist):
        """
        Add a hydrogen bond to the list of hydrogen bonds in this node.
        The bond is from PDB resdiue number resnum1 (in this node) to
        PDB residue number resnum2 in other_node with distance dist (Angstroms).

        Parameters:
           other_node - PTNode the bond is to from this node
           resnum1    - PDB residue number, must be in this node
           resnum2    - PDB residue number, must be in other_node
           dist       - (float) N..O distance of the bond (Angstroms)

        Uses data members (write):
           hydrogen_bond_list - list of (ptnode, resnum1, resnum2, dist) tuples

        Return value: None
        """
        assert(isinstance(other_node, PTNode))
        self.hydrogen_bond_list.append((other_node, resnum1, resnum2, dist))


    def get_hbond_list(self):
        """
        Return list of hydrogen bonds from this node in form of list of
        (ptnode, resnum1, resnum2, dist) tuples.

        Parameters: None
        Uses member data (readonly):
           hydrogen_bond_list - list of (other_node, resnum1, resnum2, dist)
                                tuples
        Return value:
           list of (ptnode, resnum1, resnum2,  dist) tuples.
           The actual list (ref)
           in this node, not a new copy.
        """
        return self.hydrogen_bond_list


    def get_residue_ordinal(self, pdb_resseq):
        """
        Given a PDB residue sequence number (may have insertino code) string
        for a residue in this SSE, return its ordinal position int the
        sequence of residues in this SSE (starting at 1).

        Parameters:
          pdb_resseq - string PDB residue sequence number e.g. '60' or '60A'
        Return value:
          integer >=1 ordinal number of the resiude in sequnce of residues
           for this SSE (from N to C terminal end)
        Uses data members (read/write):
           residue_ordinal_dict - dict { pdb_resseq : ordinal } mapping
              the PDB sequence number string to integer orinal position
              in sequence. Used to memoize this funtino so fast time it
              is called the dict is built, subsequent calls look up in
              dictionary.
        """
        if self.residue_ordinal_dict:
            return self.residue_ordinal_dict[pdb_resseq]
        self.residue_ordinal_dict = {}
        ordinal = 1
        for residue in self.get_residue_list():
            res_seq = biopdbresid_to_pdbresseq(residue.get_id())
            self.residue_ordinal_dict[res_seq] = ordinal
            ordinal += 1
        return self.residue_ordinal_dict[pdb_resseq]
        
        
    def get_residue_list(self):
        """
        Return the list of Bio.PDB Residue objects for the residues in this
        PTNode (strand or helix).
        Also store list in residue_list data member, and return stored
        value if already there.

        Parameters:
           None
        Uses data members:
           residue_list (read/write).
           pdb_residue_list, pdb_resid_dict (readonly)
           start_res_seq,end_res_seq,chainid (readonly)

        """
        if self.residue_list: # memoization: return if already computed
            return self.residue_list

        start_indx = self.pdb_resid_dict[(self.chainid,self.start_res_seq)]
        end_indx = self.pdb_resid_dict[(self.chainid,self.end_res_seq)]
        self.residue_list = self.pdb_residue_list[start_indx : end_indx + 1]
        return self.residue_list


    def set_sideways(self, sideways):
        """
        Label this strand with the sideways flag. See comments in __init__

        Parameters: sideways - True/False sideways flag
        Return value: None
        Uses data members (write): sideways
        """
        self.sideways = sideways

    def get_sideways(self):
        """
        Get the value of the sideways flag. See comments in __init__

        Parameters: None.
        Return value: True/False value of sideways flag
        Uses data members (readnly): sideways
        """
        return self.sideways


    def set_reversed(self, reversed):
        """
        Label this node with the reversed flag. See comments in __init__

        Parameters: reversed - True/False reversed flag
        Return value: None
        Uses data members (write): reversed
        """
        self.reversed = reversed

    def get_reversed(self):
        """
        Get the value of the reversed flag. See comments in __init__

        Parameters: None.
        Return value: True/False value of reversed flag
        Uses data members (readnly): reversed
        """
        return self.reversed




    def is_resnum_nearer_top(self, resnum):
        """
        Return True if the supplied residue sequence number is nearer
        the 'top' than the 'bottom' of this strand or helix (SSE).
        I.e. if it is nearer
        to the C- than the N- terminal end of the strand when the
        SSE is not reversed (i.e. drawn pointing 'up'). Or, if the
        SSE is reversed (drawn pointing 'down'), if it is nearer he
        N terminal than the C terminal end of the SSE.

        Parameters:
           resnum - residue sequence number, must be in this strand

        Uses data members (readonly)
            start_res_seq, end_res_seq, reversed

        Return value:
            True iff resnum is nearer C term and strand is not reversed OR
                     resnum is nearer N terman and strand is reversed
        """
        assert(self.is_in_interval(resnum))
        midpoint_resnum = ( self.get_residue_ordinal(self.start_res_seq) +
                            (self.get_residue_ordinal(self.end_res_seq)
                             - self.get_residue_ordinal(self.start_res_seq))
                            / 2 )
        if resnum > midpoint_resnum:
            # closer to C terminal end
            return not self.reversed
        else:
            # closer to N terminal end
            return self.reversed



    def add_closenode(self, other_node, dist):
        """
        Add a (node, dist) tuple
        to the list of nodes that are close to this one.

        Parameters:
           other_node - PTNode which is below threshold distance from this node.
           dist - distance to the other_node (Angstroms)

        Uses member data (write):
          closenode_list - the list of nodes close to this one.

        Return value: None
        """
        assert(isinstance(other_node, PTNode))
        self.closenode_list.append((other_node, dist))


    def get_closenode_list(self):
        """
        Return list of (node, dist) tuples for nodes
        that are below threshold distance from this one.

        Parameters: None

        Uses member data (readonly):
          closenode_list - list of close nodes

        Return value:
           List of (ptnode, dist) tuples.  The actual list (ref) in this
           node, not a new copy.
        """
        return self.closenode_list


    def set_endchain(self, endchain):
        """
        Label this node with the endchain flag. See comments in __init__

        Parameters: reversed - True/False endchain flag
        Return value: None
        Uses data members (write): endchain
        """
        self.endchain = endchain

    def get_endchain(self):
        """
        Get the value of the endchain flag. See comments in __init__

        Parameters: None.
        Return value: True/False value of endchain flag
        Uses data members (readnly): endchain
        """
        return self.endchain

    def get_degree(self):
        """
        Return the degree of this node (number of edges incident with it).
        This is equal to the number of spatial edges
           (nubmer of spatially adjacent nodes) plus the number of
           sequence edges which are implicit - all nodes have two except
           the first and last in a chain which have only one.

        Parameters:
           None
        Return value:
           Degree of the node.
        """
        deg = len(self.closenode_list)
        if self.endchain:
            deg += 1
        else:
            deg += 2
        return deg


    def axis_dihedral_angle(self, SSE1, SSE2, pdb_struct):
        """
        We fit an axis to each of the two SSEs, and this SSE, and
        compute the dihedral angle between the two planes defined by
        (1) the axis lines SSE1 and self, and (2) SSE2 and self.

        3 consecutive vectors between 4 points can define a dihedral
        angle. Call the points A,B,C,D in order. ABC is a plane and
        BCD is a plane and we can calculate the angle between those
        planes by the conventional formula (using Bio.PDB Vector
        module).

        Here, instead of having this vectors as actual bonds between
        atoms (as in hbonds_dihedral_angle()), we are using purely the
        abstraction of the SSEs as straight lines. Like
        Arun Konagurthu's method in TableauCreator (email 05Oct2007)
        we choose the points so that the vector between each of the two
        SSE axes and self SSE axis is the shortest line between
        the axes (mutual perpendicular). So we choose A and B as the
        points on SSE1 axis and self axis respectively such that AB is
        the shortest line between the SSE1 axis and this axis.
        Similarly C and D are chosen so that CD is the shortest line
        between this axis and SSE2 axis.




              SSE1       self        SSE2

                |          |          |
                |          |          |            AB is shortest line between
                |          |    v3    |            SSE1 and self, defining
                |        C *--------->* D          vector v1
                |          ^          |
                |          |          |            CD is shortest line between
                |          |          |            self and SSE2, defining
                |          | v2       |            vector v3
                |          |          |
                |         (|)theta    |            v2 is then defined by the
                |          |          |            line BC
              A *--------->* B        |
                |   v1     |          |            and the dihedral angle theta
                |          |          |            between planes ABC and BCD
                |          |          |            is given by:
                |          |          |
                                                            |v2|v1 . (v2 x v3)
                                               tan(theta) = --------------------
                                                           (v1 x v2) . (v2 x v3)


        Parameters:
            SSE1 - first SSE node (strand or helix) to test
            SSE2 - 2nd SSE node (strand or helix) to test
            pdb_struct - The Bio.PDB parsed PDB struct (atomic co-ordinates)
                          for this protein.

        Return value:
            The angle (in (-pi, pi]) between the planes formed between
            the axes fitted to the two SSEs with this SSE in common.
            or None if no common perpendicular can be found or one
            (or more) of the required axes cannot be computed.


        NOTE: this method is only for STRAND and HELIX nodes, which have
        a fit_axis() method. This method is identical in operation between
        these types of nodes so is shared, but fit_axis() works differently
        so PTNodeHelix and PTNodeStrand each have their own implemetnation
        of it.
        """
        SSE1_axis = SSE1.fit_axis(pdb_struct)
        SSE2_axis = SSE2.fit_axis(pdb_struct)
        self_axis = self.fit_axis(pdb_struct)
        if SSE1_axis == None or SSE2_axis == None or self_axis == None:
            return None

        (SSE1_dircos, SSE1_centroid) = SSE1_axis
        (SSE2_dircos, SSE2_centroid) = SSE2_axis
        (self_dircos, self_centroid) = self_axis

        # Get pa and pb, endpoints of the shortest line segment between
        # SSE1 axis line and self axis line
        # Note using Numeric.array '*' operator here for
        # element-wise multiplication
        # as Bio.PDB.Vector '*' operator is vector dot product.
        # (Note also must have Vector + int * array and NOT
        # int * array + Vector due to Python type coercion rules).
        s1self_line = LineLineIntersect(
            SSE1_centroid, SSE1_centroid
                              + ALPHA * SSE1_dircos.get_array(),
            self_centroid, self_centroid
                           + ALPHA * self_dircos.get_array() )
        if s1self_line == None:
            if verbose:
                sys.stderr.write('no common perpendicular for axes ' +
                                 SSE1.nodeid + ', ' + self.nodeid + '\n')
            return None
        else:
             (pa, pb, mua, mub) = s1self_line

        # and pc, pd similarly for self and SSE2
        # Note using Numeric.array '*' operator here for
        # element-wise multiplication
        # as Bio.PDB.Vector '*' operator is vector dot product.
        # (Note also must have Vector + int * array and NOT
        # int * array + Vector due to Python type coercion rules).
        s2self_line = LineLineIntersect(
            self_centroid, self_centroid
                           + ALPHA * self_dircos.get_array(),
            SSE2_centroid, SSE2_centroid
                              + ALPHA * SSE2_dircos.get_array() )
        if s2self_line == None:
            if verbose:
                sys.stderr.write('no common perpendicular for axes ' +
                                 SSE2.nodeid + ', ' + self.nodeid + '\n')
            return None
        else:
            (pc, pd, muc, mud) = s2self_line

#        angle = calc_dihedral(pa, pb, pc, pd)

        v1 = pb - pa
        v2 = pc - pb
        v3 = pd - pc

#        print 'xxx',self.nodeid,SSE1.nodeid,SSE2.nodeid,pa,pb,pc,pd

        # Using Bio.PDB.Vector class; v*v is dot product, v**v is cross product
        # This is the same result as calling Vector.calc_dihedral() anyway
#         theta = atan2( Vector((v2.norm() * v1.get_array())) * (v2 ** v3),
#                        (v1 ** v2) * (v2 ** v3) )

#         print 'xxxx',theta,
#         theta0 = theta
        
        # and this is the more elegant way, not using atan2() (as per Arun).
        normals_dotprod = (v1 ** v2).normalized() * (v2 ** v3).normalized()
        # handle roundoff errors
        normals_dotprod = min(normals_dotprod, 1)
        normals_dotprod = max(normals_dotprod, -1)
        theta = acos(normals_dotprod)
        normals_crossprod = (v1 ** v2).normalized() ** (v2 ** v3).normalized()
        stp = v2 * normals_crossprod
        if stp < 0:
            theta = -theta

#         print 'yyyy',theta
#         assert(abs(theta0 - theta) < EPSILON)


        if verbose:
            sys.stderr.write('axis_dihedral_angle ' + self.nodeid + ' '
                             + SSE1.nodeid + ' '
                             + SSE2.nodeid + ' '
                             + str(theta) + '\n')

        return theta



    def relative_angle(self, SSE1, pdb_struct):
        """
        We fit an axis to this SSE and the supplied other SSE (SSE1), and
        compute the relative angle (omega) between the two axes.
        This is the angle use to define the Tableau
        (See Kamat and Lesk 2007 Proteins 66:869-876 and
        Konagurthu, Stuckey and Lesk 2008 'Structural search and retrieval using
        a tableau representation of protein folding patterns' Bioinformatics
        (advance access, to be published Jan 5 2008).

        As per Arun Konagurthu's method in TableauCreator (email
        05Oct2007) we choose the points so that the vector between
        each of the two SSE axes is the shortest line between the axes
        (mutual perpendicular). So we choose A and B as the points on
        SSE1 axis and self axis respectively such that AB is the
        shortest line between the SSE1 axis and this axis.



              SSE1       self

                |          |
                |          |                       BC is shortest line between
                |          |                       SSE1 and self, defining
                |          |                       vector v2
              A *          * D
                |         /|\
                |          |                       v1 and v3 are vectors
                |          |                       defining the axes of SSE1
                |v1        | v3                    and self respectively.
                |          |
               \|/  omega #|                       The relative angle omega
              B *---(-)--->* C                     (interaxial angle) is the
                |#      v2 |                       smallest angle required to
                |          |                       reorient v1 to eclipse v3
                |          |                       (or vice versa):
                |          |
                                                            |v2|v1 . (v2 x v3)
                                               tan(omega) = --------------------
                                                           (v1 x v2) . (v2 x v3)


        Parameters:
            SSE1 - The other SSE (helix or strand) to get relative angle to self
            pdb_struct - The Bio.PDB parsed PDB struct (atomic co-ordinates)
                          for this protein.

        Return value:
            The angle (in (-pi, pi]) required to reorient self axis
            to eclipse axis of SSE1 (or vice versa) 'looking along' the
            shortest line (mutual perpendicular) between the two axes.


        NOTE: this method is only for STRAND and HELIX nodes, which have
        a fit_axis() method. This method is identical in operation between
        these types of nodes so is shared, but fit_axis() works differently
        so PTNodeHelix and PTNodeStrand each have their own implemetnation
        of it.
        """
        sse1_axis = SSE1.fit_axis(pdb_struct)
        self_axis = self.fit_axis(pdb_struct)
        if sse1_axis == None or self_axis == None:
            return None
        
        (SSE1_dircos, SSE1_centroid) = sse1_axis
        (self_dircos, self_centroid) = self_axis

        pa = SSE1_centroid + ALPHA * SSE1_dircos.get_array()
        pd = self_centroid + ALPHA * self_dircos.get_array()
        
        # Get pb and pc, endpoints of the shortest line segment between
        # SSE1 axis line and self axis line
        # Note using Numeric.array '*' operator here for
        # element-wise multiplication
        # as Bio.PDB.Vector '*' operator is vector dot product.
        # (Note also must have Vector + int * array and NOT
        # int * array + Vector due to Python type coercion rules).
        s1self_line = LineLineIntersect(SSE1_centroid, pa,
                                        self_centroid, pd)
        if s1self_line == None:
            if verbose:
                sys.stderr.write('relative_angle: ' +
                                 'no common perpendicular for axes ' +
                                 str(self) + ',' + str(SSE1) + '\n')
            return None
        else:
             (pb, pc, mub, muc) = s1self_line

#         # DEBUG
#         if verbose:
#             sys.stderr.write('mutual perpendicular ' + self.nodeid + ' '
#                              + SSE1.nodeid +
#                              ': pb = ' + str(array(list(pb))) +
#                              '; pc = ' + str(array(list(pc))) + '\n')
#         # END DEBUG

            
#        omega = calc_dihedral(pa, pb, pc, pd)

        v1 = pb - pa
        v2 = pc - pb
        v3 = pd - pc

        # Using Bio.PDB.Vector class; v*v is dot product, v**v is cross product
        # This is the same result as calling Vector.calc_dihedral() anyway
#         omega = atan2( Vector((v2.norm() * v1.get_array())) * (v2 ** v3),
#                        (v1 ** v2) * (v2 ** v3) )

#         print 'xxxx',omega,
                
        # and this is the more elegant way, not using atan2() (as per Arun).
        normals_dotprod = (v1 ** v2).normalized() * (v2 ** v3).normalized()
        # handle roundoff errors
        normals_dotprod = min(normals_dotprod, 1)
        normals_dotprod = max(normals_dotprod, -1)
        omega = acos(normals_dotprod)
        normals_crossprod = (v1 ** v2).normalized() ** (v2 ** v3).normalized()
        stp = v2 * normals_crossprod
        if stp < 0:
            omega = -omega

#         print 'yyyy',omega


#         DEBUG
#         if verbose:
#             sys.stderr.write('relative_angle ' + self.nodeid + ' '
#                              + SSE1.nodeid + ' '
#                              + str(omega) + '\n')
#         END DEBUG

        return omega


    def build_resname_sequence(self):
        """
        Build list of (3 letter) residue names in sequence for the residues
        in this node (SSE). E.g. and matching
        list of PDB residue sequence numbers 
        
        Parameters:
            None.

        Return value: None
        Uses data member (write): resname_list
                                  resid_list
        """
        residue_list = self.get_residue_list()
        self.resname_list = [residue.get_resname() for residue in residue_list]
        # id of a residue in Bio.PDB is tuple (hetatm, resseqnum, icode)
        self.resid_list = [str(residue.get_id()[1]) +
                           char_if_not_blank(residue.get_id()[2])
                           for residue in residue_list]

    def set_seq_component_id(self, seqc_id):
        """
        Label this node with a sequence connected component identifier

        Parameters:
           seqc_id - int seq connected component id
        Uses data members:
           seq_component_id (write)
        Return value: None
        """
        self.seq_component_id = seqc_id

    def get_seq_component_id(self):
        """
        Return the id of the seq connectec component to which
        this node belongs (may be None)

        Parameters: None
        Uses data members (readonly): seq_component_id
        Return value: Seq cc id set in this node (int) or None
        """
        return self.seq_component_id


class PTNodeHelix(PTNode):
    """
    The PTNodeHelix class is the type of PTNode for an alpha helix.
    """

    def __init__(self, helixtype="ALPHA", *args):
        """
        Construct PTNodeHelix with supplied nodeid and type.
        Parameters:
             helixtype - string "ALPHA" or "PI" or "310". default "ALPHA".
             +Variable parameter list: straight to PTNode constructor (q.v.).

        This extends PTNode by adding is_positioned,
        and other than that just calls PTNode constructor.
        Raises exceptions:
            TypeError if helixtype argument is invalid.
        """
        if helixtype not in ["ALPHA", "PI", "310"]:
            raise TypeError("Invalid helixtype " + helixtype + "\n")

        PTNode.__init__(self, *args)
        self.helixtype = helixtype

        #
        # Flags used in heuristic helix placement
        #
        self.is_positioned = False # set to True if helix is already positioned
                                   # (drawn). Used when old_helix_placement
                                   # (i.e. -i option NOT supplied) is in use,
                                   # sometimes as a special case we position
                                   # helices before calling write_helices_svg()
                                   # and this flag is set to mark these as
                                   # already positioned.
                                   # use get/set_is_positioned
        self.is_interstrand = False # Set to True to mark special case
                                    # of helix between strands on same
                                    # vert axis in same sheet, treated
                                    # specially when using heuristic helix
                                    # placement to
                                    # force interstand helices beside sheet.
                                    # If this is set then the reversed
                                    # flag of this helix node (in base class)
                                    # is set to reversed flag of
                                    # the first strand N-terminal to it.
                                    # (see write_insterstrand_helices_svg()).
                                    # Use get/set is_interstrand.
        self.is_last_in_seq = False # Set to True if this helix is the last
                                    # in a sequence of helices all aligned
                                    # on an axis, AND n-terminal strand
                                    # is on different axis.
                                    # Used so that if
                                    # connector into this helix is on
                                    # another axis then we know to use
                                    # the top (reversed) else bottom
                                    # port instead of what what normally
                                    # be used (reveresed is set as per
                                    # is inter_strand (see comments above).
                                    # Use get/set is_last_in_seq().
        self.cluster_id = None      # If use helix 'clustering' then this is
                                    # the intege id of the cluster to which this
                                    # helix belongs (ids assigned sequentially
                                    # starting from 1).
                                    # Use get/set cluster_id()

        self.axis_centroid = None # Bio.PDB.Vector representing the
                                  # centroid of midpoints of consecutive
                                  # c_alpha atoms of consecutive residue triples
                                  # in the helix. Set by fit_axis()

        self.axis_direction_cosines = None # Bio.PDB.Vector representing the
                                           # direction cosines of the axis line
                                           # fitted to this helix. Set by
                                           # fit_axis()
        self.axis_nterm_point = None # Bio.PDB.Vector of projection most
                                     # N-terminal point of SSE onto axis.
                                     # Set by fit_axis()
        self.axis_cterm_point = None # Bio.PDB.Vector of projection of most
                                     # C-terminal point of SSE onto axis.
                                     # Set by fit_axis()
                                     


    def __str__(self):
        """
        Return String representation of the node as
        'TYPE id [startResNum..endResNum]
        """
        return self.helixtype + " " +\
               self.nodeid + "[" + str(self.start_res_seq) + \
               ".." + str(self.end_res_seq)  + "]"

    def get_is_positioned(self):
        """
        Return True if the helix is marked already positioned
        See is_positioned in __init__()
        Parameters: None
        Return value: True if helix is marked already positioned else False
        Uses member data (readonly): is_positioned
        """
        return self.is_positioned

    def set_is_positioned(self, is_pos):
        """
        Set the is_positioned flag to the supplied boolean value.
        See is_positioned in __init__()
        Parmeters:
            is_pos - True to mark as already positioned, False to unmark
        Return value: None
        Uses member data (WRITE): is_positioned
        """
        self.is_positioned = is_pos

    def get_type(self):
        """
        Return the helix type 'ALPHA', 'PI' or '310'
        Parameters: None
        Return value: helix type as above.
        """
        return self.helixtype

    def get_is_interstrand(self):
        """
        Return True if the helix is marked already interstrand
        See is_interstrand in __init__()
        Parameters: None
        Return value: True if helix is marked already interstrand else False
        Uses member data (readonly): is_interstrand
        """
        return self.is_interstrand

    def set_is_interstrand(self, is_intstr):
        """
        Set the is_interstrand flag to the supplied boolean value.
        See is_interstrand in __init__()
        Parmeters:
            is_pos - True to mark as already interstrand, False to unmark
        Return value: None
        Uses member data (WRITE): is_interstrand
        """
        self.is_interstrand = is_intstr

    def get_is_last_in_seq(self):
        """
        Return True if the helix is marked already last_in_seq
        See is_last_in_seq in __init__()
        Parameters: None
        Return value: True if helix is marked already last_in_seq else False
        Uses member data (readonly): is_last_in_seq
        """
        return self.is_last_in_seq

    def set_is_last_in_seq(self, is_last):
        """
        Set the is_last_in_seq flag to the supplied boolean value.
        See is_last_in_seq in __init__()
        Parmeters:
            is_pos - True to mark as already last_in_seq, False to unmark
        Return value: None
        Uses member data (WRITE): is_last_in_seq
        """
        self.is_last_in_seq = is_last

    def get_cluster_id(self):
        """
        Return the cluster id of cluster to which this helix belongs.
        See cluster_id in __init__()
        Parmaters: None
        Return value: cluster id (integer from 1) or None if no cluster
        USes member data (readonly): cluster_id
        """
        return self.cluster_id

    def set_cluster_id(self, cluster_id):
        """
        Set the cluster id of the cluser to which this helix belongs.
        Parameters:
           cluster_id: integer > 0 or None
        Return value: None
        Uses member data (write): cluster_id
        """
        self.cluster_id = cluster_id


    def fit_axis(self, pdb_struct, mfile_fh = None):
        """
        Approximate this helix as a straight line by fitting a total least
        squares line through the midpoints of planes formed by
        C-alpha atoms all consecutive residue triples in the helix.

        This is the method used in Arun Konagurthu's TableauCreator program,
        this is more or less a re-implementation in Python of the C++ function
        used there.
        (See Kamat and Lesk 2007 Proteins 66:869-876 and
        Konagurthu, Stuckey and Lesk 2008 'Structural search and retrieval using
        a tableau representation of protein folding patterns' Bioinformatics
        (advance access, to be published Jan 5 2008).
        Parameters:
             pdb_struct - The Bio.PDB parsed PDB struct (atomic co-ordinates)
                          for this protein.
             mfile_fh - (Default None)
                        filehandle (open write) to write MATLAB commands
                        for plotting helix axis data to, or None for no MATLAB.

        Return value:
             tuple (direction_cosines, centroid)
             direction_coisines is Vector (a,b,c), the direction cosines
              of the axis
             centroid is Vector (x,y,z) the centroid of the midpoints of
              planes formed by C_alpha atoms of consecutive residue triples.

        Uses data members:
             axis_centroid (read/write) - Vector of centroid of triplet planes.
             axis_direction_cosines (read/write) - Vector of direction cosines
                                                   of axis line

        Note: this function is memoized - the first time it is called it
              computes the axis and stores it in data members as well as
              returning it; subsequent calls return the stored values.
              This is because we don't necessarily need this computation
              for all helices so don't want to always compute it up front
              but also may need it multiple times for some helices.
        """
        # return stored values if already computed
        if self.axis_direction_cosines != None:
            return (self.axis_direction_cosines, self.axis_centroid)

        # We use Numeric rather than newer NumPy since BioPython requires
        # Numeric anyway.

        # get list of position vectors of C_alpha atoms in this helix
#        c_alpha_veclist = [ residue['CA'].get_vector()
#                            for residue in self.get_residue_list() ]
        # need to do it this uglier way to handle exception
        c_alpha_veclist = []
        for residue in self.get_residue_list():
            try:
                c_alpha_vector = residue['CA'].get_vector()
            except KeyError: # sometimes happens e.g.  d1mj4a_.ent with pmml
                sys.stderr.write('WARNING: no CA atom for residue ' +
                                 str(residue) + ', using C atom\n')
                c_alpha_vector = residue['C'].get_vector()
            c_alpha_veclist.append(c_alpha_vector)


        if len (c_alpha_veclist) < 3:
            sys.stderr.write('WARNING: helix ' + str(self) +
                             ' has only ' + str(len(c_alpha_veclist))
                             + ' residues, cannot fit axis\n')
            return None

        # build list of position vectors of midpoints of each three
        # consecutive C_alpha atoms.
        midpoints = []
        for i in range(1, len(c_alpha_veclist)-1):
            vec1 = c_alpha_veclist[i-1] - c_alpha_veclist[i]
            vec2 = c_alpha_veclist[i+1] - c_alpha_veclist[i]
            midpoint = (vec1 + vec2) / 2
            # convert to position vector relative to C_alpha of residue i
            midpoint = c_alpha_veclist[i] + midpoint
            midpoints.append(midpoint)

        # calculate centroid of the midpoints (axis must pass through this)
        centroid = Vector([0,0,0])
        for i in range(len(midpoints)):
            centroid += midpoints[i]
        centroid /= len(midpoints)

        if len(c_alpha_veclist) >= 4:
            # build array A where each row is vector from centroid
            # to midpoint
#            midpoints.reverse() #XXX
            A = array([list(mp - centroid) for mp in midpoints])

            # compute the SVD of the array A, giving singular values on
            # diagonal of s and right singular vectors, forming
            # an orthonormal basis for the space spanned by rows of A,
            # as columns of v
            # (i.e. rows of vt, vt = transpose(v)).
            # TODO: handle exception LinAlgError and fit axis by other
            # method (as when fewer than 4 residues). Never actually seen
            # this exception occur here but it could.
            (u, s, vt) = singular_value_decomposition(A)
            # the direction cosine is the first row of vt, this is the basis
            # vector corresponding to first (largest) singular value,
            # i.e. the dimension in which variance is greatest, hence this
            # must be the major axis of the helix.
            dircos = Vector(vt[0,:])

            # get projection of most N-terminal and most C-terminal midpionts
            # onto the axis line
            self.axis_nterm_point = ProjectPointOntoLine(centroid, centroid +
                                               ALPHA*dircos.get_array(),
                                               midpoints[0])
            self.axis_cterm_point = ProjectPointOntoLine(centroid, centroid +
                                               ALPHA*dircos.get_array(),
                                               midpoints[-1])
            # The dircos gives the axis of the helix, but there is no
            # guarantee that it 'points the right way' i.e. from the N-
            # to C- terminus. We test if it does by finding the angle
            # between it and the vector pointing frmo nterm_piont to
            # cterm_point. This angle must either be 0 or pi, if 0 then the
            # axis pionts the right way, otherwise it is pointing the wrong
            # way so we reverse it.

            # FIXME: instead of using Vector.angle() should just compute
            # ((cterm_point-nterm_point) dotproduct dircos)
            #   / (cterm_point-nterm_point).norm() * dircos.norm()
            # and test for 1 or -1, avoiding arccosine computation.
            # (Actually don't need to normalize, could just check +ve or -ve).
            # In fact, we probably don't really need to poject points onto
            # line at all, could just compute angle between dircos and
            # vector through N- and C- terminal midpoints.
            angle = (self.axis_cterm_point-self.axis_nterm_point).angle(dircos)
            assert (abs(angle) < EPSILON or abs(angle - pi) < EPSILON)
            if (abs(angle - pi) < EPSILON):
                dircos = Vector(-1 * dircos.get_array())
                if verbose:
                    sys.stderr.write('reversed axis for ' + str(self) + '\n')
                    
                
#             # DEBUG
#             print 'helix axis svd ',str(self),':'
#             print 'A = '
#             print A
#             print 's = '
#             print s
#             print 'vt= '
#             print vt
#             print 'dircos='
#             print dircos
#             # END DEBUG
            
            if mfile_fh != None:
                mfile_write_helix(mfile_fh, c_alpha_veclist,
                                  midpoints, centroid, dircos,
                                  self.axis_nterm_point,
                                  self.axis_cterm_point,
                                  self.nodeid)
        else:
            # if we have fewer than 4 residues in the helix then we cannot
            # compute an axis with SVD using midpoints of consecutive
            # residue triples.
            # If here are 3 residues then we can compute two midpionts for
            # consecutive triples, an will just say the axis is the line
            # through these 2 points. Otherwise, for 1 or 2 residues,
            # we cannot sensibly compute any axis.
            assert(len(c_alpha_veclist) == 3)
            mp1 = (c_alpha_veclist[1] - c_alpha_veclist[0]) / 2 + \
                  c_alpha_veclist[0]
            mp2 = (c_alpha_veclist[2] - c_alpha_veclist[1]) / 2 + \
                  c_alpha_veclist[1]
            centroid = (mp1 + mp2) / 2
            v = mp2 - mp1
            dircos = v / v.norm()
            if mfile_fh != None:
                # get projection of most N-terminal and most
                # C-terminal C_alpha atom points onto the axis line
                self.axis_nterm_point = ProjectPointOntoLine(centroid,
                                                   centroid +
                                                   ALPHA*dircos.get_array(),
                                                   c_alpha_veclist[0])
                self.axis_cterm_point = ProjectPointOntoLine(centroid,
                                           centroid +
                                           ALPHA*dircos.get_array(),
                                           c_alpha_veclist[-1])
                mfile_write_helix(mfile_fh, c_alpha_veclist,
                                  [mp1, mp2], centroid, dircos,
                                  self.axis_nterm_point,
                                  self.axis_cterm_point,
                                  self.nodeid)

        self.axis_direction_cosines = dircos
        self.axis_centroid = centroid
        return (dircos, centroid)


class PTNodeStrand(PTNode):
    """
    The PTNodeStrand class is the type of PTNode for a beta strand.

    The additions it has are the bridge_list and methods for using it -
    beta strands can have beta bridges to other beta strands, and a beta
    sheet label.
    """
    def __init__(self, *args):
        """
        Construct a PTNodeStrand with supplied nodeid and type, and empty
        bridge bond list.

        Parameters:
             Variable parameter list: straight to PTNode constructor (q.v.).

        This extends PTNode by adding bridge_list,
        and other than that just calls PTNode constructor.

        """
        PTNode.__init__(self, *args)
        self.bridge_list = [] # list of (other_node, bdir, side) tuples
                              # bdir is 'N' or 'P' for antiparallel or
                              # parallel respectively
                              # side is '+' or '-' to indicate relative
                              # side of this strand bridge partners are
                              # on, as per Westhead et al 1999
                              # (see label_strand_bridge_sides())
                              # use add_bridge(), remove_brdige(),
                              # get_bridge_list()

        self.sheet_id = None  # id of sheet that this strand belongs to
                              # single char 'A', 'B', etc.
                              # use get/set_sheet_id()


        self.align_pos = 0    # The relative 'vertical' (assuming
                              # strands are drawn as arrows 'pointing'
                              # up or down) alignmenent position of this strand
                              # on its vertical axis.
                              # This value is in 'residues' (i.e. number
                              # of residues from 'top' of neighbour strand)
                              # to 'top' of this strand.
                              # Maybe be positive or negative (or zero).
                              # (see build_sheet_constraints() in ptgraph2.py)
                              # use get/set_align_pos()
        self.barrel_edge = False # True if this strand is one of the two
                                 # in a beta-barrel that was opened up
                                 # by breaking the bridge between them
                                 # use get/set_barrel_edge()

        self.axis_centroid = None # Bio.PDB.Vector representing the
                                  # centroid of c_alpha atoms of residues
                                  # in the strand. Set by fit_axis()

        self.axis_direction_cosines = None # Bio.PDB.Vector representing the
                                           # direction cosines of the axis line
                                           # fitted to this strand. Set by
                                           # fit_axis()
        self.axis_nterm_point = None # Bio.PDB.Vector of projection most
                                     # N-terminal point of SSE onto axis.
                                     # Set by fit_axis()
        self.axis_cterm_point = None # Bio.PDB.Vector of projection of most
                                     # C-terminal point of SSE onto axis.
                                     # Set by fit_axis()
                                     



    def __str__(self):
        """
        Return String representation of the node as
        'TYPE id [startResNum..endResNum]'
        """
        return "STRAND" + " " + self.nodeid + "[" + str(self.start_res_seq) + \
               ".." + str(self.end_res_seq)  + "]"


    def add_bridge(self, other_node, bdir):
        """
        Add a bridge to another strand to the table of bridges in this node.
        The bridge is to other_node and bdir is 'N' for antiparallel or
        'P' for parallel.

        Does not add duplicates: if there is already an edge to the supplied
        node, this new one is not added; also, ensure that these edges
        are undirected by adding the reverse edge at the same time
        (this is really the only reason we need to check for duplicates:
        it seems that stride sometimes gives a bridge partner from one
        strand to another but not back the other way, but not always).

        Bridge side label ('+' or '-') is set as '.', this is actually set
        by label_bridge_sides() to be called afterwards.

        Parameters:
           other_node - PTNode the bond is to from this node (see NOTE)
           bdir       - 'N' or 'P' for antiparrallel or parallel resp.

        Uses data members (write):
           brdige_list - list of (ptnodestrand, bdir, side) tuples

        NOTE: also modifies other_node by adding an edge in its bridge_list
        directly.

        Return value: None
        """

        assert(isinstance(other_node, PTNodeStrand))
        assert(bdir == 'N' or bdir == 'P')
        if other_node not in [ node for (node, bdir_unused, side_unused)
                               in self.bridge_list ]:
            self.bridge_list.append((other_node, bdir, '.'))
            other_node.bridge_list.append((self, bdir, '.'))

        #----- debug TESTING FOR STRANDS WITH MORE THAN 2 PARTNERS ---
        if verbose:
            if len(self.bridge_list) > 2:
                sys.stderr.write(self.nodeid + " has " \
                                 + str(len(self.bridge_list)) +\
                                 " adjacent strands\n")
        #----- end -----


    def remove_bridge(self, other_node):
        """
        Remove the bridge to other_node from the list of bridges in this node.

        Parameters:
            other_node - PTNode the bridge is to from this node (see NOTE)

        Uses data members (write):
           brdige_list - list of (ptnodestrand, bdir, side) tuples


        NOTE: also modifies other_node by adding an edge in its bridge_list
        directly.

        Return value: None

        Raises exceptions:
           KeyError if other_node not found in bridge list

        """
        found = False
        for i in range(len(self.bridge_list)):
            if self.bridge_list[i][0] == other_node:
                found = True
                break
        if found:
            self.bridge_list.pop(i)
            # now remove other node's bridge to this one
            found = False
            for i in range(len(other_node.bridge_list)):
                if other_node.bridge_list[i][0] == self:
                    found = True
                    break
            assert(found) # logic error if there wasn't a matching bridge
            other_node.bridge_list.pop(i)
        else:
            raise KeyError("node not found")



    def set_sheet_id(self, sheet_id):
        """
        Label this strand with a sheet id to which it belongs.

        Parameters:
           sheet_id - single char sheet identifier
        Uses data members:
           sheet_id - sheet id (write)
        Return value: None
        """
        self.sheet_id = sheet_id

    def get_sheet_id(self):
        """
        Return the id of the sheet to which this strand belongs.

        Parameters: None
        Uses data members (readonly): sheet_id
        Return value: Sheet id set in this node (single char)
        """
        return self.sheet_id

    def set_align_pos(self, align_pos):
        """
        Set the (integer) relative alignment position for this strand.
        See comments on align_pos in __init__

        Parmeters: align_pos (integer)
        Return value: None
        Uses data members (write): align_pos
        """
        self.align_pos = align_pos

    def get_align_pos(self):
        """
        Get the relative alignment position of this strand.
        See comments on align_pos in __init__

        Paramneters: None
        Return value: align_pos value (integer)
        Uses data members (Readonly): align_pos
        """
        return self.align_pos

    def set_barrel_edge(self, barrel_edge):
        """
        Set the barrel_edge flag in this strand to the supplied value.
        See comments on barrel_edge in __init__
        Parameters: barrel_edge: True or False
        REturn value: None
        Uses data members (write): barrel_edge
        """
        self.barrel_edge = barrel_edge

    def get_barrel_edge(self):
        """
        Return the value of the barrel_edge flag in this node
        See comments on barrel_edge in __init__
        Parameter: None
        Return value: barrel_edge flag (True or False)
        Uses data members (readonly): barrel_edge
        """
        return self.barrel_edge

    def get_bridge_list(self):
        """
        Return list of bridges from this node in form of list of
        (ptnodestand, bdir, side) tuples.

        Parameters: None
        Uses member data (readonly):
           bridge_list - list of bridge tuples
        Return value:
           list of (ptnodestrand, bdir, side) tuples. This list is pointer to
           the list in this node, not a new copy.
        """
        return self.bridge_list


    def num_neighbours(self):
        """
        Return the number of neighbouring strands that this strand has.
        This is just the number of elements in the bridge list.

        Parameters: None
        Uses member data (Readonly):
           bridge_list - list of brdige tuples
        Return value:
           number of neighbouring strands (length of bridge list)
        """
        return len(self.bridge_list)

    def is_neighbour(self, strand):
        """
        Return True iff the supplied strand is a neighbour of this one
        (i.e. is in the bridge list)

        Parameters:
            strand - strand to test for being in the bridge list of this node
        Uses data members (readonly):
            bridge_list - list of bridge tuples (ptnodestrand, bdir, side)
        Return value:
            True if strand is in the list of bridges at this node, else False
        """
        assert(isinstance(strand, PTNodeStrand))
        for (node, bdir_unused, side_unused) in self.bridge_list:
            if node == strand:
                return True
        return False


    def label_strand_bridge_sides(self, pdb_struct):
        """
        Label bridge edges with '+' or '-' indicating relative side of
        the strand the bridge partners are on, using hydrogen bond overlap.
        For one strand with two more more neighbours, if there is H bond
        overlap (i.e. one residue in the reference strand has more than
        one neighbour) then those neighbours must be on opposite sides
        of the reference strand.

        Otherwise, they may be on the same side, (but not necessarily).
        This has to be determined by geometric criteria, specifically
        the dihedral angle (NB this is not the usual phi/psi meaning
        of dihedral, but the generic meaning of an angle between planes)
        formed between two carbon alpha atoms on the reference strand
        and one on each of the other two strands.

        Every PTNodeStrand with neighbour(s) has a bridge edge (tuple
        in the bridge_list) going out to the neighbour, and that neighbour
        symmetrically has one going back (see PTNodeStrand.add_bridge()).
        So for a strand with one neighbour, the label on the bridge edge
        is '+'; if there is more than one neighbour, the edges are labelled
        '+' for those on one side and '-' on the other.

        This algorithm is described in
        Westhead et al 1999 ``Protein structural topology: Automated
        analysis and diagrammatic representation'' Protein Science 8:897-904
        (see pp.901-902 and Figures 4 and 5).

        Parameters:
             pdb_struct - The Bio.PDB parsed PDB struct (atomic co-ordinates)
                          for this protein.

        Uses data members
               (read/write):
             bridge_list - list of bridge (node, bdir, side) tuples
               (readonly):
             hydrogen_bond_list - list of H-bond (node, resnum1, resnum2, dist)
                                  tuples

        Return value: None
        """
        if len(self.bridge_list) == 0:
            return # no neighbours, nothing to do here

        # now go through all pairs of neighbour strands (if any)
        # in turn, labelling as opposite sides if H-bond nesting criteria
        # tells us they must be on different sides of this strand.
        for i in range(len(self.bridge_list)):
            (node_i, bdir_i, side_i) = self.bridge_list[i]
            for j in range(i+1, len(self.bridge_list)):
                (node_j, bdir_j, side_j) = self.bridge_list[j]
                # FIXME: really has_hbond_strand_overlap() implies
                # has_strand_extent_overlap() (or should do), so
                # should be able to remove the former.
                if (self.has_hbond_strand_overlap(node_i, node_j) or
                    self.has_strand_extent_overlap(node_i, node_j) ):
                    # node_j has overlap with node_i on this strand, so
                    # set its side to opposite of that of node_i
                    if side_i == '+':
                        self.bridge_list[j] = (node_j, bdir_j, '-')
                    elif side_i == '-':
                        self.bridge_list[j] = (node_j, bdir_j, '+')
                    else:
                        self.bridge_list[i] = (node_i, bdir_i, '+')
                        self.bridge_list[j] = (node_j, bdir_j, '-')
                    if verbose:
                        sys.stderr.write('overlap; opposite sides of '+
                                         self.nodeid + ': ' +
                                         self.bridge_list[i][0].nodeid + ','+
                                         self.bridge_list[j][0].nodeid + '\n')


        # now for all pairs where H bond overlap did not manage to set
        # a relative side, use geometry test where no overlap
        for i in range(len(self.bridge_list)):
            (node_i, bdir_i, side_i) = self.bridge_list[i]
            for j in range(i+1, len(self.bridge_list)):
                (node_j, bdir_j, side_j) = self.bridge_list[j]
                if ((side_j == '+' or side_j == '-') and
                    (side_i == '+' or side_i == '-')):
                    continue # already set, skip it
                if self.strands_on_opposite_sides(node_i, node_j, pdb_struct):
                    # node_j and node_i are on different sides of this
                    # strand, so set its side to opposite of that of
                    # node_i, if side set
                    if verbose:
                        sys.stderr.write('OPPOSITESIDES of '+self.nodeid+': '+
                              node_i.nodeid+','+node_j.nodeid+'\n')
                    if side_j != '+' and side_j != '-':
                        if side_i == '+':
                            self.bridge_list[j] = (node_j, bdir_j, '-')
                        elif side_i == '-':
                            self.bridge_list[j] = (node_j, bdir_j, '+')
                        else:
                            self.bridge_list[i] = (node_i, bdir_i, '+')
                            self.bridge_list[j] = (node_j, bdir_j, '-')
                    elif side_i != '+' and side_i != '-':
                        if side_j == '+':
                            self.bridge_list[i] = (node_i, bdir_i, '-')
                        elif side_j == '-':
                            self.bridge_list[i] = (node_i, bdir_i, '+')
                        else:
                            self.bridge_list[i] = (node_i, bdir_i, '+')
                            self.bridge_list[j] = (node_j, bdir_j, '-')
                    else:
                        # can't happen since skipped if both sides set
                        assert(False)
                else:
                    # they must be on the same side of this strand
                    if verbose:
                        sys.stderr.write('SAMESIDE of '+self.nodeid+': '+
                              node_i.nodeid+','+node_j.nodeid+'\n')
                    if side_i == '+' or side_i == '-':
                        self.bridge_list[j] = (node_j, bdir_j, side_i)
                    elif side_j == '+' or side_j == '-':
                        self.bridge_list[i] = (node_i, bdir_i, side_j)
                    else:
                        self.bridge_list[i] = (node_i, bdir_i, '+')
                        self.bridge_list[j] = (node_j, bdir_j, '+')


    def strands_on_opposite_sides(self, strand1, strand2, pdb_struct):
        """
        Test if the two strands are on opposite sides of this strand
        by using a geometric test like that described in Westhead et al 1999.

        This is determined by geometric criteria, specifically
        the dihedral angle (NB this is not the usual phi/psi meaning
        of dihedral, but the generic meaning of an angle between planes)
        between the strands.

        The dihedral angle as calculated by the Bio.PDB
        Vector module is in the range (-pi, pi] so if the absolute
        value of this angle is larger than pi/2 we say the strands
        are on opposite sides of the reference strand, else they are
        on the same side.

        If the above method fails because a common perpendicular
        cannot be found, the alternative method of averaging the dihedral
        angles for all H bonds along the SSEs is used.

        Parameters:
            strand1 - first STRAND node to test, a bridge partner of this one
            strand2 - 2nd STRAND node to test, a bridge partner of this one
            pdb_struct - The Bio.PDB parsed PDB struct (atomic co-ordinates)
                          for this protein.

        Return value:
            True if the geometric criteria between this and
            the neighouring strands indicate they (the neighbour strands)
            are on opposite sides of this strand. Else
            False if it indcates they are on the same side.
        """
        angle = self.axis_dihedral_angle(strand1, strand2, pdb_struct)
        if angle == None:
            angle = self.calc_average_hbonds_dihedral_angle(strand1, strand2,
                                                            pdb_struct)

        if abs(angle) > pi/2:
            return True
        else:
            return False




    def calc_average_hbonds_dihedral_angle(self, strand1, strand2, pdb_struct):
        """
        Return the average abosolute value of the dihedral angle
        formed between two carbon alpha atoms on the reference strand
        and one on each of the other two strands, over all such
        hydrogen bonds.

        Parameters:
            strand1 - first STRAND node to test, a bridge partner of this one
            strand2 - 2nd STRAND node to test, a bridge partner of this one
            pdb_struct - The Bio.PDB parsed PDB struct (atomic co-ordinates)
                          for this protein.

        Return value:
           The average absolute value of the dihedral angle
           between planes formed by hydrogen bonds between the
           strands, which is in the range [0, pi].
        """

        hbonds_to_strand1 = self.get_hbonds_to_neighbour(strand1)
        hbonds_to_strand2 = self.get_hbonds_to_neighbour(strand2)
        if len(hbonds_to_strand1) == 0 or len(hbonds_to_strand2) == 0:
            # sometimes the (modified) stride FA1,FA2 records will
            # indicate bridge partners but we won't have any corresponding
            # H bonds from DNR/ACC records, so we don't find hbonds
            # here. Will just return True (opposite sides) since that's
            # what the most common case is.
            sys.stderr.write(
                'WARNING: no hbonds available between neighbour strands ' +
                             strand1.nodeid + ', ' +
                             self.nodeid + ', ' +
                             strand2.nodeid +'\n')
            return True

        # compute dihedral angle between planes formed by each pair of
        # hydrogen bonds, and return average absolute value.
        sum_abs_angles = 0.0
        num_angles = 0
        for hb1 in hbonds_to_strand1:
            for hb2 in hbonds_to_strand2:
                angle = self.hbonds_dihedral_angle(hb1, hb2, pdb_struct)
                if verbose:
                    sys.stderr.write('angle '+ self.nodeid + ' '
                                     + strand1.nodeid + ' ' + str(hb1[1]) + ' '
                                     + strand2.nodeid + ' ' + str(hb2[1]) + ' '
                                     + str(angle) + '\n')
                num_angles += 1
                sum_abs_angles += abs(angle)

        avg_abs_angle = sum_abs_angles / num_angles
        if verbose:
            sys.stderr.write('avg abs angle ' + self.nodeid + ' '
                             + strand1.nodeid + ' '
                             + strand2.nodeid + ' '
                             + str(avg_abs_angle) + '\n')
        return avg_abs_angle


    def hbonds_dihedral_angle(self, hbond_to_strand1, hbond_to_strand2,
                              pdb_struct):
        """
        Calculate the dihedral angle formed with hbonds from this
        strand to two others.

        3 consecutive bonds between 4 atoms can define a dihedral
        angle. Call the atoms A,B,C,D in order. ABC is a plane and
        BCD is a plane and we can calculate the angle between those
        planes by the conventional formula (using Bio.PDB Vector
        module). We choose B and C as the carbon alpha atoms of
        residues on the reference strand (self), of residues that
        have H-bonds to residues on the first and second test
        strands respectively, choosing their carbon alpha atoms to
        be A and D.

        Parameters:
            bhond_to_strand1 - hbond tuple (node, resnum1, resnum2, dist)
                               of H bond to first STRAND node to test,
                               a bridge partner of this one
            hbond_to_strand2 - hbond tuple (node, resnum1, resnum2, dist)
                               of H bond to 2nd STRAND node to test,
                               a bridge partner of this one
            pdb_struct - The Bio.PDB parsed PDB struct (atomic co-ordinates)
                          for this protein.

        Return value:
            The angle (in (-pi, pi]) between the planes formed between two
            C_alpha atoms on the reference (this) strand and respectively
            a C_alpha atom on strand1 and strand2 is (absolute value)
        """

        chainid_self = stride_chainid_to_pdb_chainid(self.chainid)
        chainid_strand1 = stride_chainid_to_pdb_chainid(
                                             hbond_to_strand1[0].get_chainid())
        chainid_strand2 = stride_chainid_to_pdb_chainid(
                                             hbond_to_strand2[0].get_chainid())

        pdb_model = pdb_struct[0] # TODO always using model 0 for now

        # hbond is tuple (node, resnum1, resnum2, dist)
        # resnum1 is the residue seqnum in self, resnum2 in node
        CA_B = pdb_model[chainid_self][hbond_to_strand1[1]]['CA']
        CA_C = pdb_model[chainid_self][hbond_to_strand2[1]]['CA']
        CA_A = pdb_model[chainid_strand1][hbond_to_strand1[2]]['CA']
        CA_D = pdb_model[chainid_strand2][hbond_to_strand2[2]]['CA']

        angle = calc_dihedral(CA_A.get_vector(),
                              CA_B.get_vector(), CA_C.get_vector(),
                              CA_D.get_vector())
        return angle


    def fit_axis(self, pdb_struct, mfile_fh = None):
        """
        Approximate this strand as a straight line by fitting a total least
        squares line through the midpoints of consecutive
        C-alpha atoms of its residues.

        Parameters:
             pdb_struct - The Bio.PDB parsed PDB struct (atomic co-ordinates)
                          for this protein.
             mfile_fh - (Default None)
                        filehandle (open write) to write MATLAB commands
                        for plotting strand data to, or None for no MATLAB.

        Return value:
             tuple (direction_cosines, centroid)
             direction_coisines is Vector (a,b,c), the direction cosines
              of the axis
             centroid is Vector (x,y,z) the centroid of the midpoints of
               consectuvie c_alpha atoms.

        Uses data members:
             axis_centroid (read/write) - Vector of centroid of c_alpha atoms
             axis_direction_cosines (read/write) - Vector of direction cosines
                                                   of axis line

        Note: this function is memoized - the first time it is called it
              computes the axis and stores it in data members as well as
              returning it; subsequent calls return the stored values.
              This is because we don't necessarily need this computation
              for all strands so don't want to always compute it up front
              but also may need it multiple times for some strands.
        """
        # return stored values if already computed
        if self.axis_direction_cosines != None:
            return (self.axis_direction_cosines, self.axis_centroid)

        # We use Numeric rather than newer NumPy since BioPython requires
        # Numeric anyway.

        c_alpha_veclist = []
        for residue in self.get_residue_list():
            try:
                c_alpha_vector = residue['CA'].get_vector()
            except KeyError: # sometimes happens e.g. residue 10 in d1cm3a_
                sys.stderr.write('WARNING: no CA atom for residue ' +
                                 str(residue) + ', using C atom\n')
                try:
                    c_alpha_vector = residue['C'].get_vector()
                except KeyError: # e.g. residue 33 in d1qdla_  with pmml
                    sys.stderr.write('WARNING: no CA or C atom for residue ' +
                                    str(residue) + ', skipping this residue\n')
                    continue
            
            c_alpha_veclist.append(c_alpha_vector)


        # calculate centroid of the points (axis must pass through this)
        centroid = Vector([0,0,0])
        for i in range(len(c_alpha_veclist)):
            centroid += c_alpha_veclist[i]
        centroid /= len(c_alpha_veclist)

        # find midpoint of each pair of consecutive c_alpha atoms in
        # order to reduce error in fitting line due to pleat of beta strand
        # (Cohen et al 1981, J. Mol. Biol. 148(3):253-272)
        # and build list of vectors from centroid to each midpoint
        if len(c_alpha_veclist) > 3:
            centroid_mp_veclist = []
            for i in range(len(c_alpha_veclist)-1):
                midpoint = ((c_alpha_veclist[i+1] - c_alpha_veclist[i])/2) + \
                           c_alpha_veclist[i]
                centroid_mp_veclist.append(midpoint - centroid)

            # build array A where each row is vector from centroid
            # to midpoint
            A = array([list(vect) for vect in centroid_mp_veclist])

            # compute the SVD of the array A, givng singular values on
            # diagonal of s and right singular vectors as columns of v
            # (i.e. rows of vt, vt = transpose(v)).
            # TODO: handle exception LinAlgError and fit axis by other
            # method (as when fewer than 3 residues). Never actually seen
            # this exception occur here but it could.
            (u, s, vt) = singular_value_decomposition(A)
            # the direction cosine is the first row of vt
            dircos = Vector(vt[0,:])

            # get projection of most N-terminal and most C-terminal 
            # C_alpha astoms onto the axis line
            self.axis_nterm_point = ProjectPointOntoLine(centroid, centroid +
                                               ALPHA*dircos.get_array(),
                                               c_alpha_veclist[0])
            self.axis_cterm_point = ProjectPointOntoLine(centroid, centroid +
                                               ALPHA*dircos.get_array(),
                                               c_alpha_veclist[-1])
            # The dircos gives the axis of the strand, but there is no
            # guarantee that it 'points the right way' i.e. from the N-
            # to C- terminus. We test if it does by finding the angle
            # between it and the vector pointing frmo nterm_piont to
            # cterm_point. This angle must either be 0 or pi, if 0 then the
            # axis pionts the right way, otherwise it is pointing the wrong
            # way so we reverse it.

            # FIXME: instead of using Vector.angle() should just compute
            # ((cterm_point-nterm_point) dotproduct dircos)
            #   / (cterm_point-nterm_point).norm() * dircos.norm()
            # and test for 1 or -1, avoiding arccosine computation.
            angle = (self.axis_cterm_point-self.axis_nterm_point).angle(dircos)
            assert (abs(angle) < EPSILON or abs(angle - pi) < EPSILON)
            if (abs(angle - pi) < EPSILON):
                dircos = Vector(-1 * dircos.get_array())
                if verbose:
                    sys.stderr.write('reversed axis for ' + str(self) + '\n')
            
            if mfile_fh != None:
                mfile_write_strand(mfile_fh, c_alpha_veclist,
                                   centroid_mp_veclist, centroid, dircos,
                                   self.axis_nterm_point, self.axis_cterm_point,
                                   str(self.seqnum))
            # DEBUG
            #coords=array([list(vec) for vec in c_alpha_veclist])
            #print str(self),coords,centroid,A,vt,dircos
            # END DEBUG

            self.axis_direction_cosines = dircos
            self.axis_centroid = centroid
            return (dircos, centroid)
        else:
            # When we have fewer than 3 c-alpha atoms, just use vector
            # between the two (if only two) or vector between the two
            # midpoints (if three).
            if len(c_alpha_veclist) == 3:
                mp1 = (c_alpha_veclist[1] - c_alpha_veclist[0]) / 2 + \
                      c_alpha_veclist[0]
                mp2 = (c_alpha_veclist[2] - c_alpha_veclist[1]) / 2 + \
                      c_alpha_veclist[1]
                v = mp2 - mp1
            elif len(c_alpha_veclist) == 2:
                v = c_alpha_veclist[1] - c_alpha_veclist[0]
            else:
                # only 1 residue - cannot fix an axis to this strand
                sys.stderr.write('WARNING: strand ' + str(self) +
                                 ' has only ' + str(len(c_alpha_veclist))
                                 + ' residues, cannot fit axis\n')
                return None
                
            dircos = v / v.norm()

            # get projection of most N-terminal and most C-terminal C_alphas
            # onto the axis line
            self.axis_nterm_point = ProjectPointOntoLine(centroid, centroid +
                                               ALPHA*dircos.get_array(),
                                               c_alpha_veclist[0])
            self.axis_cterm_point = ProjectPointOntoLine(centroid, centroid +
                                               ALPHA*dircos.get_array(),
                                               c_alpha_veclist[-1])

            # DEBUG
            #coords=array([list(vec) for vec in c_alpha_veclist])
            #print str(self),'(short)',coords,centroid,dircos
            # END DEBUG
            if mfile_fh != None:
                mfile_write_strand(mfile_fh, c_alpha_veclist,
                                   None ,centroid, dircos,
                                   self.axis_nterm_point, self.axis_cterm_point,
                                   str(self.seqnum))
            self.axis_direction_cosines = dircos
            self.axis_centroid = centroid
            return (dircos, centroid)


    def get_endpoint_projections(self, pdb_struct):
        """
        Return the projections of the most N-terminal and most C-terminal
        C_alpha midpoints onto the strand axis.

        Parameters:
             pdb_struct - The Bio.PDB parsed PDB struct (atomic co-ordinates)
                          for this protein.

        Return value: tuple (axis_nterm_point, axis_cterm_point) 
                      Bio.PDB Vector of projection of respectively most N-term
                      and most C-term endpoints onto the axis line.
                      
        Uses data members (read/write):
                axis_direction_cosines
                axis_centroid
                axis_nterm_point
                axis_cterm_point
          (basically just calls fit_axis() if necessary to compute the points
          otherwise just returns them, all the above dat members are
          memoizations for fit_axis()).
        """
        if self.axis_nterm_point == None:
            self.fit_axis(pdb_struct)
        return (self.axis_nterm_point, self.axis_cterm_point)
    

    def has_hbond_strand_overlap(self, strand1, strand2):
        """
        For two strands neighbouring (i.e. having bridge edges
        from&to) this one, determine if the H-bonds from those strands
        to this strand are interleaved or nested, indicating that the
        strands must be on opposite sides of this strand.

        Parameters:
            strand1 - a strand with a bridge edge to this strand (i.e.
                      node is in self.bridge_list)
            strand2 - another strand, not strand1, with a bridge to this strand

        Uses data members (readonly):
            hydrogen_bond_list - list of H-bond (node, resnum1, resnum2, dist)
                                  tuples

        Return value:
            True if H bonds from strand1 and strand2 to this strand are
            interleaved. Otherwise False (all bonds from strand1 are
            before (or after) all bonds to strand2 on this strand).

        """
        try:
            min_resnum_to_strand1 = self.min_resnum_hbond_to_neighbour(strand1)[0]
            min_resnum_to_strand2 = self.min_resnum_hbond_to_neighbour(strand2)[0]
            max_resnum_to_strand1 = self.max_resnum_hbond_to_neighbour(strand1)[0]
            max_resnum_to_strand2 = self.max_resnum_hbond_to_neighbour(strand2)[0]
        except ValueError:
            # sometimes get no H bonds even to a neighbour strand (at least
            # with STRIDE), resulting in ValueError min/max arg is an empty
            # sequence.
            # We'll just have to return False and let
            # geometric criteria try to sort it out
            return False

        if verbose:
            sys.stderr.write('test hbond overlap ' + self.nodeid + ': ' +
                             strand1.nodeid + ' ' + str(min_resnum_to_strand1) +
                             ',' + str(max_resnum_to_strand1) + ' ' +
                             strand2.nodeid + ' ' + str(min_resnum_to_strand2) +
                             ',' + str(max_resnum_to_strand2) + '\n')


        if (max_resnum_to_strand1 >= min_resnum_to_strand2 and
            max_resnum_to_strand2 >= min_resnum_to_strand1):
            return True
        else:
            return False


    def min_resnum_hbond_to_neighbour(self, strand):
        """
        Return lowest residue number in this strand that has a hydrogen bond
        to a resdiue in the supplied neighbour (bridge partner) strand.

        Parameters:
            strand - neighbouring strand to find min resnum of the H bonds to

        Return value:
            tuple (resnum_this, resnum_neighbour) where resnum_this is
             smallest residue number in this strand with an H-bond to a residue
             in the supplied neighbour strand.

        Uses data members (readonly):
            hydrogen_bond_list - list of H-bond (node, resnum1, resnum2, dist)

        """
        assert(self.is_neighbour(strand))
        return min([(self.get_residue_ordinal(resnum1),
                     strand.get_residue_ordinal(resnum2))
                    for (node, resnum1, resnum2, dist) in
                    self.get_hbonds_to_neighbour(strand)])


    def max_resnum_hbond_to_neighbour(self, strand):
        """
        Return highest residue number in this strand that has a hydrogen bond
        to a resdiue in the supplied neighbour (bridge partner) strand.

        Parameters:
            strand - neighbouring strand to find min resnum of the H bonds to

        Return value:
            tuple (resnum_this, resnum_neighbour) where resnum_this is
            largest residue number in this strand with an H-bond to a residue
             in the supplied neighbour strand.

        Uses data members (readonly):
            hydrogen_bond_list - list of H-bond (node, resnum1, resnum2, dist)

        """
        assert(self.is_neighbour(strand))
        return max([(self.get_residue_ordinal(resnum1),
                     strand.get_residue_ordinal(resnum2))
                    for (node, resnum1, resnum2, dist) in
                    self.get_hbonds_to_neighbour(strand)])


    def calc_strand_neighbour_occupation(self, strand):
        """
        Return the first and last residues on this strand that are
        'occupied' by the supplied neighbouring strand, i.e. the residues
        that have residues in the neighbour strand beside them when they
        are aligned according to H bonds.

        Parameters:
            strand - neighbour strand to find occupancy of

        Return value:
             tuple (first_resnu, last_resnum) of the lowest and highest
             residue sequence numbers in this strand 'occupied' by the
             neighbour strand

        This is calculated as follows: if residue P is the residue in this
        strand with maximum sequence number with an H bond to a residue Q
        in the neighbour strand, then let a = Q - strand_min and
        b = strand_max - Q where strand_min and strand_max are the min
        and max residue sequence numbers respectively in the neighbouring
        strand. Then the 'occupied' residues in this strand (those with
        neighbours in the neighbouring strand) are those from
        P - b through to P + a (inclusive).


                                  neighbour
                         self      strand
                           .
                          /#\
              P + a -------#--------------------
                           #          #       |
                           #          #       a
                           #P         #Q      |
                    -------*==========*---------
                           #  H bond  #       |
                           #  w/ max  #       |
                           #  res num #       b
                           #  in self #       |
                           #  strand \./      |
              P - b -------#--------------------
                           #
                           #
                           #
                           #

        Note: this does not take into account beta-bulges at all, we are
        considering the represnetatino of strands as perfectly straight and
        lined up beside each other perfectly in the cartoon.

        Note2: The call to max_resnum_hbond_to_neighbour() may generate
        a ValueError exception sometimes when using STRIDE which idnicates
        bridge partners but no H bonds sometimes.

        """
        (P, Q) = self.max_resnum_hbond_to_neighbour(strand)
        a = Q - strand.get_residue_ordinal(strand.get_start_res_seq())
        b = strand.get_residue_ordinal(strand.get_end_res_seq()) - Q
        first_resnum = P - b
        last_resnum = P + a
        assert(last_resnum - first_resnum + 1 == strand.get_span())
        return (first_resnum, last_resnum)


    def has_strand_extent_overlap(self, strand1, strand2):
        """
        For two strands neighbouring (i.e. having bridge edges
        from&to) this one, determine if the strands would overlap
        if they are both aligned according to the H bonds on the same
        side of this strand

        Parameters:
            strand1 - a strand with a bridge edge to this strand (i.e.
                      node is in self.bridge_list)
            strand2 - another strand, not strand1, with a bridge to this strand

        Return value:
            True if the interval of residues on this strand that is
            'occupied' by the residues on strand1 (i.e. are adjacent if
            aligned according to H bonds) overlaps that of strand2,
            else False (the intervals are disjoint and the two strands
            could possibly be drawn on the same side of this strand
            without overlapping).

        """
        try:
            (min_resnum_to_strand1, max_resnum_to_strand1) = \
                    self.calc_strand_neighbour_occupation(strand1)
        except ValueError:
            # sometimes get no H bonds even to a neighbour strand (at least
            # with STRIDE), resulting in ValueError min/max arg is an empty
            # sequence.
            # We'll just have to return False and let
            # geometric criteria try to sort it out
            sys.stderr.write(
                'WARNING: no hbonds available between neighbour strands ' +
                             self.nodeid + ', ' +
                             strand1.nodeid +'\n')
            return False
        try:
            (min_resnum_to_strand2, max_resnum_to_strand2) = \
                    self.calc_strand_neighbour_occupation(strand2)
        except ValueError:
            # sometimes get no H bonds even to a neighbour strand (at least
            # with STRIDE), resulting in ValueError min/max arg is an empty
            # sequence.
            # We'll just have to return False and let
            # geometric criteria try to sort it out
            sys.stderr.write(
                'WARNING: no hbonds available between neighbour strands ' +
                             self.nodeid + ', ' +
                             strand2.nodeid +'\n')
            return False

        if verbose:
            sys.stderr.write('test strand extent overlap '+self.nodeid + ': ' +
                             strand1.nodeid + ' ' + str(min_resnum_to_strand1) +
                             ',' + str(max_resnum_to_strand1) + ' ' +
                             strand2.nodeid + ' ' + str(min_resnum_to_strand2) +
                             ',' + str(max_resnum_to_strand2) + '\n')


        if (max_resnum_to_strand1 >= min_resnum_to_strand2 and
            max_resnum_to_strand2 >= min_resnum_to_strand1):
            return True
        else:
            return False

    def get_hbonds_to_neighbour(self, strand):
        """
        Return list of hydrogen bond tuples (node, resnum1, resnum2, dist) from
        this strand to the the supplied neighbouring strand (i.e. one
        that is in the bridge list of this strand). Because it is a bridge
        partner of this strand it should have H bond(s) to it (or from it).

        Parameters:
            strand - neighbouring strand to find the H bonds to

        Return value:
            list of (node, resnum1, resnum2, dist) tuples of hbonds
            where resnum1 is in self and resnum2 is in strand
            (parameter).

            NB we are storing donor H bonds, not necessarily having
            the symmetical edge from strand back to self, so it
            can happen that there is no hbond from self to strand
            even if there is a bridge between them (though there must
            be a bond from strand to self) - in that case we
            get the bond from strand to self and swap the residue numbers
            in the tuple so resnum1 is the residue number in self
            and resnum2 is the reisude number in the parameter strand.

        Uses data members (readonly):
            hydrogen_bond_list - list of H-bond (node, resnum1, resnum2, dist)

        """
        assert(isinstance(strand, PTNodeStrand))

        hblist = []
        for hbond in self.hydrogen_bond_list:
            if hbond[0] == strand:
                hblist.append(hbond)

        # also add bond(s) from other strand to this one
        for (node, resnum1, resnum2, dist) in strand.get_hbond_list():
            if node == self:
                # note resnum1,resnum1 swapped
                hblist.append((strand, resnum2, resnum1, dist))

        return hblist



    def get_side_of_neighbouring_strand(self, strand):
        """
        Using the labels in the bridge_list set by
        label_strand_bridge_sides() (q.v.), return the same/other (+/-)
        side label for the supplied strand, which is expected to be
        in the bridge list for this node.


        Parameters:
            strand   - strand to find side label for in the bridge list

        Uses data members (readonly):
             bridge_list - list of (node, bdir, side) tuples

        Return value:
            same/other ('+' or '-') side label for the strand if found
            in the bridge list or None if not found.
        """
        assert(isinstance(strand, PTNodeStrand))
        for (node, bdir, side) in self.bridge_list:
            if node == strand:
                return side
        return None


    def is_parallel(self, strand):
        """
        Return True iff the supplied strand is a parallel bridge partner
        of this one.

        Parameters:
            strand - strand to test for being in the bridge list of this node

        Uses data members (readonly):
            bridge_list - list of bridge tuples (ptnodestrand, bdir, side)

        Return value:
            True if strand marked as parallel in the bridge list, else False.

        Raises exceptions:
             KeyError if strand not found in bridge list.
        """
        assert(isinstance(strand, PTNodeStrand))
        for (node, bdir, side_unused) in self.bridge_list:
            if node == strand:
                return (bdir == 'P')
        raise KeyError('is_parallel(): neighbour strand not found')


    def get_is_positioned(self):
        """
        Strands are counted as always positioned, unlike helices,
        so just always return True (This is used only for distance
        matrix placement).
        """
        return True


class PTNodeTerminus(PTNode):
    """
    The PTNodeTerminus class is the type of PTNode for a (N or C) terminus.
    """

    def __init__(self, termtype, pseudo, *args):
        """
        Construct PTNodeTerminus with supplied nodeid and type.
        Parameters:
             termtype - string "N" or "C"
             pseudo - Boolean True for pseudo-terminus (domain boundary),
                       else False (actually N or C terminal of chain).
             +Variable parameter list: straight to PTNode constructor (q.v.).

        This extends PTNode by adding the termtype (note nodeid is more
        informative, includes information about break for domain etc.
        but termtype just N or C is useful for checking in code)
        and other than that just calls PTNode constructor.
        Raises exceptions:
            TypeError if termtype argument is invalid.

        """
        if termtype not in ['N', 'C']:
            raise TypeError("PTNodeTerminus bad termtype " + termtype + '\n')
        PTNode.__init__(self, *args)
        self.termtype = termtype

        self.is_positioned = False # set to True if helix is already positioned
                                   # (drawn). Used when old_helix_placement
                                   # (i.e. -i option NOT supplied) is in use,
                                   # sometimes as a special case we position
                                   # helices before calling write_helices_svg()
                                   # and this flag is set to mark these as
                                   # already positioned.
                                   # use get/set_is_positioned

        self.pseudo = pseudo       # set to True if the terminus is a
                                   # pseudo-terminus to mark domain boundary
                                   # Else False.
                                   # Use get/set_pseudo

        self.adjnode = None        # for pseudo nodes, the node that is
                                   # immediately adjacent, i.e. the most
                                   # N-terminal SSE for pseudo-N-Terminus
                                   # and the most C-terminal SSE for
                                   # pseudo-N-terminus. Else None.
                                   # use get/set_adjnode

        


    def __str__(self):
        """
        Return String representation of the node as 'TYPE id [resnum]'
        """
        return "TERMINUS" + " " + self.nodeid + "[" + str(self.start_res_seq)+"]"

    def get_termtype(self):
        """
        Just return the type of this terminus node, N or C
        Parmeters: None
        Return value: 'N' or 'C' for N or C terminus resp.
        """
        return self.termtype

    def get_is_positioned(self):
        """
        Return True if the node is marked already positioned
        See is_positioned in __init__()
        Parameters: None
        Return value: True if node is marked already positioned else False
        Uses member data (readonly): is_positioned
        """
        return self.is_positioned

    def set_is_positioned(self, is_pos):
        """
        Set the is_positioned flag to the supplied boolean value.
        See is_positioned in __init__()
        Parmeters:
            is_pos - True to mark as already positioned, False to unmark
        Return value: None
        Uses member data (WRITE): is_positioned
        """
        self.is_positioned = is_pos

    def get_pseudo(self):
        """
        Return True if the node is marked as pseudo-terminus (domain boundary)
        See pseudo in __init__()
        Parameters: None
        Return value: True if node is marked as pseudoterminus else False
        Uses member data (readonly): pseudo
        """
        return self.pseudo

    def set_pseudo(self, pseudo):
        """
        Set the psuedo flag to the supplied boolean value.
        See pseudo  in __init__()
        Parmeters:
            pseudo - True to mark as pseudo-terminus, False to unmark
        Return value: None
        Uses member data (WRITE): pseudo
        """
        self.pseudo = pseudo

    def get_adjnode(self):
        """
        Return the adjnode. Only for pseudo nodes.
        See adjnode in __init__()
        Parameters: None
        Return value: The adjacent node to this pseudo node
        Uses member data (readonly): adjnode
        """
        return self.adjnode

    def set_adjnode(self, adjnode):
        """
        Set the adjnode to the supplied PTNode.Only for pseudo nodes.
        See adjnode  in __init__()
        Parmeters:
            adjnode  - The PTNode to set as the adjnode.
        Return value: None
        Uses member data (WRITE): adjnode
        Raises exceptions:
           TypeError if adjnode is not a PTNode instance.
        """
        if not isinstance(adjnode, PTNode):
            raise TypeError('bad adjnode parameter')
        self.adjnode = adjnode

class PTNodeLoop(PTNode):
    """
    The PTNodeLoop class representes loops (coils) rather than SSEs as such.
    Used for domain decomposition (domainfc.py) so that all residues
    are represented in some node.
    """
    def __str__(self):
        """
        Return String representation of the node as
        'TYPE id [startResNum..endResNum]
        """
        return "LOOP" + " " +\
               self.nodeid + "[" + str(self.start_res_seq) + \
               ".." + str(self.end_res_seq)  + "]"


#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

def dfs_strands_from(start_strand, visited, dfs_list, from_node):
    """
    Make a depth-first search traversal of STRAND nodes
    using bridge (not sequence)
    edges starting at the specfieid strand,
    returning list of (node,from_node) tuples in DFS traversal order
    where from_node is the node from which node is reached.

    Parameters:
       start_strand - STRAND node to start at
       visited - (in/out) dictionary of {ptnode:True} visited nodes
       dfs_list - (in/out) list of (node, from_node) visited in dfs order
       from_node - node from which we are being (recursively) called

    Recursive function. call initially as
        dfslist = []
        dfs_strands_from(startnode, {}, dfslist, None)

    Return value:
        None. (output is dfs_list parameter)

    """
    visited[start_strand] = True
    dfs_list.append((start_strand,from_node))
    for (node, bdir_unused, side_unused) in start_strand.get_bridge_list():
        if node not in visited:
            dfs_strands_from(node, visited, dfs_list, start_strand)




def compute_align_positions(node, from_node):
    """
    set relative vertical position based on offsets of
    maximum (or minimum, for reversed) from last in strand
    (or first, for reversed) residue sequence numbers in
    strand of H bonds to previous neighbour strand (starting
    at 0 for first node)
    Called by build_sheet_constraints() in ptgraph2.py, and may be
    called again to recompute align positions if we reverse the order of
    strands in a sheet. (TODO: should be a more efficient way of just
    recalcuating these without calling this again, but since we need the dfs
    order anyway, it does not really matter much).

    Parameters:
        node - not to set align positin in
        from_node - node we reach this from in DFS order
    Return value:
        None.

    Sets the align_pos in node.
    """
    # set relative vertical position based on offsets of
    # maximum (or minimum, for reversed) from last in strand
    # (or first, for reversed) residue sequence numbers in
    # strand of H bonds to previous neighbour strand (starting
    # at 0 for first node)

    # TODO: introduction of get_residue_ordinal() has made his even more
    # overly complicated than it was - since first residue is now always 1
    # using ordinal per-strand numbers rather than pdb sequence numbers,
    # could simplify this eg by not having to get the
    # ordinal for start_res_seq since it is always 1 anyway.;
    # and in fact after changing again so that we use the pdb_resid_dict
    # to get index in full sequence list don't even need get_residue_ordinal()
    # any more... really need to clean all the code up.
    
    if node.get_reversed():
        try:
            (this_bond_resnum, neighbour_bond_resnum) = \
                          node.min_resnum_hbond_to_neighbour(from_node)
        except ValueError:
            # get 'min() arg is an empty sequence' when no H bonds
            # somtimes happens e.g. 1CSE (ptnode.py gives a warning)
            sys.stderr.write('WARNING: no H bonds found between ' +
                             from_node.nodeid + ' and ' +
                             node.nodeid + ', position may be wrong\n')
            this_bond_resnum = node.get_residue_ordinal(node.get_start_res_seq())
            if node.is_parallel(from_node):
                neighbour_bond_resnum = from_node.get_residue_ordinal(from_node.get_start_res_seq())
            else:
                neighbour_bond_resnum = from_node.get_residue_ordinal(from_node.get_end_res_seq())

        this_bond_offset = this_bond_resnum -\
                           node.get_residue_ordinal(node.get_start_res_seq())
    else:
        try:
            (this_bond_resnum, neighbour_bond_resnum) = \
                           node.max_resnum_hbond_to_neighbour(from_node)
        except ValueError:
            sys.stderr.write('WARNING: no H bonds found between ' +
                             from_node.nodeid + ' and ' +
                             node.nodeid + ', position may be wrong\n')
            this_bond_resnum = node.get_residue_ordinal(node.get_end_res_seq())
            neighbour_bond_resnum = from_node.get_residue_ordinal(from_node.get_end_res_seq())
            if node.is_parallel(from_node):
                neighbour_bond_resnum = from_node.get_residue_ordinal(from_node.get_end_res_seq())
            else:
                neighbour_bond_resnum = from_node.get_residue_ordinal(from_node.get_start_res_seq())

        this_bond_offset = node.get_residue_ordinal(node.get_end_res_seq()) -\
                           this_bond_resnum
    if from_node.get_reversed():
        neighbour_bond_offset = neighbour_bond_resnum - \
                                from_node.get_residue_ordinal(from_node.get_start_res_seq())
    else:
        neighbour_bond_offset = from_node.get_residue_ordinal(from_node.get_end_res_seq()) - \
                                neighbour_bond_resnum
    # node_align_pos is relative to neigbour (from_node)
    node_align_pos = neighbour_bond_offset - this_bond_offset
#            print 'rrrrrr',str(node),node_align_pos
    # set the align_pos to pos relative to first (leftmost) node
    # which is to say the cumulative align_pos
    node.set_align_pos(node_align_pos + from_node.get_align_pos())
            

def ptnode_set_verbose(verb):
    """
    set the module global verbose flag in this module to supplied value
    Parameters: verb - True (for verbose output) or False
    Return value: None
    Uses globals: verbose (in this module)
    """
    global verbose
    verbose = verb
