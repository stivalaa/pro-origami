###############################################################################
#
# ptdistmatrix.py - Protein distance matrices
#
# File:    ptdistmatrix.py
# Author:  Alex Stivala
# Created: October 2007
#
# $Id$
#
#
###############################################################################

from sets import Set  # note not using builtin set, so we can use python 2.3

import Bio.PDB
import numpy.oldnumeric as Numeric

from ptnode import *
from ptutils import biopdbresid_to_pdbresseq

# TODO: quite a lot of wasted computation in building the whole
# residue distance matrix then building the SSE distance matrix from
# it - we don't really need the whole residue distance matrix, really
# only the residues that are part of SSEs. So maybe should compute
# only those and save a fair bit of computation and memory.  
# having to do all the mapping from Bio.PDB Residue objects to array
# indices is wasteful/overly complex as well, now that we maintain
# a sequential residue number system in ptgraph2.py (pdb_resid_dict)
# anyway.

#-----------------------------------------------------------------------------
#
# Module globals
#
#-----------------------------------------------------------------------------

#
# global variables
#

# dict of { Residue : bool } just to prevent excessive repeated
# error message in calc_residue_dist()
residue_errmsg_dict = {}

#-----------------------------------------------------------------------------
#
# Class definitions 
#
#-----------------------------------------------------------------------------

class PTDistMatrix:
    """

    PTDistMatrix gives residue and secondary structure element distance
    matrices for a given list of residues and SSEs.
    
    See Peter Cock's Python programming pages for how easy it is to compute
    contact maps with with Bio.PDB and Numeric:

    http://www2.warwick.ac.uk/fac/sci/moac/currentstudents/peter_cock/python/protein_contact_map/

    We don't actually compute contact maps here, but different kinds
    of distance matrices (for residues, and for SSEs), and provide
    some methods on them, specifically to find the closet SSE to a given SSE.

    A third kind of distance map is also used, a sheet/sse distance map,
    where the elements are either sheets or helices. This is also a matrix
    but it is a bit odd in that since each element is the distance between
    (possibly) two sheets or a sheet and a helix, the idea that running
    along the matrix index follows sequence no longer applies (since a sheet
    consists of several strands, which may be all over the place in terms
    of the sequence numbers of the residues that make them up).
    So it is a bit unusual, but we still use a matrix for convenenience here,
    and it is actually used in ptgraph2 to find the closest sheet or helix
    to a given sheet or helix. It is built here from the sse distance
    map, which in turn is built from the residue distance map.
    """

    def __init__(self, residue_list, ptnode_list, sheet_dict, pdb_struct):
        """
        Construct the PTDistMatrix with the supplied list of residues.
        This creates a residue distance matrix and an SSE (PTNode) distance
        matrix. Various methods can then be used to get information from
        these matrices.

        Parameters:
           residue_list - List of Bio.PBD Residue objects
           (ie iterable over Resdiue), from the pdb_struct

           ptnode_list - list of PTNode objects (ie iterable of PTNode)
                         representing the SSEs (helices,strands) we want
                         a distance matrix for.

           sheet_dict - dictionary of {sheet_id : nodelist} representing
                        each sheet as list of strand ptnodes.
                        May be None for no sheet/sse distance map.

           pdb_struct - The Bio.PDB parsed PDB struct (atomic co-ordinates)
                          for this protein.


        """
        self.pdb_struct = pdb_struct
        
        # square symmetric
        # matrix (Numeric.array), of dimensions len(ersidue_list) x
        # len(residue_list) where each element (Numeric.Float) is
        # distance between C-alpha atoms of residues i,j.
        self.dist_matrix = None # build by calc_dist_matrix

        

        #  dict of { residue : array_index } mapping a Bio.PDB
        # residue object to index in dist_matrix
        self.index_map = {}     # built by calc_dist_matrix

        # list of Bio.PDB residue objects mapping
        # index back to Bio.PDB (i.e. reverse of index_map)
        self.reverse_index_map = []  # built by calc_dist_matrix


        # square symmetric Numeric.array matrix of dimensions
        # len(ptnode_list) x len(ptnode_list) where each element
        # is distance between the two SSEs represented by the
        # ptnodes, as defined by calc_sse_dist() (min residue distance)
        self.sse_dist_matrix = None # but by calc_sse_dist_matrix

        # dict of { ptnode : array_index } mapping a PTNode SSE object to
        # index in sse_dist_matrix
        self.sse_index_map = {} # build by calc_sse_dist_matrix

        # list of PTNode objects mapping index back to PTNode i.e.
        # reverseof sse_index_map
        self.reverse_sse_index_map = [] # built by calc_sse_dist-matrix

        #  dict of {(ptnode1, ptnode2) : (residue1, residue2)}
        #  which for every pair of sses gives the residue
        #  in each which are closest (used in the distance
        #  matrix). Note both (ptnode1,ptnode2) and
        #  (ptnode2,ptnode1) are stored, with residues
        #  swapped appropriately.
        #                     
        self.sse_residue_map = {} # built by calc_sse_dist_matrix


        # square symmetric matrix (Numeric.array) of dimensions
        # n x n where n is the number of sheets + number of helices
        # each element is Numeric.Float and is the distance between
        # the two sheets (or two helices, or sheet and helix).
        self.sheet_matrix = None # build by calc_sheet_dist_matrix

        # dict of { id : array_index } where id is either a sheet identifier
        # (a single char 'A' or 'B' etc.) or a helix identifier string
        # as used in ptgraph2 (string like "HELIX_A_1") etc.
        # FIXME: This seems a bit hacky/dangerous, maybe sheet should be
        # a proper object like a PTNode
        self.sheet_index_map = {} #built by calc_sheet_dist_matrix

        # list of identifiers where an identifier is either sheet id e.g.
        # 'A' or helix id e.g. "HELIX_A_1" mapping index back to these
        # identifiers i.e. the inverse of sheet_index_map
        self.reverse_sheet_index_map = [] # build by calc_sheet_dist_matrix

        # dict of {(sheet_id1, id2) : strand_node}
        # which for every sheet id and (for id2) object id (sheet
        # id or ptnodehelix) gives the PTNodeStrand in the sheet
        # that was the closest (use in the sheet_matrix).
        self.sheet_strand_map = {} # built by calc_sheet_dist_matrix

        
        self.calc_dist_matrix(residue_list)
        self.calc_sse_dist_matrix(ptnode_list)
        if sheet_dict != None:
            self.calc_sheet_dist_matrix(sheet_dict, ptnode_list)

    

    def calc_dist_matrix(self, residue_list):
        """
        Compute the matrix of C-alpha distances between residues

        Parameters:
           residue_list - List of Bio.PBD Residue objects

        Return value:
           None. (sets data members)
           

        Uses data members (WRITE):

          dist_matrix - square symmetric
              matrix (Numeric.array), of dimensions len(ersidue_list) x
              len(residue_list) where each element (Numeric.Float) is
              distance between C-alpha atoms of residues i,j.

          index_map - dict of { residue : array_index } mapping a Bio.PDB
               residue object to index in dist_matrix

          reverse_index_map - list of Bio.PDB residue objects mapping
               index back to Bio.PDB (i.e. reverse of index_map)
        
        """
        self.dist_matrix = Numeric.zeros((len(residue_list), len(residue_list)),
                                         Numeric.Float)

        self.reverse_index_map = len(residue_list) * [ -1 ] # will in 0..len-1
        index_maplist = list(enumerate(residue_list))
        for i in range(len(index_maplist)):
            row, residue_one = index_maplist[i]
            self.index_map[residue_one] = row
            self.reverse_index_map[row] = residue_one
            for j in range(i+1, len(index_maplist)):
                col, residue_two = index_maplist[j]
                dist = calc_residue_dist(residue_one, residue_two)
                self.dist_matrix[row, col] = dist
                self.dist_matrix[col, row] = dist

#        print self.dist_matrix
        

    def get_distance(self, residue_one, residue_two):
        """
        Get the distance from residue_one to residue_two, from the
        data members already computed by calc_dist_matrix()
        Note: calc_residue_dist() actually calculates this distance,
        this function retrieves the previously calculated value from
        the matrix.

        Parameters:
           residue_one - Bio.PDB Residue object
           residue_two - Bio.PDB Residue object

        Uses data members (readonly):
           index_map
           dist_matrix

        Return value:
           distance (Angstroms) between residue_one and residue_two C-alphas
        """
        try:
            row = self.index_map[residue_one]
        except KeyError:
            # this happens when domain decomposition has broken an SSE e.g.
            # 1CTN with DDomain and DSSP.
            if not residue_errmsg_dict.has_key(residue_one):
                sys.stderr.write('WARNING: Residue ' +
                                 str(residue_one) +
                                 ' not found,\n  probably due to domain '
                                 'decomposition breaking an SSE.'
                                 '\n  Distance set to infinity\n')
            residue_errmsg_dict[residue_one] = True
            return float("inf")
        try:
            col = self.index_map[residue_two]
        except KeyError:
            if not residue_errmsg_dict.has_key(residue_two):
                sys.stderr.write('WARNING: Residue ' +
                                 str(residue_two) +
                                 ' not found,\n  probably due to domain '
                                 'decomposition breaking an SSE.'
                                 '\n  Distance set to infinity\n')
            residue_errmsg_dict[residue_one] = True
            return float("inf")
            
        dist = self.dist_matrix[row, col]
        return dist
    

    def get_max_distance_residue(self, residue):
        """
        Get the residue with maxmum distance from supplied residue,
        from the data members already computed by calc_dist_matrix()

        Paremeters:
           residue = Bio.PDB residue to get min distance to

        Uses data members (readonly):
           index_map
           reverse_index_map
           dist_matrix

        Return value:
           Bio.PDB residue that has max distance from supplied residue
           
        """
        row = self.index_map[residue]
        maxdist_index = Numeric.argmax(self.dist_matrix[row])
        maxdist_residue = self.reverse_index_map[maxdist_index]
        return maxdist_residue
    

    def calc_sse_dist(self, sse1, sse2):
        """
        Calculate the distance between two SSEs (helices or strands,
        represented by PTNode objects).
        This distance is defined as the smallest distance between any two
        residues, one in each of the SSEs, i.e. the distance betwee the
        two parts of the SSEs tha are closest.

        This is calculated from the residue distance matrix, i.e.
        calc_dist_matrix is assumed to have been already called.

        Parameters:
            sse1 - PTNode representing one SSE
            sse2 - PTNode representing the other SSE

        Uses data members (readonly):
            pdb_struct - The Bio.PDB parsed PDB struct (atomic co-ordinates)
                          for this protein.


        Return value:
            tuple (distance, residue1, residue2) where
              distance is the distance (Angstroms) between the
              closest residues, one from each
              of the two SSEs.
              residue1 is the residue in sse1 used in this min distance
              residue2 is the residue in sse2 used in this min distance
        """
        assert (isinstance(sse1, PTNodeHelix) or isinstance(sse1, PTNodeStrand)
                or isinstance(sse1, PTNodeLoop))
        assert (isinstance(sse2, PTNodeHelix) or isinstance(sse2, PTNodeStrand)
                or isinstance(sse2, PTNodeLoop))
        
        min_dist = float("inf")
        min_res1 = None
        min_res2 = None
        sse1_residues = sse1.get_residue_list()
        sse2_residues = sse2.get_residue_list()
        for res1 in sse1_residues:
            for res2 in sse2_residues:
                dist = self.get_distance(res1, res2)
                if dist < min_dist:
                    min_dist = dist
                    min_res1 = res1
                    min_res2 = res2
        return (min_dist, min_res1, min_res2)
    

                
    def calc_sse_dist_matrix(self, ptnode_list):
        """
        Build the matrix of SSE distance, i.e. min dist between residues
        in the SSEs.

        NOTE: the self-distance (i.e. matrix elements [i,i]) are set to
        infinity rather than 0, so we can efficiently use argmin
        in get_sse_min_distance() to find SSE (not same one) with min
        distance - we are wanting to find minimum distances
        in ptgraph2, not maximum distances.

        Parameters:
          ptnode_list - iterable over PTNode objects represneting SSEs


        Return value:
           None. (sets data members)
           

        Uses data members (WRITE):

          sse_dist_matrix - square symmetric
                            Numeric.array matrix of dimensions
                            len(ptnode_list) x len(ptnode_list) where each
                            elementis distance between the two SSEs
                            represented by the ptnodes,
                            as defined by calc_sse_dist() (min residue distance)

          sse_index_map - dict of { ptnode : array_index } mapping a PTNode
                          object to index in sse_dist_matrix

          reverse_sse_index_map - list of PTNode objects mapping
                                 index back to PTNode
                                 (i.e. reverse of index_map)

          sse_residue_map - dict of {(ptnode1, ptnode2) : (residue1, residue2)}
                             which for every pair of sses gives the residue
                             in each which are closest (used in the distance
                             matrix). Note both (ptnode1,ptnode2) and
                             (ptnode2,ptnode1) are stored, with residues
                             swapped appropriately.
                              
        
        """
        self.sse_dist_matrix =Numeric.zeros((len(ptnode_list),len(ptnode_list)),
                                            Numeric.Float)

        # set the self-distances to infinity (see comments above and in
        # get_sse_min_distance()
        # TODO: maybe if we used NaN instead of inf, this would allow
        # both min/max and argmin/argmax rather than just min/argmin
        # (as we actualy use) to be useful. I tried it with Python 2.5.1
        # on Linux and it worked (ie NaN is neither max nor min) but
        # not really sure how reliable that behaviour is... so sticking
        # with inf for now since we only need min/argmin anyway.
        for i in range(0, Numeric.size(self.sse_dist_matrix,0)):
            self.sse_dist_matrix[i,i] = float("inf")

        self.reverse_sse_index_map = len(ptnode_list) * [ -1 ] #will in 0..len-1
        index_maplist = list(enumerate(ptnode_list))
        for i in range(len(index_maplist)):
            row, sse_one = index_maplist[i]
            self.sse_index_map[sse_one] = row
            self.reverse_sse_index_map[row] = sse_one
            for j in range(i+1, len(index_maplist)):
                col, sse_two = index_maplist[j]
                (dist, res_one, res_two) = self.calc_sse_dist(sse_one, sse_two)
                self.sse_dist_matrix[row, col] = dist
                self.sse_dist_matrix[col, row] = dist
                self.sse_residue_map[sse_one, sse_two] = (res_one, res_two)
                self.sse_residue_map[sse_two, sse_one] = (res_two, res_one)

#        print self.sse_dist_matrix

        
    def get_sse_distance(self, sse1, sse2):
        """
        Get the distance from sse1 to sse2, from the
        data members already computed by calc_sse_dist_matrix()
        Note: calc_sse_dist() actually calculates this distance,
        this function retrieves the previously calculated value from
        the matrix.

        Parameters:
           sse1 - PTNode for an SSE (helix or strand)
           sse2 - PTNode for an SSE (helix or strand)

        Uses data members (readonly):
           sse_index_map
           sse_dist_matrix

        Return value:
           distance (Angstroms) between sse1 and sse1 as defined
           by calc_sse_distance()
        """
        row = self.sse_index_map[sse1]
        col = self.sse_index_map[sse2]
        dist = self.sse_dist_matrix[row, col]
        return dist

    def get_min_distance_sse(self, ptnode):
        """
        Get the SSE with minimum distance from supplied SSE,
        from the data members already computed by calc_sse_dist_matrix()

        Optionally, set the element that was found to infinity so that
        this routine can be used iteratively to find only elements that
        have not already been found.

        Paremeters:
           ptnode - PTNode representing an SSE (helix,strand)

        Uses data members (readonly):
           sse_index_map
           reverse_sse_index_map
           sse_dist_matrix

        Return value:
           PTNode representing SSE that has min distance from supplied SSE
           
        """
        # NB: use of argmin depends on having set diagonal (self distance)
        # elements to inf instead of 0 in calc_sse_dist_matrix().
        row = self.sse_index_map[ptnode]
        mindist_index = Numeric.argmin(self.sse_dist_matrix[row])
        mindist_ptnode = self.reverse_sse_index_map[mindist_index]
        return mindist_ptnode


    def calc_sse_sheet_dist(self, sse, sheet_node_list):
        """
        Calculate the distance between an SSE and a sheet, where a sheet
        is defined by the suplied list of PTNodes representing strands.
        This distance is defined as the smallest distance between the
        supplied sse and any of the strands in the sheet,
        i.e. the distance between a sheet and some other SSE is the distance
        between the SSE and the closest strand to it in the sheet.
        
        This is calculated from the SSE distance matrix, i.e.
        calc_sse_dist_matrix is assumed to have been already called.

        Parameters:
            sse - PTNode representing one SSE
            sheet_node_list - list of PTNodes of strands in the sheet.

        Return value:
            tuple (dist, strand) where
            dist is the
            distance (Angstroms) between the closest residues, one from each
            of the two SSEs and
            strand is PTNodeStrand of the strand in the sheet used for this
            minimum distance
        """
        assert isinstance(sse, PTNodeHelix) or isinstance(sse, PTNodeStrand)
        
        min_dist = float("inf")
        min_dist_strand = None
        for strand in sheet_node_list:
            assert isinstance(strand, PTNodeStrand)
            dist = self.get_sse_distance(sse, strand)
            if dist < min_dist:
                min_dist = dist
                min_dist_strand = strand
        return (min_dist, min_dist_strand)
    
    def calc_sheet_sheet_dist(self, sheet1, sheet2):
        """
        Calculate the distance between two sheets, where the sheets
        are represnetd by lists of PTNodes represneting the strands in them.
        This distance is defined as the distance between the two closest
        strands, one from each sheet.

        This is calculated from the SSE distance matrix, i.e.
        calc_sse_dist_matrix is assume already to have been called.

        Parameters:
           sheet1 - list of PTNodes of strands in the sheet
           sheet2 - list of PTNodes of strands in the other sheet

        Return value:
           tuple (dist, strand1, strand2) where dist is the
           distance (Angstroms) between the two sheets, as defined above
           and strand1 is the strand in sheet1 that was used in this
           minimum distance calculation, and strand2 is that in sheet2
        """
        min_dist = float("inf")
        min_dist_strand1 = None
        min_dist_strand2 = None
        for strand in sheet1:
            assert isinstance(strand, PTNodeStrand)
            (dist, strand2) = self.calc_sse_sheet_dist(strand, sheet2)
            if dist < min_dist:
                min_dist = dist
                min_dist_strand1 = strand
                min_dist_strand2 = strand2
        return (min_dist, min_dist_strand1, min_dist_strand2)


    def calc_sheet_dist_matrix(self, sheet_dict, ptnode_list):
        """
        Build the matrix of distances between all the helices and sheets.

        NOTE: the self-distance (i.e. matrix elements [i,i]) are set to
        infinity rather than 0, so we can efficiently use argmin
        in get_sse_min_distance() to find SSE (not same one) with min
        distance - we are wanting to find minimum distances
        in ptgraph2, not maximum distances.

        Parameters:
          ptnode_list - iterable over PTNode objects represneting SSEs,
                        we only use the helices though.

          sheet_dict - dictionary of {sheet_id : nodelist} representing
                        each sheet as list of strand ptnodes.

        Uses data members (WRITE):
            sheet_matrix -
               square symmetric matrix (Numeric.array) of dimensions
               n x n where n is the number of sheets + number of helices
               each element is Numeric.Float and is the distance between
               the two sheets (or two helices, or sheet and helix).
        
            sheet_index_map - 
               dict of { id : array_index } where id is either a sheet id
               (a single char 'A' or 'B' etc.) or a helix identifier string
               as used in ptgraph2 (string like "HELIX_A_1") etc.
               FIXME: This seems a bit hacky/dangerous, maybe sheet should be
               a proper object like a PTNode

            reverse_sheet_index_map - 
              list of identifiers where an identifier is either sheet id e.g.
              'A' or helix id e.g. "HELIX_A_1" mapping index back to these
              identifiers i.e. the inverse of sheet_index_map

            sheet_strand_map - 
              dict of {(sheet_id1, id2) : strand_node}
              which for every sheet id and (for id2) object id (sheet
              id or ptnodehelix) gives the PTNodeStrand in the sheet
              that was the closest (use in the sheet_matrix).


        Return value:
           None. (sets data members)
        
        """
        helix_list = [ ptnode for ptnode in ptnode_list if
                       isinstance(ptnode, PTNodeHelix) ]
        sheet_id_list = sheet_dict.keys()
        objlist = helix_list + sheet_id_list
        n = len(objlist) # n = number of helices + number of sheets
        self.sheet_dist_matrix = Numeric.zeros((n, n), Numeric.Float)

        # set the self-distances to infinity (see comments above and in
        # get_sse_min_distance()
        # TODO: maybe if we used NaN instead of inf, this would allow
        # both min/max and argmin/argmax rather than just min/argmin
        # (as we actualy use) to be useful. I tried it with Python 2.5.1
        # on Linux and it worked (ie NaN is neither max nor min) but
        # not really sure how reliable that behaviour is... so sticking
        # with inf for now since we only need min/argmin anyway.
        for i in range(0, n):
            self.sheet_dist_matrix[i,i] = float("inf")

        self.reverse_sheet_index_map = n * [ -1 ] # will be in 0..n-1
        index_maplist = list(enumerate(objlist))
        for i in range(len(index_maplist)):
            row, obj1 = index_maplist[i]
            if isinstance(obj1, PTNode):
                obj1_id = obj1.nodeid
            else:
                obj1_id = obj1 # it is a sheet id e.g. 'A'
                assert(obj1.isalpha())
            self.sheet_index_map[obj1_id] = row
            self.reverse_sheet_index_map[row] = obj1_id
            for j in range(i+1, len(index_maplist)):
                col, obj2 = index_maplist[j]
                if isinstance(obj2, PTNode):
                    obj2_id = obj2.nodeid
                else:
                    obj2_id = obj2 # it is a sheet id e.g. 'A'
                    assert(obj2.isalpha())
                if isinstance(obj1, PTNode) and isinstance(obj2, PTNode):
                    # both are helices
                    dist = self.get_sse_distance(obj1, obj2)
                elif isinstance(obj1, PTNode):
                    # obj1 is a helix, obj2 is a sheet
                    (dist, strand) = \
                           self.calc_sse_sheet_dist(obj1, sheet_dict[obj2_id])
                    self.sheet_strand_map[(obj2_id, obj1_id)] = strand
                elif isinstance(obj2, PTNode):
                     # obj1 is a sheet, obj2 is a helix
                    (dist, strand) = \
                           self.calc_sse_sheet_dist(obj2, sheet_dict[obj1_id])
                    self.sheet_strand_map[(obj1_id, obj2_id)] = strand
                else:
                    # both are sheets
                    (dist, strand1, strand2) = \
                           self.calc_sheet_sheet_dist(sheet_dict[obj1_id],
                                                      sheet_dict[obj2_id])
                    self.sheet_strand_map[(obj1_id, obj2_id)] = strand1
                    self.sheet_strand_map[(obj2_id, obj1_id)] = strand2
                    
                self.sheet_dist_matrix[row, col] = dist
                self.sheet_dist_matrix[col, row] = dist

#        print self.sheet_strand_map
#        print self.reverse_sheet_index_map
#        print self.sheet_dist_matrix


    def get_min_distance_objid(self, objid, not_objid_set, sheets_only=False):
        """
        Get the sheet or helix with minimum distance from the supplied
        sheet or helix, specified by id (e.g. 'A' for a sheet or
        'HELIX_A_10' for a helix).

        Optionally, set the element that was found to infinity so that
        this routine can be used iteratively to find only elements that
        have not already been found.

        Paremeters:
           objid  - sheet id (e.g. 'A') or helix id (e.g. 'HELIX_A_10')
                    of the object to find the id of the closest object for.
           not_objid_set - set of objids that we do NOT want to find.
                    Used so we can find the nearest element to an
                    already positioned element that is not itself
                    an already positioned element.
           sheets_only - (Default False) only find sheets, not helices.


        Uses data members (readonly):
           sheet_index_map
           reverse_sheet_index_map
           sheet_dist_matrix  

        Return value:
           tuple (id, dist) where
           id (as per the objid paramter) of the closest sheet or helix
            to the speicfied one and
            dist is that smallest distance, and it is not in the
            not_objid_set.
           
        """

        row = self.sheet_index_map[objid]
        mindist_index = Numeric.argmin(self.sheet_dist_matrix[row])

        # get 1d array of object ids sorted (ascending) by their distance
        # from the target objid in the sheet dist matrix
        # NB: use of argsort depends on having set diagonal (self distance)
        # elements to inf instead of 0 in calc_sse_dist_matrix().
        objids_sorted_by_dist = Numeric.argsort(self.sheet_dist_matrix[row])

        # find the first (i.e. smallest distance) id that is not in
        # the not_objid_set
        mindist_index = None
        for mindist_index in objids_sorted_by_dist:
            mindist_objid = self.reverse_sheet_index_map[mindist_index]
            if ( (mindist_objid not in not_objid_set) and
                 (not sheets_only or len(mindist_objid) == 1) ): 
                dist = self.sheet_dist_matrix[row, mindist_index]
                break
        return (mindist_objid, dist)


    def get_strand_nearest_element(self, sheet_id, element):
        """
        Get the PTNodeStrand in the supplied sheet that is nearest
        to the supplied object (specified by object id). Uses the
        sheet_strand_map built by calc_sheet_dist_matrix() to do this;
        the idea is that nearest objects are found with
        get_min_distance_objid(), and in the case that an element is
        nearest to a sheet, this function is then called to find
        the paritcular strand in that sheet that was used as the nearest
        element.

        Parameters:
           sheet_id - id of sheet to find strand nearest objid
           element - element (sheet id or PTNodeHelix) to find the
                      strand in the sheet nearest to.

        Return value:
        tuple (strand1, strand2) where strand 1 is the
           PTNodeStrand for the strand in the sheet that is nearest to element
           and strand2 is the PTNodeStrand for the strand in element that
           was closest, if element is a sheet id, or None if element is a helix.
           
        Uses data members (readonly):
            sheet_strand_map - 
              dict of {(sheet_id1, id2) : strand_node}
              which for every sheet id and (for id2) object id (sheet
              id or ptnodehelix) gives the PTNodeStrand in the sheet
              that was the closest (use in the sheet_matrix).
        """
        if isinstance(element, PTNodeHelix):
            element_objid = element.nodeid
        else:
            element_objid = element # objid of sheet id is just sheet id
        strand1 = self.sheet_strand_map[(sheet_id, element_objid)]
        assert isinstance(strand1, PTNodeStrand)
        assert strand1.get_sheet_id() == sheet_id
        if isinstance(element, PTNodeHelix):
            strand2 = None
        else:
            strand2 = self.sheet_strand_map[(element_objid, sheet_id)]
            assert isinstance(strand2, PTNodeStrand)
            assert strand2.get_sheet_id() == element_objid
#        print 'zzz',sheet_id, str(element),str(strand1),str(strand2)
        return (strand1, strand2)
    

    def get_nearest_sse_residues(self, sse1, sse2):
        """
        Find the residue in each of the two SSEs that are nearest to each
        other and were used in building the SSE distance matrix.
        Uses the sse_residue_map built by calc_sse_dist_matrix() to do this;
        the idea is that nearest SSEs are found with get_min_distance_sse()
        or other functions using the SSE distance matrix, then if required
        this functino is used to retrieve the particular residues that
        were used in calculating the min distance between SSEs.
        
        Parameters:
           sse1 - PTNode for helix/strand 1 
           sse2 - PTNode for helix/strand 2

        Return value:
           tuple (res_seq_num_1, res_seq_num_2) where res_seq_num_1 and
           res_seq_num_2 are the residue sequence numbers in sse1 and sse2
           respectively that have min distance to each other (of all
           residues in sse1 and sse2)

        Uses data members:
            sse_residue_map - 
              dict of {(ptnode1, ptnode2) : (residue1, residue2)}
              which for every pair of sses gives the residue
              in each which are closest (used in the distance
              matrix). Note both (ptnode1,ptnode2) and
              (ptnode2,ptnode1) are stored, with residues
              swapped appropriately.
                             
        """
        (residue1, residue2) = self.sse_residue_map[sse1, sse2]
        # id of a residue in Bio.PDB is tuple (hetatm, resseqnum, icode)
        res_seq_num_1 = biopdbresid_to_pdbresseq(residue1.get_id())
        res_seq_num_2 = biopdbresid_to_pdbresseq(residue2.get_id())
        return (res_seq_num_1, res_seq_num_2)
    
    
    def find_nearest_sses_in_sets(self, ref_set, test_set):
        """
        Find the nearest sse (helix or strand)  in test_set to any
        of the sses in the supplied ref_set.

        Parameters:
           ref_set - set of PTNodes to find the nearest
                         element to any of them, that is in the test_test
           test_set - set of PTNodes to find neareset
                      from.
        Uses data members (read):
           distmatrix - The PTDistMatrix that has been built already
        Return value:
           tuple (set_element, close_element, dist)
           where set_element is an element in the supplied ref_set and
           close_element is the
           PTNode element
           in the test_test which has minimum distance
           to set_element
           and dist is the distance between the two.
           
        """
        min_dist = float("inf")
        close_element = None
        set_element = None

        for ref_node in ref_set:
            for test_node in test_set:
                this_dist = self.get_sse_distance(ref_node, test_node)
                if this_dist < min_dist:
                    set_element = ref_node
                    close_element = test_node
                    min_dist = this_dist

        return (set_element, close_element, min_dist)


    def find_nearby_sses_in_sets(self, ref_set, test_set, dist_threshold):
        """
        Find all sses (helices or strands) in test_set whose distance
        from any of the sses in the supplied ref_set is below a
        threshold.

        Parameters:
           ref_set - set of PTNodes to find nearby
                         elements to any of them, that is in the test_test
           test_set - set of PTNodes to find nearby sses
                      from.
           dist_threshold - threshold (Angstroms) below which SSEs are 'nearby'
        Uses data members (read):
           distmatrix - The PTDistMatrix that has been built already
        Return value:
           List of tuples (dist, sse) of SSEs in test_set that are less
           than dist_treshold from some SSE in ref_set, and the distance
           (dist) that each is from its closest SSE in the ref_set
           (Note dist is before sse in tuple to make sorting by dist easy).
           
        """
        close_dict = {} # dict of {node : dist} for min dist to test node
        for ref_node in ref_set:
            for test_node in test_set:
                this_dist = self.get_sse_distance(ref_node, test_node)
                if this_dist < dist_threshold:
                    if (not close_dict.has_key(test_node) or
                        this_dist < close_dict[test_node]):
                        close_dict[test_node] = this_dist
        # items() converts dict to list of (node,dist) tuples then
        # swap each tuple so we have list of (dist,node) tuples for ease
        # of sorting
        return [ (dist, node) for (node, dist) in close_dict.items() ]



    
#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues

    Paramters
       residue_one - Bio.PDB Residue object
       residue_two - Bio.PDB Residue object

    Return value:
       distance in Angstroms between Carbon alpha atoms of residue_one and
       residue_two

    Uses globals (read/write):
       residue_errmsg_dict - map Residue to bool flagging error msg issued

    Based on Peter Cock's Python programming pages:

    http://www2.warwick.ac.uk/fac/sci/moac/currentstudents/peter_cock/python/protein_contact_map/

    but of course life is never quite that simple with PDB files...
    note this function is almost entirely error handling, only two
    lines actually do the calculation, everything else handles exceptions.
    """
    try:
        res1_ca_coord = residue_one["CA"].coord
    except KeyError: # this happens ocassionaly on some PDB files e.g. 1BRD
        if not residue_errmsg_dict.has_key(residue_one):
            sys.stderr.write('WARNING: no Carbon-alpha atom for residue ' +
                             str(residue_one) +
                             '\nDistance in matrix set to infinity\n')
            residue_errmsg_dict[residue_one] = True
        return float('inf')
    try:
        res2_ca_coord = residue_two["CA"].coord
    except KeyError:
        if not residue_errmsg_dict.has_key(residue_two):
            sys.stderr.write('WARNING: no Carbon-alpha atom for residue ' +
                             str(residue_two) +
                             '. Distance in matrix set to infinity\n')
            residue_errmsg_dict[residue_two] = True
        return float('inf')

    diff_vector = res1_ca_coord - res2_ca_coord
    return Numeric.sqrt(Numeric.sum(diff_vector * diff_vector))


def calc_point_residue_dist(residue_one, point) :
    """Returns the distance between the C_alpha of a residue and
       a point.

    Paramters
       residue_one - Bio.PDB Residue object
       point - Bio.PDB Vector representatin of a point in 3d space

    Return value:
       distance in Angstroms between Carbon alpha atom of residue_one and
       point

    Uses globals (read/write):
       residue_errmsg_dict - map Residue to bool flagging error msg issued
    """
    try:
        res1_ca_coord = residue_one["CA"].coord
    except KeyError: # this happens ocassionaly on some PDB files e.g. 1BRD
        if not residue_errmsg_dict.has_key(residue_one):
            sys.stderr.write('WARNING: no Carbon-alpha atom for residue ' +
                             str(residue_one) +
                             '\nDistance in matrix set to infinity\n')
            residue_errmsg_dict[residue_one] = True
        return float('inf')
    diff_vector = res1_ca_coord - point.get_array()
#        print 'debug pointresdist',res1_ca_coord,point.get_array(),diff_vector
    return Numeric.sqrt(Numeric.sum(diff_vector * diff_vector))


def calc_sse_sse_dist(sse1, sse2, pdb_struct):
    """
    Calculate the distance between two SSEs (helices or strands,
    represented by PTNode objects).
    This distance is defined as the smallest distance between any two
    residues, one in each of the SSEs, i.e. the distance betwee the
    two parts of the SSEs tha are closest.

    Parameters:
        sse1 - PTNode representing one SSE
        sse2 - PTNode representing the other SSE
        pdb_struct - The Bio.PDB parsed PDB struct (atomic co-ordinates)
                      for this protein.


    Return value:
        tuple (distance, residue1, residue2) where
          distance is the distance (Angstroms) between the
          closest residues, one from each
          of the two SSEs.
          residue1 is the residue in sse1 used in this min distance
          residue2 is the residue in sse2 used in this min distance
    """
    assert (isinstance(sse1, PTNodeHelix) or isinstance(sse1, PTNodeStrand)
            or isinstance(sse1, PTNodeLoop))
    assert (isinstance(sse2, PTNodeHelix) or isinstance(sse2, PTNodeStrand)
            or isinstance(sse2, PTNodeLoop))

    min_dist = float("inf")
    min_res1 = None
    min_res2 = None
    sse1_residues = sse1.get_residue_list()
    sse2_residues = sse2.get_residue_list()
    for res1 in sse1_residues:
        for res2 in sse2_residues:
            dist = calc_residue_dist(res1, res2)
            if dist < min_dist:
                min_dist = dist
                min_res1 = res1
                min_res2 = res2
    return (min_dist, min_res1, min_res2)


def calc_sse_sse_midpoint_dist(sse1, sse2, pdb_struct):
    """
    Calculate the midpoint distance between two SSEs (helices or strands,
    represented by PTNode objects).
    This distance is defined as the distance between the midpoints of
    the two SSE axes (as calculated by fit_axis() methods).

    Parameters:
        sse1 - PTNode representing one SSE
        sse2 - PTNode representing the other SSE
        pdb_struct - The Bio.PDB parsed PDB struct (atomic co-ordinates)
                      for this protein.


    Return value:
          distance (Angstroms) between the midpoints of the axes fitted
          to each SSE; or None if no axis could be found.
    """
    sse1_axis = sse1.fit_axis(pdb_struct)
    sse2_axis = sse2.fit_axis(pdb_struct)
    if sse1_axis == None or sse2_axis == None:
        return None

    (sse1_dircos, sse1_centroid) = sse1_axis
    (sse2_dircos, sse2_centroid) = sse2_axis

    diff_vector = sse1_centroid - sse2_centroid
    return Numeric.sqrt(Numeric.sum(diff_vector * diff_vector))


def compute_sse_midpoint_dist_matrix(ptnode_list, pdb_structure):
    """
    Return a distance matrix between midpoints of each SSE
    by computing axis midpoint distances between all SSEs in the ptnode_list

    Parameters:
        ptnode_list - list of PTNode objects (ie iterable of PTNode)
                         representing the SSEs (helices,strands) the
                         SSE midpoint distance matrix is for.
        pdb_structure - parsed Bio.PDB structure

    Return value:
       Numeric.array square symmetric (order length of ptnode_list) where
       each entry is distance (Angstroms) between midpoints of SSE axes.
       Main diagonal entries set to SSE type (0 strand, 1 alpha helix,
       2 pi helix, 3 3_10 helix).
    
    """
    n = len(ptnode_list)
    dist_array = Numeric.zeros((n, n), Numeric.Float)
    for i in range(n):
        for j in range(i+1, n):
            dist = calc_sse_sse_midpoint_dist(ptnode_list[i], ptnode_list[j],
                                              pdb_structure)
            if dist == None:
                dist_array[i, j] = float('NaN')
            else:
                dist_array[i, j] = dist
            dist_array[j, i] = dist_array[i, j]

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
        dist_array[i,i] = v

    return dist_array


