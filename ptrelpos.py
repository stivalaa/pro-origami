###############################################################################
#
# ptrelpos.py - find relative positions to place protein cartoon elements
#
# File:    ptrelpos.py
# Author:  Alex Stivala
# Created: October 2007
#
# $Id$
#
#
###############################################################################

import Bio.PDB

from ptnode import *
from ptdistmatrix import PTDistMatrix, calc_residue_dist
from pttableau import PTTableau


# TODO: have use_longest_for_orientation to choose to use longest rather
# than nearest strand in sheet for orientation. Currently using nearest
# (using longest makes 2PEE-3 different from other serpins (1QLP, 1MTP, etc.)
# for example, due to very bent longest strand in large sheet).
# Maybe should use strand with best fitting axis for oriention instead?

#-----------------------------------------------------------------------------
#
# Module globals
#
#-----------------------------------------------------------------------------


# constants

RELPOS_ABOVE = 0
RELPOS_BELOW = 1
RELPOS_LEFT  = 2
RELPOS_RIGHT = 3

# global variables

verbose = False



#-----------------------------------------------------------------------------
#
# Class definitions 
#
#-----------------------------------------------------------------------------

class PTRelativePosition:
    """
    PTRelativePosition is a class for finding the relative position
    of SSEs to each other for laying them out in the cartoon, using information
    from the PDB structure and the distance matrices and sheet (strand
    position) information that has already been determined.
    """

    def __init__(self, pdb_struct, distmatrix, sheet_strandlists_dict, tableau,
                 chain_dict, sheet_dict):
        """

        Parameters:
           pdb_struct - The Bio.PDB parsed PDB struct (atomic co-ordinates)
                          for this protein.
           distmatrix - The PTDistMatrix distance matrices for this protein.
           sheet_strandlists_dict  -
                dictionary of { sheet_id : list of list of nodes }
                where the list of list of nodes is
                described in build_sheet_constraints()
           tableau - the PTTableau which has been built for this protein
           chain_dict -   Each value of the chain_dict is a
                          List of nodes in order from N to C terminus
                          so chain_dict is { chainid : node_list }
           sheet_dict - dict of {sheet_id : ptnode_list} represneting sheets
        """
        self.pdb_struct = pdb_struct
        self.distmatrix = distmatrix
        self.sheet_strandlists_dict = sheet_strandlists_dict
        self.tableau = tableau
        self.chain_dict = chain_dict
        self.sheet_dict = sheet_dict


    def get_strand_posnum(self, strand, sheet_strandlists_dict = None):
        """
        Return the index of the supplied strand in its sheet
        in the outermost ('horizontal') list i.e. the number of strands
        it is from the 'leftmost' strand.

        Parameters:
           strand   - PTNode strand to find position number of
           sheet_strandlists_dict - the sheet strandlists dict to use
                      for this strand. Default None. If None, use
                      the data member sheet_strandlists_dict
                      (This is to enable this function to be used
                      for other domains, not the one this object is for).
                      The strand has to belong to the same domain as
                      the sheet_strandlists_dict, otherwise this makes no
                      sense.
           
        Uses data members (read):
             sheet_strandlists_dict  -
                dictionary of { sheet_id : list of list of nodes }
                where the list of list of nodes is
                described in build_sheet_constraints()

        Return value - index in outermost list of entry for this sheet id
                       that the strand is in.
        """
        assert(isinstance(strand, PTNodeStrand))
        sheet_id = strand.get_sheet_id()
        if sheet_strandlists_dict != None:
            ssd = sheet_strandlists_dict
        else:
            ssd = self.sheet_strandlists_dict
        horiz_order_list = ssd[sheet_id]
        for posnum in range(len(horiz_order_list)):
            if strand in horiz_order_list[posnum]:
                return posnum
        assert(False) # strand must be in its own sheet somewhere
            

    def any_strands_before_or_after_strand(self, strand1, strandlist):
        """
        Return True if any strand in strandlist
        immediately follows or precedes strand1 in sequence,
        i.e. is some strand in strandlist
        is the SSE immeidately C-terminal or N-terminal of strand1
        in the same chain.

        Parameters:
           strand1 - PTNodeStrand of strand to check if any strand is after
           strandlist - list of PTNodeStrand to check if any of them
                         immediately follow strand1 in sequence

        Return value:
           True if some strand in strandlist
           is immediately C-terminal or N-terminal of strand1 in chain
           else False

        Uses data members (Readonly):
           chain_dict

        Note index() raises ValueError exception if strand is not
        found in the list of nodes for its chain, which should never
        happen (ie if this exception is raised there is some internal
        inconsistency in the chain dict or strand structure).
        """
        assert(isinstance(strand1, PTNodeStrand))
        chainid = strand1.get_chainid()
        nodelist = self.chain_dict[chainid]
        # FIXME index() is probably a linear search, should
        # maybe build some dictionary to do this faster, but doesn't
        # really matter that much (probably)
        strand1_index = nodelist.index(strand1)
        next_index = strand1_index + 1
        prev_index = strand1_index - 1
        if next_index >= len(nodelist) and prev_index < 0:
            return False
        if next_index < len(nodelist):
            nextnode = nodelist[next_index]
        else:
            nextnode = None
        if prev_index >= 0:
            prevnode = nodelist[prev_index]
        else:
            prevnode = None
        if (not isinstance(nextnode, PTNodeStrand) and
            not isinstance(prevnode, PTNodeStrand)):
            return False
        for strand2 in strandlist:
            if strand2.get_chainid() != chainid:
                continue
            if nextnode == strand2 or prevnode == strand2:
                return True
        return False


    def get_longest_strand(self, horiz_order_list):
        """
        Return the strand and its length (as number of residues)
        of the longest starnd in the sheet specified
        by its list of list of strands (horizontal outer list, each
        elment list aligned vertically).

        Parameters:
           horiz_order_list - the sheet strand list for the sheet as built by
                               build_sheet_constraints
        Return value:
           tuple (ptnodestrand, length)
           where ptnodestrand is PTNodeStrand of longest strand and length is
           number of residues in longest strand in the sheet.

        Uses no data members.
        """
        longest_strand_length = 0
        longest_strand = None
        for vert_list in horiz_order_list:
            if len(vert_list) < 1:
                continue # should not happen anyway
            strand = vert_list[0] # NOTE: assumes longest is single
            # FIXME: this assumption is not always good, e.g. 1W81
            # where 16 and 12 are on same vert axis both neighbours
            # of 11, and 16 is about same length as 11 so 16 and 12 
            # together definitely longer than 11 causing overlap on figure
            if strand.get_span() > longest_strand_length:
                longest_strand_length = strand.get_span()
                longest_strand = strand
        return (longest_strand, longest_strand_length)


    def flip_all_strands_in_sheet(self, sheet_id):
        """
        Turn the sheet 'upside-down'.
        Flip the reverse flag in each strand of the sheet i.e. set if not
        set, unset if set. Initially (in build_sheet_constraints(), these
        flags are set based on the first ('leftmost') strand being set as
        not-reversed, but after we find orientations we may actually want
        that strand the other way, so we just flip all the reversed flags.

        Not only do we flip the reverse flag, we also have to
        shift the align position as it was calculated (in
        bulid_sheet_constraints()) with the reverse flag as it was before
        (obviously). So now the offset is changed to be relative to the
        other (i.e. after reversing) end of the strand, and no special
        case is needed for reversed sideways strands when laying out the sheet
        for the diagram.
        
        Parmeters:
           sheet_id - id of the sheet to flip reverse flags in

        Return value: None
        
        Uses data members (read/write):
            sheet_strandlists_dict  -
                dictionary of { sheet_id : list of list of nodes }
                where the list of list of nodes is
                described in build_sheet_constraints()
            
        """
        # align positions were relative to the (original) top of this strand
        first_strand_len= self.sheet_strandlists_dict[sheet_id][0][0].get_span()
        
        for strandlist in self.sheet_strandlists_dict[sheet_id]:
            for strand in strandlist:
                strand.set_reversed(not strand.get_reversed())
                # now make the align position relative to the other end
                strand.set_align_pos(first_strand_len - strand.get_align_pos()
                                     - strand.get_span())


    def reverse_strand_order_in_sheet(self, sheet_id, sheet_strandlists_dict):
        """
        Flip the sheet left-to-right.
        Reverse the order of the strands in the sheet.
        Not only do we need to rervese the horiz order list of strands,
        we also have to adjust the align positions as calculaated in
        build_sheet_constraints() accordingly. These offsets were relative
        to the first strand in the list (which itself is offset 0), now
        that strand is the last so we need to adjust them all so new
        first is offset 0 and others relative to it.
        This is not as easy as going through the horiz_order_list because
        of bifurcations and the order of offsets being added is the dfs
        order used in the original build_sheet_constraints(), so we
        recompute the align positions from scratch.
        (TODO: should be a more efficient way of just
        recalcuating these without calling
        compute_align_positions() again to do it from scratch,
        but since we need the dfs
        order anyway, it does not really matter much).

        Parameters:
           sheet_id - id of sheet to reverse
           sheet_strandlists_dict -  IN/OUT
                                    the sheet strandlists dict that contains
                                    the sheet identified by sheet_id
        Return value:
           None

        Uses data members:
            None

        Note strand nodes are also modified (the align_pos value), only
        nodes that are in the sheet are referenced.
        """
        # first recompute the relative align positions
        start_node = sheet_strandlists_dict[sheet_id][-1][0] # start at end
        dfs_list = []
        dfs_strands_from(start_node, {}, dfs_list, None)
        assert(start_node == dfs_list[0][0] and dfs_list[0][1] == None)
        start_node.set_align_pos(0)
        for (node, from_node) in dfs_list[1:]:
            compute_align_positions(node, from_node)

        # now reverse the list
        sheet_strandlists_dict[sheet_id].reverse()

        
    def set_all_sheet_strands_sideways(self, sheet_id):
        """
        Set the sideways flag in every strand of a sheet.

        Parameters:
           sheet_id - id of the sheet to set sideways flags in

        Return value: None

        Uses data members (read/write):
            sheet_strandlists_dict  -
                dictionary of { sheet_id : list of list of nodes }
                where the list of list of nodes is
                described in build_sheet_constraints()
        """
        for strandlist in self.sheet_strandlists_dict[sheet_id]:
            for strand in strandlist:
                strand.set_sideways(True)
          

    def get_relative_position(self, reference_element, test_element):
        """
        Find the relative position of test_element relative to
        reference_element.

        Parameters:
            reference_element - an element (either sheet id e.g. 'A' or
                   helix (PTNodeHelix object) to find position relative to
            test_element - and element (as per reference_element) to find
                   position of relative to reference_element

           NOTE: reversed and sideways flags in test_element may be set
           by this function.
           
        Uses data members (read):
            distmatrix - the distance matrix
            sheet_strandlists_dict  -
                dictionary of { sheet_id : list of list of nodes }
                where the list of list of nodes is
                described in build_sheet_constraints()
            
        Return value:
            tuple (relpos, ref_strand, test_strand) where relpos is
            RELPOS_ABOVE, RELPOS_BELOW, RELPOS_LEFT or RELPOS_RIGHT
            ref_strand is strand in reference sheet it is relative to
            or None if reference element is not a sheet
            test_strand is strand in test sheet that is relative to
            reference_element or None if test element is not a sheet
        """
        assert(isinstance(reference_element, PTNodeHelix) or
               len(reference_element) == 1) # sheet id
        assert(isinstance(test_element, PTNodeHelix) or
               len(test_element) == 1) # sheet id

        if isinstance(reference_element, PTNodeHelix):
            ref_strand = None
            (relpos, ref_sse, test_strand) = \
                     self.get_relpos_to_helix(reference_element, test_element)
        else:
            (relpos, ref_sse, test_strand) = \
             self.get_relpos_to_sheet(reference_element, test_element)


        return (relpos, ref_sse, test_strand)



    def get_relpos_helix_to_helix(self, reference_helix, test_element,
                                  nearest_ref_resnum, nearest_test_resnum):
        """
        Find the relative position of a helix to a helix.

        Parameters:
           reference_helix - helix to place relative to
           test_element - helix to place relative to reference helix
           nearest_ref_resnum - residue number in reference_helix that test
                               helix is closest to
           nearest_test_resnum - residue number in test_helix that is closest
                                to reference helix

        Return value:
           relpos of test to reerence helix

        Uses no data members
        """
        if reference_helix.is_resnum_nearer_top(nearest_ref_resnum) \
               and not test_element.is_resnum_nearer_top(nearest_test_resnum):
            if reference_helix.get_sideways():
                relpos = RELPOS_LEFT
            else:
                relpos = RELPOS_ABOVE
        else:
            if reference_helix.get_sideways():
                relpos = RELPOS_RIGHT
            else:
                relpos = RELPOS_BELOW
        return relpos


    def get_relpos_sheet_to_helix(self, reference_helix, closest_test_strand,
                                  nearest_ref_resnum, nearest_test_resnum,
                                  test_sheet_strandlists_dict):
        """
        Find the relative position of a sheet to a helix.

        Parameters:
          reference_helix - helix to find relpos of sheet to
          closest_test_strand - strand in sheet closest to reference helix
          nearest_ref_resnum - residue number in reference_helix that test
                               strand is closest to
           nearest_test_resnum - residue number in closest_test_strand
                                that is closest
                                to reference helix

          test_sheet_strandlists_dict - strandlists dict of test sheet

        Return value:
           relpos of sheet to helix

        """
        test_strand_posnum = self.get_strand_posnum(closest_test_strand,
                                                    test_sheet_strandlists_dict)

        if test_strand_posnum == 0:
            if reference_helix.get_sideways():
                relpos = RELPOS_BELOW
            else:
                relpos = RELPOS_RIGHT
        elif test_strand_posnum == \
       len(test_sheet_strandlists_dict[closest_test_strand.get_sheet_id()]) - 1:
            if reference_helix.get_sideways():
                relpos = RELPOS_ABOVE
            else:
                relpos = RELPOS_LEFT
        else:
            # need to decide ABOVE/BELOW
            # decide based on nearby residues at top/bottom
            if reference_helix.is_resnum_nearer_top(nearest_ref_resnum) \
               and not closest_test_strand.is_resnum_nearer_top(
                                                      nearest_test_resnum):
                if reference_helix.get_sideways():
                    relpos = RELPOS_LEFT
                else:
                    relpos = RELPOS_ABOVE
            else:
                if reference_helix.get_sideways():
                    relpos = RELPOS_RIGHT
                else:
                    relpos = RELPOS_BELOW

        return relpos
    
        
    def get_relpos_to_helix(self, reference_helix, test_element,
                            use_longest_for_orientation = False):
        """
        Find the relative position of test_element relative to
        the helix reference_helix

        Parameters:
            reference_helix - PTNodeHelix to find position relative to
            test_element - and element (sheet id or helix) to find
                   position of relative to reference_helix
            use_longest_for_orientation - (default True) if True, use
                   the longest strand in each sheet to determine the
                   relative orientations using tableau, otherwise uses
                   the closest strands (the ones used to determine relative
                   position).

           NOTE: reversed and sideways flags in test_element may be set
           by this function.

        Uses data members (read):
            distmatrix - the distance matrix
            sheet_strandlists_dict -
                dictionary of { sheet_id : list of list of nodes }
                where the list of list of nodes is
                described in build_sheet_constraints()
        Return value:
            tuple (relpos, reference_helix, test_strand) where relpos is
            RELPOS_ABOVE, RELPOS_BELOW, RELPOS_LEFT or RELPOS_RIGHT
            and test_strand is strand in test_element it is relative
            to or None if test element is not a sheet
            reference_helix is just the parameter 
        """
        assert(isinstance(reference_helix, PTNodeHelix))
        assert(isinstance(test_element, PTNodeHelix) or
               len(test_element) == 1) # sheet id

        if isinstance(test_element, PTNodeHelix):
            closest_test_strand = None
            # orientation needs to be taken into account (tableau)
            if self.tableau != None:
                try:
                    tabcode = self.tableau[(reference_helix, test_element)]
                    if verbose:
                        sys.stderr.write('  orientation ' +
                                         str(reference_helix) + ', ' +
                                         str(test_element) + ': ' +
                                         tabcode +
                                         '\n')
                except:
                    sys.stderr.write('WARNING: no tableau entry for ' +
                                     str(reference_helix) + ',' + 
                                     str(test_element) + '.' +
                                     'Using PE (parallel).\n')
                    tabcode = 'PE'
            else:
                tabcode = 'PE'

            # if test helix is crossing- Left or Right of reference,
            # and reference is not sideways, set test helix sideways,
            # and for crossing-Right set reversed flag if referce is not
            # reversed (and for xing-Left set reversedflag if reference IS
            # reversed).
            # Otherwise, if helices are antiparallel then set the
            # reversed flag in the test helix to the opposite value of
            # that in the reference helix.
            # FIXME: should clean this up and use resolve_orientation()
            if ( (tabcode[0] == 'L' or tabcode[0] == 'R')
                 and not reference_helix.get_sideways() ):
                test_element.set_sideways(True)
                if ( (tabcode[0] == 'R' and not reference_helix.get_reversed())
                     or
                     (tabcode[0] == 'L' and reference_helix.get_reversed()) ):
                    test_element.set_reversed(True)
            elif ( tabcode[0] == 'O' ):
                test_element.set_reversed(not reference_helix.get_reversed())
            
            # decide on placement based on nearest residues in the helices

            # FIXME: this is really no good, need to take account of
            # orientation and find some way of deciding if helices are really
            # 'beside' each other (esp when antiparalel for example)..
            # Note 'helix clustering' (partially) solves this problem for
            # the special case of being near a sheet to use as reference
            # for positioning, see get_helixcluster_relative_position().

            (nearest_ref_resnum, nearest_test_resnum) = \
                   self.distmatrix.get_nearest_sse_residues(reference_helix,
                                                            test_element)
            relpos = self.get_relpos_helix_to_helix(reference_helix,
                                                    test_element,
                                                    nearest_ref_resnum,
                                                    nearest_test_resnum)
        else:
            # the test element is a sheet. Place above or below, aligning
            # strand with helix, or if strand is on edge of sheet possibly
            # left/right of helix.
            (closest_test_strand, unused) = \
                  self.distmatrix.get_strand_nearest_element(test_element,
                                                             reference_helix)
            if verbose:
                sys.stderr.write('  relpos to helix: test is ' +
                                 str(closest_test_strand) + ' in sheet ' +
                                 test_element + '\n')

            # orientation needs to be taken into account (tableau)
            if self.tableau != None:
                if use_longest_for_orientation:
                    (orientation_test_strand, unused_length2) = \
                            self.get_longest_strand(
                               self.sheet_strandlists_dict[test_element])
                else:
                    orientation_test_strand= closest_test_strand
                try:
                    tabcode = self.tableau[(reference_helix,
                                            orientation_test_strand)]
                    if verbose:
                        sys.stderr.write('  orientation ' +
                                         str(reference_helix) + ', ' +
                                         str(orientation_test_strand) + ': ' +
                                         tabcode +
                                         '\n')
                except:
                    sys.stderr.write('WARNING: no tableau entry for ' +
                                     str(reference_helix) + ',' + 
                                     str(orientation_test_strand) + '.' +
                                     'Using PE (parallel).\n')
                    tabcode = 'PE'
            else:
                tabcode = 'PE'

            
            # if ref helix and test strand are antiparallel but flagged as
            # same direction in nodes, or parallel but flagged as different
            # direction in nodes, then flip them all strands in the
            # test sheet.
            if (((tabcode[0] == 'O') and   
                    reference_helix.get_reversed() ==     
                    closest_test_strand.get_reversed())  or  
                ((tabcode[0] == 'P') and  
                    reference_helix.get_reversed() !=     
                   closest_test_strand.get_reversed())):
                self.flip_all_strands_in_sheet(test_element)
            # if test strand is crossing- Left or Right of reference,
            # and reference is not sideways, set test sheet sideways
            # FIXME:  should clean this up and use resolve_orientation()
            elif ( (tabcode[0] == 'L' or tabcode[0] == 'R') and
                 not reference_helix.get_sideways() ):
                self.set_all_sheet_strands_sideways(test_element)
                if verbose:
                    sys.stderr.write('  sheet ' + test_element +
                                     ' is sideways (' + tabcode[0] + ')\n')
                # un-reversed ('up') when sideways is left-pointing
                if tabcode[0] == 'R':
                    self.flip_all_strands_in_sheet(test_element)

            (nearest_ref_resnum, nearest_test_resnum) = \
                   self.distmatrix.get_nearest_sse_residues(reference_helix,
                                                            closest_test_strand)
            relpos = self.get_relpos_sheet_to_helix(reference_helix,
                                                    closest_test_strand,
                                                    nearest_ref_resnum,
                                                    nearest_test_resnum,
                                                    self.sheet_strandlists_dict)

        if verbose:
            sys.stderr.write('  relpos to helix: test is ' +
                             ptrelpos_to_str(relpos) + ' reference.\n')


        return (relpos, reference_helix, closest_test_strand)




    def get_helixcluster_relative_position(self, reference_helix, test_helix,
                                           ref_strand):
        """
        Find the relative position of test_element helix relative to
        the helix reference_helix in a helix cluster, in which the
        first helix in the cluster is algined on the seq_strand axis.

        Parameters:
            reference_helix - PTNodeHelix to find position relative to
            test_helix - and element (PTNodeHelix) to find
                   position of relative to reference_helix
            ref_strand - The PTNodeStrand that we are deeming to be sharing
                    an axis with the reference_helix, used to align that helix.
                    For the first helix in the cluster, this is the strand
                    that the helix is immediately C-terminal of. For subsequent
                    helices, it is returned from this subroutine as the
                    strand we have decided it will be aligned with based
                    on dihedral angle (same/other side) calculation.

           NOTE: reversed and sideways flags in test_element may be set
           by this function.

        Uses data members (read):
            distmatrix - the distance matrix
            sheet_strandlists_dict -
                dictionary of { sheet_id : list of list of nodes }
                where the list of list of nodes is
                described in build_sheet_constraints()
        Return value:
            tuple (relpos, test_strand) where relpos is one of
            RELPOS_ABOVE, RELPOS_BELOW, RELPOS_LEFT or RELPOS_RIGHT
            and test_strand is the strand that we have decided the test-helix
            is on the same side of the ref_strand as.
        """
        assert(isinstance(reference_helix, PTNodeHelix))
        assert(isinstance(test_helix, PTNodeHelix))

        # orientation needs to be taken into account (tableau)
        if self.tableau != None:
            try:
                tabcode = self.tableau[(reference_helix, test_helix)]
                if verbose:
                    sys.stderr.write('  orientation ' +
                                     str(reference_helix) + ', ' +
                                     str(test_helix) + ': ' +
                                     tabcode +
                                     '\n')
            except:
                sys.stderr.write('WARNING: no tableau entry for ' +
                                 str(reference_helix) + ',' + 
                                 str(test_helix) + '.' +
                                 'Using PE (parallel).\n')
                tabcode = 'PE'
        else:
            tabcode = 'PE'

        # if test helix is crossing- Left or Right of reference,
        # and reference is not sideways, set test helix sideways,
        # and for crossing-Right set reversed flag if referce is not
        # reversed (and for xing-Left set reversedflag if reference IS
        # reversed).
        # Otherwise, if helices are antiparallel then set the
        # reversed flag in the test helix to the opposite value of
        # that in the reference helix.
        # FIXME: should clean thi sup and use resolve_orientation()
        if ( (tabcode[0] == 'L' or tabcode[0] == 'R')
             and not reference_helix.get_sideways() ):
            test_helix.set_sideways(True)
            if ( (tabcode[0] == 'R' and not reference_helix.get_reversed())
                 or
                 (tabcode[0] == 'L' and reference_helix.get_reversed()) ):
                test_helix.set_reversed(True)
        elif ( tabcode[0] == 'O' ):
            test_helix.set_reversed(not reference_helix.get_reversed())

        # decide on placement of test helix relative to reference helix
        # using the sheet containting the seq_strand as reference.
        # We will do this by using the dihedral angle calculation similar
        # to that used in deciding relative sides of strads in a sheet
        # (see strands_on_opposite_sides() in ptnode.py).
        # The already placed (reference) helix will be assumed to be
        # aligned on the axis of some strand
        # in the nearby sheet (this is the ref_strand parameter).
        # Then we compute the dihedral angle between
        # the planes formed by the axes of the test_helix and a neighbour
        # of that strand, with the reference strand
        # in common. If the absolute value of this angle is < pi/2 then
        # the test helix is on the same side of the reference strand
        # as the neighbour strand (i.e. we will say it is on the same
        # axis as the neihgbour strand), otherwise on the other side.
        #
        # TODO: for angles close to 0, should align on same as referene
        # helix, ie.. above/below it not left/right.

        # if one strand is sideways, all are
        sheet_is_sideways = ref_strand.get_sideways()
        
        ref_strand_posnum = self.get_strand_posnum(ref_strand)
        if (ref_strand_posnum ==
            len(self.sheet_strandlists_dict[ref_strand.get_sheet_id()])-1):
            neighbour_strand_posnum = ref_strand_posnum - 1
            other_side_strand_posnum = None # on the rightmost side of sheet
            if sheet_is_sideways:
                neighbour_relpos = RELPOS_ABOVE # XXX check this
                other_relpos = RELPOS_BELOW
            else:
                neighbour_relpos = RELPOS_LEFT
                other_relpos = RELPOS_RIGHT
        else:
            neighbour_strand_posnum = ref_strand_posnum + 1
            if ref_strand_posnum > 0:
                other_side_strand_posnum = ref_strand_posnum - 1
            else:
                other_side_strand_posnum = None # leftmost side of sheet
            if sheet_is_sideways:
                neighbour_relpos = RELPOS_BELOW # XXX check this
                other_relpos = RELPOS_ABOVE
            else:
                neighbour_relpos = RELPOS_RIGHT
                other_relpos = RELPOS_LEFT
        neighbour_strand = \
                   self.sheet_strandlists_dict[ref_strand.get_sheet_id()]\
                   [neighbour_strand_posnum][0]
        if other_side_strand_posnum != None:
            other_side_strand = \
                  self.sheet_strandlists_dict[ref_strand.get_sheet_id()]\
                  [other_side_strand_posnum][0]
        else:
            other_side_strand = None
        angle = ref_strand.axis_dihedral_angle(neighbour_strand, test_helix,
                                               self.pdb_struct)
        # FIXME: arbitrarily choosing 'same side' if angle cannot be calculated
        if angle == None or abs(angle) < pi/2:  # same side
            test_strand = neighbour_strand
            relpos = neighbour_relpos
        else:
            test_strand = other_side_strand
            relpos = other_relpos

        
        if verbose:
            sys.stderr.write('  helixcluster relpos helix: test is ' +
                             ptrelpos_to_str(relpos) + ' reference.\n')

        # FIXME: if ref strand is last in sheet and we end up on other
        # side from neighbour, there is no new test_strand to return,
        # need to do something else in this case.
        if test_strand == None:
            sys.stderr.write('WARNING: (helix clustering) '
                             'no reference strand for helix ' +
                             str(test_helix) + '\n')
            test_strand = ref_strand # FIXME: just use end strand for now
            

        return (relpos, test_strand)



    def get_relpos_helix_to_sheet(self, closest_ref_strand,
                                  nearest_ref_resnum):
        """
        Find the relative position of helix relative to sheet

        Parameters:
           closest_ref_strand - strand in sheet closest to the test helix
           nearest_ref_resnum - residue number in the closest_ref_strand that
                                is nearest the test helix.

        Note that the test helix itself is not needed in this funtion, it
        just uses the position of the nearest_ref_resnum to determine
        the relative position

        Return value:
          relpos (ABOVE/LEFT/etc.) to the ref strand

        Uses data members (readonly):
          distmatrix - the distance matrix
            sheet_strandlists_dict  - 
                dictionary of { sheet_id : list of list of nodes }
                where the list of list of nodes is
                described in build_sheet_constraints()

        """
        reference_sheetid = closest_ref_strand.get_sheet_id()
        ref_strand_posnum = self.get_strand_posnum(closest_ref_strand)
        if ref_strand_posnum == 0 or \
           ref_strand_posnum == \
           len(self.sheet_strandlists_dict[reference_sheetid]) - 1:
            # strand is on edge of sheet so could place helix beside
            # it if appropriate
            if ref_strand_posnum == 0:
                if closest_ref_strand.get_sideways():
                    relpos = RELPOS_ABOVE
                else:
                    relpos = RELPOS_LEFT
            else:
                if closest_ref_strand.get_sideways():
                    relpos = RELPOS_BELOW
                else:
                    relpos = RELPOS_RIGHT
        else:
            # not near strand on edge of sheet, place above/below
            if closest_ref_strand.is_resnum_nearer_top(nearest_ref_resnum):
                if closest_ref_strand.get_sideways():
                    relpos = RELPOS_LEFT
                else:
                    relpos = RELPOS_ABOVE
            else:
                if closest_ref_strand.get_sideways():
                    relpos = RELPOS_RIGHT
                else:
                    relpos = RELPOS_BELOW
        return relpos


    def get_relpos_sheet_to_sheet(self, closest_ref_strand,closest_test_strand,
                                  test_sheet_strandlists_dict,
                                  tabcode,
                                  enable_changes=False):
        """
        Find the relative position of a sheet relative to sheet

        Parameters:
           closest_ref_strand - strand in ref sheet closest to the test sheet
           closest_test_strand - strand in test sheet closest to ref sheet
           test_sheet_strandlists_dict - The sheet_strandlists_dict for the
                              test sheet
           tabcode - two character tableau code for relative orinetation
                     between the two sheets
           enable_changes - (default False) If True, the function can
                            change reverse/sideways flags in strands of test
                            sheet, otherwise does not change them.
        
        Return value:
          relpos (ABOVE/LEFT/etc.) to the ref strand

        Uses data members:
            distmatrix - the distance matrix
            sheet_strandlists_dict  - (read, write (only if enable_changes=True)
                dictionary of { sheet_id : list of list of nodes }
                where the list of list of nodes is
                described in build_sheet_constraints()
                NB writing to this refers to changing the orientation
                (reversed/sideways)
                flags in the strand nodes if necessary, or to reversing
                the order of the (outermost) node list if necessary.
                This is only done if parameter enable_changes=True
        
        """
        # if test and ref strands are both first/last in sheet,
        # and sheets are parallel or antiparallel (not crossing),
        # then place the sheets side-by-side. if test strand is
        # on 'right' end of sheet and ref is also on 'right' end of sheet
        # (or similarly for left), then we will 'flip' the test sheet
        # over by reverseing order of strands so the ref and test
        # strands are beside each other.
        test_strand_posnum = self.get_strand_posnum(closest_test_strand,
                                                    test_sheet_strandlists_dict)
        ref_strand_posnum = self.get_strand_posnum(closest_ref_strand)
        ref_left_edge = (ref_strand_posnum == 0)
        reference_sheetid = closest_ref_strand.get_sheet_id()
        test_sheetid = closest_test_strand.get_sheet_id()
        ref_right_edge = (ref_strand_posnum == 
                  len(self.sheet_strandlists_dict[reference_sheetid]) - 1)
        test_left_edge = (test_strand_posnum == 0)
        test_right_edge = (test_strand_posnum == 
                  len(test_sheet_strandlists_dict[test_sheetid]) - 1)
        crossing = (tabcode[0] == 'R' or tabcode[0] == 'L') #not par/antipar
        if ( not crossing and ( (ref_left_edge or ref_right_edge) and
                                (test_left_edge or test_right_edge) ) ):
            if ( (ref_left_edge and test_left_edge) or
                 (ref_right_edge and test_right_edge) ):
                if enable_changes:
                    self.reverse_strand_order_in_sheet(test_sheetid,
                                                  test_sheet_strandlists_dict)
                    if verbose:
                        sys.stderr.write('  reversed strand order for sheet ' +
                                         test_sheetid +
                                         ' so test strand is near ref strand\n')
            if ref_left_edge:
                if closest_ref_strand.get_sideways():
                    relpos = RELPOS_ABOVE
                else:
                    relpos = RELPOS_LEFT
            else:
                assert(ref_right_edge)
                if closest_ref_strand.get_sideways():
                    relpos = RELPOS_BELOW
                else:
                    relpos = RELPOS_RIGHT

        else:
            # need to decide ABOVE/BELOW
            # decide based on nearby residues at top/bottom
            # For now, let's try a kind of dodgy method of taking
            # the 'central' strand in the test and the longest
            # in the reference, and finding the test residue to the
            # 'top' and 'bottom' residues in the reference strand.
            # If it is nearer top, position above else below.
            # FIXME: should have something more principled here, e.g.
            # actually using sheet centroids and projecting onto
            # plane of largest sheet or something.
            # Did try using projectino of c and n term points onto
            # axis, but results are even more inconsistent, on serpins
            # esp. since longest strand has high curvature and irregularity
            # (see notebook 10/2/08 - maybe should try using strand
            # with best fitted axis for this (and orientation)).
            # So now to stop small differences in relative 3d positions
            # causing different representation of topologically similar
            # structures, we have a fudge factor and assume BELOW
            # unless 'significantly' closer to other end.

#                 (nearest_ref_resnum, nearest_test_resnum) = \
#                    self.distmatrix.get_nearest_sse_residues(closest_ref_strand,
#                                                             closest_test_strand)
#                 if closest_ref_strand.is_resnum_nearer_top(nearest_ref_resnum) \
#                    and not closest_test_strand.is_resnum_nearer_top(
#                                                           nearest_test_resnum):
            test_central_strand = \
                  test_sheet_strandlists_dict[test_sheetid] \
                     [len(test_sheet_strandlists_dict[test_sheetid])/2][0]
            (ref_longest_strand, length_unused) = \
                  self.get_longest_strand(\
                           self. sheet_strandlists_dict[reference_sheetid])
            # residue lists are ordered in N to C direction
            test_residue_list = \
                    test_central_strand.get_residue_list()
            ref_residue_list = \
                    ref_longest_strand.get_residue_list()
            test_central_residue = \
                                test_residue_list[len(test_residue_list)/2]
            if enable_changes:
                dist_to_ref_nterm = \
                      self.distmatrix.get_distance(test_central_residue,
                                                   ref_residue_list[0])
                dist_to_ref_cterm = \
                      self.distmatrix.get_distance(test_central_residue,
                                                   ref_residue_list[-1])
            else:
                # something of an abuse of the name of this variable; it
                # is only True when using an external (outside this domain)
                # element as the test element, so in such a case we cannot
                # use the self.distmatrix so explicitly calculate the
                # distances instead
                dist_to_ref_nterm = \
                      calc_residue_dist(test_central_residue,
                                        ref_residue_list[0])
                dist_to_ref_cterm = \
                      calc_residue_dist(test_central_residue,
                                        ref_residue_list[-1])
            near_cterm = (dist_to_ref_cterm < dist_to_ref_nterm)
            FUDGE = 0.15 # difference must be more than 15% of min dist
            is_signficant = (abs(dist_to_ref_nterm - dist_to_ref_cterm)
                          > FUDGE*min(dist_to_ref_nterm,dist_to_ref_cterm))
            if verbose:
                sys.stderr.write('  sheet relpos test strand ' +
                                 str(test_central_strand) + ' ref strand ' +
                                 str(ref_longest_strand) + '\n')
                sys.stderr.write('    cterm dist = ' +str(dist_to_ref_cterm)
                                 + ' nterm dist = ' +str(dist_to_ref_nterm))
                sys.stderr.write('; is_signifcant = ' +
                                 str(is_signficant) + '\n')
            if is_signficant:
                if ref_longest_strand.get_reversed():
#                    print 'zzzzzzzz reversed'
                    near_top = not near_cterm
                else:
                    near_top = near_cterm

#                print 'zzz',near_cterm,near_top

                if near_top:
                    if closest_ref_strand.get_sideways():
                        relpos = RELPOS_LEFT
                    else:
                        relpos = RELPOS_ABOVE
                else:
                    if closest_ref_strand.get_sideways():
                        relpos = RELPOS_RIGHT
                    else:
                        relpos = RELPOS_BELOW
            else:
                if closest_ref_strand.get_sideways():
                    relpos = RELPOS_RIGHT
                else:
                    relpos = RELPOS_BELOW


            # make sure the closest test strand is drawn close to the
            # ref sheet. Since sheet is always drawn by strands in order
            # they are in the list of list of strands in the
            # sheet_strandlists_dict for the sheet, we may need to
            # reverse this list.
            if ( enable_changes and
                 ((ref_strand_posnum <
                   (len(self.sheet_strandlists_dict[reference_sheetid]))/2
                   and
                   test_strand_posnum >
                   (len(self.sheet_strandlists_dict[test_sheetid])-1) / 2 )
                  or
                  (ref_strand_posnum >
                   (len(self.sheet_strandlists_dict[reference_sheetid]))/2
                   and
                   test_strand_posnum <
                   (len(test_sheet_strandlists_dict[test_sheetid])) / 2) ) ):
                self.reverse_strand_order_in_sheet(test_sheetid,
                                                   self.sheet_strandlists_dict)
                if verbose:
                    sys.stderr.write('  reversed strand order for sheet ' +
                                     test_sheetid +
                                     ' so test strand is near ref strand\n')

        return relpos

        
    
    def get_relpos_to_sheet(self, reference_sheetid, test_element,
                            use_longest_for_orientation=True):
        """
        Find the relative position of test_element relative to
        the sheet reference_sheetid

        Parameters:
            reference sheet_id - sheet id of sheet to find position
                  of test_element relative to
            test_element - and element (sheet id or helix) to find
                   position of relative to reference sheet
            use_longest_for_orientation - (default True) if True, use
                   the longest strand in each sheet to determine the
                   relative orientations using tableau, otherwise uses
                   the closest strands (the ones used to determine relative
                   position).

           NOTE: reversed and sideways flags in test_element may be set
           by this function. It may also reverse the horiz_order_list
           for the test sheet.

        Uses data members (read):
            distmatrix - the distance matrix
            sheet_dict - dict of {sheet_id : nodelist }
            sheet_strandlists_dict  - (read/write)
                dictionary of { sheet_id : list of list of nodes }
                where the list of list of nodes is
                described in build_sheet_constraints()
                NB writing to this refers to changing the orientation
                (reversed/sideways)
                flags in the strand nodes if necessary, or to reversing
                the order of the (outermost) node list if necessary.
                
        Return value:
            tuple (relpos, ref_strand, test_strand) where relpos is
            RELPOS_ABOVE, RELPOS_BELOW, RELPOS_LEFT or RELPOS_RIGHT
            and ref_strand is strand in reference sheet to which it is
            relative and
            test_strand is strand in test element relative to it or
            None if test element is not a strand
        """
        assert(len(reference_sheetid) == 1)
        assert(isinstance(test_element, PTNodeHelix) or
               len(test_element) == 1) # sheet id

        # find the strand in reference sheet the test obj is closest to
        (closest_ref_strand, closest_test_strand)= \
                  self.distmatrix.get_strand_nearest_element(reference_sheetid,
                                                             test_element)
        if verbose:
            sys.stderr.write('  relpos to sheet: reference is ' +
                             str(closest_ref_strand) + ' in sheet ' +
                             reference_sheetid + '\n')
        if isinstance(test_element, PTNodeHelix):
            # place helix close to (above, below, or ,for strands on edge
            # sheet, left/right) the reference strand
            (nearest_ref_resnum, nearest_test_resnum) = \
                   self.distmatrix.get_nearest_sse_residues(closest_ref_strand,
                                                            test_element)
            relpos = self.get_relpos_helix_to_sheet(closest_ref_strand,
                                                    nearest_ref_resnum)
            # orientation needs to be taken into account (tableau)
            if use_longest_for_orientation:
                (orientation_ref_strand, unused_length) = \
                        self.get_longest_strand(
                           self.sheet_strandlists_dict[reference_sheetid])
            else:
                orientation_ref_strand = closest_ref_strand
            if self.tableau != None:
                try:
                    tabcode=self.tableau[(orientation_ref_strand, test_element)]
                                          
                    if verbose:
                        sys.stderr.write('  orientation ' +
                                         str(orientation_ref_strand) + ', ' +
                                         str(test_element) + ': ' + tabcode +
                                         '\n')
                except:
                    sys.stderr.write('WARNING: no tableau entry for ' +
                                     str(orientation_ref_strand) + ',' + 
                                     str(test_element) + '.' +
                                     'Using PE (parallel).\n')
                    tabcode = 'PE'
            else:
                tabcode = 'PE'

            # if ref strand and helix are not crossing, but ref sheet is
            # sideways and helix isn't, then set helix sideways if
            # ref sheet is (Note sheets start out not sideways)
            # FIXME: should clean this up and use resolve_orientation()
            crossing = (tabcode[0] == 'R' or tabcode[0] == 'L') #not par/antipar
            if ( not crossing and closest_ref_strand.get_sideways() ):
                test_element.set_sideways(True)
                
            # if test helix is crossing- Left or Right of reference,
            # and reference is not sideways, set test helix sideways
            if ( (tabcode[0] == 'L' or tabcode[0] == 'R')
                 and not orientation_ref_strand.get_sideways() ):
                test_element.set_sideways(True)

        else:
            # the test element is a sheet. Place above or below, aligning
            # strands, or, if ref and test strand are both on edge
            # of sheet, left or right.
            if verbose:
                sys.stderr.write('  relpos to sheet: test is ' +
                                 str(closest_test_strand) + ' in sheet ' +
                                 test_element + '\n')

            # orientation needs to be taken into account (tableau)
            if use_longest_for_orientation:
                (orientation_ref_strand, unused_length) = \
                        self.get_longest_strand(
                           self.sheet_strandlists_dict[reference_sheetid])
                (orientation_test_strand, unused_length2) = \
                        self.get_longest_strand(
                           self.sheet_strandlists_dict[test_element])
            else:
                orientation_ref_strand = closest_ref_strand
                orientation_test_strand= closest_test_strand
            if self.tableau != None:
                try:
                    tabcode = self.tableau[(orientation_ref_strand,
                                            orientation_test_strand)]
                    if verbose:
                        sys.stderr.write('  orientation ' +
                                         str(orientation_ref_strand) + ', ' +
                                         str(orientation_test_strand) + ': ' +
                                         tabcode +
                                         '\n')
                except:
                    sys.stderr.write('WARNING: no tableau entry for ' +
                                     str(orientation_ref_strand) + ',' + 
                                     str(orientation_test_strand) + '.' +
                                     'Using OS (antiparallel).\n')
                    tabcode = 'OS' #2nd char is arbitrary
            else:
                tabcode = 'OS'


            crossing = (tabcode[0] == 'R' or tabcode[0] == 'L') #not par/antipar

            # heuristic test for 'folded over' sheets (like sandwiches)
            # where if the tabcode is  antiparallel, we
            # actually want to reverse it (antipar->par)
            # so that it is as if we have 'unfolded' the sheets along
            # the 'hinge' formed by the coils between adjacent in
            # sequence strands. (see notes 26/2/08-26/2/08 (FIXME:
            # should describe this better here rather than referecning
            # handwritten notes!))
            # This test is that if at least two strands
            # in the test sheet immediately follow strands in ref sheet
            # or vice versa
            # and orientation is antiparallel then convert it to parallel.
            #  NB the HH and KK codes are only supposed to be used
            # for strands in the same sheet (not between strands in different
            # sheets) and that is now (10June2008) what is implemented.
            # so we never check KK or HH here since we are dealing with
            # strands in different sheets, instead always stick to P or O.
            if tabcode[0] == 'P':
                ADJSTRAND_COUNT_THRESHOLD = 2 # at least this many to reverse
                adjstrand_count = 0
                for strand1 in self.sheet_dict[reference_sheetid]:
                    if self.any_strands_before_or_after_strand(
                        strand1, self.sheet_dict[test_element]):
                        adjstrand_count += 1
#                print 'xxxx',reference_sheetid,test_element,adjstrand_count
                if adjstrand_count >= ADJSTRAND_COUNT_THRESHOLD:
                    if verbose:
                        sys.stderr.write('  sheet ' + reference_sheetid +
                                         ' and  sheet ' + 
                                         test_element +
                                         ' folded over: reversing orientation\n')
                    tabcode = 'OS' # 2nd char is arbitrary.

            # if ref and test strands are not crossing, but one sheet is
            # sideways and other isn't, then set test sheet sideways if
            # ref sheet is (Note sheets start out not sideways)
            if ( not crossing and closest_ref_strand.get_sideways() ):
                self.set_all_sheet_strands_sideways(test_element)
                
            # if ref and test strands are antiparallel but flagged as
            # same direction in nodes, or parallel but flagged as different
            # direction in nodes, then flip them all strands in the
            # test sheet.
#            print 'qqqqq',tabcode,orientation_test_strand,orientation_test_strand.get_reversed(),orientation_ref_strand,orientation_ref_strand.get_reversed()
            if (((tabcode[0] == 'O') and   
                    orientation_ref_strand.get_reversed() ==     
                    orientation_test_strand.get_reversed())  or  
                ((tabcode[0] == 'P') and  
                    orientation_ref_strand.get_reversed() !=     
                   orientation_test_strand.get_reversed())):
#                print 'zzzzzzzzzzzzzzzzzzzzzzzzzz',orientation_test_strand,orientation_ref_strand
                self.flip_all_strands_in_sheet(test_element)

            # if test strand is crossing- Left or Right of reference,
            # and reference is not sideways, set test sheet sideways
            elif ( crossing and not closest_ref_strand.get_sideways() ):
                self.set_all_sheet_strands_sideways(test_element)
                if verbose:
                    sys.stderr.write('  sheet ' + test_element +
                                     ' is sideways (' + tabcode[0] + ')\n')
                # un-reversed ('up') when sideways is left-pointing
                if ( (tabcode[0] == 'R' and
                      not orientation_test_strand.get_reversed()) or
                     (tabcode[0] == 'L' and
                      orientation_test_strand.get_reversed()) ):
                    self.flip_all_strands_in_sheet(test_element)

            relpos = self.get_relpos_sheet_to_sheet(closest_ref_strand,
                                  closest_test_strand,
                                  self.sheet_strandlists_dict,
                                  tabcode,
                                  enable_changes=True)

        if verbose:
            sys.stderr.write('  relpos to sheet: test is ' +
                             ptrelpos_to_str(relpos) + ' reference.\n')
        
        
        return (relpos, closest_ref_strand, closest_test_strand)


    def get_external_relpos(self, reference_element, test_element,
                            closest_ref_strand,
                            closest_test_strand,
                            nearest_ref_resnum,
                            nearest_test_resnum,
                            tabcode,
                            test_sheet_strandlists_dict):
        """
        Find the relative position of test_element relative to
        reference_element, where test_element is not an element in this
        domain. Used for releative placement of domains.

        Parameters:
            reference_element - an element (either sheet id e.g. 'A' or
                   helix (PTNodeHelix object) to find position relative to,
                   the element is in this domain
            test_element - and element (as per reference_element) to find
                   position of relative to reference_element,
                   the element is not in this domain (cannot use member
                   data sheet_strandlists_dict etc. for infomratin on
                   this element)
            closest_ref_strand - strand in reference sheet if reference_element
                                is a sheet, else None.
           closest_test_strand - strand in test sheet if test_element is a sheet
                                else None.
           nearest_ref_resnum - residue number in reference SSE that test
                               element is closest to
           nearest_test_resnum - residue number in test SSE that is closest
                                to reference element
           tabcode - two char tableau code for relative orientation of the
                                external domain with this domain.
           test_sheet_strandlists_dict - sheet strandlists dict for test
                                element when it is a sheet (else None).
                                Note this is neeed as the test element is
                                 not part of this domain. (For the reference
                                 element in this domain, the data member
                                 sheet_strandlists_dict can be used).

               
        Uses data members (read):
            distmatrix - the distance matrix
            sheet_strandlists_dict  -
                dictionary of { sheet_id : list of list of nodes }
                where the list of list of nodes is
                described in build_sheet_constraints()
            
        Return value:
            relpos where relpos is
            RELPOS_ABOVE, RELPOS_BELOW, RELPOS_LEFT or RELPOS_RIGHT.
        """
        assert(isinstance(reference_element, PTNodeHelix) or
               len(reference_element) == 1) # sheet id
        assert(isinstance(test_element, PTNodeHelix) or
               len(test_element) == 1) # sheet id
        
        if isinstance(reference_element, PTNodeHelix):
            if isinstance(test_element, PTNodeHelix):
                relpos = self.get_relpos_helix_to_helix(reference_element,
                                                   test_element,
                                                   nearest_ref_resnum,
                                                   nearest_test_resnum)
            else:
                relpos = self.get_relpos_sheet_to_helix(reference_element,
                                                   closest_test_strand,
                                                   nearest_ref_resnum,
                                                   nearest_test_resnum,
                                                   test_sheet_strandlists_dict)
        else:
            # reference element is a sheet
            if isinstance(test_element, PTNodeHelix):
                relpos = self.get_relpos_helix_to_sheet(closest_ref_strand,
                                                   nearest_ref_resnum)
            else:
                relpos = self.get_relpos_sheet_to_sheet(closest_ref_strand,
                                                   closest_test_strand,
                                                   test_sheet_strandlists_dict,
                                                   tabcode,
                                                   enable_changes = False)
        return relpos
        
    ##########################################################################
    
    
#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

def ptrelpos_set_verbose(verb):
    """
    set the module global verbose flag in this module to supplied value
    Parameters: verb - True (for verbose output) or False
    Return value: None
    Uses globals: verbose (in this module)
    """
    global verbose
    verbose = verb
    

def ptrelpos_to_str(relpos):
    """
    Return string representation of relative position RELPOS_ABOVE etc.
    for verbose output/debugging.
    Parameters: relpos - RELPOS_ABOVE, etc.
    Return value: string corresponding to relpos
    """
    if relpos == RELPOS_ABOVE:
        s = "ABOVE"
    elif relpos == RELPOS_BELOW:
        s = "BELOW"
    elif relpos == RELPOS_LEFT:
        s = "LEFT of"
    elif relpos == RELPOS_RIGHT:
        s = "RIGHT of"
    else:
        s = "*BAD RELPOS (" + str(relpos) + ") *"
    return s


def resolve_orientation(tabcode, ref_sse, test_sse):
    """
    Resolve the orientation encoded in tableau code (see pttableau.py)
    between ref_sse and test_sse into a (sideways, reversed) tuple.

    Parameters:
       tabdoe - two charater tableau code for orienmtation between ref_sse and
                test_sse
       ref_sse - PTNode (strand or helix) as reference (sideways and reversed
                 taken to be already fixed in this node)
       test_sse - PTNode (strand or helix) to return sideways/reversed flags
                 for, relative to ref_sse, using tabcode

    Return value:
       tuple (sideways, reversed) where sideways and reversed are Boolean
       describing if test_sse needs to be sideways or reversed to have
       correct relationship to ref_sse according to tabcode
    """
    crossing = (tabcode[0] == 'R' or tabcode[0] == 'L')
    if ( (crossing and not ref_sse.get_sideways()) or
         (not crossing and ref_sse.get_sideways()) ):
        sideways = True
    else:
        sideways = False

    parallel = (tabcode[0] == 'P' or tabcode[0] == 'K')
    if ( (parallel and ref_sse.get_reversed()) or
         (not parallel and not ref_sse.get_reversed()) ):
        reversed = True
    else:
        reversed = False

    return (sideways, reversed)


