###############################################################################
#
# ptsvgconstraints.py - represent layout constraints for writing as Dunnart SVG
#
# File:    ptsvgconstraints.py
# Author:  Alex Stivala
# Created: February 2008
#
# $Id$
#
# A PTVSVGConstraint may be one of several types for separation constraints,
# distribution constraints, alignment constraints etc.
#
###############################################################################

import sys

from ptsvgnode import PTSVGNode,PTGRAPH_NS,get_residue_strings,PTSVGNodeTerminus
from ptsecstruct import stride_chainid_to_pdb_chainid
from ptutils import get_int_icode,char_if_not_blank,biopdbresid_to_pdbresseq
from Bio.PDB import *

#-----------------------------------------------------------------------------
#
# Constants
#
#-----------------------------------------------------------------------------

DUNNART_ALIGN_TOP       = 0      # dunnart constants (guideline.h)
DUNNART_ALIGN_MIDDLE    = 1
DUNNART_ALIGN_BOTTOM    = 2
DUNNART_ALIGN_LEFT      = 3
DUNNART_ALIGN_CENTER    = 4
DUNNART_ALIGN_RIGHT     = 5
DUNNART_GUIDE_TYPE_VERT = 100
DUNNART_GUIDE_TYPE_HORI = 101

# list of connector colors. Use only dark colors (and full
# opacity) since connectors are drawn as thin lines so light
# colors will be hard to see.
DUNNART_DEFAULT_LINE_COLOR = "000000ff"              # black
DUNNART_LINE_COLORS = [ DUNNART_DEFAULT_LINE_COLOR,
                        "0000ffff",                  # blue
                        "551a8bff",                  # purple4
                        "8b008bff",                  # magenta4
                        "8b0000ff",                  # darkred
                        "8b2323ff",                  # brown4
                        "8fbc8fff",                  # dark sea green
                        "191970ff",                  # midnight blue
                        "9400d3ff"                   # dark violet
                      ]

#-----------------------------------------------------------------------------
#
# Class definitions
#
#-----------------------------------------------------------------------------


class PTSVGIndGuide:
    """
    PTSVGIndGuide represents a guideline for aligning shapes on.
    """
    def __init__(self, xmlid, pos, direction):
        """
        Create a PTSVGIndGuide given unique XML id, position, and direction.

        Parameters:
            xmlid - unique XML id for this guide
            pos - position (x or y depending on direction)
            direction - DUNNART_GUIDE_TYPE_VERT or _HORI

        Raises Exceptions:
           ValueError for bad direction
        """
        if direction not in [DUNNART_GUIDE_TYPE_HORI, DUNNART_GUIDE_TYPE_VERT]:
            raise ValueError('bad direction '  + direction)
        self.xmlid = xmlid
        self.pos = pos
        self.direction = direction

    def write_svg(self, fh):
        """
        Write this indguide to the SVG file
                
        Parameters:
             fh - open filehandle to write SVG XML text to

        Return value:
            None

        """
        fh.write('  <dunnart:node dunnart:type="indGuide" ' +
                 'dunnart:position="' + str(self.pos) + '" ' +
                 'dunnart:direction="' + str(self.direction) + '" ' + 
                 'id="' + str(self.xmlid) + '"/>\n')

    def translate(self, xshift, yshift):
        """
        Move theindguide left/right by xshift (-/+) and/or up/down by
        yshift (-/+)..

        Parameters:
           xshift - amount to move left(-) or right(+) by
           yshift - amount to mvoe up(-) or down(+) by

        Return value:
           None

        Modifies data members:
           pos
        """
        if self.direction == DUNNART_GUIDE_TYPE_HORI:
            self.pos += yshift
        else: # VERT
            self.pos += xshift
        
        
        
class PTSVGDistribution:
    """
    PTSVGDistribution represents the distribution handle for distribution
    constraint.
    """

    def __init__(self, xmlid, direction, sepdistance, position):
        """
        Construct a PTSVGDistribution, given direction, separation
        distance, position  and id.

        Parameters:
            xmlid - unique XML identifier for this object
            direction - DUNNART_GUID_TYPE_HORI or _VERT
            sepdistance - separation distance
            position - position (x or y depending on HORI or VERT) for handle
            
        Raises Exceptions:
           ValueError for bad direction
        """
        if direction not in [DUNNART_GUIDE_TYPE_HORI, DUNNART_GUIDE_TYPE_VERT]:
            raise ValueError('bad direction '  + direction)

        self.xmlid = xmlid
        self.direction = direction
        self.sepdistance = sepdistance
        self.position = position

    def write_svg(self, fh):
        """
        Write this distribution to the SVG file
                
        Parameters:
             fh - open filehandle to write SVG XML text to

        Return value:
            None

        """
        fh.write('  <dunnart:node dunnart:type="distribution" ' +
                 'dunnart:direction="' + str(self.direction) + '" ' +
                 'dunnart:sepDistance="' +
                 str(self.sepdistance) + '" ' +
                'dunnart:position="' + str(self.position) + '" ' +
                 'id="' + str(self.xmlid) + '"/>\n')

    def translate(self, xshift, yshift):
        """
        Move the handle left/right by xshift (-/+) and/or up/down by
        yshift (-/+)..

        Parameters:
           xshift - amount to move left(-) or right(+) by
           yshift - amount to mvoe up(-) or down(+) by

        Return value:
           None

        Modifies data members:
           position
        """
        if self.direction == DUNNART_GUIDE_TYPE_HORI:
            self.position += xshift
        else: # VERT
            self.position += yshift
        


class PTSVGDistroConstraint:
    """
    PTSVGDistroConstraint represents a distribution constraint, that is
    a common separation between a set of indguides. It refers
    back to its PTSVGDistribution which represents the distribution handle.
    """

    def __init__(self, indguide1, indguide2, distro):
        """
        Construct a PTSVGDistroConstraint given two indguides and a
        distribution.

        Parameters:
            indguide1 - PTSVGIndGuide for one object to align
            indguide2 - PTSVGIndGuide for the other object to align
            distro    - PTSVGDistribution to align induides with
        """
        self.indguide1 = indguide1
        self.indguide2 = indguide2
        self.distro = distro

    def write_svg(self, fh):
        """
        Write this constraint to the SVG file
                
        Parameters:
             fh - open filehandle to write SVG XML text to

        Return value:
            None

        """
        fh.write('  <dunnart:node dunnart:type="constraint" ' +
                 'isMultiway="1" relType="distribution" ' +
                 'constraintID="' +str(self.distro.xmlid)+'" '+
                 'objOneID="'+str(self.indguide1.xmlid)+'" '+
                 'objTwoID="' +str(self.indguide2.xmlid)+'"/>\n')

    def translate(self, xshift, yshift):
        """
        Translation has no menaing for an actual constraint
        """
        return
    

class PTSVGAlignmentConstraint:
    """
    PTSVGAlignmentConstraint represents an alignment constraint, that is
    a shape is aligned on an indguide.
    """
    def __init__(self, indguide, svgnode, alignpos):
        """
        Construct a PTSVGAlignmentConstraint, given a PTSVGIndGuide object,
        PTSVGNode derived object (for a shape) and an alignment position.

        Parameters:
            indguide - PTSVGIndGuide to align on
            svgnode  - PTSVGNode to align
            alignpos - DUNNART_ALIGN_LEFT/CENTER/RIGHT/TOP/MIDDLE/BOTTOM
        Raises Exceptions:
             ValueError for invalid alignpos
             TypeError for wrong type of indguide or svgnode
        """
        if not isinstance(svgnode, PTSVGNode):
            raise TypeError('wrong type for svgnode')
        if not isinstance(indguide, PTSVGIndGuide):
            raise TypeError('wront type for indguide')
        if (alignpos not in [DUNNART_ALIGN_LEFT, DUNNART_ALIGN_CENTER,
                            DUNNART_ALIGN_RIGHT, DUNNART_ALIGN_TOP,
                            DUNNART_ALIGN_MIDDLE,DUNNART_ALIGN_BOTTOM]):
            raise ValueError('invalid alignpos ' + alignpos)
        self.indguide = indguide
        self.svgnode = svgnode
        self.alignpos = alignpos

    def write_svg(self, fh):
        """
        Write this constraint to the SVG file
                
        Parameters:
             fh - open filehandle to write SVG XML text to

        Return value:
            None

        """
        fh.write('  <dunnart:node dunnart:type="constraint" ' +
                 'isMultiway="1" relType="alignment" ' +
                 'constraintID="' + str(self.indguide.xmlid) +
                 '" ' +
                'objOneID="' + str(self.svgnode.xmlid) + '" ' +
                'alignmentPos="' + str(self.alignpos) + '" />\n')

    def translate(self, xshift, yshift):
        """
        Translation has no menaing for an actual constraint
        """
        return


# This is not a constraint but might a well go here
class PTSVGConnector:
    """
    PTSVGConnector represents a connector (line) between shapes
    """
    def __init__(self, xmlid, src, dest, srcFlags, dstFlags, color,
                 directed = False):
        """
        Create a PTSVGConnector.

        Parameters:
           xmlid - unique XML id for this connector
           src - PTSVGNode connector is FROM
           dest - PTSVGNode connector is TO
           srcFlags  - connector flags (for ports, DUNNART_DEFAULT_PORT, etc.)
           destFlags - connector flags
           color - line color hex RGB string
           directed - (default False) If True, puts arrowhead on dest end
        """
        self.xmlid = xmlid
        self.src = src
        self.dest = dest
        self.srcFlags = srcFlags
        self.dstFlags = dstFlags
        self.color = color
        self.directed = directed

        # The following are built by build_resname_sequence():
        
        self.issued_discontinuous_warning = False
        self.residue_list = [] # list of Bio.PDB Residues in the coil region
        self.resname_list = [] # list of 3 letter residue names in this coild
                               # region

        self.resid_list = []   # list of string PDB residue sequence numbers

        self.nterm_resname = None # 3 letter residue name of residue immediately
                                  # N-terminal of this coil region, or ""
                                  # for N-terminus
                                  
        self.nterm_resid = None   # PDB residue sequence number of residue
                                  # immediately N-terminal of this coil region,
                                  # or None for N-terminus.
                                  
        self.cterm_resname = None # 3 letter residue name of residue immediately
                                  # C-terminal of this coil region, or ""
                                  # for C-terminus
                                  
        self.cterm_resid = None   # PDB residue sequence number of residue
                                  # immediately C-terminal of this coil region,
                                  # or None for C-terminus.


    
        
    def get_residue_list(self, pdb_residue_list, pdb_resid_dict):
        """
        Return the list of Bio.PDB Residue objects for the residues in the
        coil region represneted by this connector.

        Parameters:
          pdb_residue_list - list of all residues (for all chains) in the protein
          pdb_resid_dict -  dict of { {chainid,pdb_resseq) : seqindx }
                               where chainid and pdb_resseq make up
                               the PDB residue identifier, the pdb_resseq
                               being string resnum+icode if any e.g.
                               '60' or '60A', seqindx is the indiex
                               into sequential list of all residues
                                pdb_residue_list.
        

        Return value:
          tuple (nterm_residue, cterm_residue, residue_list)
          where
          nterm_residue is the Bio.PDB Residue object of residue immediately
            nterminal of this coil region, or None if N terminus
          cterm_residue is the Bio.PDB Residue object of residue immediately
            cterminal of this coil regino, of None if C terminus
          residue_list is list of Bio.PDB Residue objects of residues in
            this coil region.

        Uses data members (read/write):
           residue_list - used to memoize building of the residue_list
        """
        if isinstance(self.src, PTSVGNodeTerminus):
            nterm_residue = None
            start_indx = 0
        else:
            nterm_residue = self.src.get_residue_list()[-1]
            start_indx = pdb_resid_dict[(self.src.chainid,
                                         self.src.get_end_res_seq())] + 1
            
        if isinstance(self.dest, PTSVGNodeTerminus):
            cterm_residue = None
            end_indx  = len(pdb_residue_list) - 1
        else:
            cterm_residue = self.dest.get_residue_list()[0]
            end_indx = pdb_resid_dict[(self.dest.chainid,
                                       self.dest.get_start_res_seq())] - 1

        if self.residue_list: # memoization: use if already computed
            resisdue_list = self.residue_list

        residue_list = pdb_residue_list[start_indx : end_indx + 1]
        self.residue_list = residue_list

        return (nterm_residue, cterm_residue, residue_list)

                               
    def build_resname_sequence(self, pdb_residue_list, pdb_resid_dict):
        """
        Build list of (3 letter) residue names in sequence for the residues
        in this node (SSE). E.g. and matching
        list of PDB residue sequence numbers 
        
        Parameters:
          pdb_residue_list - list of all residues (for all chains) in the protein
          pdb_resid_dict -  dict of { {chainid,pdb_resseq) : seqindx }
                               where chainid and pdb_resseq make up
                               the PDB residue identifier, the pdb_resseq
                               being string resnum+icode if any e.g.
                               '60' or '60A', seqindx is the indiex
                               into sequential list of all residues
                                pdb_residue_list.

        Return value: None
        Uses data member (write): resname_list
                                  resid_list
        """
        (nterm_residue, cterm_residue, residue_list) =\
             self.get_residue_list(pdb_residue_list, pdb_resid_dict)
        self.resname_list = [residue.get_resname() for residue in residue_list]
        # id of a residue in Bio.PDB is tuple (hetatm, resseqnum, icode)        
        self.resid_list = [str(residue.get_id()[1]) +
                           char_if_not_blank(residue.get_id()[2])
                           for residue in residue_list]
        if nterm_residue:
            self.nterm_resid = nterm_residue.get_id()[1]
            self.nterm_resname = nterm_residue.get_resname()
        else:
            self.nterm_resid = None
            self.nterm_resname = ""
        if cterm_residue:
            self.cterm_resid = cterm_residue.get_id()[1]
            self.cterm_resname = cterm_residue.get_resname()
        else:
            self.cterm_resid = None
            self.cterm_resname = ""


    def write_svg(self, fh):
        """
        Write this connector to the SVG file
                
        Parameters:
             fh - open filehandle to write SVG XML text to

        Return value:
            None

        """
        if self.directed:
            directed_str = "1"
        else:
            directed_str = "0"
        #  put first and last residue sequence numbers at left and right
        #  of hovertext string, with residues in sequence in between e.g.
        #  "134 ASP LYS ARG 136". For only single residue it will
        # be just like single-residue hovertext in shapes e.g. "ASP 134"
        # and for no residues (connector between two
        # adajcnet SSEs with no coil regino in between) we will put it like
        #  "(134-135)" indicating the two residue sequence numbers it joins.
        # TODO: have per-residue hovertext like helices and strands.
        (residue_names, residue_ids) = get_residue_strings(self.resname_list,
                                                           self.resid_list)
        if len(self.resname_list) == 0:
            if self.nterm_resid and self.cterm_resid:
                hovertext = '(' + str(self.nterm_resid) + '-' +\
                            str(self.cterm_resid) + ')'
            elif self.cterm_resid:
                hovertext = '(N-' + str(self.cterm_resid) + ')'
            else:
                hovertext = '(' + str(self.nterm_resid) + '-C)'
                
        elif len(self.resname_list) == 1:
            hovertext = self.resname_list[0] + " " + str(self.resid_list[0])
        else:
            hovertext = str(self.resid_list[0]) + " " + residue_names + " " +\
                        str(self.resid_list[-1])
        fh.write('  <dunnart:node id="' + str(self.xmlid) + '" ' +
                 'dunnart:srcID="' + str(self.src.xmlid) + '" ' +
                 'dunnart:dstID="' + str(self.dest.xmlid) + '" ' +
                 'dunnart:srcFlags="' + str(self.srcFlags) + '" ' +
                 'dunnart:dstFlags="' + str(self.dstFlags) + '" ' +
                 'dunnart:directed="' + directed_str + '" ' +
                 'dunnart:lineColour="' + self.color + '" ' +
                 PTGRAPH_NS + ':' + 'residueNames="' +
                 residue_names + '" ' +
                 PTGRAPH_NS + ':' + 'residueSeqNums="' +
                 residue_ids +
                 '" ' +
                 PTGRAPH_NS + ':' + 'hovertext="' + hovertext + '" '
                 'dunnart:type="connAvoidPoly"/>\n')
        



class PTSVGSeparation:
    """
    PTSVGSeparation represents information about separatino between indguides
    to be used with PTSVGSeparationConstraint.

    """

    def __init__(self, xmlid, direction, sepdistance, position):
        """
        Construct a PTSVGSeparation, given direction, separation
        distance, position  and id.

        Parameters:
            xmlid - unique XML identifier for this object
            direction - DUNNART_GUID_TYPE_HORI or _VERT
            sepdistance - separation distance
            position - position (x or y depending on HORI or VERT) for handle
        """
        self.xmlid = xmlid
        self.direction = direction
        self.sepdistance = sepdistance
        self.position = position

    def write_svg(self, fh):
        """
        Write this separation to the SVG file
                
        Parameters:
             fh - open filehandle to write SVG XML text to

        Return value:
            None

        """
        fh.write('  <dunnart:node dunnart:type="separation" ' +
                 'dunnart:direction="' + str(self.direction) + '" ' +
                 'dunnart:sepDistance="' +
                 str(self.sepdistance) + '" ' +
                'dunnart:position="' + str(self.position) + '" ' +
                 'id="' + str(self.xmlid) + '" ' +
                 'dunnart:equality="1" ' +
                 '/>\n')





##############################################################################
#
# Obsolete code below here. Separation (equality) constraints are no longer
# used in Dunnart (as of 0.15), now distribution constraints can be used
# instead.
#
##############################################################################
        
class PTSVGSeparationConstraint:
    """
    PTSVGSeparationConstraint represents a separation constraint, that is
    a separation between a pair of indguides.  It refers back to its
    PTSVGSeparation which specifies the separation distance, direction etc.
    """

    def __init__(self, indguide1, indguide2, sep):
        """
        Construct a PTSVGDistroConstraint given two indguides and a
        distribution.

        Parameters:
            indguide1 - PTSVGIndGuide for one object
            indguide2 - PTSVGIndGuide for the other object
            sep  - PTSVGSeparation object specifying the separation
        """
        self.indguide1 = indguide1
        self.indguide2 = indguide2
        self.sep = sep

    def write_svg(self, fh):
        """
        Write this constraint to the SVG file
                
        Parameters:
             fh - open filehandle to write SVG XML text to

        Return value:
            None

        """
        fh.write('  <dunnart:node dunnart:type="constraint" isMultiway="1" '\
                 'relType="separation" constraintID="' + str(self.sep.xmlid)
                 +'" '\
                 'objOneID="' + str(self.indguide1.xmlid) + '" ' \
                 'objTwoID="' + str(self.indguide2.xmlid) + '" />\n')


