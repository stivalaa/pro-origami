###############################################################################
#
# psvgnode.py - represents a node as a shape for writing as SVG (for Dunnart)
#
# File:    ptsvgnode.py
# Author:  Alex Stivala
# Created: February 2008
#
# $Id$
#
# an SVGNode contains information such as height, width, x,y position, color,
# etc. needed for writing a shape as SVG (XML) for the diagram editor.
#
###############################################################################

import sys

from ptnode import *
from color import rgb_tuple_to_hex_str,hex_str_to_rgb_tuple


#-----------------------------------------------------------------------------
#
# Constants
#
#-----------------------------------------------------------------------------

DUNNART_SPAN_LENGTH_FACTOR = 21  # multiply number of residues by this 
DUNNART_HELIX_SPAN_FACTOR  = 10  # for helices (above if for strands)
                                 # NB may need to adjust value for x/y offset
                                 # for SVG/firefox shape hack 
                                 # for hover text residue index in 
                                 # webserver/javascript/pohovertext.js if
                                 # this value is changed.
DUNNART_HELIX_WIDTH        = 30  # width of helices
DUNNART_STRAND_WIDTH       = 40  # width of strands
DUNNART_SHEET_GAP_SIZE     = 90 # space to leave between sheets
DUNNART_TERMINUS_SIZE      = 32  # width and height of terminus nodes
DUNNART_DOMAIN_GAP_SIZE    = 50  # space to leave between domain bounding boxes

DUNNART_TOP_PORT     = 49  # Dunnart input XML SVG code on flags for port
DUNNART_BOTTOM_PORT  = 50
DUNNART_LEFT_PORT    = 52
DUNNART_RIGHT_PORT   = 56
DUNNART_DEFAULT_PORT = 544 # default flags for no specified port

SHAPE_ALPHA = 1.0 # alpha channel (opacity) value for shapes

PTGRAPH_NS = "proorigami" # Namespace prefix for pro-origami properties
                          # TODO: should have DUNNART_NS instead of hardcoded
                          
#-----------------------------------------------------------------------------
#
# Class definitions
#
#-----------------------------------------------------------------------------

class PTSVGNode:
    """
    A PTSVGNode contains information such as height, width, x,y
    position, color, etc. needed for writing a shape as SVG (XML) for the
    diagram editor.
    
    We will end up using ineritance to inherit both from this and from the
    relevant PTNode (PTNodeHelix etc.) to get a class that has all the
    protein structural information (from the PTNode derived class) and
    augmented with graph layout information from this class.
    """

    def __init__(self):
        """
        Construct a PTSVGNode, intially all data is set to None.
        """
        self.color  = None         # Color for drawing shape for this node.
                                   # (r, g, b) tuple. Use get/set_color()
        self.color_hex = None      # color for drawing shape for this node.
                                   # hex string 'rrggbb'. Use get/set_color()
                                   # or get/set_color_hex().
                                   # set_color() and set_color_hex()
                                   # both will set both color and
                                   # color_hex.

        # set the following with set_svginfo (in derived classes) but just
        # just read them directly, sick of pointless get/set methods
        # (color only has one as I copied&pasted from PTNode before I created
        # this class and got sick of unnecessary get/set methods bloating code).

        self.xmlid  = None         # XML identifier for this shape.
        self.xpos   = None         # initial x co-ordinate for layout of diagram
        self.ypos   = None         # initial y co-ordinate for layout of diagram
        self.width  = None         # width of shape on diagram
        self.height = None         # height of shape on diagram
        self.label  = None         # label to put on shape
        self.indguide = None       # PTSVGIndGuide this node is aligned on
        self.nterm_port = None     # DUNNART_x_PORT for incoming connector
        self.cterm_port = None     # DUNNART_x_PORT for outgoing connector
        self.sseseqnum = None      # SSE sequential number
        self.headLabel = ''        # label for 'head' (C-terminal end) of shape
        self.tailLabel = ''        # label for 'tail' (N-terminal end) of shape

    def set_color(self, color):
        """
        Label this node with color (r,g,b) tuple. See comments in __init__

        Parameters: color - (r,g,b) color tuple
        Return value: None
        Uses data members (write): color, color_hex
        """
        self.color = color
        self.color_hex = rgb_tuple_to_hex_str(color)

    def get_color(self):
        """
        Get the value of the color tuple. See comments in __init__

        Parameters: None.
        Return value: (r,g,b) color tuple
        Uses data members (readnly): color
        """
        return self.color

    def set_color_hex(self, hexstr):
        """
        Label this node with color rrggbb hex string. See comments in __init__

        Parameters:  hexstr - rrggbb hex color string
        Return value: None
        Uses data members (write): color, color_hex
        """
        self.color = hex_str_to_rgb_tuple(hexstr)
        self.color_hex = hexstr
        
    def get_color_hex(self):
        """
        Get the hex color string. See comments in __init__

        Parameters: None
        Return value: hex color string 'rrggbb'
        Uses data members (readonly): color_hex
        """
        return self.color_hex
    

    def get_empty_port(self):
        """
        For a node that has only one port used, return the flags
        (top/bottom/left/right) for the port to use for the other connector
        (i.e. the opposite one to that used, so top if bottom used, etc.)

        Must only be used when the node has exactly one port already set,
        will raise exception otherwise.

        Parameters:
           None
        Retrun value:
           DUNNART_x_PORT value for the port to put second connector on.

        Uses data members (readonly):
           nterm_port, cterm_port
        Raises exceptions:
            ValueError if no or both ports are set, or either is
            DUNNART_DEFAULT_PORT (or invalid port value)
        """
        new_port = None
        if self.nterm_port != None and self.cterm_port != None:
            raise ValueError('both ports set')
        if self.nterm_port != None:
            used_port = self.nterm_port
        elif self.cterm_port != None:
            used_port = self.cterm_port
        else:
            raise ValueError('no port set')
        if used_port == DUNNART_TOP_PORT:
            new_port = DUNNART_BOTTOM_PORT
        elif used_port == DUNNART_BOTTOM_PORT:
            new_port = DUNNART_TOP_PORT
        elif used_port == DUNNART_LEFT_PORT:
            new_port = DUNNART_RIGHT_PORT
        elif used_port == DUNNART_RIGHT_PORT:
            new_port = DUNNART_LEFT_PORT
        else:
            raise ValueError('bad port value ' + str(used_port))
        return new_port


class PTSVGNodeHelix(PTNodeHelix, PTSVGNode):
    """
    A Helix PTNode augmneted with SVG information
    """
    def __init__(self, *args):
        """
        Construct PTSVGNodeHelix with supplied nodeid and type.
        Parameters:
        Variable parameter list: straight to PTNodeHelix constructor (q.v.).

        Raises exceptions:
            TypeError if helixtype argument is invalid.
        """
        PTNodeHelix.__init__(self, *args)
        PTSVGNode.__init__(self)


    def set_svginfo(self, xmlid, xpos, ypos, label, sseseqnum):
        """
        Set the SVG information in this node.
        
        Parameters:
             xmlid   -  the current XML id for this node
             xpos    - x coorindate to write helix
             ypos    - y coordinate to write helix
             label   - label to write on this node
             sseseqnum - SSE sequential number (from 1 to n, same as for
                         'sequential' labelling scheme).
        Return value:
             None
        Writes data members:
             xmlid
             xpos
             ypos
             label
             height
             width
             sseseqnum
        """
        self.xmlid = xmlid
        self.xpos = xpos
        self.ypos = ypos
        self.label = label
        self.height = self.get_span() * DUNNART_HELIX_SPAN_FACTOR
        # make sure we don't get helices drawn wrong way due to helix width
        # being greater than height
        if self.height < DUNNART_HELIX_WIDTH + 2: 
            self.height = DUNNART_HELIX_WIDTH + 2
        self.width = DUNNART_HELIX_WIDTH
        if self.get_sideways():
            temp = self.height
            self.height = self.width
            self.width = temp
        self.sseseqnum = sseseqnum
        

    def write_svg(self, fh):
        """
        Write this helix to the SVG file.

        Mostly there are Dunnart namespace attributes for Dunnart to use
        in layout out diagram, but also some proorigami namespace attributes
        that are ignored by Dunnart but retained in its finished diagram
        SVG output to be used by interactive SVG for hover text showing
        residues and so on.
        
        Parameters:
             fh - open filehandle to write SVG XML text to

        Return value:
            None
        """
        shape = 'bioHelix'
        # NB reversed used by Dunnart only for head and tail labels
        # but is also needed in interactive
        # SVG in webserver (pohovertext.js) for getting residue position
        # for hover text; so we use the proorigami namespace and the
        # dunnart namespace
        if self.get_reversed():
            reversed_str = "1"
        else:
            reversed_str = "0"
        clr = rgb_tuple_to_hex_str(self.get_color(), SHAPE_ALPHA)
        style = 'fill:#' + rgb_tuple_to_hex_str(self.get_color()) + ';'
        (residue_names, residue_ids) = get_residue_strings(self.resname_list,
                                                           self.resid_list)
        if self.get_cluster_id():
            clusterid = str(self.get_cluster_id())
        else:
            clusterid = ''
        if self.domain_id:
            domainid = self.domain_id
        else:
            domainid = ''
        fh.write('  <dunnart:node id="'+str(self.xmlid)+'" ' +
                 'dunnart:label="' + self.label + '" ' +
                 'dunnart:width="' + str(self.width) + '" ' +
                 'dunnart:height="' + str(self.height) +'" ' +
                 'dunnart:xPos="' + str(self.xpos) + '" ' +
                 'dunnart:yPos="' + str(self.ypos) + '" ' +
                 PTGRAPH_NS + ':reversed="' + reversed_str + '" ' +
                 'dunnart:reversed="' + reversed_str + '" ' +
                 'dunnart:fillColour="' + clr + '" ' + # color 
                 'style="' + style + '" ' +    # color
                 'dunnart:type="' + shape + '" ' +
                 PTGRAPH_NS + ':' + 'residueNames="' +
                 residue_names + '" ' +
                 PTGRAPH_NS + ':' + 'residueSeqNums="' +
                 residue_ids +
                 '" ' +
                 PTGRAPH_NS + ':' + 'helixType="' +
                 self.get_type().lower() + '" '  +
                 PTGRAPH_NS + ':' + 'sseseqnum="' +
                 self.sseseqnum + '" ' +
                 PTGRAPH_NS + ':' + 'chainId="' +
                 self.get_chainid() + '" ' +
                 PTGRAPH_NS + ':' + 'clusterId="' +
                 clusterid + '" ' +
                 PTGRAPH_NS + ':' + 'domainId="' +
                 domainid + '" ' +
                 # hovertext is updated by Javascript in interactive SVG
                 # as is selected flag
                 PTGRAPH_NS + ':' + 'hovertext="' + self.nodeid +'" ' +
                 PTGRAPH_NS + ':' + 'selected="0" ' +
                 'onmousemove="updateHoverText(evt)" ' +
                 'onclick="handleClickEvent(evt)"' + ' ' +
                 'dunnart:headLabel="' + self.headLabel + '" ' +
                 'dunnart:tailLabel="' + self.tailLabel + '" ' 
                 ' />\n')


class PTSVGNodeTerminus(PTNodeTerminus, PTSVGNode):
    """
    A Terminus PTNode augmented with SVG information
    """
    def __init__(self, *args):
        """
        Construct PTSVGNodeTerminus with supplied nodeid and type

        Parameters:
           variable parameter list: straight to PTNodeTerminus construtor (q.v.)
        """
        PTNodeTerminus.__init__(self, *args)
        PTSVGNode.__init__(self)
        

    def set_svginfo(self, xmlid, xpos, ypos, label):
        """
        Set the SVG information in this node.
        
        Parameters:
             xmlid   -  the current XML id for this node
             xpos    - x coorindate to write helix
             ypos    - y coordinate to write helix
             label   - label to write on this node
        Return value:
             None
        Writes data members:
             xmlid
             xpos
             ypos
             label
             height
             width
        """
        self.xmlid = xmlid
        self.xpos = xpos
        self.ypos = ypos
        self.label = label
        self.height = DUNNART_TERMINUS_SIZE
        self.width = DUNNART_TERMINUS_SIZE
        self.label = label


    def write_svg(self, fh):
        """
        Write this termimus to the SVG file.
        
        Parameters:
             fh - open filehandle to write SVG XML text to

        Return value:
             None
        """
        if self.domain_id:
            domainid = self.domain_id
        else:
            domainid = ''
        #shape='flowEndOfProc'
        shape='rect' # FIXME errors with flowEndOfProc, even without ports
        clr = rgb_tuple_to_hex_str(self.get_color(), SHAPE_ALPHA)
        style = 'fill:#' + rgb_tuple_to_hex_str(self.get_color()) + ';'
        fh.write('  <dunnart:node id="'+str(self.xmlid)+'" ' +
                 'dunnart:label="' + self.label + '" ' +
                 'dunnart:width="' + str(self.width) + '" ' +
                 'dunnart:height="' + str(self.height) +'" ' +
                 'dunnart:xPos="' + str(self.xpos) + '" ' +
                 'dunnart:yPos="' + str(self.ypos) + '" ' +
                 'dunnart:fillColour="' + clr + '" ' + # color 
                 'style="' + style + '" ' +    # color
                 'dunnart:type="' + shape + '" ' +
                 PTGRAPH_NS + ':' + 'chainId="' +
                 self.get_chainid() + '" ' +
                 PTGRAPH_NS + ':' + 'domainId="' +
                 domainid + '" ' +
                 PTGRAPH_NS + ':' + 'hovertext="chain ' + self.chainid + '" ' +
                 'onclick="handleClickEvent(evt)"' +
                 '/>\n')
        

class PTSVGNodeStrand(PTNodeStrand, PTSVGNode):
    """
    A Strand PTNode augmented with SVG information
    """
    def __init__(self, *args):
        """
        Construct PTSVGNodeStrand with supplied noeid.
        Parameters:
           Variable parameter list: straight to PTNodeStrand constructor (q.v.).
        """
        PTNodeStrand.__init__(self, *args)
        PTSVGNode.__init__(self)

    def set_svginfo(self, xmlid, xpos, ypos, label, indguide, sseseqnum):
        """
        Set the SVG information in this node.
        
        Parameters:
             xmlid   -  the current XML id for this node
             xpos    - x coorindate to write helix
             ypos    - y coordinate to write helix
             label   - label to write on this node
             indguide - PTSVGIndGuide this strand is aligned on
             sseseqnum - SSE sequential number (from 1 to n, same as for
                         'sequential' labelling scheme).
             
        Return value:
             None
        Writes data members:
             xmlid
             xpos
             ypos
             label
             height
             width
        """
        self.xmlid = xmlid
        self.xpos = xpos
        self.ypos = ypos
        self.label = label
        self.height = self.get_span() * DUNNART_SPAN_LENGTH_FACTOR
        # make sure we don't get strands drawn wrong way due to strand width
        # being greater than height
        if self.height < DUNNART_STRAND_WIDTH + 2: 
            self.height = DUNNART_STRAND_WIDTH + 2
        self.width = DUNNART_STRAND_WIDTH
        if self.get_sideways():
            temp = self.height
            self.height = self.width
            self.width = temp
        self.indguide = indguide
        self.sseseqnum = sseseqnum
        


    def write_svg(self, fh):
        """
        Write this strand to the SVG file.

        Mostly there are Dunnart namespace attributes for Dunnart to use
        in layout out diagram, but also some proorigami namespace attributes
        that are ignored by Dunnart but retained in its finished diagram
        SVG output to be used by interactive SVG for hover text showing
        residues and so on.
        
        Parameters:
             fh - open filehandle to write SVG XML text to

        Return value:
             None
        """
        shape = 'bioStrand'

        if self.get_reversed():
            reversed_str = "1"
        else:
            reversed_str = "0"
        clr = rgb_tuple_to_hex_str(self.get_color(), SHAPE_ALPHA)
        style = 'fill:#' + rgb_tuple_to_hex_str(self.get_color()) + ';'
        (residue_names, residue_ids) = get_residue_strings(self.resname_list,
                                                           self.resid_list)
        if self.get_sheet_id():
            sheetid = self.get_sheet_id()
        else:
            sheetid = ''
        if self.domain_id:
            domainid = self.domain_id
        else:
            domainid = ''
        fh.write('  <dunnart:node id="'+str(self.xmlid)+'" ' +
                 'dunnart:label="' + self.label + '" ' +
                 'dunnart:width="' + str(self.width) + '" ' +
                 'dunnart:height="' + str(self.height) +'" ' +
                 'dunnart:xPos="' + str(self.xpos) + '" ' +
                 'dunnart:yPos="' + str(self.ypos) + '" ' +
                 'dunnart:reversed="' + reversed_str + '" ' +
                 'dunnart:fillColour="' + clr + '" ' + # color 
                 'style="' + style + '" ' +    # color
                 'dunnart:type="' + shape + '" ' + 
                 PTGRAPH_NS + ':' + 'residueNames="' +
                 residue_names + '" ' +
                 PTGRAPH_NS + ':' + 'residueSeqNums="' +
                 residue_ids +
                 '" ' +
                 PTGRAPH_NS + ':' + 'sseseqnum="' +
                 self.sseseqnum + '" ' +
                 PTGRAPH_NS + ':' + 'chainId="' +
                 self.get_chainid() + '" ' +
                 PTGRAPH_NS + ':' + 'sheetId="' +
                 sheetid + '" ' +
                 PTGRAPH_NS + ':' + 'domainId="' +
                 domainid + '" ' + 
                 # hovertext is updated by Javascript in interactive SVG
                 PTGRAPH_NS + ':' + 'hovertext="' + self.nodeid + '" ' +
                 PTGRAPH_NS + ':' + 'selected="0" ' +
                 'dunnart:headLabel="' + self.headLabel + '" ' +
                 'dunnart:tailLabel="' + self.tailLabel + '" ' 
                 'onmousemove="updateHoverText(evt)" ' +
                 'onclick="handleClickEvent(evt)"' +
                 ' />\n')


        
        

#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

def get_residue_strings(resname_list, resid_list):
    """
    Give a list of residue names and list of residue PDB sequence identifiers,
    return string represnetation of the lists to put in SVG for use e.g
    in hovertext on interactive SVG.

    Parameters:
       resname_list - list of 3-letter residue names
       resid_list - list of residue PDB sequence numbers

    Return value:
       tuple (residue_names, residue_list) where
       residue_names is string with all residue names space-separated
       residue_list is string with all sequene numbers space-separarated
    """
    # can have empty residue name and id lists for eg connectors
    # that are between a terminus and an SSE node, and no coil region
    # between terminus and that SSE
    if len(resname_list) > 1:
        residue_names = reduce(lambda a,b : a + ' ' + b, resname_list)
        residue_ids = reduce(lambda a,b : str(a) + ' ' + str(b), resid_list)
    elif len(resname_list) == 1:
        residue_names = str(resname_list[0])
        residue_ids = str(resid_list[0])
    else:
        residue_names = ""
        residue_ids = ""
    return (residue_names, residue_ids)
