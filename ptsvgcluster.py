###############################################################################
#
# ptsvgcluster.py - represent a cluster, being a set of nodes grouped together
#
# File:    ptsvgcluster.py
# Author:  Alex Stivala
# Created: February 2008
#
# $Id$
#
# A PTSVGCluster contains a list of PTSVGNodes, and a shading color.
#
###############################################################################

from color import color_gradient, rgb_tuple_to_hex_str,DUNNART_DEFAULT_CLUSTER_FILL_COLOR
from ptsvgnode import PTSVGNode



#-----------------------------------------------------------------------------
#
# Class definitions
#
#-----------------------------------------------------------------------------

class PTSVGCluster:
    """
    A PTSVGCluster represents a group of PTSVGNodes to be clustered together
    with a shaded background to indicating the clustering.
    """

    def __init__(self, svgnodelist, xmlid,
                 clusterid,color_num = 0,
                 color=DUNNART_DEFAULT_CLUSTER_FILL_COLOR):
        """
        Construct the PTSVGCluster from list of PTSVGNodes (shapes).

        Parmaeters:
           svgnodelist - list of PTSVGNode objects
           xmlid - unique XML identifier for this cluster
           clusterid - some identifiying string to label cluster with
           color_num - integer 'color' marking distinct groups of clusters,
                       then translated to an actual color in the color member.
           color - color of background shading (defaults to default color),
                   as RGB hex string.
        """
        self.svgnodelist = svgnodelist
        self.xmlid = xmlid
        self.clusterid = clusterid
        self.color_num = color_num
        self.color = color

    def __str__(self):
        """
        Return string representation of node for debugging etc.
        """
        return "PTSVGCluster " + self.clusterid + \
               ' [' + str(self.color_num) + ']'

    def write_svg(self, fh):
        """
        Write Dunnart constraints in SVG format for 'clusters' of helices.
        Called by write_helix_cluster().

        Parameters:
             fh - open filehandle to write SVG XML text to

        Uses data members (readonly):
             svgnodelist, xmlid
        """
        i = 0
        cluster_member_xmlids = ""
        for node in self.svgnodelist:
            assert(isinstance(node, PTSVGNode))
            if cluster_member_xmlids != "":
                cluster_member_xmlids += " "
            cluster_member_xmlids += str(node.xmlid)
        cluster_fill_color = 'dunnart:fillColour="' +\
                             self.color +\
                             '"'
        fh.write('  <path id="' + str(self.xmlid) + '" ' +
                 'class="cluster" ' +
                 'dunnart:contains="' + cluster_member_xmlids + '" ' +
                 'dunnart:type="cluster" ' +
                 'dunnart:rectangular="1" ' +
                 cluster_fill_color + '/>\n')
        
