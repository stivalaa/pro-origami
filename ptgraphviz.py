###############################################################################
#
# ptgraphviz.py - write H-bond or topology graphs with GraphViz
#
# File:    ptgraphviz.py
# Author:  Alex Stivala
# Created: July 2007
#
# $Id$
#
# Output H-bond or topology graph using GraphViz (for the -h, -n and -d options
# of ptgraph2).
#
#  pydot: python interface to GraphViz dot language (0.9.10)
#  http://dkbza.org/pydot.html
#
#  which in turn requires pyparsing
#  http://pyparsing.sourceforge.net/
#
#  For the -n and -d options, GraphViz itself is also required (2.12)
#  http://www.research.att.com/sw/tools/graphviz
#
###############################################################################

import pydot
###from ptgraph2 import PTGraph2
from ptnode import *

#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

def make_graphviz_graph(ptgraph):
    """
    Build a pydot (interface to GraphViz) representation of
    this PTGraph2 object. This can then be used to create the image
    representation with write_gif() etc. methods of the pydot Dot object
    returned.

    Parameters: ptgraph - PTGraph2 object to draw graph from

       Uses PTGraph2 member data: (readonly)
           nodelist - list of nodes
           stride_struct - StrideStruct representing the STRIDE output

            Note also uses (Readonly) member data in each node

    Return value:
        pydot Dot object representing the graph for GraphViz


    Precondition: node_list is sorted (by start res seq ascending);
                  this is done by build_graph_from_stride() before calling.


    """
    g = pydot.Dot()
    g.set_label(ptgraph.secstruct.pdb_header)

    g.set_overlap('scale') #overlap removal

    for nodelist in ptgraph.iter_chains():
        # add every node and an edge from each node to next in sequence
        prevnode = None
        for node in nodelist:
            dot_node = pydot.Node(node.nodeid)
            dot_node.set_fixedsize(True) # don't expand width to fit label
            if isinstance(node, PTNodeStrand):
                dot_node.set_shape('rect')
                dot_node.set_width(float(node.get_span()) / 5.0) # inches
                if node.get_sheet_id() != None:
                    dot_node.set_label(node.nodeid +
                                       " (sheet " + node.get_sheet_id()+")")
            elif isinstance(node, PTNodeHelix):
                dot_node.set_shape('ellipse')
                dot_node.set_width(float(node.get_span()) / 5.0) # inches
            g.add_node(dot_node)
            if prevnode != None:
                dot_edge = pydot.Edge(prevnode.nodeid, node.nodeid)
                g.add_edge(dot_edge)
            prevnode = node

        if ptgraph.use_hbonds:
            # add an edge for every hydrogen bond between structural elements
            for node in nodelist:
                for (other_node, rn1_unused, rn2_unused, dist) in  \
                 node.get_hbond_list():
                    g.add_edge(pydot.Edge(node.nodeid, other_node.nodeid,
                                          style='dotted',
                                          label=str(dist),
                                          len=dist/6.0,  # inches
                                          weight=100.0
                                          )
                               )
        else:
            # add an edge for bridges between beta strands
            for node in [ node for node in nodelist \
                          if isinstance(node, PTNodeStrand) ]:
                for (other_node, bdir, side) in node.get_bridge_list():
                    g.add_edge(pydot.Edge(node.nodeid, other_node.nodeid,
                                          style='dotted',
                                          label=bdir + ' ' + side,
                                          )
                               )



    return g
