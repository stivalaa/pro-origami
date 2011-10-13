###############################################################################
#
# dographviz.py - draw protein topology graphs with GraphViz
#
# File:    dographviz.py
# Author:  Alex Stivala
# Created: July 2007
#
# $Id$
#
# Output topology graph using GraphViz (for -n or -d options of domainfc.py)
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
from ptnode import *

#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

def make_graphviz_graph(pdb_header, chain_iter):
    """
    Build a pydot (interface to GraphViz) representation of
    this PTGraphD object. This can then be used to create the image
    representation with write_gif() etc. methods of the pydot Dot object
    returned.

    Parameters:
         pdb_header - PDB header for graph label
         chain_iter - Iterable over all chains (each yields a list of PTNodes)

    Return value:
        pydot Dot object representing the graph for GraphViz


    Precondition: node_list is sorted (by start res seq ascending);
                  this is done by build_graph_from_stride() before calling.


    """
    g = pydot.Dot()
    g.set_type('graph') # undirected (other option is 'digraph')
    g.set_label(pdb_header)

    g.set_overlap('scale') #overlap removal

    for nodelist in chain_iter:
        # add every node and an edge from each node to next in sequence
        prevnode = None
        for node in nodelist:
            dot_node = pydot.Node(node.nodeid)
            dot_node.set_label(str(node.seqnum) + ': ' + node.nodeid)
            dot_node.set_fixedsize(True) # don't expand width to fit label
            if isinstance(node, PTNodeStrand):
                dot_node.set_shape('rect')
                dot_node.set_width(float(node.get_span()) / 5.0) # inches
            elif isinstance(node, PTNodeHelix):
                dot_node.set_shape('ellipse')
                dot_node.set_width(float(node.get_span()) / 5.0) # inches
            g.add_node(dot_node)
            if prevnode != None:
                dot_edge = pydot.Edge(prevnode.nodeid, node.nodeid,
                                      style = 'bold')
                g.add_edge(dot_edge)
            prevnode = node

        # add an edge for every node pair with dist less than threshold
        for node in nodelist:
            for (other_node, dist) in node.get_closenode_list():
                g.add_edge(pydot.Edge(node.nodeid, other_node.nodeid,
                                      )
                           )
    return g
