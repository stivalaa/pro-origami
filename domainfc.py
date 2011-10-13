#!/usr/bin/env python
###############################################################################
#
# domainfc.py - Domain decomposition using FastCommunity or MCL on SSE graph.
#
# File:    domainfc.py
# Author:  Alex Stivala
# Created: December 2007
#
# $Id$
#
# This program is an experimental implementation of domain decomposition
# using the FastCommunity algorithm on an undirected graph representing
# protein topology where nodes are secondary structure elements.
# This is similar to the graph used in DomainICA
# (Emmert-Streib and Mushegian 2007 "A topological algorithm for the
# identification of structural domains of proteins"
# BMC Bioinformatics 8:237).
#
# See docstring in main() or usage etc.
#
# FastCommunity is available from
#
# http://cs.unm.edu/~aaron/research/fastmodularity.htm
#
# where the usage and file formats are also described. Note that the
# fastcom.py script (in this directory) must be used instead of FastCommunityMH
# directly.
#
# The relevant publication to cite is
#
# A. Clauset, M.E.J. Newman and C. Moore, "Finding community structure
# in very large networks."  Phys. Rev. E 70, 066111 (2004).
#
#
# Instad of FastCommunity, the MCL clustering program can be used instead.
# MCL is available from
#
# http://micans.org/mcl/
#
# and the citation is (see website above for more citations)
#
# Stijn van Dongen. A cluster algorithm for graphs. Technical Report
# INS-R0010, National Research Institute for Mathematics and Computer
# Science in the Netherlands, Amsterdam, May 2000
# http://www.cwi.nl/ftp/CWIreports/INS/INS-R0010.ps.Z
# 
# MCL version 06-58 was used.
#
###############################################################################

import warnings # so we can suppress the annoying tempnam 'security' warning
import sys,os
import getopt
import re
from time import strftime,localtime
import resource # for getrusage()
from math import sqrt

from numpy.oldnumeric import *
from Bio.PDB import *

from ptnode import *
from ptdistmatrix import PTDistMatrix
import ptsecstruct
from ptdomain import *
import dographviz
from ptutils import cleanup_tmpdir,get_int_icode,biopdbresid_to_pdbresseq
from domeval import *
import getdomains

#-----------------------------------------------------------------------------
#
# Constants
#
#-----------------------------------------------------------------------------

DEFAULT_DISTANCE_THRESHOLD = 15.0 # Angstroms

PEPTIDE_BOND_LENGTH = 1.33  # Angstroms

FASTCOM  = "fastcom.py" # FastCommunity wrapper script (in PATH)
MCL_PROG = "mcl"        # MCL binary (in PATH)

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

#
# Real classes
#
class PTGraphD:
    """
    The topology graph consists of a sequence of structure (helix, strand)
    nodes with sequence edges in and out of them in sequence from N-terminus to
    C-terminus and adjacency edges for SSEs less than a threshold distance
    apart.

    Note there may be multiple such sequences (one for each
    chain).

    Also the nodes are all labelled with start and end residue sequence
    numbers, and node types etc. but this is not used at all in the
    actual community/clustering procedure on the graph, it is only
    used in calculating the distance matrix (which needs to get residues
    in the node to calculate their distances using PDB data) and
    for annotation in the graph reprsentation, and because this
    code was reused from another program (ptraph2.py) which does
    require the node labelling.
    """

    #
    # member functions
    #

    def __init__(self, pdb_structure, pdbid, dist_threshold,
                 include_310_helices = False, include_pi_helices = False,
                 add_loop_nodes = False):
        """
        Construct empty PTGraphD. To build the graph call
        build_graph_from_secstruct().

        Parameters:
            pdb_structure - parsed PDB structure from Bio.PDB
            pdbid - PDB identifier
            dist_threshold - distance threshold (Ansgtroms) for adding
                              spatial edges.
            include_310_helices - include 3_10 helices in the graph if True
            include_pi_helices - include pi_helices in the graph if True
            add_loop_nodes - include nodes for loop regions between SSEs if True

        """
        self.pdb_struct = pdb_structure
        self.pdbid = pdbid
        self.pdb_resid_dict = None # dict of { {chainid,pdb_resseq) : seqindx }
                                   # where chainid and pdb_resseq make up
                                   # the PDB residue identifier, the pdb_resseq
                                   # being string resnum+icode if any e.g.
                                   # '60' or '60A', seqindx is the indiex
                                   # into sequential list of all residues
                                   # residue_list.
        self.residue_list = None   # list of all residues (for all chains)
                                   # in sequence, built by get_residue_list()
        self.dist_threshold = dist_threshold
        self.chain_dict = None  # Each value of the chain_dict is a
                                # List of nodes in order from N to C terminus
                                # so chain_dict is { chainid : node_list }
        self.seqnum_dict = {}   # dictionary of { node_sequence_number : node }
                                # mapping the sequential node ids required
                                # by FastCommunity back to PTNode objects.
        self.distmatrix = None  # PTDistMatrix built in build_dist_matrix
        self.include_310_helices = include_310_helices
        self.include_pi_helices = include_pi_helices
        self.add_loop_nodes = add_loop_nodes
        self.num_edges = 0      # total number of edges (sequence + spatial)
        self.adjmatrix = None   # adjacency matrix built by build_adjmatrix()
        self.edge_index_list = None # list of (i,j) tuples of nonzero entries
                                   # in adjmatrix built by build_adjmatrix()
        # weight for sequence edges used when
        # weighted edges are available (MCL)
        self.sequence_edgeweight = PEPTIDE_BOND_LENGTH
                                        
        
        
    def iter_chains(self):
        """
        This generator function iterates over all chains in this PTGraph.
        A chain is just a list of nodes so it yields a node list for each
        chain.

        Parameters: Nonde.
        Return value: YIELDs a node list.
        Uses data members (readony):
            chain_dict - dict of {chainid:node_list}
        """
        # FIXME: can we just 'return self.chain_dict.itervalues()' here?
        for nodelist in self.chain_dict.itervalues():
            yield nodelist
            

    def iter_nodes(self):
        """
        This generator function iterates over all the node in this PTGraph.

        Parameters: None
        Return Value: YIELDs a node.
        Uses data members: (readonly):
             chain_dict - dict of {chainid_node_list}
        """
        for nodelist in self.iter_chains():
            for ptnode in nodelist:
                yield ptnode
                

    def build_graph_from_secstruct(self, secstruct, chainid, domain=None):
        """
        Build the list of nodes from the the supplied PTSecStruct
        object. Add edges for hydrogen bonds between structural
        elements also.

        Parameters:
            secstruct - PTSecStruct (ptsecstruct.py) object to build from
            chainid - If not None, use only this chain.
            domain - (default None). If not None, use only SSEs and residues
                     in the domain (PTDomain object)

        Uses member data (write):
            chain_dict - dict of { chainid : node_list } where node_list is
                          list of nodes in order, built in this function
            secstruct - keeps a pointer to the supplied secstruct
            seqnum_dict - dict of { node_seq_num : node } mapping the sequential
                          node numbers back to PTNode objects.
            num_edges - total number of edges in graph
            chainid - the chainid this graph is for (or None)

          (readonly):
            pdb_struct - The Bio.PDB parsed PDB struct (atomic co-ordinates)
                         for this protein.
            include_310_helices, include_pi_helices - if true, include
                         these kinds of helices.
            add_loop_nodes - If True, add nodes for
                           loop regions (ie between SSEs)



        Raises exceptions:
           NoSSE_Exception if no helices or strands found
        
        Return value:
            None.
            
        """

        self.secstruct = secstruct
        self.chainid = chainid
        
        num_helices = 0
        num_strands = 0
        
        nodenum = 0  # arbitrary enumeration starting 0, used by FastCommunity
        # Node these nodenum values are used as the node.seqnum and
        # they are used differently from ptgraph2 - these node numbers
        # are sequential and unique and are required to be so in this
        # use - for instance they are used as indices of the adjancecy
        # matrix.


        #
        # Build dictionary mapping (chainid, pdb_resid) to index in residue_list
        # for ALL residues, not just those in this domain.
        #
        self.residue_list = self.get_residue_list(self.pdb_struct,
                                                  PTDomain(None, None))
        self.pdb_resid_dict = {}
        seq_indx = 0
        while seq_indx < len(self.residue_list):
            residue = self.residue_list[seq_indx]
            self.pdb_resid_dict[( ptsecstruct.pdb_chainid_to_stride_chainid(
                                                residue.get_full_id()[2]), 
                                  biopdbresid_to_pdbresseq(
                                              residue.get_id()) )] = seq_indx
            seq_indx += 1

        
        self.chain_dict = {} # dict of {chainid : node_list}

        for (start_chainid, start_resnum, end_chainid, end_resnum, helixtype) \
              in secstruct.helix_list:
            assert(start_chainid == end_chainid) #helix must be same chain
            if (chainid != None and chainid != start_chainid):
                continue
            if (domain != None and
                not domain.is_in_domain(start_chainid,
                                        get_int_icode(start_resnum)[0])):
                continue
            if helixtype == "H":
                idprefix = "ALPHAHELIX_"
                htype = "ALPHA"
            elif helixtype == "I":
                idprefix = "PIHELIX_"
                htype = "PI"
                if not self.include_pi_helices:
                    continue
            elif helixtype == "G":
                idprefix = "310HELIX_"
                htype = "310"
                if not self.include_310_helices:
                    continue
            else: # shouldn't happen
                sys.stderr.write("ERROR: bad helix type " + helixtype+"\n")
            ah_node = PTNodeHelix(htype,
                                  idprefix + start_chainid+"_" +\
                                  str(num_helices),
                                  nodenum,
                                  start_resnum, end_resnum, start_chainid,
                                  self.residue_list, self.pdb_resid_dict)
            nodenum += 1
            num_helices += 1
            if not self.chain_dict.has_key(start_chainid):
                self.chain_dict[start_chainid] = []
            self.chain_dict[start_chainid].append(ah_node)

        for (start_chainid, start_resnum, end_chainid, end_resnum) \
                in secstruct.strand_list:
            assert(start_chainid == end_chainid) # must be in same chain
            if (chainid != None and chainid != start_chainid):
                continue
            if (domain != None and
                not domain.is_in_domain(start_chainid,
                                        get_int_icode(start_resnum)[0])):
                continue
            bs_node = PTNodeStrand("STRAND_"+start_chainid +"_"+\
                                   str(num_strands),
                                   nodenum,
                                   start_resnum, end_resnum, start_chainid,
                                   self.residue_list, self.pdb_resid_dict)
            num_strands += 1
            nodenum += 1
            if not self.chain_dict.has_key(start_chainid):
                self.chain_dict[start_chainid] = []

            self.chain_dict[start_chainid].append(bs_node)


        # raise an exception if there are no SSEs at all 
        if num_helices == 0 and num_strands == 0:
            raise NoSSE_Exception


        delete_chainid_list = [] # list of chainids to delete from chain_dict
        for (chainid, nodelist) in self.chain_dict.iteritems():
            # sort in order of start residue id ascending (all must be disjoint)
            nodelist.sort()

            if len(nodelist) < 1:
                # There are no SSEs in this chain, get rid of it.
                sys.stderr.write('WARNING: no SSEs in chain ' + chainid +
                                 '; chain ignored\n')
                delete_chainid_list.append(chainid) # don't delete while in loop
                continue
            else:
                # Check for chain with only SSEs that will not be used
                # (i.e. pi or 310 helices), and delete those too
                found_useful_node = False
                for ptnode in nodelist:
                    if isinstance(ptnode, PTNodeStrand):
                        found_useful_node = True
                        break
                    elif isinstance(ptnode, PTNodeHelix):
                        if ptnode.get_type() == "ALPHA":
                            found_useful_node = True
                            break
                        elif ((ptnode.get_type() == "310" and
                                 self.include_310_helices) or
                                (ptnode.get_type() == "PI" and
                                 self.include_pi_helices)):
                            found_useful_node = True
                            break
                if not found_useful_node:
                    sys.stderr.write('WARNING: only pi or 310 helices in chain '
                                     + chainid +
                                     '; chain ignored\n')
                    delete_chainid_list.append(chainid)
                    continue
                
            # (note PDB residue numbers don't necessarily start at 1,
            #  may be higher)
            lowest_res_seq = get_int_icode(nodelist[0].get_start_res_seq())[0]
            highest_res_seq = get_int_icode(nodelist[-1].get_end_res_seq())[0]


        # delete chains from chain_dict that were marked earlier for deletion
        for chainid in delete_chainid_list:
            self.chain_dict.pop(chainid)

        # add loop nodes (coil regions between SSEs) if requested
        if self.add_loop_nodes:
            self.add_loop_nodes_to_graph(nodenum)
                
        # buld the dictionary mapping sequential numbers to PTNode objects
        for node in self.iter_nodes():
            self.seqnum_dict[node.seqnum] = node

        # flag the first and last nodes in each chain
        for nodelist in self.iter_chains():
            nodelist[0].set_endchain(True)
            nodelist[-1].set_endchain(True)
            self.num_edges += len(nodelist) - 1


    def add_loop_nodes_to_graph(self, start_nodenum):
        """
        Add in nodes reprsenting the loop (aka coil) regions in between
        SSEs (helices and strands).

        Parameters:
           start_nodenum - the node number to start at
        Return value:
           New next node number
        Uses data members: (read/write)
            chain_dict - dict of { chainid : node_list } where node_list is
                          list of nodes in order (modifies node lists)
        Preconditions:
            For each chain the nodelist is sorted by residue sequence number
            ascending.
        """
        num_loops = 0
        nodenum = start_nodenum
        for nodelist in self.iter_chains():
            nodelist = list(nodelist) # copy as we modify the real one
            # nodelist is already sorted by res seq num ascending
            prevnode = None
            for node in nodelist:
                chainid = node.get_chainid()
                if prevnode != None:
                    loop_start = prevnode.get_end_res_seq() + 1
                    loop_end = node.get_start_res_seq() - 1
                    if loop_end >= loop_start: # may be no coil (adjacent SSEs)
                        loop_node = PTNodeLoop("LOOP_"+chainid +"_"+\
                                               str(num_loops),
                                               nodenum,
                                               loop_start,
                                               loop_end,
                                               chainid)
                        self.chain_dict[chainid].append(loop_node)
                        nodenum += 1
                        num_loops += 1
                prevnode = node
        # need to sort the nodelists again now that nodes are added.
        for nodelist in self.iter_chains():
            nodelist.sort()
        return nodenum
            
    def get_residue_list(self, pdb_struct, domain, getchainid = None):
        """
        Return list of Bio.PDB Residue objects in this domain, and optionally
        in the specified chain.,

        Parameters:
             pdb_struct - Bio.PDB parsed PDB struct for the protein
             domain -  PTDomain (ptdomain.py) object listing the segment(s)
                         that make up this domain (only one domain processed at a
                         time).
             getchainid - chain identifier to get residues in (default None -
                       all chains).

        Return value:
             list of Bio.PDB Residue objects in the domain (and optionally chain).
        Raises exceptions:
           NoSSE_Exception for empty structure (happens eg on d1oayi_.ent)

        """
        residue_list = []
        try:
            pdb_model = self.pdb_struct[0] # TODO always using model 0 for now
        except KeyError:
            raise NoSSE_Exception
        for chain in pdb_model:
            chainid = ptsecstruct.pdb_chainid_to_stride_chainid(chain.get_id())
            if getchainid and getchainid != chainid:
                continue # this is not the chain we want

            # Build a list of Bio.PDB Residue objects that are in this
            # domain.
            # id of a residue in Bio.PDB is tuple (hetatm, resseqnum, icode)
            # so we choose those where residue PDB number
            # (in the current chain) is in the domain.
            # TODO: maybe should use polypeptide builder for this instead
            # (and indeed should probably use it right from the beginning) -
            residue_list += [ residue for residue in chain.get_unpacked_list()
                              if is_aa(residue) and
                              domain.is_in_domain(chainid, residue.get_id()[1])
                            ]
            if getchainid:
                break # if getchainid specified, we now have it so can quit
        return residue_list


    def build_dist_matrix(self, domain=None):
        """
        Build the PTDistMatrix member ptdistmat containing the
        residue and SSE distance maps. This is built using information
        from the Bio.PDB Residue objects contained in the
        member pdb_struct for residues in the supplied domain which
        we are working with, and (for SSEs) the secondary structures
        defined by node lists under the chain_dict member built by
        build_graph_from_sec_struct().

        Parameters:
           domain - PTDomain object listing segments(s) that make up
                    this domain, only these are considered.
                    Default None (all residues consdiered).
           
        Uses member data (write):
           distmatrix - the PTDistanceMatrix class containing residue and
                       SSE distance maps.
          (readonly):
            chain_dict - dict of { chainid : node_list } where node_list is
                          list of nodes in order
            secstruct - keeps a pointer to the supplied secstruct
            pdb_struct - The Bio.PDB parsed PDB struct (atomic co-ordinates)
                         for this protein.

        Return value:
            None.
   
        """
        residue_list = []
        pdb_model = self.pdb_struct[0] # TODO always using model 0 for now
        for chain in pdb_model:
            chainid = ptsecstruct.pdb_chainid_to_stride_chainid(chain.get_id())
            # Build a list of Bio.PDB Residue objects that are in this
            # protein.
            # id of a residue in Bio.PDB is tuple (hetatm, resseqnum, icode)
            # so we choose those where HETATM is blank and is an amino acid
            # and is in the domain.
            # FIXME: should use Bio.PDB Polypetpide PPBuilder for this instead
            # (and indeed should useit and operate on Polypeptide objects
            # instead right from the start)
            residue_list += [ residue for residue in chain.get_list()
                              if residue.get_id()[0] == ' ' and
                              is_aa(residue) and
                              (domain == None or
                              domain.is_in_domain(chainid, residue.get_id()[1]))
                            ]

        # Also build list of all PTNodes
        ptnode_list = []
        for nodelist in self.iter_chains():
            for node in nodelist:
                ptnode_list.append(node)
            
        self.distmatrix = PTDistMatrix(residue_list, ptnode_list,
                                       None,
                                       self.pdb_struct)

    def add_spatial_edges(self, dist_threshold):
        """
        Add edges between distinct nodes that are not also adjacent in sequence,
        that are less than the dist_threshold
        apart. Uses the distmatrix to compute this.
        Note these edges are cnosidered undirected and we only add the one
        from from lower enumerated node to higher enumerated node, and not
        vice versa.

        Precondition:
           build_dist_matrix() has been previously called.

        Parameters:
           dist_threshold - (float) distance below or equal to which
                            two nodes have an edge added (Ansgtroms).

        Uses data members (readonly):
           distmatrix - the distance matrix
           nodelist - the list of nodes (via iter_chains())
           (write):
           num_edges - total number of edges in graph     
           NB writes to PTNodes  (adding edges via add_closenode())

        """
        nodelist = list(self.iter_nodes())
        for i in range(len(nodelist)-2):
            if nodelist[i+1].get_chainid() == nodelist[i].get_chainid():
                j0 = i + 2  # +2 so don't include next in chain
            else:
                j0 = i + 1
            for j in range(j0, len(nodelist)):
                dist = self.distmatrix.get_sse_distance(nodelist[i],
                                                        nodelist[j])
                if (dist < dist_threshold):
                    nodelist[i].add_closenode(nodelist[j], dist)
                    self.num_edges += 1


    def build_adjmatrix(self):
        """
        Build the adjacency matrix for the graph. The adjacency matrix
        is a symmetric binary nxn matrix
        (note binary: we do not use edge weights at
        all in the matrix) where entry (i,j) is 1 if node i and j are
        connected else 0. i and j are the node.seqnum values in [0, n-1]
        i.e. row i of the adjmatrix is the entries for the neighbours of
        the node with seqnum i.

        Also build a list of indices ((i,j) tuples) of nonzero entries
        in the adjacency matrix. This is a little redundant since
        the adjacency list information (the closenode_list in each node)
        also stores (most of) it - but not the sequence neighbour edges.
        So it is more convenient to build this list with all nonzero
        entry indices in it once and for all. Note that in this list
        we only store one of (i,j) and (j,i), since the matrix is
        symmetric it actually has twice this many nonzero entries.
        
        Paramaters:
           None
           
        Return value:
           None
           
        Uses data members:
          (write):
             adjmatrix - the adjacency matrix which is created
             edge_index_list - list of indices in adjmatrix
                         with nonzero entries
          (readonly):
             nodelist - list of nodes (via iter_nodes())
               and data in each node

        Preconditions:
           The graph hasalready been built with build_graph_from_secstruct()
            and add_spatial_edges().
        """
        nodelist = list(self.iter_nodes())
        n = len(nodelist)
        self.adjmatrix = zeros((n, n), Int)
        edgecount = 0
        self.edge_index_list = []
        for node_i in nodelist:
            for (node_j, dist) in node_i.get_closenode_list():
                self.adjmatrix[node_i.seqnum, node_j.seqnum] = 1
                self.adjmatrix[node_j.seqnum, node_i.seqnum] = 1
                edgecount += 1
                self.edge_index_list.append((node_i.seqnum, node_j.seqnum))

        # that just added the spatial edges. We also need to add
        # the implicit sequence edges.
        for nodelist in self.iter_chains():
            prevnode = None
            for node in nodelist:
                if prevnode != None:
                    self.adjmatrix[node.seqnum, prevnode.seqnum] = 1
                    self.adjmatrix[prevnode.seqnum, node.seqnum] = 1
                    edgecount += 1
                    self.edge_index_list.append((node.seqnum, prevnode.seqnum))
                prevnode = node

        assert(edgecount == self.num_edges)
        assert(len(self.edge_index_list) == edgecount)


    def build_matrices(self, do_build_adjmatrix, domain=None):
        """
        Build the matrices in this PTGraphD object.
        Just calls othermember functions ot do this
        Parameters:
           do_build_adjmatrix - only build adjacency matrix if True
           domain - PTDomain for domain to build matrice for, or None
               to build for all residues. Default None.
          
        Return value:
           None.
        Uses data members (write):
           distmatrix - distance matrix
           adjmatrix - adjacency matrix
           NOTE also writes to nodes in nodelist (adding edges)
           (readonly):
           use_modularity_cut - only build adjacency matrix if True
           dist_threshold - threshold for spatial edges (Angstroms)

        """
        # build the PTDistMatrix distance maps
        self.build_dist_matrix(domain)

        # and add the endges between nodes dist less than threshold
        self.add_spatial_edges(self.dist_threshold)

        # build the graph adjacency matrix
        if do_build_adjmatrix:
            self.build_adjmatrix()

        
        
    def write_pairsfile(self, fh, use_weights=False):
        """
        Write the edges in the .pairs file as used by FastCommunity,
        or MCL with the --abc option.
        The nodes are numbered from 0 and each line is an undirected
        edge, with two nodes spearated by a tab.
        Optionally (for MCL only, FastCommunity must be unweighted) if
        the use_weights parameter is True, write a third field which is
        the edge weight (sequence_edgeweight data member for sequence
        edges, else distance between SSEs for spatial edges).

        Parameters:
           fh - filehandle open for write to write .pairs output to
           weights - If True, add weights to edges (default False)
        Return value: None
        Uses data members:
           chain dict and node list(via iter_nodes())
           sequence_edgeweight - weight to add to sequence edges
           
        """
        # first write sequence edges
        for nodelist in self.iter_chains():
            prevnode = None
            for node in nodelist:
                if prevnode != None:
                    fh.write(str(prevnode.seqnum) + '\t' + str(node.seqnum))
                    if use_weights:
                        fh.write('\t' + str(self.sequence_edgeweight))
                    fh.write('\n')
                prevnode = node

        # now spatial edges
        for node1 in self.iter_nodes():
            for (node2, dist) in node1.get_closenode_list():
                fh.write(str(node1.seqnum) + '\t' + str(node2.seqnum))
                if use_weights:
                    fh.write('\t' + str(dist))
                fh.write('\n')

    def parse_fastcom_groups(self, fh):
        """
        Read the output of the FastCommunity program run over the pairs
        file, and convert the groups of sequential node numbers to
        a list of lists representing the clustering, each list being
        a list of PTNodes.

        Parameters:
           fh - filehandle open for reading of the FastCommunity output

        Return value:
           list of list of PTNodes, each list being a group according to fastcom

        Uses data members (readonly):
           seqnum_dict - Dict { node_seqnum : node } mapping the sequntial
                         node numbers to PTNode objects.
        """
        # FastCommunity output looks like:
        #
        # GROUP[ 6 ][ 9 ]
        # 6
        # 10
        #...
        # etc. Each GROUP starts a new group with each  id on a new
        # line. The first num in brackets is first id in group, secnod
        # is number of ids in group.
        #
        grouplist = []
        nodelist = []
        for line in fh:
            if line[:5] == "GROUP":
                if nodelist != []:
                    grouplist.append(nodelist)
                    nodelist = []
            else:
                line = line.rstrip('\n')
                if line.isdigit():
                    seqnum = int(line)
                    nodelist.append(self.seqnum_dict[seqnum])
                else:
                    sys.stderr.write('bad FastCommunity output line: ' + line+'\n')
        if nodelist != []:
            grouplist.append(nodelist)
        return grouplist


    def parse_mcl_clusters(self, fh):
        """
        Read the output of MCL 
        a list of lists representing the clustering, each list being
        a list of PTNodes.

        Parameters:
           fh - filehandle open for reading of the FastCommunity output

        Return value:
           list of list of PTNodes, each list being a group according to fastcom

        Uses data members (readonly):
           seqnum_dict - Dict { node_seqnum : node } mapping the sequntial
                         node numbers to PTNode objects.
        """
        # MCL output is just one line per cluster, with each line
        # consisting of tab-delimited labels (the sequential ids) in the cluster
        clusterlist = []
        nodelist = []
        for line in fh:
            for seqnum in line.split():
                nodelist.append(self.seqnum_dict[int(seqnum)])
            clusterlist.append(nodelist)
            nodelist = []
        return clusterlist


    def build_chaindict(self, nodelist):
        """
        Build a local chaindict { chainid : nodelist } for each different
        chain in the input nodelist paraemeter

        Parameters:
          nodelist - list of nodes to build chain dictionary for
        Return value:
          dict of { chainid : nodelist }
        Uses no data members.
        """
        chaindict = {}
        for node in nodelist:
            chainid = node.get_chainid()
            if chaindict.has_key(chainid):
                chaindict[chainid].append(node)
            else:
                chaindict[chainid] = [node]
        return chaindict


    def coalesce_nodelist(self, nodelist):
        """
        we are only having SSEs (not loops) as nodes
        so there are lots of segments with gaps in domain -
        coalesce all the segments together from adjacent nodes

        Parameters:
           nodelist - list of PTNodes

        Uses data members:
           chain_dict - dict of { chainid : node_list } where node_list is
                          list of nodes in order

        Return value
           List of PTSegments with start/end res seq nums consisting of the
           start and end of maximal subsequences of adjacent PTNodes in the
           input nodelist
        """
        # Build a local chaindict { chainid : nodelist } for each different
        # chain in the input nodelist paraemeter
        chaindict = self.build_chaindict(nodelist)

        # Now process each chain separately.
        # sort the list of nodes in the chain (sorts by res seq num)
        # and for each adjacent node, merge them as part of one segment
        # if there are no other nodes between them in entire
        # nodelist (ie. the member data nodelist)
        segment_list = []
        for chainid in chaindict.keys():
            in_nodes = chaindict[chainid]
            in_nodes.sort()
            all_nodes = self.chain_dict[chainid] # already sorted
            i = 0
            i0 = 0
            while (i < len(in_nodes) - 1):
                # FIXME: using index() is very inefficient, should use dict
                if (all_nodes.index(in_nodes[i]) + 1 ==
                    all_nodes.index(in_nodes[i+1])):
                     i += 1
                else:
                    # found a break, make this a segment now
                    segment = PTSegment(chainid,
                                        get_int_icode(in_nodes[i0].get_start_res_seq())[0],
                                        get_int_icode(in_nodes[i].get_end_res_seq())[0])
                    segment_list.append(segment)
                    i0 = i + 1
                    i += 1
            # make the remainder (possibly everything) a segment
            segment = PTSegment(chainid,
                                get_int_icode(in_nodes[i0].get_start_res_seq())[0],
                                get_int_icode(in_nodes[i].get_end_res_seq())[0])
            segment_list.append(segment)
        return segment_list


    def nodelist_to_segmentlist(self, nodelist):
        """
        Convert list of nodes representing a domain into a list of
        PTSegment objects. This is used when we are including loop regions
        as nodes so every residue is represented in (exactly one) node.
        
        Parameters:
           nodelist - list of PTNodes

        Uses data members:
           chain_dict - dict of { chainid : node_list } where node_list is
                          list of nodes in order

        Return value
           List of PTSegments with start/end res seq nums consisting of the
           start and end of maximal subsequences of adjacent PTNodes in the
           input nodelist
        """
#        print str([str(node) for node in nodelist])
        
        # Build a local chaindict { chainid : nodelist } for each different
        # chain in the input nodelist paraemeter
        chaindict = self.build_chaindict(nodelist)

        # Now process each chain separately.
        # sort the list of nodes in the chain (sorts by res seq num)
        # and for each adjacent node, merge them as part of one segment
        # if the end res seq of the first is immediately precending the
        # start res seq of the second.
        segment_list = []
        for chainid in chaindict.keys():
            in_nodes = chaindict[chainid]
            if len(in_nodes) < 2:
                segment = PTSegment(chainid,
                                    get_int_icode(in_nodes[0].get_start_res_seq())[0],
                                    get_int_icode(in_nodes[0].get_end_res_seq())[0])
                segment_list.append(segment)
                continue
                
            in_nodes.sort()
            i = 1
            i0 = 0
            while (i < len(in_nodes)):
                if (in_nodes[i-1].get_end_res_seq() + 1 ==
                    in_nodes[i].get_start_res_seq()):
                    i += 1
                else:
                    # found a break, make this a segment now
                    segment = PTSegment(chainid,
                                        get_int_icode(in_nodes[i0].get_start_res_seq())[0],
                                        get_int_icode(in_nodes[i-1].get_end_res_seq())[0])
                    segment_list.append(segment)
                    i0 = i
                    i += 1
            # make the remainder (possibly everything) a segment
            segment = PTSegment(chainid,
                                get_int_icode(in_nodes[i0].get_start_res_seq())[0],
                                get_int_icode(in_nodes[i-1].get_end_res_seq())[0])
            segment_list.append(segment)
        return segment_list


    def calculate_modularity(self, node_group_dict, adjmatrix = None):
        """
        Calculate the modularity Q as defined by
        
         A. Clauset, M.E.J. Newman and C. Moore, 'Finding community structure
         in very large networks.'  Phys. Rev. E 70, 066111 (2004).

         which is

         Q = \frac{1}{2m} \sum{vw}(A_{vw} - \frac{k_v k_w}{2m}) \delta(c_v,c_w)

        where m is the number od edges, v,w are nodes, k is node degree
        and c_v and c_w and the communities to which v and w belong resp.
        and \delta is the Kronecker delta function.

        for the community (domain) decomposition represented by the
        node_group_dict which maps each PTNode to a group number (0, 1, ...)

        Parameters:
           node_group_dict - dict of { PTnode : groupnumber } mapping
             each PTNode to a community (or domain or group or cluster...)
             numbered from 0.
           adjmatrix - (default None). If not None, the adjacency matrix
              to use, otherwise the adjmatrix member is used.
           
        Return value:        return Q
           The modularity value Q in [0,1] resulting from the 2 domain
           decomposition by cutting there
        Uses data members:
           nodelist - the list of nodes (via iter_chains())
           adjmatrix - the adjacency matrix
           num_edges - total number of edges in graph
        """
        if adjmatrix != None:
            A = adjmatrix
        else:
            A = self.adjmatrix
            
        two_m = self.num_edges * 2
        Q = 0.0
        for v in self.iter_nodes():
            for w in self.iter_nodes():
                if node_group_dict[v] == node_group_dict[w]:
                    Q += (float(A[v.seqnum, w.seqnum]) -
                          float(v.get_degree() * w.get_degree()) / float(two_m))
        Q /= float(two_m)
        return Q


    def jackknife_modularity_stddev(self, node_group_dict):
        """
        Estimate the standard deviation of modularity by jackknife procedure.
        One edge at a time is removed from the graph and the modularity
        calculated for the graph without that edge. Then an estimate
        of the variance can be obtained from

        sigma**2 = ((n-1)/n) \sum{i=1}{n} (Q_{(i)} - Q_{(\cdot)})

        where Q_{(i)} is the value of Q with edge i removed
        and Q_{(\cdot)} = (1/n) \sum{i=1}{n} Q_{(i)}, the mean value of
          modularity Q over all single-edge removals.

        Parameters:
           node_group_dict - dict of { PTnode : groupnumber } mapping
             each PTNode to a community (or domain or group or cluster...)
             numbered from 0.
           
        Return value:  
           The jacknnife estimate of variance of Q
           
        Uses data members (readonly):
           adjmatrix - the adjacency matrix
           num_edges - total number of edges in graph
           nodelist - list of nodes (via iter_nodes())
           edge_index_list - list of (i,j) indicies of nonzero adjmatrix values

        """
        A = array(self.adjmatrix)  # copy of adjmatrix, we can modify it
        assert(shape(A)[0] == shape(A)[1])
        num_nodes = shape(A)[0]
        n = self.num_edges
        Qvec = zeros(n, Float) # Qvec[i] is value of Q with edge i removed
        i = 0
        for (k,l) in self.edge_index_list:
            saved_Aij = A[k,l]
            assert(A[l,k]==saved_Aij)
            A[k,l] = 0
            A[l,k] = 0
            Qvec[i] = self.calculate_modularity(node_group_dict, A)
            A[k,l] = saved_Aij
            A[l,k] = saved_Aij
            i += 1
        assert(i == n)
        Qmean = add.reduce(Qvec) / n
        sigma_squared = (float(n-1)/n) * add.reduce((Qvec - Qmean)**2)
        return sqrt(sigma_squared)
        
        
    def calculate_modularity_1cut(self):
        """
        Compute a list of modularity values for the two-community (one cut)
        decomposotion, between each one of the adjacent in sequence nodes.

        Note that for multichain proteins each chain is in N to c terminus
        order, one chain after another so cuts can occur between chains
        in that case.
        
        Parameters:
           None
        Return value:
           list of Q values, for cut position just before each node in
           iter_nodes() (sequence N to C terminus) order. Note this
           means the last value in the list is actually the value for no cut.
    
        """
#         ## DEBUG ##
#         n = len(list(self.iter_nodes()))
#         edgecount = 0
#         for i in range(n):
#             for j in range(n):
#                 if self.adjmatrix[i, j] > 0:
#                     edgecount += 1
#         print edgecount,self.num_edges,'ccccc'
#         assert(edgecount/2 == self.num_edges)
#         ## END DEBUG ##
                
        # TODO we could no doubt do this more efficiently rather than
        # recomputing Q for each different cut along the backbone.

        # build dict mapping each node to group 0
        node_group_dict = dict([(node, 0) for node in self.iter_nodes()])

        # now calculate for 1 cut at each position along backbone
        # (note the very last value is in fact the value for 0 cuts).
        Qlist = []
        for node in self.iter_nodes():
            node_group_dict[node] = 1
            Qlist.append( self.calculate_modularity(node_group_dict) )

        return Qlist


    def calculate_RatioCut(self, node_group_dict, k):
        """
        Calculate the RatioCut

            RatioCut(A_0,...A_n) = \sum{i=1}{k} \frac{ cut(A_i, \not A_i)}
                                                     { |A_i| }

        as descirbed in von Luxburg 2006 'A Tutorial on Spectral Clustering'
        Tech. Rep. TR-149 Max-Planck-Institut fur Biologische Kybernetik
        (citing Hagen and Kahng 1992)

        cut(A,B) = \sum_{i \in A, j \in B} w_{ij}
        i.e. the sum of edge weights between clusters A and B.

        for the community (domain) decomposition represented by the
        node_group_dict which maps each PTNode to a group number (0, 1, ...)

        Parameters:
           node_group_dict - dict of { PTnode : groupnumber } mapping
             each PTNode to a community (or domain or group or cluster...)
             numbered from 0.
           k - the number of clusters (group numbers in node_group_dict
               numebered 0,...,k-1)
           
        Return value:       
           The value of RatioCut for the breakdown of nodes into clustesr
           described by node_group_dict

        Uses data members:
           nodelist - the list of nodes (via iter_chains())
           adjmatrix - the adjacency matrix

        """
        A = self.adjmatrix
        nodelist = list(self.iter_nodes())
        sumcutval = 0.0
        for i in range(k):
            groupsize = 0
            cutval = 0.0
            for v in nodelist:
                if node_group_dict[v] == i:
                    groupsize += 1
                for w in nodelist:
                    if (node_group_dict[v] == i and node_group_dict[w] != i):
                       cutval += A[v.seqnum, w.seqnum]
            cutval /= groupsize # the normalization by size of cluster
            sumcutval += cutval
        return sumcutval


    def calculate_RatioCut_1cut(self):
        """
        Compute a list of RatioCut values for the two-community (one cut)
        decomposotion, between each one of the adjacent in sequence nodes.

        Note that for multichain proteins each chain is in N to c terminus
        order, one chain after another so cuts can occur between chains
        in that case.

        In general clustering, minimizing RatioCut is an NP-hard problem and
        is approximated using eigenvectors of graph Laplactian (spectral
        clustering). But in these small instances and only cutting along
        backbone we are computing all values and minimizing exactly.
        
        Parameters:
           None
        Return value:
           list of RatioCut values, for cut position just before each node in
           iter_nodes() (sequence N to C terminus) order. Note this
           means the last value in the list is actually the value for no cut.
    
        """
        # TODO we could no doubt do this more efficiently rather than
        # recomputing RatioCut for each different cut along the backbone.

        # build dict mapping each node to group 0
        node_group_dict = dict([(node, 0) for node in self.iter_nodes()])

        # now calculate for 1 cut at each position along backbone
        # (note the very last value is in fact the value for 0 cuts,
        # and RatioCut cannot be computed there as then would div by zero
        # as one cluster has size 0).
        nodelist = list(self.iter_nodes())
        RatioCut_list = []
        for i in range(len(nodelist)-1):
            node = nodelist[i]
            node_group_dict[node] = 1
            RatioCut_list.append( self.calculate_RatioCut(node_group_dict, 2) )

        return RatioCut_list
            
#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------



def make_graph(pdb_filename,
               pdb_struct,
               secstruct_program,
               dist_threshold,
               chainid=None,
               use_dot=False, use_neato=False,
               include_310_helices = False,
               include_pi_helices = False,
               add_loop_nodes = False,
               do_build_adjmatrix = False):
    """
    For the supplied filemame, read PDB format data from that file
    and create graph for that structre.


    Paramteters:
       pdb_filename - filename of PDB file to read
       pdb_struct - Bio.PDB parsed PDB structure
       secstruct_program - secondary structure definition program
                       ('stride' or 'dssp' or 'pdb') to use.
       dist_threshold - float (Angstroms) two nodes have an edge added between
                        them if distance between is leq this.
       chainid - (default None). If not None, use only this chain.
       use_dot    - If True use dot from GraphViz to make PostScript output
       use_neato  - If True use neato from GraphViz to make PostScript output
       include_310_helices - if True, include 3_10 helices in the graph
       include_pi_helices - if True, include pi helices in the graph
       add_loop_nodes - if True, include loop (coil regions in the graph)
       do_build_adjmatrix - If True build adjacency matrix for graph.

    For use_dot or use_neato,
    writes .ps or files named as per description in main() below
    WARNING: overwrites these .ps files

    Return value: The PTGraphD object that has been constructed.
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

    ptg = PTGraphD(pdb_struct, pdbid, dist_threshold,
                   include_310_helices, include_pi_helices,
                   add_loop_nodes)

    # build PTGraphD from secondary structure
    try:
        ptg.build_graph_from_secstruct(secstruct, chainid)
    except NoSSE_Exception:
        if chainid == None:
            chainname = ''
        else:
            chainname = ' chain ' + chainid
        sys.stderr.write('WARNING: No helices or strands found in ' +
                         pdbid + chainname +
                         ': quitting\n')
        return None


    if verbose:
        for nodelist in ptg.iter_chains():
            for node in nodelist:
                sys.stderr.write(str(node) + '\n')
        
    ptg.build_matrices(do_build_adjmatrix)
    
    if (use_dot or use_neato):
        # use GraphViz to make PostScript output file
        dg = dographviz.make_graphviz_graph(ptg.secstruct.pdb_header,
                                            ptg.iter_chains())
        outfilename = pdbid + '.ps'
        if use_dot:
            progname = "dot"
        else:
            progname = "neato"
        sys.stdout.write('writing file ' + outfilename + '\n')
        dg.write_ps(outfilename,prog=progname) # or use ps2 not ps to make pdf later
    return ptg


def get_fastcom_domains(ptg):
    """
    For the supplied PTGraphD protein SSE graph,
    get groups from FastCommunity and
    convert this to domain decomposition in residues.


    Paramteters:
       ptg - The PTGraph2 protein SSE graph.
    Return value:
       list of PTDomain
    """
    # to avoid problems with blocking on pipes, we'll use a tempfile for input
    # to fastcom.py rather than writing to its stdin with popen2()
    tempfile_name = os.tempnam() # too bad we can't specify suffix only prefix
    tempfile_name += '.pairs' # as we need to dodgily do this for fastcom
    ptg.write_pairsfile(open(tempfile_name, 'w'))
    try:
        command = FASTCOM + " " + tempfile_name
        if verbose:
            sys.stderr.write('running ' + command + '...\n')
        fc_stdout = os.popen(command)
        # now parse the output back into a list of list of PTNodes
        grouplist = ptg.parse_fastcom_groups(fc_stdout)
        # and convert them to PTDomain objects and output them from from
        domainlist = nodegroups_to_domains(grouplist, ptg)
    finally:
        os.unlink(tempfile_name)

    return domainlist


def get_qcut_domains(ptg):
    """
    for the supplied PTGraphD protein SSE graph, get groups from the
    modularity cut algorithm.

    This algorithm is to compute the modularity (Q) value for a single
    cut at each point between SSEs along the backbone, and choose
    to make a cut at the max Q value position if it is higher than
    the Q value for no cut.
    Then this is done recursivly for  th resulting domains.

    Parameters:
       ptg - the PTGraphD protein SSE graph.
    Return value:
       list of PTDomain
    Uses globals:
           qfile_fh - if not None, open fh for write to write Q values to

    """
    global qfile_fh
    
    nodelist = list(ptg.iter_nodes())
        
    if (len(nodelist) < 2):  # can't further decompose with only 1 SSE
        return nodegroups_to_domains([nodelist], ptg) # single domain

    Qlist = ptg.calculate_modularity_1cut()

    if qfile_fh != None:
        for q in Qlist:
            qfile_fh.write(str(q) + '\n')
        qfile_fh.close()
        qfile_fh = None

    Q_nocut = Qlist[-1] # last val is Q for no cut ('cut' after last node)
    Qmax = 0.0
    Qmax_index = 0
    for i in range(len(Qlist)):
        if Qlist[i] > Qmax:
            Qmax = Qlist[i]
            Qmax_index = i
    maxQnode = nodelist[Qmax_index]
    
    if verbose:
        sys.stderr.write('q_nocut is ' + str(Q_nocut) + '\n')
        sys.stderr.write('max Q is ' + str(Qmax) + ' at node ' +  str(maxQnode)
                         + '\n')
#    if Qmax >  Q_nocut:

#     Q_mean = sum(Qlist) / len(Qlist)
#     Q_stddev = sqrt(sum([(Qi - Q_mean)**2 for Qi in Qlist]) / len(Qlist))
#     SD_MULTIPLIER = 2
#     if Qmax > Q_mean + SD_MULTIPLIER *  Q_stddev:


    # mean filter the Q list for 'smoothing'
    smoothQlist = list(Qlist)
    for i in range(3,len(Qlist)-3):  # note 3,len-3 due to python range operator
        smoothQlist[i] = (Qlist[i-3] + Qlist[i-2] + Qlist[i-1] + Qlist[i] + Qlist[i+1] + Qlist[i+2] + Qlist[i+3])/7

    # find local maxima
    local_maxima = []  # list of (index, value) tuples
    for i in range(1,len(smoothQlist)-1):
        if (smoothQlist[i] > smoothQlist[i-1] and smoothQlist[i] > smoothQlist[i+1]):
            local_maxima.append((i, smoothQlist[i]))

    # find highest local maximum (need not be global maximum because of ends)
    if len(local_maxima) > 0:
        Qmax = 0.0
        Qmax_index = 0
        for i in range(len(local_maxima)):
            if local_maxima[i][1] > Qmax:
                Qmax = local_maxima[i][1]
                Qmax_index = local_maxima[i][0]
        maxQnode = nodelist[Qmax_index]
    else:
        Qmax_index = 0    
        
    if verbose:
        sys.stderr.write('local maxima:\n')
        for (indx, qval) in local_maxima:
            sys.stderr.write(str(indx) + ': ' + str(qval) + '\n')

    if Qmax_index > 0:
        domain_left = nodegroups_to_domains([nodelist[:Qmax_index+1]], ptg)[0]
        ptg_left = PTGraphD(ptg.pdb_struct, ptg.pdbid, ptg.dist_threshold,
                            ptg.include_310_helices, ptg.include_pi_helices,
                            ptg.add_loop_nodes)
        ptg_left.build_graph_from_secstruct(ptg.secstruct, ptg.chainid,
                                            domain_left)
        ptg_left.build_matrices(True, domain_left)

        domain_right = nodegroups_to_domains([nodelist[Qmax_index+1:]], ptg)[0]

        ptg_right = PTGraphD(ptg.pdb_struct, ptg.pdbid, ptg.dist_threshold,
                             ptg.include_310_helices, ptg.include_pi_helices,
                             ptg.add_loop_nodes)
        ptg_right.build_graph_from_secstruct(ptg.secstruct, ptg.chainid,
                                             domain_right)
        ptg_right.build_matrices(True, domain_right)

        return get_qcut_domains(ptg_left) + get_qcut_domains(ptg_right)
    else:
        return nodegroups_to_domains([nodelist], ptg) # single domain
                            
    

def get_RatioCut_domains(ptg):
    """
    for the supplied PTGraphD protein SSE graph, get groups from the
    RatioCut algorithm.

    This algorithm is to compute the RatioCut value for a
    cut at each point between SSEs along the backbone, and choose
    to make a cut at the min value position if it is less than
    the RatioCut value for no cut.
    Then this is done recursivly for  th resulting domains.

    Parameters:
       ptg - the PTGraphD protein SSE graph.
    Return value:
       list of PTDomain
    Uses globals:
           qfile_fh - if not None, open fh for write to write RatioCut values to

    """
    global qfile_fh
    
    nodelist = list(ptg.iter_nodes())
        
    if (len(nodelist) < 2):  # can't further decompose with only 1 SSE
        return nodegroups_to_domains([nodelist], ptg) # single domain

    RatioCut_list = ptg.calculate_RatioCut_1cut()

    if qfile_fh != None:
        for rcval in RatioCut_list:
            qfile_fh.write(str(rcval) + '\n')
        qfile_fh.close()
        qfile_fh = None

    RatioCutmin = float('inf')
    RatioCutmin_index = 0
    for i in range(len(RatioCut_list)):
        if RatioCut_list[i] < RatioCutmin:
            RatioCutmin = RatioCut_list[i]
            RatioCutmin_index = i
    minRatioCutnode = nodelist[RatioCutmin_index]
    
    if verbose:
        sys.stderr.write('min RatioCut is ' + str(RatioCutmin) + ' at node ' +  str(minRatioCutnode)
                         + '\n')

    if (RatioCutmin_index > 0 and
        RatioCutmin < 5.0): # FIXME: this is bogus, need some way to decide this
        domain_left = nodegroups_to_domains([nodelist[:RatioCutmin_index+1]], ptg)[0]
        ptg_left = PTGraphD(ptg.pdb_struct, ptg.pdbid, ptg.dist_threshold,
                            ptg.include_310_helices, ptg.include_pi_helices,
                            ptg.add_loop_nodes)
        ptg_left.build_graph_from_secstruct(ptg.secstruct, ptg.chainid,
                                            domain_left)
        ptg_left.build_matrices(True, domain_left)

        domain_right = nodegroups_to_domains([nodelist[RatioCutmin_index+1:]], ptg)[0]

        ptg_right = PTGraphD(ptg.pdb_struct, ptg.pdbid, ptg.dist_threshold,
                             ptg.include_310_helices, ptg.include_pi_helices,
                             ptg.add_loop_nodes)
        ptg_right.build_graph_from_secstruct(ptg.secstruct, ptg.chainid,
                                             domain_right)
        ptg_right.build_matrices(True, domain_right)

        return get_RatioCut_domains(ptg_left) + get_RatioCut_domains(ptg_right)
    else:
        return nodegroups_to_domains([nodelist], ptg) # single domain
                            
    



def get_mcl_domains(ptg, use_edgeweights=False):
    """
    For the supplied PTGraphD protein SSE graph,
    get cluster from MCL and
    convert this to domain decomposition in residues and write to fh.


    Paramteters:
       ptg - The PTGraph2 protein SSE graph.
       use_edgeweights - (default False) If True use weighted edges.
    Return value:
       list of PTDomain objects.
    """
    # Note we don't really have to use the sequential ids for MCL,
    # we could just write the node ids (which are also unique) directly
    # but it's easier to just be consistent with what we had to do for
    # FastCommunity.
    
    # to avoid problems with blocking on pipes, we'll use a tempfile for input
    # to MCL rather than writing to its stdin with popen2()
    tempfile_name = os.tempnam()
    ptg.write_pairsfile(open(tempfile_name, 'w'), use_edgeweights)
    if verbose:
        sys.stderr.writelines(open(tempfile_name).readlines())
    try:
        # tell MCL to run with output to stdout and input in label pairs format
        command = MCL_PROG + " " + tempfile_name + " --abc -o -"
        if verbose:
            sys.stderr.write('running ' + command + '...\n')
        else:
            command += " 2>/dev/null"
        fc_stdout = os.popen(command)
        # now parse the output back into a list of list of PTNodes
        grouplist = ptg.parse_mcl_clusters(fc_stdout)
        # and convert them to PTDomain objects and output them from from
        domainlist = nodegroups_to_domains(grouplist, ptg)
    finally:
        os.unlink(tempfile_name)

    return domainlist


def nodegroups_to_domains(grouplist, ptg):
    """
    Convert the supplied list of node lists (each list of nodes being a
    domain) into a list of PTDomain ojbects (see ptdomain.py) representing
    the domains.

    Parameters:
       grouplist - list of list of PTNodes, eacst)
        h list representing a domain
       ptg - the PTGraphD object the gruoplist was clustered from

    Return value:
       list of PTDomain objects representing the domain decomposition
    """
#    debug_write_groups(sys.stdout, grouplist)
    domainid = 0
    domain_list = []
    for nodelist in grouplist:
        if ptg.add_loop_nodes:
            segment_list = ptg.nodelist_to_segmentlist(nodelist)
        else:
            # we are only having SSEs (not loops) as nodes
            # so there are lots of segments with gaps in domain -
            # so coalesce all the segments together from adjacent nodes
            segment_list = ptg.coalesce_nodelist(nodelist)
        domain = PTDomain(str(domainid), segment_list)
        domain_list.append(domain)
        domainid += 1
    return domain_list
            



def debug_write_groups(fh, grouplist):
    """
    Write to fh the list of nodes in each group in the supplied grouplist.
    Used for debugging.
    
    Parameters:
       fh - open for write filehandle to write decomposition to
       grouplist - list of list of PTNodes, each list representing a domain

    Return value:
       None
    """
    fh.write(str(len(grouplist)) + ' domains:\n')
    for i in range(len(grouplist)):
        fh.write(str(i) + ': ')
        j = 0
        for node in grouplist[i]:
            if j > 0:
                fh.write(', ')
            fh.write(str(node))
            j += 1
        fh.write('\n')




def get_domainfc_domains(pdbid, pdb_filename, pdb_struct,
                         chainid,
                         dist_threshold,
                         use_mcl = True,
                         secstruct_program='pdb',
                         use_dot=False, use_neato=False,
                         include_310_helices = False,
                         include_pi_helices = False,
                         use_edgeweights = False,
                         add_loop_nodes = False,
                         use_modularity_cut = False,
                         use_ratio_cut = False):

    """
    Build SSE graph and get domains using clustering (MCL or FastCommunity).
    This function is suitable for passing to the domeval.py
    run_on_pdomains_list() higher-order function for benchmarking
    (the pdbid, pdb_filename, pdb_struct parameters are used there,
    others must be passed as additional args).

    Parameters:
       pdbid - PDB identifier
       pdb_filename - PDB file name
       pdb_struct - Bio.PDB parsed PDB structure
       chainid - if not None, the chain to use in order to restrict to one chain
       dist_threshold - float (Angstroms) two nodes have an edge added between
                        them if distance between is leq this.
       use_mcl - (Default True). If True use MCL else use FastCommunity
       secstruct_program - secondary structure definition program
                       ('stride' or 'dssp' or 'pdb') to use. Default 'pdb'.
       use_dot    - If True use dot from GraphViz to make PostScript output
       use_neato  - If True use neato from GraphViz to make PostScript output
       include_310_helices - if True, include 3_10 helices in the diagram
       include_pi_helices - if True, include pi helices in the diagram
       use_edgeweights - (default False) Use weighted edges. Only valid
                         if using MCL.
       add_loop_nodes - add nodes for loop regions
       use_modularity_cut - use the modularity cut algorithm not MCL or
                            FastCommunity
       use_ratio_cut - use ratio cut algorithm not MCL or FastCommunity or
                        modularity cut.

    Return value:
       list of PTDomain objects representing the domain decomposition.
    """
    # make graph(s) from PDB file
    ptg = make_graph(pdb_filename,
                     pdb_struct,
                     secstruct_program,
                     dist_threshold,
                     chainid,
                     use_dot, use_neato,
                     include_310_helices, include_pi_helices,
                     add_loop_nodes,
                     do_build_adjmatrix = (use_modularity_cut or use_ratio_cut))
    if ptg == None:
        return None
    if use_modularity_cut:
        domainlist = get_qcut_domains(ptg)
    elif use_mcl:
        domainlist = get_mcl_domains(ptg, use_edgeweights)
    elif use_ratio_cut:
        domainlist = get_RatioCut_domains(ptg)
    else:
        domainlist = get_fastcom_domains(ptg)
    return domainlist

#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------

# TODO: too many options now; instead of using -q,-m,-r etc. change to
# have a -m <method> or something (like -t struct_prog and -e domain_prog)
# Also change name of -Q option since used for other values as well.

def usage(progname):
    """
    Print usage message and exit
    """
    sys.stderr.write("Usage: " + progname +
            " [-35dlnvmqr] [-a disance_threshold] [-t struct_prog] [-Q filename] "
            "[-e domain_prog] <PDBfile [chainid] | -b pDomains_file PDBroot>\n")
    sys.stderr.write("  -a <distance_threshold> in Ansgstroms (default " +
                     str(DEFAULT_DISTANCE_THRESHOLD) + ")\n")
    sys.stderr.write("  -3 include 3_10 helices\n")
    sys.stderr.write("  -5 include pi helices\n")
    sys.stderr.write("  -d use GraphViz dot to draw graph\n")
    sys.stderr.write("  -e evaluate against domain program/db\n"
                     "     valid values are none (default), "
                     "ddomain, cath:cdffile, pdomains:pdomainsfile\n")
    sys.stderr.write("  -l add nodes for loop regions\n")
    sys.stderr.write("  -m use MCL instead of FastCommunity\n")
    sys.stderr.write("  -q use modularity cut algorithm instead of MCL or FastCommunity\n")
    sys.stderr.write("  -r use ratio cut instead of MCL or FastCommunity or modularity cut\n")
    sys.stderr.write("  -Q <qfilename> write modularity (Q) values to <qfilename?> when using -q\n")
    sys.stderr.write("  -n use GraphViz neato to draw graph\n")
    sys.stderr.write("  -t struct_prog : use struct_prog define " \
                     "secondary structure\n")
    sys.stderr.write("       supported is 'pdb' (default) or 'stride' or 'dssp'\n")
    sys.stderr.write("  -v print verbose debugging messages to stderr\n")
    sys.stderr.write("  -b is an alternate mode that runs over all chains in \n"
                     "     the specified pDomains file with PDBroot as the \n"
                     "     root of the divided PDB hierarchy\n")
    sys.exit(1)


def main():
    """
    main for domainfc.py

    Usage: domainfc [-35dlnvmqr][-t structprog] [-a distance_threshold]
                    [-Q qfilenmae ]
                    [-e domainprog] <PDBfile [chainid] | -b pdomains_file PDBroot>

    Output for -n or -d is a postscript file named pdb.ps where pdb is
    the pdb identifier from the PDB file supplied.

    chainid if specified causes only the specified chain to be procssed.
    by default all chains are processed and domains can contain segments
    from different chains, but many other methods and databases do not
    allow this, for for interoperation with them, this forces only
    one chain to be considered.

    -b Instead of specifying PDBfile and optional chainid, if the -b mode
       is used the process is run over all chains in the specified pDomains
       benchmark file using the PDBroot as the root of the PDB divided
       hierarchy to get PDB files.

    -3 specifies to include 3_10 helices in the diagram. Default is only
       alpha helices.

    -5 specifies to include pi helices in the diagram. Defaul is only
       alpha helices.

    -a distance_threshold specifies the threshold distance between SSEs
       to connect SSE nodes (if they are less than the threshold).
       In Angstroms, default is 15.0.
       
    -d specifies to use GraphViz dot to draw the graph.
       This outputs a PostScript file.

    -e evaluate the performance against another method/database.
       Valid values are 'none' (default), 'ddomain', 'cath:cdf_filename'.

    -l add nodes for loop (coil) regions i.e. those in between SSEs.
    
    -m use MCL instead of FastCommunity.

    -q use the modularity cut algorithm instead of MCL or FastCommunity.

    -r use ration cut algorithm instead of MCL or FastCommunity or
       modularity cut.

    -Q <filename> when using -q, write list of Q values at each cut position
       to <filename>. WRANING: <filename> is overwritten if it exists.
    
    -n specifies to use GraphViz neato to draw the graph.
       This outputs a PostScript file.

    -t specifies the secondary structure assignment program to use.
       Currently suppoed is 'pdb' and 'dfh,ssp' and 'stride'. Default 'pdb'.

    -v specifies verbose mode: debugging output is written to stderr.

    -w specifies to use weighted edges (MCL only). Edge weights are the
       distances between the SSEs.

    """
    global verbose
    global qfile_fh
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "35b:de:lmqrQ:na:t:vw?")
    except getopt.GetoptError:
        usage(os.path.basename(sys.argv[0]))

    valid_secstruct_programs = ["dssp", "stride", "pdb"]
    valid_domain_programs = getdomains.valid_domain_programs + [r"none"]
    valid_domain_programs_re = [ re.compile(re_str) for re_str in
                                 valid_domain_programs ]

    use_dot = False
    use_neato = False
    verbose = False # global (python globals are only 'global' to module though)
    qfile_fh = None # global filehand to write Q values to for -Q 
    secstruct_program = "pdb"
    include_310_helices = False
    include_pi_helices = False
    dist_threshold = DEFAULT_DISTANCE_THRESHOLD
    use_mcl = False
    eval_domain_program = "none"
    benchmark_mode = False
    use_edgeweights = False
    add_loop_nodes = False
    use_modularity_cut = False
    use_ratio_cut = False

    for opt,arg in opts:
        if opt == "-b": # use all chains in pDomains bechmark file
            benchmark_mode = True
            pdomains_filename = arg
        elif opt == "-3": # include 3_10 helices
            include_310_helices = True
        elif opt == "-5": # include pi helices
            include_pi_helices = True
        elif opt == "-a": # distance threshold
            dist_threshold = float(arg)
        elif opt == "-d": # use dot
            use_dot = True
        elif opt == "-e": # domain parsing program/db to evaluate against
            eval_domain_program = None
            for valid_domarg_re in valid_domain_programs_re:
                if valid_domarg_re.match(arg):
                    eval_domain_program = arg
                    break
            if eval_domain_program == None:
                sys.stderr.write("valid values for -e are: " +
                                 str(valid_domain_programs) + "\n")
                usage(sys.argv[0])
        elif opt == "-l":  # add nodes for loop regions
            add_loop_nodes = True
        elif opt == "-m": # use MCL not FastCommunity
            use_mcl = True
        elif opt == "-q": # use the modularity cut algorithm not FC or MCL
            use_modularity_cut = True
        elif opt == "-r": # use ratio cut not modularity or FC or MCL
            use_ratio_cut = True
        elif opt == "-Q": # write Q values for -q to the filename in arg
            qfile_fh = open(arg, "w")
        elif opt == "-n": # use neato
            use_neato = True
        elif opt == "-t":
            if arg not in valid_secstruct_programs:
                sys.stderr.write("valid values for -t are: " +
                                 str(valid_secstruct_programs) + "\n")
                usage(sys.argv[0])
            secstruct_program = arg
        elif opt == "-v": # verbose
            verbose = True # this module only
            ptnode_set_verbose(True) # ptnode module
            ptsecstruct.ptsecstruct_set_verbose(True) # ptsecstruct module
            ptdomain_set_verbose(True) # ptdomain module
        elif opt == "-w": # use weighted edges
            use_edgeweights = True
        else:
            usage(sys.argv[0])

    if use_dot and use_neato:
        sys.stderr.write(
            '-d (dot) and -n (neato) options are mutually exclusive\n')
        usage(sys.argv[0])

    if use_edgeweights and not use_mcl:
        sys.stderr.write('Edge weights (-w) can only be used with MCL (-m)\n')
        usage(os.path.basename(sys.argv[0]))

    xoptcount = 0
    if use_mcl:
        xoptcount += 1
    if use_modularity_cut:
        xoptcount += 1
    if use_ratio_cut:
        xoptcount += 1
    if xoptcount > 1:
        sys.stderr.write('Modularity cut (-q) and MCL (-m) and ratio cut (-r) are all pairwise mutually exclusive\n')
        usage(os.path.basename(sys.argv[0]))
        
    if benchmark_mode:
        if len(args) != 1:
            sys.stderr.write('benchmark (-b) mode requires pDomains file '
                             'and PDB divided hierarchy root.\n')
            usage(sys.argv[0])
        pdbroot = args[0]
        timestamp = strftime("%d%b%Y %H:%M:%S", localtime())
        sys.stdout.write("# Generated on " + timestamp
                         + " by domainfc.py $Revision$\n"
                         + "# " + " ".join(sys.argv) + "\n")

        (num_undercut, num_overcut, num_correct, avgscore,
         num_correct_assign, numdomains_dict, num_processed) = \
                       run_on_pdomains_file(pdbroot, pdomains_filename,
                                            True,
                                            get_domainfc_domains,
                                            dist_threshold, use_mcl,
                                            secstruct_program,
                                            use_dot, use_neato,
                                            include_310_helices,
                                            include_pi_helices,
                                            use_edgeweights,
                                            add_loop_nodes,
                                            use_modularity_cut,
                                            use_ratio_cut)
        # time.clock() is too problematic (CLOCKS_PER_SEC etc.) so we will
        # use getrusage() instead
        # (Although note this time includes time spent parsing files,
        # doing overlap calculations, etc.)
        user_cpu_time = (resource.getrusage(resource.RUSAGE_SELF).ru_utime + 
                         resource.getrusage(resource.RUSAGE_CHILDREN).ru_utime)
        sys.stdout.write('\nUser cpu seconds: ' + str(user_cpu_time) + '\n\n')
        numdomains_list = list(numdomains_dict.iteritems())
        numdomains_list.sort()
        for (num_domains, (freq,total,avg,undercut,overcut,domain_num_correct,
                           domain_assign_correct)) in numdomains_list:
             sys.stdout.write(str(num_domains) + ' domain proteins:\n')
             print_scores(freq, undercut, overcut, domain_num_correct,
                          domain_assign_correct, avg, indent=2)
             sys.stdout.write('\n')
        sys.stdout.write("\nTotals:\n")
        print_scores(num_processed,num_undercut, num_overcut, num_correct,
                     num_correct_assign, avgscore, indent=0)
    else:
        if len(args) == 2:
            chainid = args[1].upper()
            if len(chainid) != 1:
                sys.stderr.write('invalid chain identifier\n')
                usage(os.path.basename(sys.argv[0]))
        elif len(args) == 1:
            chainid = None
        else:
            usage(os.path.basename(sys.argv[0]))

        pdb_filename = args[0]

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
            domainlist = get_domainfc_domains(pdbid, our_pdb_filename,
                                              pdb_struct, chainid,
                                              dist_threshold,
                                              use_mcl, secstruct_program,
                                              use_dot, use_neato,
                                              include_310_helices,
                                              include_pi_helices,
                                              use_edgeweights,
                                              add_loop_nodes,
                                              use_modularity_cut,
                                              use_ratio_cut)
            if domainlist != None:
                write_domains(sys.stdout, domainlist)
                if eval_domain_program != "none":
                    print evaluate_domains(domainlist, eval_domain_program,
                                           pdbid,
                                           our_pdb_filename, pdb_struct,
                                           chainid)
        finally:
            if used_tmp_file:
                cleanup_tmpdir(TMPDIR)
    

if __name__ == "__main__":
    warnings.filterwarnings('ignore', 'tempnam', RuntimeWarning) 
    main()
