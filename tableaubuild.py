###############################################################################
#
# tableaubuild.py - class for building protein tableaux in Python
#
# File:    tableaubuild.py
# Author:  Alex Stivala
# Created: May 2008 (moved from pytableaucreate.py)
#
# $Id$
#
###############################################################################

"""
Build a protein tableau.
The implemntation is actually in pttableau.py which is used by ptgraph2.py
(Pro-Origami),

Also used to create SSE midpoint distance matrix.

Tableaux are described by Kamat and Lesk 2007
'Contact Patterns Between Helices and Strands of Sheet Define Protein
 Folding Patterns' Proteins 66:869-876
and Lesk 2003 'From Electrons to Proteins and Back Again'
Int. J. Quant. Chem. 95:678-682
and Lesk 1995 'Systematic representation of folding patterns'
J. Mol. Graph. 13:159-164.

The implementation is based on Arun Konagurthu's TableauCreator program, see
Konagurthu, Stuckey and Lesk 2008 'Structural search and retrieval using
a tableau representation of protein folding patterns' Bioinformatics
(advance access, to be published Jan 5 2008).

Filenames may be either in the format above or the pdbq1lp.pdb format.
Compressed pdb files are supported (gzip) (e.g. pdb1qlp.ent.gz).

It is written in Python and depends on some Python libraries:

. BioPython (including Bio.PDB)
  http://www.biopython.org

  Reference for Bio.PDB is:
  Hamelryck and Manderick 2003 "PDB parser and structure class implemented
  in Python" Bioinformatics 19:2308-2310

  which in turn depends on Numeric
  http://sourceforge.net/projects/numpy


Developed on Linux 2.6.9 (x86_64) with Python 2.5.1
and BioPython 1.43 with Numeric 24.2

"""

import warnings # so we can suppress the annoying tempnam 'security' warning
import sys,os
import getopt
import numpy.oldnumeric as Numeric
from Bio.PDB import *

from ptnode import *
import pttableau
import ptsecstruct
from ptdomain import *
from ptutils import cleanup_tmpdir,get_int_icode,biopdbresid_to_pdbresseq
from domeval import *
import getdomains
from ptdistmatrix import compute_sse_midpoint_dist_matrix

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
class TableauBuild:
    """
    The protein representation consists of a sequence of structure
    (helix, strand) nodes with sequence edges in and out of them in
    sequence from N-terminus to C-terminus and adjacency edges for
    SSEs less than a threshold distance apart.

    Note there may be multiple such sequences (one for each
    chain).

    Also the nodes are all labelled with start and end residue
    sequence numbers, and node types etc. but this is not used at all
    in the here, it is only included because this code was reused from
    another program (ptraph2.py) which does require the node
    labelling.
    """

    #
    # member functions
    #

    def __init__(self, pdb_structure, pdbid,
                 include_310_helices = False, include_pi_helices = False,
                 add_loop_nodes = False):
        """
        Construct empty TableauBuild. To build the structure call
        build_graph_from_secstruct().

        Parameters:
            pdb_structure - parsed PDB structure from Bio.PDB
            pdbid - PDB identifier
            include_310_helices - include 3_10 helices in the graph if True
            include_pi_helices - include pi_helices in the graph if True
            add_loop_nodes - include nodes for loop regions between SSEs if True

        """
        self.pdb_struct = pdb_structure
        self.pdbid = pdbid
        self.chain_dict = None  # Each value of the chain_dict is a
                                # List of nodes in order from N to C terminus
                                # so chain_dict is { chainid : node_list }
        self.seqnum2node = {}   # dictionary of {  seqnum : PTNode }
                                # maping int sequence numbers to PTNode objs
        self.tableau = None     # PTTableau build in build_tableau
        self.include_310_helices = include_310_helices
        self.include_pi_helices = include_pi_helices
        self.pdb_resid_dict = None # dict of { {chainid,pdb_resseq) : seqindx }
                                   # where chainid and pdb_resseq make up
                                   # the PDB residue identifier, the pdb_resseq
                                   # being string resnum+icode if any e.g.
                                   # '60' or '60A', seqindx is the indiex
                                   # into sequential list of all residues
                                   # residue_list.
        self.residue_list = None   # list of all residues (for all chains)
                                   # in sequence, built by get_residue_list()


        
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



    def build_graph_from_secstruct(self, secstruct, domain, chainid=None,
                                   ignore_insertion_codes=False):
        """
        Build the list of nodes from the the supplied PTSecStruct
        object. 


        Parameters:
            secstruct - PTSecStruct (ptsecstruct.py) object to build from
            domain - PTDomain (ptdomain.py) object listing the segment(s)
                     that make up this domain (only one domain processed at a
                     time).
                     (in/out) NOTE: may be modified by having a segment
                     added if SSE is only partly in domain.
            chainid - chain identifier to build graph for only this chain,
                      or None for all chains (default)
            ignore_insertion_codes - If True, a hack to make it work with
                      PMML (only) which does not report insertion codes
                      unlike DSSP and STRIDE

        Uses member data (write):
            chain_dict - dict of { chainid : node_list } where node_list is
                          list of nodes in order, built in this function
            secstruct - keeps a pointer to the supplied secstruct

          (readonly):
            pdb_struct - The Bio.PDB parsed PDB struct (atomic co-ordinates)
                         for this protein.
            include_310_helices, include_pi_helices - if true, include
                         these kinds of helices.

        Raises exceptions:
           NoSSE_Exception if no helices or strands found
        
        Return value:
            None.
            
        """

        self.secstruct = secstruct

        helix_num = 1 
        strand_num = 1

        num_helices_in_domain = 0
        num_strands_in_domain = 0

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
                                              residue.get_id(),
                                              ignore_insertion_codes)
                                                               )] = seq_indx
            seq_indx += 1
        
        # Note that now we are only adding elements in the supplied domain,
        # so the so-called 'chains' may really be segments, i.e. subsequences
        # of chains (rest of chain may be in other domain(s)

        self.chain_dict = {} # dict of {chainid : node_list}

        for (start_chainid, start_resnum, end_chainid, end_resnum, helixtype) \
              in secstruct.helix_list:
            assert(start_chainid == end_chainid) #helix must be same chain
            if chainid and chainid != start_chainid:
                continue # chainid specified, skip ones not in that chain
            # will consider structures in domain if first residue is in domain
            if domain.is_in_domain(start_chainid,
                                   get_int_icode(start_resnum)[0]):
                num_helices_in_domain += 1
                if helixtype == "H":
                    idprefix = "ALPHAHELIX_"
                    htype = "ALPHA"
                    this_helix_num = helix_num
                    helix_num += 1
                elif helixtype == "I":
                    if not self.include_pi_helices:
                        continue
                    idprefix = "PIHELIX_"
                    htype = "PI"
                    this_helix_num = helix_num
                    helix_num += 1
                elif helixtype == "G":
                    if not self.include_310_helices:
                        continue
                    idprefix = "310HELIX_"
                    htype = "310"
                    this_helix_num = helix_num
                    helix_num += 1
                else: # shouldn't happen
                    sys.stderr.write("ERROR: bad helix type " + helixtype+"\n")
                ah_node = PTNodeHelix(htype,
                                      idprefix + start_chainid+"_" +\
                                      str(this_helix_num),
                                      this_helix_num,
                                      start_resnum, end_resnum, start_chainid,
                                      domain.domainid,
                                      self.residue_list, self.pdb_resid_dict)
                if not self.chain_dict.has_key(start_chainid):
                    self.chain_dict[start_chainid] = []
                self.chain_dict[start_chainid].append(ah_node)

                # we must already have handled the case of SSEs that cross
                # domain boundaries (by moving whole SSE to one of the domains)
                assert( domain.is_in_domain(end_chainid,  get_int_icode(end_resnum)[0]) )

        for (start_chainid, start_resnum, end_chainid, end_resnum) \
                in secstruct.strand_list:
            assert(start_chainid == end_chainid) # must be in same chain
            if chainid and chainid != start_chainid:
                continue # chainid specified, skip ones not in that chain
            if domain.is_in_domain(start_chainid,
                                   get_int_icode(start_resnum)[0]):
                num_strands_in_domain += 1
                bs_node = PTNodeStrand("STRAND_"+start_chainid +"_"+\
                                       str(strand_num),
                                       strand_num,
                                       start_resnum, end_resnum, start_chainid,
                                       domain.domainid,
                                       self.residue_list,
                                       self.pdb_resid_dict)
                strand_num += 1
                if not self.chain_dict.has_key(start_chainid):
                    self.chain_dict[start_chainid] = []

                # we must already have handled the case of SSEs that cross
                # domain boundaries (by moving whole SSE to one of the domains)
                assert( domain.is_in_domain(end_chainid,  get_int_icode(end_resnum)[0]) )
                self.chain_dict[start_chainid].append(bs_node)


        # raise an exception if there are no SSEs at all in this domain
        if num_helices_in_domain == 0 and num_strands_in_domain == 0:
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
                # Check for chain with only SSEs that will not be drawn
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
                

        # delete chains from chain_dict that were marked earlier for deletion
        for chainid in delete_chainid_list:
            self.chain_dict.pop(chainid)

        # -------------------------------------------------------------------

        # This is needed only for labelling sheets for HH and KK codes
        # (see dfs_strands() etc. below)
        
        # add edges for hydrogen bonds
        # uses secstruct and chainid member data
        # these are used for determining which side bridge partners are
        # on (and also for drawing a hydrogen bond graph if requested)
        self.add_hbond_edges_from_secstruct()
        
        # add edges for bridge partners
        # uses secstruct and chainid member data
        self.add_bridge_edges_from_secstruct()

        #---------------------------------------------------------------------


        # for sequential numbering, we'll build this dictionary mapping
        # sequential number (note NOT restarting for each chain)
        # to PTNode
        # so that sequential numbers as used in ptgraph2 -b sequential
        # option.
        # this is a dictionary of { seqnum : PTNode }
        self.seqnum2node = {}
        for (seqnum, node) in \
            enumerate([node for node in self.iter_nodes() if \
                       not ( (isinstance(node, PTNodeTerminus)) or
                              (isinstance(node, PTNodeHelix) and
                               ( (node.get_type() == "310" and
                                  not self.include_310_helices) or
                                 (node.get_type() == "PI" and
                                  not self.include_pi_helices) ) ) ) ]):
                self.seqnum2node[seqnum+1] = node # start at 1 not 0

    # ------------------------------------------------------------------------

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


    # iter_strands(), dfs_strands(),
    # find_connected_components() and label_sheets() are needed to assign
    # strands to sheets in order for HH and KK codes to be used for strands
    # only when they are in the same sheet.
    # TODO this code is cut&pasted from ptgraph2.py, probably should
    # have a base class that this at PTGraph2 both inherit from or somethibng
    # rather than all this duplication.

    def add_hbond_edges_from_secstruct(self):
        """
        Add edges between structural elements for hydrogen bonds between
        those nodes. Called by build_graph_from_secstruct().

        NB: adds bonds between STRANDs only, not between HELIXes (helices).
        
        Parameters: None.
        Return value: None.
        Uses data members:
           readonly:
              secstruct - PTSecStruct object to get hbonds from
              chainid - chainid of chain in PTSecStruct to use
           read/write:
              chain_dict - dict by chainid of
                           list of nodes (changes node data, not list as such)

        Precondition: each nodelist in chain_dict
                      is sorted (by start res seq ascending);
                      this is done by build_graph_from_secstruct()
                      before calling.

        """
        hbond_list = self.secstruct.hbond_list
        # TODO: do this more efficiently using presorting (ie how it used to
        # be done when only one chain)
        for (chainid1, resnum1, chainid2, resnum2, dist) in hbond_list:
            for ptnode in self.iter_strands():
                if chainid1 == ptnode.get_chainid() and \
                   resnum1 >= ptnode.get_start_res_seq() and \
                   resnum1 <= ptnode.get_end_res_seq():
                    dest_node = self.find_node_containing_seqnum(resnum2,
                                                                 chainid2)
                    if dest_node != None and \
                           isinstance(dest_node, PTNodeStrand): # only STRANDs
                        ptnode.add_hbond(dest_node, resnum1, resnum2, dist)

    def add_bridge_edges_from_secstruct(self):
        """
        Add edges between strand nodes representing beta brdiges between
        those nodes (add just one edge between any two strands).
        Called by build_graph_from_secstruct().

        NB: adds bonds between STRANDs only, not between HELIXes (helices).

        Parameters: None.
        Return value: None.
        Uses data members:
           readonly:
              secstruct - PTSecStruct object to get hbonds from
              chainid - chainid of chain in PTSecStruct to use
           read/write:
              chain_dict - dict by chainid of
                            list of nodes (changes node data, not list as such)

        """

        bridge_list =  self.secstruct.bridgeres_list
        #         (chainid1, resnum1, chainid2, resnum2, bdir)
 
        # TODO: do this more efficiently using presorting (ie how it used to
        # be done when only one chain)

        for ptnode in self.iter_strands():
            for (chainid1, resnum1, chainid2, resnum2, bdir) in bridge_list:
                if chainid1 == ptnode.get_chainid() and \
                   resnum1 >= ptnode.get_start_res_seq() and \
                   resnum1 <= ptnode.get_end_res_seq():
                    try:
                        dest_node = self.find_node_containing_seqnum(resnum2,
                                                                     chainid2)
                    except KeyError:
                        dest_node = None
                        sys.stderr.write('WARNING: chain ' + chainid2 + \
                                         ' involved in beta bridge not found.'+\
                                         '\n  Probably due to domain parsing' +\
                                         ' breaking a beta sheet.\n')
                    if dest_node != None and \
                           isinstance(dest_node, PTNodeStrand): # only STRANDs
                        if ptnode == dest_node:
                            sys.stderr.write('WARNING: ignoring self-bridge ' +
                                             ptnode.nodeid + '\n')
                        else:
                            ptnode.add_bridge(dest_node, bdir)

    
    def iter_strands(self):
        """
        This generator function iterates over all strands in this PTGraph
        object. I.e. it yields a strand for each strand in the 
        node lists.

        Parameters: None.
        Return value: YIELDs a strand.
        Uses data members (readonly):
           self.chain_dict - dict of { chainid : list of nodes }
        """
        for nodelist in self.iter_chains():
            for ptnode in nodelist:
                if isinstance(ptnode, PTNodeStrand):
                    yield ptnode

    def find_node_containing_seqnum(self, res_seqnum, chainid):
        """
        Find and return node in node list for chain chainid
        containing supplied PDB residue
        sequence number.

        Parameters:
           res_seqnum - PDB residue sequence number to find node for
           chainid - chain identifier to find node in

        Return value:
           PTNode pointer of PTNode containing the supplied residue seq num
           in supplied chainid
           or None if the residue is not in a structural element PTNode

        Uses data members (readonly):
           chain_dict - chainid dict of list of PTNodes
        """
        # TODO: since node_list is sorted should use binary search here
        #       (maybe try the Python bisect module)
        if not self.chain_dict.has_key(chainid):
            return None # no such chain, can happen due to domain parsing
        for ptnode in self.chain_dict[chainid]:
            if ptnode.is_in_interval(res_seqnum):
                return ptnode
        return None

    def dfs_strands(self, start_strand, visited, dfs_list, from_node,
                    back_edge_list,
                    sheet_id=None):
        """
        Make a depth-first search traversal of STRAND nodes
        using bridge (not sequence)
        edges starting at the specfied strand.

        Parameters:
           start_strand - STRAND node to start at
           visited - (in/out) dictionary of {ptnode:True} visited nodes
           dfs_list - (in/out) list of ptnodes visited in dfs order
           from_node - node from which we are being (recursively) called
           back_edge_list - list of (node, node) tuples representing an
                             edge between the two nodes, which is a back
                             edge, i.e. from a node to an ancestor of that
                              node in the spanning tree. The back edge
                              means there is a cycle of which the back
                              edge forms a part.
           sheet_id - identifier of this sheet (connected component) to mark
                      each strand in it with, or None to not mark at all
                      (default).
                      

        Recursive function. call initially as
            dfslist = []
            back_edge_list = []
            dfs_strands(startnode, {}, dfslist, None, back_edge_list)

        Return value:
            None. (Output is dfs_list, back_edge_list parameters)
        
        Uses members (readonly):
           chain_dict - dict by chainid of list of PTNodes
           
        """
        visited[start_strand] = True
        if sheet_id != None:
            start_strand.set_sheet_id(sheet_id)
            #print 'xxx',str(start_strand),sheet_id
        dfs_list.append(start_strand)
        for (node, bdir_unused, side_unused) in start_strand.get_bridge_list():
            if node not in visited:
                self.dfs_strands(node, visited, dfs_list, start_strand,
                                 back_edge_list, sheet_id)
            elif node != from_node: #not parent of start_strand in spanning tree
                # don't add duplicate back edges
                # ((node1,node2) is same as (node2,node1))
                duplicate = False
                for (a,b) in back_edge_list:
                    if ((start_strand == a and node == b) or
                        (node == a and start_strand == b)):
                        duplicate = True
                        break
                if not duplicate:
                    if verbose:
                        sys.stderr.write('dfs_strands back edge from ' +
                                         str(start_strand) + ' to ' +
                                         str(node) +
                                         '\n')
                    back_edge_list.append((start_strand, node))



    def find_connected_components(self):
        """
        Find the connected components (considering only STRAND nodes
        and bridge [not sequence] edges in the graph).

        This is done by a DFS traversal at every node in the graph
        (skipping already visited ones), giving us the partition of
        the graph into connected components.
        
        Parameters: None

        Uses member data: 
            chain_dict - dict by chainid of list
                        of PTNodes in the graph (modifies PTNodes not list)

            (WRITE):

             sheet_dict - 
              dictionary of { sheet_id : ptnode_list } where sheet_id is 'A',
              'B', etc. and ptnode_list is a list of PTNodeStrand instances
               in that connected component (sheet).

               self.sheet_backedges_dict  -
                 dict of {sheet_id : ((node1,node2))}
                 listing 'back edges' i.e. edges
                 to an ancestor in DFS spanning tree
                 in the connected component (sheet).
                 note (node1,node2) and (node2,node1)
                 are the same (undirected graph) and
                 only one of the two is present in the


        Labels each strand node with the sheet id it belongs to as it goes.
        """

        sheet_id = 'A' # sheet id is single alpha char A, B, etc.
                       # (will be a problem for more than 26 sheets... eg
                       # this actually happens on 2J28), wrap to lowercase
                       
        visited = {}   # dictionary of {ptnode : True} visited nodes
        back_edge_list = []  # list of (ptnode, ptnode) tuples for back edges
        self.sheet_dict = {} # dictionary of {sheet_id : nodelist}
        self.sheet_backedges_dict = {} # dict of {sheet_id : ((node1,node2))}
                                       # listing 'back edges' i.e. edges
                                       # to an ancestor in DFS spanning tree
                                       # in the connected component (sheet).
                                       # note (node1,node2) and (node2,node1)
                                       # are the same (undirected graph) and
                                       # only one of the two is present in the
                                       # list.
        for node in self.iter_strands():
            if node not in visited:
                connected_node_list = []
                back_edge_list = []
                self.dfs_strands(node, visited, connected_node_list, None,
                                 back_edge_list,
                                 sheet_id)
                self.sheet_dict[sheet_id] = list(connected_node_list)
                self.sheet_backedges_dict[sheet_id] = list(back_edge_list)
                sheet_id = chr(ord(sheet_id)+1)
                if sheet_id == '[':
                    sheet_id = 'a' # if go past Z, wrap to lowercase

    
    def label_sheets(self):
        """
        Label strands with sheet id to which each belongs by finding
        connected components; strands in a connected componenent of
        the graph (considering nonly STRAND nodes and bridge edges)
        form a sheet.

        Parameters: None

        Uses member data:
           node_list - list of nodes. Modifies nodes by labelling them.

        Return value:
           Returns the sheet dictionary (dictionary of
           { sheet_id : ptnode_list }) from find_connected_components.
        """
        # ACtually don't do anything except call find_connected_components()
        # which does the labeling itself (more efficient since it knows
        # as each one is added which sheet it is added to)
        return self.find_connected_components()
            

    # -------------------------------------------------------------------------
    
    def build_tableau(self, pdbid, domain, ptnode_list = None,
                      use_hk = True):
        """
        Build the tableau data member (see PTTableau in pttableau.py)
        by calling function in  pttableau.py.

        Parameters:
           pdbid - PDB identifier of the strucutre
           domain - The PTDomain object for our current domain
           ptnode_list - list of PTNodes (in sequence order, but not
                         necessarily continguous) to build the tableau for,
                         or None to use all nodes in domain.
                         Default None.
        use_hk - If True, use the HH and KK codes for respectively
                 antiparallel and parallel strands. Default True.

        Return value: None
        Uses data members (WRITE):
           tableau - created by this function
           (readonly):
           chain_dict - dict { chainid : ptnode_list } of nodes in chains
           pdb_structure - Bio.PDB parsed PDB structure
        """
        if ptnode_list == None:
            # Build list of all helix and strand PTNodes
            ptnode_list = []
            for nodelist in self.iter_chains():
                for node in nodelist: # these nodes are only those in our domain
                    if (not isinstance(node, PTNodeTerminus)): # not terminii
                        ptnode_list.append(node)

        self.tableau = pttableau.compute_tableau(ptnode_list, self.pdb_struct,
                                                 use_hk)



    def build_omega_matrix(self, pdbid, domain, ptnode_list = None):
        """
        Return the relative angles matrix by calling function in pttableau.py

        Parameters:
           pdbid - PDB identifier of the strucutre
           domain - The PTDomain object for our current domain
           ptnode_list - list of PTNodes (in sequence order, but not
                         necessarily continguous) to build the tableau for,
                         or None to use all nodes in domain.
                         Default None.
        Return value: Numeric.array Omega matrix.
        Uses data members:
           (readonly):
           chain_dict - dict { chainid : ptnode_list } of nodes in chains
           pdb_structure - Bio.PDB parsed PDB structure
        """
        if ptnode_list == None:
            # Build list of all helix and strand PTNodes
            ptnode_list = []
            for nodelist in self.iter_chains():
                for node in nodelist: # these nodes are only those in our domain
                    if (not isinstance(node, PTNodeTerminus)): # not terminii
                        ptnode_list.append(node)

        return pttableau.compute_omega_matrix(ptnode_list, self.pdb_struct)


    def build_sse_dist_matrix(self, pdbid, domain, ptnode_list = None):
        """
        Return SSE axis midpoint distance matrix by calling function
        in ptdistmatrix.py

        Parameters:
           pdbid - PDB identifier of the strucutre
           domain - The PTDomain object for our current domain
           ptnode_list - list of PTNodes (in sequence order, but not
                         necessarily continguous) to build the matrix for,
                         or None to use all nodes in domain.
                         Default None.
        Return value: Numeric.array SSE midpoint distance matrix.
        Uses data members:
           (readonly):
           chain_dict - dict { chainid : ptnode_list } of nodes in chains
           pdb_structure - Bio.PDB parsed PDB structure
        """
        if ptnode_list == None:
            # Build list of all helix and strand PTNodes
            ptnode_list = []
            for nodelist in self.iter_chains():
                for node in nodelist: # these nodes are only those in our domain
                    if (not isinstance(node, PTNodeTerminus)): # not terminii
                        ptnode_list.append(node)

        return compute_sse_midpoint_dist_matrix(ptnode_list, self.pdb_struct)



            
#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------



def make_tableaux(pdb_filename,
                pdb_struct,
                secstruct_program,
                domain_program,
                include_310_helices = False,
                include_pi_helices = False,
                use_numeric = False,
                sse_id_list = None,
                use_hk = False,
                min_sse_len = None,
                build_dist_matrix = False,
                chainid = None,
                domainid = None):
    """
    For the supplied filemame, read PDB format data from that file
    and create tableaux or SSE distance matrix for that structre.
    This function is called by get_tableaux() (below), which handles
    compressed files etc.

    Note: for multi-domains, will be multiple output tableaux, one for
    each domain.

    Paramteters:
       pdb_filename - filename of PDB file to read
       pdb_struct - Bio.PDB parsed PDB structure
       secstruct_program - secondary structure definition program
                       ('stride' or 'dssp' or 'pdb' or 'pmml') to use.
       domain_progam - domain decompositino method ('ddomain','cath', etc.)
       include_310_helices - if True, include 3_10 helices in the graph
       include_pi_helices - if True, include pi helices in the graph
       use_numeric - If True, use numeric matrix rather than tableau
       sse_id_list - list of ints representing SSE sequential id numbers
                     to include in tableau. Default None.
                     When None, all SSEs are included.
       use_hk - If True, use HH and KK codes for strands.
       min_sse_len - if not None, the minimum length of SSE to include
                    in tableau.
       build_dist_matrix - If True, build SSE midpoint distance matrix
                   instead of tableau.
       chainid - If not None, only build tableau for that chain id.
       domainid - If note None, only build tableau for that domain id.

    Return value: tuple (tableaux_list, sse_string_list)
                   where tableaux_list is
                   list of tableaux (only one in list unless domain decomp
                   is used and finds multiple domains);
                   or list of omega matrices (Numeric.array) if use_numeric
                   is True
                   or list of SSE axis midpiont distance matrices
                   (Numeric.array) if build_dist_matrix is True
                   and 
                   sse_string_list is SSE string description e.g. 'EEHHE' etc.
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
        elif secstruct_program == "pmml":
            secstruct = ptsecstruct.read_secstruct_from_pmml(pdb_filename)
        else:
            assert(False)

    if domain_program != None:
        domain_list = getdomains.get_domains(domain_program,
                                             pdbid, pdb_filename, pdb_struct)
    else:
        domain_list = [PTDomain(None, None)] # one-domain protein, no further info


    # for SSEs that cross domain boundaries, move whole SSE to one of the domains
    fixup_crossdomain_sses(secstruct, domain_list)

    tableaux_list = [] # NB may be list of PTTableau or list of Numeric.array
    sse_str_list = []
    for domain in domain_list:
        if domainid and domain.domainid != domainid:
            if verbose:
                sys.stderr.write("skipped domainid " + domainid + "\n")
            continue
        
        ptg = TableauBuild(pdb_struct, pdbid,
                           include_310_helices, include_pi_helices)
        # build tableaubuild object from secondary structure
        if secstruct_program == "pmml":
            ignore_insertion_codes = True #  PMML does not report them
        else:
            ignore_insertion_codes = False
        try:
            ptg.build_graph_from_secstruct(secstruct, domain, chainid,
                                           ignore_insertion_codes)
        except NoSSE_Exception:
            if chainid:
                sys.stderr.write('WARNING: No helices or strands found in ' +
                                 pdbid +
                                 ' chain ' + chainid +
                                 ': skipping\n')
            else:
                sys.stderr.write('WARNING: No helices or strands found in ' +
                                 pdbid +
                                 ': skipping\n')
            continue
        except ValueError:
            # PTNode.__init__ raises ValueEror if end_res_seq < start_res_seq
            # This happens on e.g. d1dkia_.ent when using PMML (only)
            # as PMML does not handle insertion codes so end up with
            # residue numbers all wrong. Nothign we can do with this much of
            # a mess, jsut give up.
            sys.stderr.write('ERROR: out of order residue numbers in ' +
                             pdbid + ', probably due to insertino codes ' +
                             'not workin in pmml: skipping this structure\n')
            continue
                

        if use_hk: # only need to know sheets if using HH and KK codes
            if secstruct_program == "pmml":
                sys.stderr.write('WARNING: PMML has no hbond or bridge information, cannot label sheets for HH and KK codes\n')
            ptg.label_sheets()

        if verbose:
            for nodelist in ptg.iter_chains():
                for node in nodelist:
                    sys.stderr.write(str(node) + '\n')

        # if list of int SSE sequential ids supplied, convert to list of
        # PTNode objects
        if sse_id_list:
            try:
                ptnode_list = [ptg.seqnum2node[sse_id] for sse_id in sse_id_list]
            except KeyError,k:
                sys.stderr.write("SSE sequential id " + str(k)
                                 + " does not exist\n")
                sys.exit(1)
        else:
            ptnode_list = None

        if not ptnode_list:
            # Build list of all helix and strand PTNodes with len >= min_sse_len
            ptnode_list = []
            for nodelist in ptg.iter_chains():
                for node in nodelist: # these nodes are only those in our domain
                    if (not isinstance(node, PTNodeTerminus)): # not terminii
                        ptnode_list.append(node)

        if min_sse_len:
            ptnode_list = [node for node in ptnode_list
                           if node.get_span() >= min_sse_len]
            
        if build_dist_matrix:
            dist_matrix = ptg.build_sse_dist_matrix(pdbid, domain, ptnode_list)
            tableaux_list.append(dist_matrix)
        elif use_numeric:
            Omega = ptg.build_omega_matrix(pdbid, domain, ptnode_list)
            tableaux_list.append(Omega)
        else:
            ptg.build_tableau(pdbid, domain, ptnode_list, use_hk)
            tableaux_list.append(ptg.tableau)

        sse_str = ""
        for node in ptnode_list:
            if isinstance(node, PTNodeStrand):
                sse_str += 'E'
            elif isinstance(node, PTNodeHelix):
                sse_str += 'H'
            else:
                raise ValueError('bad node type ' + str(node))
        sse_str_list.append(sse_str)
        
        
    return (tableaux_list, sse_str_list)



def get_tableaux(pdb_filename,
                 secstruct_program = 'dssp',
                 domain_program = 'none',
                 include_310_helices = True,
                 include_pi_helices = True,
                 sse_id_list = None,
                 min_sse_len = None,
                 use_numeric = False, 
                 use_hk = False,
                 build_dist_matrix = False):
    
    """
    Get a tableau for a single PDB or ASTRAL pdb-style file
    (compressed files e.g. pdb1qlp.ent.gz) or uncompressed
    or the ASTRAL pdb-style hierarchy
    (uncompressed files e.g. d1qlpa_.ent).

    Parameters:
       pdb_filename - filename of PDB or ASTRAL pdb-style file, as above.
       secstruct_program - secondary structure definition program
                       ('stride' or 'dssp' or 'pdb') to use.
       domain_progam - domain decompositino method ('ddomain','cath', etc.)
       include_310_helices - if True, include 3_10 helices in the graph
       include_pi_helices - if True, include pi helices in the graph
       sse_id_list - list of ints representing SSE sequential id numbers
                     to include in tableau. Default None.
                     When None, all SSEs are included.
       min_sse_len - min number of residues in SSE to be ncluded.
                      Default None (no min length).
       use_numeric - if True build Numeric.array Omega matrix not PTTableau
       use_hk - If True build tableaux with HH and KK codes for strands in
                same sheet. default False.
       build_dist_matrix - If True, build SSE midpoint distance matrices
                   instead of tableaux.


    Return value:
        tuple (pdbid, tableaux_list, sse_string_list)
         from the pdb file, only one in lists unless
         domain decomposition is used and finds multidomains in input.
         tableaux_list is list of tableaux or omega matrices
         sse_string_list is SSE string description e.g. 'EEHHE' etc.
    """
    tableaux_list = []
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
        pdbid = name
        if len(pdbid) >= 6 and pdbid[:3].upper() == "PDB":
            pdbid = pdbid[3:7].upper()
        # parse PDB file
        pdb_parser = PDBParser()
        pdb_struct = pdb_parser.get_structure(pdbid, our_pdb_filename)
        # create the Tableaux and output them
        (tableaux_list, sse_string_list) = make_tableaux(our_pdb_filename,
                    pdb_struct,
                    secstruct_program,
                    domain_program,
                    include_310_helices,
                    include_pi_helices,
                    use_numeric,
                    sse_id_list,
                    use_hk,
                    min_sse_len,
                    build_dist_matrix)
                                   
    finally:
        if used_tmp_file:
            cleanup_tmpdir(TMPDIR)
    return (pdbid, tableaux_list, sse_string_list)


