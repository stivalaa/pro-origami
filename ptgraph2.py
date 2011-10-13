#!/usr/bin/env python
###############################################################################
#
# ptgraph2.py - generate protein topology graphs from PDB atomic co-ordinates
# plus STRIDE or DSSP information.
#
# File:    ptgraph2.py
# Author:  Alex Stivala
# Created: July 2007
#
# $Id$
#
# ptgraph2 is a program to draw two-dimensional graph representations
# of protein topology. It reads information from STRIDE -h or DSSP
# output and writes a graph representation of the topological
# structure.
#
# Note, when using STRIDE,  the 'private' stride options -\$ -i are also used,
# and in fact a modified version of stride is required (see
# ptsecstruct.py) which is supplied in the stride subdirectory.
#
# It is written in Python and depends on some Python libraries:
#
# . BioPython (including Bio.PDB)
#   http://www.biopython.org
#
#   Reference for Bio.PDB is:
#   Hamelryck and Manderick 2003 "PDB parser and structure class implemented
#   in Python" Bioinformatics 19:2308-2310
#
#   which in turn depends on Numeric
#   http://sourceforge.net/projects/numpy
#
# . (for the -h, -n and -d options; see ptgraphviz.py):
#     pydot: python interface to GraphViz dot language (0.9.10)
#     http://dkbza.org/pydot.html
#
#     which in turn requires pyparsing
#     http://pyparsing.sourceforge.net/
#
#     For the -h, -n and -d options, GraphViz itself is also required (2.12)
#     http://www.research.att.com/sw/tools/graphviz
#
# Instead of using GraphViz (dot/neato), the default output 
# is now SVG for use with the Dunnart constraint-based
# diagram editor:
#
#  http://www.csse.monash.edu.au/~mwybrow/dunnart/
#
# Dunnart (0.14+SVN) (revision 1469 {2007-08-16})
# has been modified to work with this by including
# strand and helix shapes, amongst other things (by Michael Wybrow).
# (Currently not generally available).
#
# Developed on Linux 2.6.9 (x86_64) with Python 2.5.1
# and BioPython 1.43 with Numeric 24.2
#
# Example usage:
#
#  ptgraph2.py 1QLP.pdb
#
# Filenames may be either in the format above or the pdbq1lp.pdb format.
# Compressed pdb files are supported (gzip) (e.g. pdb1qlp.ent.gz).
# SCOP/ASTRAL domains in the format like d1qp3a_.ent are also supported,
# then the output filename has the same basename (eg d1qp3a_.svg).
#
###############################################################################

import warnings # so we can suppress the annoying tempnam 'security' warning
import sys,os
import getopt
from time import strftime,localtime
import re
from sets import Set # note not using builtin set, so we can use python 2.3

from Bio.PDB import *

from genseqid import GenSeqId
from ptsvgnode import *
from ptsvgcluster import *
from ptsvgconstraints import *
from ptdomain import *
from ptdistmatrix import PTDistMatrix, calc_sse_sse_dist
from ptrelpos import *
from ptmfile import *
import ptsecstruct
import pttableau
try:
    import ptgraphviz
except:
    sys.stderr.write('WARNING could not import graphviz interface: -h,-n,-d will not work\n')
    # TODO: disable these options in this case

from color import color_gradient,rgb_tuple_to_hex_str,get_color_list,get_cluster_fill_colors,DUNNART_CLUSTER_ALPHA_HEX,get_glasbey_colors_rgb
from ptutils import cleanup_tmpdir,get_int_icode,biopdbresid_to_pdbresseq
from ptposmap import *
from ptversion import get_version
from getdomains import verify_domain_disjoint,build_domain_chaindict


#-----------------------------------------------------------------------------
#
# Module constants
#
#-----------------------------------------------------------------------------

# constants used in build_helices_svg() etc.
SEQPOS_AFTER  = 1
SEQPOS_BEFORE = 2

# amount to multiply all co-ordinates by for uniform scaling
SCALE_FACTOR = 1.6


DEFAULT_SHAPE_COLOR_HEXSTR = 'f0f0d2' # beige;  Dunnart default shape color


DEFAULT_DUNNART_STRAND_SEPARATION  = 55  # space between strands in sheet
DEFAULT_DUNNART_MIN_GAP_SIZE       = 55  # minimum space to leave between things

#-----------------------------------------------------------------------------
#
# Module globals
#
#-----------------------------------------------------------------------------

DUNNART_STRAND_SEPARATION  = DEFAULT_DUNNART_STRAND_SEPARATION  # space between strands in sheet
DUNNART_MIN_GAP_SIZE       = DEFAULT_DUNNART_MIN_GAP_SIZE # minimum space to leave between things


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

# class ChainPPBuilder(_PPBuilder):
#     """
#     This class inherits from the Bio.PDB _PPBuilder base class to
#     build Polypeptide object using no distance criteria at all
#     (Polypepetide.py provides classes for Ca-Ca or C-N distances).
#     This is so we can just build sequential
#     list of residues and ignore chain breaks: each chainid will have a
#     single polypetide.
#     """
#     def __init__(self):
#         __PPBuilder.__init__(self, radius=None)

#     def _is_connected(self, prev, next):
#         return 1 # always connected.
        

class PTGraph2:
    """
    The topology graph consists of a sequence of structure (helix, strand)
    nodes with edges in and out of them in sequence from N-terminus to
    C-terminus. Note there may be multiple such sequences (one for each
    chain).

    Edges are also added for hydrogen bonds or bridge relationships
    between strands.
    """

    #
    # member functions
    #

    def __init__(self, pdb_structure, use_hbonds=False,
                 include_310_helices = False, include_pi_helices = False):
        """
        Construct empty PTGraph2. To build the graph call
        build_graph_from_secstruct().

        Parameters:
            pdb_structure - parsed PDB structure from Bio.PDB 
            use_hbonds - If True make hydrogen bond graph instead of using
                        (modified) bridge partner information to make
                        sheets from strands.
            include_310_helices - include 3_10 helices in the diagram if True
            include_pi_helices - include pi_helices in the diagram if True
            write_mfiles - write MATLAB m-files to plot strands

        """
        self.pdb_struct = pdb_structure
        self.chain_dict = None  # Each value of the chain_dict is a
                                # List of nodes in order from N to C terminus
                                # so chain_dict is { chainid : node_list }
        self.use_hbonds = use_hbonds # hbond mode
        self.distmatrix = None  # PTDistMatrix built in build_dist_matrix
        self.ptrelpos = None    # PTRelativePosition built in build_constraints
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
    
    def iter_helices(self):
        """
        This generator function iterates over all helices in this PTGraph
        object. I.e. it yields a PTNodeHelix for each helix in the 
        node lists.

        NOTE: it excludes PI and/or 310 helices if the relevant option(s)
        are set to exlclude them.
        
        Parameters: None.
        Return value: YIELDs a PTNodeHelix.
        Uses data members (readonly):
           self.chain_dict - dict of { chainid : list of nodes }
           include_310_helices, include_pi_helices - flags to omit these
        """
        for nodelist in self.iter_chains():
            for ptnode in nodelist:
                if isinstance(ptnode, PTNodeHelix):
                    if (not ( (ptnode.get_type() == "310" and
                               not self.include_310_helices) or
                              (ptnode.get_type() == "PI" and
                               not self.include_pi_helices) ) ):
                        yield ptnode


    def get_num_helices(self):
        """
        Return the total number of helices
        NB: See NOTE on iter_helices() above re pi and 310 helices.
        Parameters: None
        Uses data members (readonly):
           chain_dict - dictionary of {chainid : nodelist}
        """
        num_helices = len(list(self.iter_helices()))
        return num_helices
        
        
    def largest_helix(self):
        """
        Return the largest helix in this PTGraph, defined as the one
        that spans the most residues.

        Parameters: None.
        Return value: PTNodeHelix of helix with most residues, or None.
        Uses data members (readonly):
           self.chain_dict - dict of { chainid : list of nodes }
        """
        max_span = 0
        max_span_helix = None
        for helix in self.iter_helices():
            span = helix.get_span()
            if span > max_span:
                max_span = span
                max_span_helix = helix
        return max_span_helix
    

    def get_node_by_id(self, nodeid):
        """
        Return the node in with the supplied id

        Parameters:
            nodeid - id of the node to find
        Return value:
            PTNode in the node list that has the supplied nodeid.
        Uses data members:
            (readonly) chain_list - the list of list of nodes
        Raises exceptions:
             KeyError if node with supplied nodeid is not found in the list
        """
        # FIXME: this is just a linear search, should have a dictionary
        # by id to make it efficient (or change code to not need this at all)
        for nodelist in self.iter_chains():
            for node in nodelist:
                if node.nodeid == nodeid:
                    return node
        raise KeyError('get_node_by_id(): ' + nodeid + ' not found')
        
        

    def build_graph_from_secstruct(self, secstruct, domain):
        """
        Build the list of nodes from the the supplied PTSecStruct
        object. Add edges for hydrogen bonds between structural
        elements also.

        Parameters:
            secstruct - PTSecStruct (ptsecstruct.py) object to build from
            domain - PTDomain (ptdomain.py) object listing the segment(s)
                     that make up this domain (only one domain processed at a
                     time).

        Uses member data (write):
            chain_dict - dict of { chainid : node_list } where node_list is
                          list of nodes in order, built in this function
            secstruct - keeps a pointer to the supplied secstruct
            domainid - domain identification string

          (readonly):
            use_hbonds - If True make hydrogen bond graph instead of using
                    bridge partner information to make
                    sheets from strands.
            pdb_struct - The Bio.PDB parsed PDB struct (atomic co-ordinates)
                         for this protein.
            include_310_helices, include_pi_helices - if true, include
                         these kinds of helices so inluce them in helix
                         sequential numbering. Otherwise we still need to
                         put them in graph to ensure lines up with
                         tableaucreator etc. but keep them with a separate
                         numbering as they won't be drawn at the end.


        Raises exceptions:
           NoSSE_Exception if no helices or strands found
        
        Return value:
            None.
            
        """

        self.secstruct = secstruct
        self.domainid = domain.domainid
        
        helix_num = 1 
        strand_num = 1
        hidden_num = 1 # for pi and/or 310 helices outside visible numbering

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
                                              residue.get_id()) )] = seq_indx
            seq_indx += 1
        
        # Note that now we are only adding elements in the supplied domain,
        # so the so-called 'chains' may really be segments, i.e. subsequences
        # of chains (rest of chain may be in other domain(s)
        
        self.chain_dict = {} # dict of {chainid : node_list}

        for (start_chainid, start_resnum, end_chainid, end_resnum, helixtype) \
              in secstruct.helix_list:
            assert(start_chainid == end_chainid) #helix must be same chain
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
                    idprefix = "PIHELIX_"
                    htype = "PI"
                    if self.include_pi_helices:
                        this_helix_num = helix_num
                        helix_num += 1
                    else:
                        this_helix_num = hidden_num
                        hidden_num += 1
                elif helixtype == "G":
                    idprefix = "310HELIX_"
                    htype = "310"
                    if self.include_310_helices:
                        this_helix_num = helix_num
                        helix_num += 1
                    else:
                        this_helix_num = hidden_num
                        hidden_num += 1
                else: # shouldn't happen
                    sys.stderr.write("ERROR: bad helix type " + helixtype+"\n")
                ah_node = PTSVGNodeHelix(htype,
                                      idprefix + start_chainid+"_" +\
                                      str(this_helix_num),
                                      this_helix_num,
                                      start_resnum, end_resnum, start_chainid,
                                         self.domainid,
                                         self.residue_list,
                                         self.pdb_resid_dict)
                ah_node.build_resname_sequence()

                if not self.chain_dict.has_key(start_chainid):
                    self.chain_dict[start_chainid] = []
                self.chain_dict[start_chainid].append(ah_node)

                # we must already have handled the case of SSEs that cross
                # domain boundaries (by moving whole SSE to one of the domains)
                assert( domain.is_in_domain(end_chainid,  get_int_icode(end_resnum)[0]) )

        for (start_chainid, start_resnum, end_chainid, end_resnum) \
                in secstruct.strand_list:
            assert(start_chainid == end_chainid) # must be in same chain
            if domain.is_in_domain(start_chainid,
                                   get_int_icode(start_resnum)[0]):
                num_strands_in_domain += 1
                bs_node = PTSVGNodeStrand("STRAND_"+start_chainid +"_"+\
                                       str(strand_num),
                                       strand_num,
                                       start_resnum, end_resnum, start_chainid,
                                          self.domainid,
                                          self.residue_list,
                                          self.pdb_resid_dict)
                strand_num += 1
                bs_node.build_resname_sequence()
                
                if not self.chain_dict.has_key(start_chainid):
                    self.chain_dict[start_chainid] = []

                # we must already have handled the case of SSEs that cross
                # domain boundaries (by moving whole SSE to one of the domains)
                assert( domain.is_in_domain(end_chainid,  get_int_icode(end_resnum)[0]) )

                self.chain_dict[start_chainid].append(bs_node)


        # raise an exception if there are no SSEs at all in this domain
        if num_helices_in_domain == 0 and num_strands_in_domain == 0:
            raise NoSSE_Exception

        # build a dictionary of {chainid: (min_res_seq, max_res_seq)}
        # which we use to determine where chain has been broken by
        # domain parsing so we can label the pseudo-terminus appropriately
        chain_minmax_dict = self.build_chain_minmax_dict(self.pdb_struct)

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
                
            # (note PDB residue numbers don't necessarily start at 1,
            #  may be higher)

            # Now that we are labelling connector objects for coil regions
            # with residue names and PDB sequence numbers, we need
            # terminus nodes to have the lowest (N) and highest (C)
            # PDB residue numbers in the chain, not the lowest/highest
            # in any SSEs found.
            residue_list = self.get_residue_list(self.pdb_struct,
                                                 domain, chainid)
            # id of a residue in Bio.PDB is tuple (hetatm, resseqnum, icode)
            lowest_res_seq = residue_list[0].get_id()[1]
            highest_res_seq =residue_list[-1].get_id()[1]

            # but we still need the lowest and highest sequence numbers 
            #  in any SSEs for pseudo-terminus nodes, useed in converting
            # them to real terminus nodes or else removing to make interdomain
            # connectors later.
            sse_lowest_res_seq=get_int_icode(nodelist[0].get_start_res_seq())[0]
            sse_highest_res_seq=get_int_icode(nodelist[-1].get_end_res_seq())[0]

            # get min and max res seq num of segments of this chain in domain
            # which we use to determine where chain has been broken by
            # domain parsing so we can label the pseudo-terminus appropriately
            (domain_min_resnum, domain_max_resnum) = \
                                domain.get_minmax_res_seq_in_chain(chainid)
            # and the lowest and highest residue number in this entire chain
            min_resnum = chain_minmax_dict[chainid][0]
            max_resnum = chain_minmax_dict[chainid][1]
            if min_resnum < 0: # sometimes this happens, eg 2H85 
                min_resnum = 0
#            print 'qqq',chainid,min_resnum,max_resnum
#            print 'rrr',chainid,domain_min_resnum,domain_max_resnum
            # TODO: this will not work if a domain can have multiple
            # segments of the same chain in it with segments in between
            # in other domains. This does not currently happen with
            # DDOMAIN and CATH but could in general happen with other
            # domain decomposition methods in future.
            
            # add N terminus node at beginning
            if domain.is_single() or domain_min_resnum == min_resnum:
                if self.num_chains() > 1:
                    label = "N" + str(chainid).lower()
                else:
                    label = "N"
                pseudo = False
                lowres = lowest_res_seq
            else: # pseudo-n-terminus
                label = "n" + str(domain.domainid) + str(chainid).lower()
                pseudo = True
                lowres = sse_lowest_res_seq
            n_terminal_node = PTSVGNodeTerminus('N', pseudo,
                                                label,
                                                0, str(lowres - 1),
                                                str(lowres - 1),
                                                chainid,
                                                self.domainid,
                                                self.residue_list,
                                                self.pdb_resid_dict,
                                                True) # fake_resids
            if pseudo: # record the most N-terminal SSE in the pseudo-terminus
                n_terminal_node.set_adjnode(nodelist[0])
            nodelist.insert(0, n_terminal_node)


            # add C terminus node on end
            if domain.is_single() or domain_max_resnum == max_resnum:
                if self.num_chains() > 1:
                    label = "C" + str(chainid).lower()
                else:
                    label = "C"
                pseudo = False
                highres = highest_res_seq
            else: # pseudo-c-terminus
                label = "c" + str(domain.domainid) + str(chainid).lower()
                pseudo = True
                highres = sse_highest_res_seq
            c_terminal_node = PTSVGNodeTerminus('C', pseudo,
                                                label,
                                                1, str(highres + 1),
                                                str(highres + 1),
                                                chainid,
                                                self.domainid,
                                                self.residue_list,
                                                self.pdb_resid_dict,
                                                True) # fake_resids
            if pseudo: # record the most C-terminal SSE in the pseudo-terminus
                c_terminal_node.set_adjnode(nodelist[-1])
            nodelist.append(c_terminal_node)

        # delete chains from chain_dict that were marked earlier for deletion
        for chainid in delete_chainid_list:
            self.chain_dict.pop(chainid)
            
        # add edges for hydrogen bonds
        # uses secstruct and chainid member data
        # these are used for determining which side bridge partners are
        # on (and also for drawing a hydrogen bond graph if requested)
        self.add_hbond_edges_from_secstruct()
        
        # add edges for bridge partners
        # uses secstruct and chainid member data
        self.add_bridge_edges_from_secstruct()


    def build_chain_minmax_dict(self, pdb_struct):
        """
        build a dictionary of {chainid: (min_res_seq, max_res_seq)}
        which we use to determine where chain has been broken by
        domain parsing so we can label the pseudo-terminus appropriately

        Parameters:
            pdb_struct - The Bio.PDB parsed PDB struct (atomic co-ordinates)
                         for this protein.

        Return value:
             dictionary of {chainid : (min_res_seq, max_res_seq)}
             mapping chainidentifier to minimum and maximum residue
             sequence numbers in that chain, as determined by
             Bio.PDB PDB parser (pdb_struct parameter)

        Note: uses no data members (need not be a member function)
        """
        chain_minmax_dict = {}
        pdb_model = pdb_struct[0] # TODO always using model 0 for now
        for chain in pdb_model:
            # id of a residue in Bio.PDB is tuple (hetatm, resseqnum, icode)
            res_seqnum_list = [ res.get_id()[1] for res in chain.get_list()
                                if res.get_id()[0] == ' ' ]
            if len(res_seqnum_list) > 0:
                min_res_seq = min(res_seqnum_list)
                max_res_seq = max(res_seqnum_list)
                chain_minmax_dict[ptsecstruct.pdb_chainid_to_stride_chainid( \
                                chain.get_id())] = (min_res_seq, max_res_seq)
        return chain_minmax_dict
    

    def num_chains(self):
        """
        Return the number of chains read from the PDB in this object.

        Parameters: Nonde
        Return value: Number of chains represented in this graph
        Uses data members (readonly):
           chain_dict - dictionary of { chainid : nodelist }
        """
        return len(self.chain_dict)

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
                if ( chainid1 == ptnode.get_chainid() and
                     ptnode.is_in_interval(resnum1) ):
                    try:
                        dest_node = self.find_node_containing_seqnum(resnum2,
                                                                     chainid2)
                    except KeyError:
                        # it seems STRIDE sometimes gives h bond to
                        # residue that does not exist - maybe because the
                        # chainid is wrong (puts same chainid for both
                        # sides even when bond is between chains? - could
                        # recover from this if so (TODO). e.g. 1EAI, 1ABI
                        sys.stderr.write('WARNING: external program' 
                                         ' reported H bond to nonexistent '
                                         'residue ' + resnum2 +
                                         ' (chain ' + chainid2 + ')\n')
                        continue
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
                if ( chainid1 == ptnode.get_chainid() and
                     ptnode.is_in_interval(resnum1) ):
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


    def build_dist_matrix(self, domain):
        """
        Build the PTDistMatrix member ptdistmat containing the
        residue and SSE distance maps. This is built using information
        from the Bio.PDB Residue objects contained in the
        member pdb_struct for residues in the supplied domain which
        we are working with, and (for SSEs) the secondary structures
        defined by node lists under the chain_dict member built by
        build_graph_from_sec_struct().

        Parameters:
            domain - PTDomain (ptdomain.py) object listing the segment(s)
                     that make up this domain (only one domain processed at a
                     time).


        Uses member data (write):
           distmatrix - the PTDistanceMatrix class containing residue and
                       SSE distance maps.
          (readonly):
            chain_dict - dict of { chainid : node_list } where node_list is
                          list of nodes in order, built in this function
            secstruct - keeps a pointer to the supplied secstruct
            pdb_struct - The Bio.PDB parsed PDB struct (atomic co-ordinates)
                         for this protein.
            sheet_dict - dict of {sheet_id : ptnode_list} represneting sheets

        Return value:
            None.
   
        """
        residue_list = self.get_residue_list(self.pdb_struct, domain)
        # Also build list of all PTNodes
        ptnode_list = []
        for nodelist in self.iter_chains():
            for node in nodelist:
                if (not isinstance(node, PTNodeTerminus)): # not terminii
                    ptnode_list.append(node)
            

        self.distmatrix = PTDistMatrix(residue_list, ptnode_list,
                                       self.sheet_dict,
                                       self.pdb_struct)


    def build_tableau(self, pdbid, domain, use_tableaucreator):
        """
        Build the tableau data member (see PTTableau in pttableau.py)
        by calling function in  pttableau.py.

        Parameters:
           pdbid - PDB identifier of the strucutre
           domain - The PTDomain object for our current domain
           use_tableaucreator - if True, use external TableauCreator program
        Return value: None
        Uses data members (WRITE):
           tableau - created by this function
           (readonly):
           chain_dict - dict { chainid : ptnode_list } of nodes in chains
           pdb_structure - Bio.PDB parsed PDB structure
        """
        # Build list of all helix and strand PTNodes
        ptnode_list = []
        for nodelist in self.iter_chains():
            for node in nodelist: # these nodes are only those in our domain
                if (not isinstance(node, PTNodeTerminus)): # not terminii
                    ptnode_list.append(node)

        if use_tableaucreator:
            self.tableau = pttableau.get_tableau_from_pdbstruct(
                pdbid, domain, self.pdb_struct, ptnode_list)
        else:
            self.tableau = pttableau.compute_tableau(ptnode_list,
                                                     self.pdb_struct)


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
            

          

    def label_bridge_sides(self):
        """
        Label bridge edges with '+' or '-' indicating relative side of
        the strand the bridge partners are on. This method just calls
        label_strand_bridge_sides() (PTNodeStrand method, see documetation
        there for details) for each strand to do this.
        
        Parameters: None
        
        Uses member data:
             chain_dict - dict by chainid of
                         list of nodes in the graph. Not iself modified,
                         but edges in the PTNodeStrands (i.e. bridge_list)
                         in those is modified by label_strand_bridge_sides()

             pdb_struct - The Bio.PDB parsed PDB struct (atomic co-ordinates)
                          for this protein.

        Return value: None
        """
        for node in self.iter_strands():
            node.label_strand_bridge_sides(self.pdb_struct)
                


    def build_constraints(self):
        """
        Build constraints for all sheets. The constraints consist of strands
        in a sheet being in a cluster, and strand ordering constraints
        for each sheet computed by build_sheet_constraints()
        
        The idea is that these constraints can then be used as input
        to a constraint-based graph layout system (Dunnart - see
        write_dunnart_svg()).


        Return value: None; member data sheet_strandlists_dict is built.

        Uses data members:
            sheet_dict - the sheet dictionary (dictionary of
                         { sheet_id : ptnode_list }) from
                         find_connected_components.

            sheet_strandlists_dict (WRITE) -
                dictionary of { sheet_id : list of list of nodes }
                where the list of list of nodes is
                described in build_sheet_constraints()

            ptrelpos (WRITE) -
                     PTRelativePosition instance created to find relative
                       positions of elements later, using dist matrix
                        and sheet strand lists.

            pdb_stucrt, distmatrix, tableau - passed to PTRelativePosition which
                                    keeps a pointer to them itself.
        """
        
        # label bridge edges + or - for relative sides of strand, as per
        # Westhead et al 1999.
        self.label_bridge_sides()

        # now that all bridge sides are labelled, process each sheet 
        self.sheet_strandlists_dict = {}
        for (sheet_id, ptnode_list) in self.sheet_dict.iteritems():
            self.sheet_strandlists_dict[sheet_id] =  \
                                      self.build_sheet_constraints(ptnode_list)

        self.ptrelpos = PTRelativePosition(self.pdb_struct,
                                           self.distmatrix,
                                           self.sheet_strandlists_dict,
                                           self.tableau,
                                           self.chain_dict,
                                           self.sheet_dict)
        

            
                                  
    def build_sheet_constraints(self, ptnode_list):
        """
        Build constraints, for a single sheet,
        on the layout of the strand nodes based on their
        sheet membership (determined by label_sheets(), which must
        be called before this subroutine, and the return value from it
        passed to this one).

        These constraints specify that strands in a sheet are positioned
        in order along a (horizontal) line. Bifurcations
        in the sheet also need to be taken into account.
        Strands are ordered according to their neighbour (bridge partner)
        relationships, and a bifurcation results in having two (or more,
        conceivably) strands both neighbours on the same side of one strand.
        These can then have further neighbours of their own, indpenedently.

        In addition to the data structure returned by this subroutine
        defining constraints, the strand nodes themselves are annotated
        with further information, such as the orientation (up or down)
        of the arrown to be drawn (based on the parallel/antiparallel
        relationships between adjacent strands as previously labelled
        by label_bridge_sides()), and the vertical position relative to the
        neighbouring strand determined by the positions (residues) on the
        strands where the H-bonds forming the beta-bridges.
        
        Parameters:
            ptnode_list - list of PTNodeStrand elements for the sheet.

        Return value:
                list of list of PTNodeStrand.
                Each element in that list is itself a list of
                strands. The outermost list for horizontal (left to right)
                position in the layout. The list of nodes at a given positino
                in that outermost list is a list of nodes at the same
                horizontal position, and each is labelled with a relative
                vertical position for the layout.
                
                This is similar to what is done in TOPS
                as per Westhead et al 1999 (see Fig 5B, p. 902).

        Uses data members (readonly):
                sheet_backedges_dict - dict of { sheet_id : (node1, node2) }
                                       of back edges from DFS in connected
                                       component (sheet)
                                       
        Note: also set properties in PTNodeStrand objects (reversed, align_pos)
        """

        # note number of undirected edges is 1/2 number of bridge edges
        # in our connected component as we have one for each direction.
        num_dir_edges= sum([ node.num_neighbours() for node in ptnode_list])
        assert(num_dir_edges % 2 == 0)
        num_edges = num_dir_edges / 2
        # because a sheet is (by definition the way we have found it)
        # is a connected component, if |E| = |V| - 1 it is acyclic
        # otherwise it must be cyclic, and in fact:
        num_cycles = num_edges - len(ptnode_list) + 1
        assert(num_cycles >= 0)

        #------ debug output
        if verbose:
            sys.stderr.write('num nodes = ' + str(len(ptnode_list)) + '\n')
            sys.stderr.write('num edges = ' + str(num_edges) + '\n')
            sys.stderr.write('hence num cycles = ' +str(num_cycles) + '\n')
        #------ end

        if len(ptnode_list) < 2:
            sys.stderr.write(
              'WARNING: sheet has fewer than 2 strands (' +
              str(ptnode_list[0]) + ') \n')
            return [ptnode_list]

        # now that we have the bridge edges labelled with relative sides,
        # we can (for each sheet) use a DFS, starting at a leaf node
        # (one with only one neighbour) to get the positioning of strands
        # as per Westhead et al 1999 (see Fig 5, p. 902).

        open_node1 = open_node2 = None
        if num_cycles > 0:
            # for a connected component, the number of back edges must be
            # equal to the number of edges - number of nodes + 1
            # (= number of (basic) cycles)
            sheet_id = ptnode_list[0].get_sheet_id()
            assert(num_cycles == len(self.sheet_backedges_dict[sheet_id]))
            # There is at least one cycle (beta barrel) so we won't be able to
            # find a strand with only one neighbour.
            # We will 'open' the barrel into a sheet.
            # At the moment we will 'break' a bridge arbitrarily
            # (by removing each 'back edge' found in the DFS to break
            # each fundamental (aka basic or basis) cycle).
            # This is consistent and and guaranteed to break the cycle but
            # it is rather arbitrary from a physical/chemical/biological
            # point of view.
            # TODO: maybe should have some more physically relevant
            # criteria for choosing which
            # bridge to break e.g. Hutchinson and Thornton 1990 (HERA)
            # open at position which minimizes number of H-bonds broken.
            sys.stdout.write("detected " + str(num_cycles) +
                             " beta barrel(s)\n")
            for (open_node1, open_node2) in self.sheet_backedges_dict[sheet_id]:
                sys.stdout.write("opened barrel by removing bridge from " +
                                 open_node1.nodeid + " to " +
                                 open_node2.nodeid +"\n")
                open_node1.set_barrel_edge(True)
                open_node2.set_barrel_edge(True)
                open_node1.remove_bridge(open_node2)


        # get node to start at (a node with only one neighbour)
        # If a barrel was broken, start at the node where the brdige was
        # removed, if that node has only one neighbour (otherwise use
        # the frist node we find that does only have on neighbour)
        if open_node1 != None and open_node1.num_neighbours() == 1:
            start_node = open_node1
        elif open_node2 != None and open_node2.num_neighbours() == 1:
            start_node = open_node2
        else:
            start_node = None
            for node in ptnode_list:
                if node.num_neighbours() == 1:
                    start_node = node
                    break
        assert(start_node.num_neighbours() == 1)

        dfs_list = []
        dfs_strands_from(start_node, {}, dfs_list, None) # in ptnode.py

        #------ debug output
        if verbose:
            for (node, from_node) in list(dfs_list):
                if from_node == None:
                    fm_nodeid = '<none>'
                else:
                    fm_nodeid = from_node.nodeid
                    sys.stderr.write(node.nodeid+" from "+ fm_nodeid+", ")
            sys.stderr.write('\n')
        #------ end

        nodelist = list(dfs_list) # nodelist is in DFS order
        horiz_order_list = []
        for i in range(len(nodelist)): # may be longer than needed
            horiz_order_list.append([])
        node_index_dict = {} # dictionary of { nodeid : horiz_order_list_index }
        firstnode = nodelist[0][0]
        firstnode.set_reversed(False)
        firstnode.set_align_pos(0)
        horiz_order_list[0] = [firstnode]
        node_index_dict[firstnode.nodeid] = 0
        prev_index = 0
        prev_side = None
        prev_node = None
        for node_index in range(1, len(nodelist)):
            (node, from_node) = nodelist[node_index]

            # put node in correct position in outer list according
            # to constraints on side of node it was reached from
            # relative to other node adjacent to that node
            found_constraint = False
            fromnode_index = node_index_dict[from_node.nodeid]
            cur_side = from_node.get_side_of_neighbouring_strand(node)
            if cur_side != '.': # side contraints are '+' or '-'; '.' means none
                prev_side = None
                for (prev_node,bdir_unused,side) in from_node.get_bridge_list():
                    if prev_node != node:
                        prev_index = None
                        for k in range(fromnode_index + 2): #upto AFTER fromnode
                            if prev_node in horiz_order_list[k]:
                                prev_index = k
                                prev_side = side
                                if prev_side != '.': # found a constraint
                                    found_constraint = True
                                    break
                    if found_constraint:
                        break
            if found_constraint:
                if cur_side == prev_side:
                    horiz_order_list_index = prev_index
                else:
                    if prev_index < fromnode_index:
                        horiz_order_list_index = fromnode_index + 1
                    else:
                        horiz_order_list_index = fromnode_index - 1
            else:
                horiz_order_list_index = fromnode_index + 1
            node_index_dict[node.nodeid] = horiz_order_list_index
            vert_list = horiz_order_list[horiz_order_list_index]
            vert_list.append(node)

            # set reversed flag depending on parallel/antiparallel relationship
            if node.is_parallel(from_node):
                node.set_reversed(from_node.get_reversed())
            else:
                node.set_reversed(not from_node.get_reversed()) # antiparallel

            # set relative vertical positions of strands based on H bonds
            compute_align_positions(node, from_node)
            

        # since outer list may have been longer than needed, remove empty lists
        for i in range(len(horiz_order_list)-1, -1, -1):
            if horiz_order_list[i] == []:
                horiz_order_list.pop(i)
                
        #------ debug output
        if verbose:
            for vert_list in horiz_order_list:
                for node in vert_list:
                    dirstr = "up"
                    if node.get_reversed():
                        dirstr = "down"
                    sys.stderr.write(node.nodeid + '(' + dirstr + ' ' +
                                     str(node.get_align_pos()) + ') ')
                                     
                sys.stderr.write('\n')
        #------ end

        return horiz_order_list



    def sheet_size(self, sheet_id):
        """
        Calculate the 'size' of a sheet, which is defined as the number
        of residues in the sheet, i.e. the sum of the number of residues
        in each strand of the sheet.
        
        Parameters:
            sheet_id - id of the sheet to find size of

        Uses member data:
            sheet_dict - dict of {sheet_id : ptnode_list} represneting sheets

        Return value:
             Number of residues in the sheet (sum of number of residues in
             each strand of the sheet)
        """
        return sum([strand.get_span() for strand in self.sheet_dict[sheet_id]])


    def largest_sheet(self):
        """
        Return the id of the 'largest' sheet int the sheet_dict, where
        size of a sheet is as defined by the sheet_size() function.

        Parameters: None.
        
        Uses member data:
            sheet_dict - dict of {sheet_id : ptnode_list} represneting sheets

        Return value:
            sheet id of the largest sheet, or None
        """
        max_size = 0
        max_size_sheet_id = None
        for sheet_id in self.sheet_dict.iterkeys():
            size = self.sheet_size(sheet_id)
            if size > max_size:
                max_size = size
                max_size_sheet_id = sheet_id
        return max_size_sheet_id


    def get_orientation_sse(self):
        """
        Return the SSE (strand or helix) to be used for relative orientation
        of this whole PTgraph2 (domain/protein). This is the longest
        strand of the largest sheet (NB not longest strand overall), or
        longest helix if there is no sheet.

        Parameters:
            None

        Return value:
            PTNode which is either a PTNodeStrand or PTNodeHelix, as described
            above.

        Uses data members (Readonly):
            sheet_dict, nodelist, ptrelpos, sheet_strandlists_dict
        """
        sheet_id = self.largest_sheet()
        if sheet_id != None:
            (sse, length) = self.ptrelpos.get_longest_strand(
                                          self.sheet_strandlists_dict[sheet_id])
        else:
            sse = self.largest_helix()
        return sse
    
        
    def write_sheet_mfiles(self, pdbid, domainid):
        """
        Write a MATLAB m-file for each sheet. The m-file contains
        commands to plot in 3D the carbon-alpha trace of each strand
        and the axis fitted to it.
        Note we can't use the same convention as our .ps and .svg filenames
        since MATLAB can't have '-' characters  in m-file filenames,
        or digits as the first character (since m-files are used as commands).
        So we will have filenames that always have the domain identifier
        (even if single domain), and will add the sheet id on the end,
        and an M on the front
        e.g. 'M1QLP1A.m' (1QLP domain 1 sheet 2).
        Then you can run 'matlab -r M1QLP1A' to plot it
        (or just enter M1QLP1A at the MATLAB prompt).

        Parameters:
            pdbid - PDB identfier
            domainid - domain identifier
            
        Return value: None

        Uses member data: (readonly)
            nodelist - list of nodes; NOTE however that as the fit_axis()
                       method is called this will compute axes and store
                       in PTNodeStrand nodes - when not writing m-files
                       only axes are computed as needed, this forces them
                       all to be computed.
            secstruct - PTSecStruct representing the secondary structure
            sheet_dict - dict of {sheet_id : ptnode_list} represneting sheets

        WARNING: overwrites the output files
        """
        for (sheet_id, ptnode_list) in self.sheet_dict.iteritems():
            filename = 'M' + pdbid + domainid + sheet_id + '.m'
            sys.stdout.write('writing file ' + filename + '\n')
            fh = open(filename, 'w')
            mfile_write_prelude(fh)
            for node in ptnode_list:
                node.fit_axis(self.pdb_struct, fh) # writes to fh itself
            mfile_write_conclusion(fh)
            fh.close()

    def write_helix_mfiles(self, pdbid, domainid):
        """
        Write a (single) MATLAB m-file for all helices. The m-file contains
        commands to plot in 3D the carbon-alpha trace of each helix
        and the axis fitted to it.
        Note we can't use the same convention as our .ps and .svg filenames
        since MATLAB can't have '-' characters  in m-file filenames,
        or digits as the first character (since m-files are used as commands).
        So we will have filenames that always have the domain identifier
        (even if single domain), and will add the sheet id on the end,
        and an MH on the front
        e.g. 'MH1QLP1.m' (1QLP domain 1).
        Then you can run 'matlab -r MH1QLP1' to plot it
        (or just enter MH1QLP1 at the MATLAB prompt).

        Parameters:
            pdbid - PDB identfier
            domainid - domain identifier
            
        Return value: None

        Uses member data: (readonly)
            nodelist - list of nodes; NOTE however that as the fit_axis()
                       method is called this will compute axes and store
                       in PTNodeHelix nodes - when not writing m-files
                       only axes are computed as needed, this forces them
                       all to be computed.
            secstruct - PTSecStruct representing the secondary structure

        WARNING: overwrites the output files
        """
        filename = 'MH' + pdbid + domainid + '.m'
        sys.stdout.write('writing file ' + filename + '\n')
        fh = open(filename, 'w')
        mfile_write_prelude(fh)
        for node in self.iter_helices():
            node.fit_axis(self.pdb_struct, fh) # writes to fh itself
        mfile_write_conclusion(fh)
        fh.close()

                   
    def build_dunnart_svg(self,
                          sse_label_scheme = 'separate',
                          use_connector_arrowheads=False,
                          heuristic_helix_placement=False,
                          sheet_shading_colors = None,
                          enable_sheet_gap_rule = False,
                          use_helix_clustering = False,
                          helix_cluster_shading_color = None,
                          connector_color_scheme = 'all',
                          color_scheme = 'none',
                          helix_proximity_shading_colors = None,
                          initial_xmlid = 1,
                          main_sideways = False,
                          main_reversed = False,
                          interdomain_connectors = False,
                          label_residue_numbers = False):
        """
        Build SVG for input to the Dunnart interactive constraint-based
        layout diagramming program.
        
        http://www.csse.monash.edu.au/~mwybrow/dunnart/
        
        Dunnart (0.14+SVN) has been modified to work with this by including
        strand and helix shapes, amongst other things (by Michael Wybrow).
        (Currently not generally available).

        The greedy algorithm for laying out the cartoon is basically:

        1. Place the largest element somewhere near the middle of canvas
        2. Use the distance map to find closest element to any already
           placed element
        3. position that closest element relative to the chosen already
           placed one, using distance/position maps and tableaux to determine
           relative position and orientation
        4. repeat from 2 until all elements placed

        Sheets are single elements in themsleves, with the relative positioning
        and layout of their strands performed by build_sheet_svg() based
        on the constraints already built in build_sheet_constraints().
        
        Parameters:
           sse_label_scheme - 'none', 'sequential' or 'separate'.
                   Note that this does not affect internal number, it is
                   just for the labels on Dunnart shapes.
           use_connector_arrowheads - If true put arrowheads on connectors
                         indicating sequence direction from N- to C- terminus.
                         Default False.
           heuristic_helix_placement - instead of using the greedy algorithm
                         with distance matrix information for helices,
                         use the old heuristic algorithm to place them
                         aligned neatly with strands nearby in sequence.
                         Still use greedy distance matrix algorithm for sheets.
           sheet_shading_colors - None (use default shade for all) or
                                  'auto' (use color gradient to shade each
                                  differently) or list of colors.
           enable_sheet_gap_rule - If True and using heuristc helix
                         placement, don't put 'too long' helices between
                         neighbouring sheets.
           use_helix_clustering - If True and using heuristic helix placement,
                        cluster sequential helices and place in cluster
                        with tableau and distance matrix rather than
                        aligning them all on strand axis.
           helix_cluster_shading_color - color to shade helix clusters
           connector_color_scheme - 'all','chain','domain','crossing' (see main)
           color_scheme: 'none', 'simple', 'gradient', 'sheet', 'fold'
                          (specc. in main)
           use_helix_proximity_color - if True and using helix clustering
                        shade nearby helix clusters the same color.

           helix_proximity_shading_colors - If not None & using helix clustering,
                                          shade nearby helix clusters the same
                                          color: 'auto' (use color gradient
                                          to shade each differently),
                                          or list of colors.
                        
           initial_xmlid - Initial XML sequential identifier. Default 1.
           main_sideways - if True, the 'main' part of the domain (largest sheet
                        or longest helix) is drawn sideways instead of vertical
           main_reversed - If True, the main part (as above) is drawn reversed
                       (down not up, or right not left when sideways)
           interdomain_connectors - If True, do NOT make pseduo-terminus nodes
                       at domain boundaries. Instead the domain boundary
                       SSEs are left ot hve connectors to other domain
                       added later. Default False.
           label_residue_numbers - If True put start and end residue numbers
                                   of SSE on head and tail of shape as labels

        Return value: None


        Uses member data: (readonly)
            secstruct - PTSecStruct representing the secondary structure
            sheet_strandlists_dict - the dictionary (by sheet id)
                               of list of lists built by
                               build_constraints()
            ptrelpos -
                     PTRelativePosition instance to find relative positions

            include_p_helices, include_310_helices
            pdb_struct - need to get residue name/id informatino for connectors

            (read/write):
              sheet_cluster_list - list of SVGCluster objects for sheets
              helix_cluster_list - list of SVGCluster objects for helix clusters
              svg_constraint_list - list of PTSVGConstraint derived objects
              svg_connector_list - list of PTSVGConnector objects
              Note also uses  member data in each node.


        NOTE: The output file is overwritten if it exists.

        Precondition: node_list is sorted (by start res seq ascending);
                      this is done by build_graph_from_secstruct()
                      before calling.


        """
        # everything in the SVG will have a unique id, sequentially
        self.xmlid_generator = GenSeqId(initial_xmlid)

        # if using sequential numbering, we'll build this dictionary mapping
        # nodeid to sequential number (note NOT restarting for each chain)
        # since we will not be
        # iterating through node_lists when writing the SVG, so will need
        # to look up the ordinal position (sequence number) for a node.
        # this is a dictionary of { nodeid : seqnum }
        seqnum_dict = {}
        for (seqnum, nodeid) in \
            enumerate([node.nodeid for node in self.iter_nodes() if \
                       not ( (isinstance(node, PTNodeTerminus)) or
                              (isinstance(node, PTNodeHelix) and
                               ( (node.get_type() == "310" and
                                  not self.include_310_helices) or
                                 (node.get_type() == "PI" and
                                  not self.include_pi_helices) ) ) ) ]):
                seqnum_dict[nodeid] = seqnum + 1 # start at 1 not 0

        init_xpos = 300  # initial x position
        init_ypos = 200  # initial y posision


        xpos = init_xpos
        ypos = init_ypos
        
        self.sheet_cluster_list = [] # list of PTSVGCluster objects for sheets
        self.helix_cluster_list = [] # list of PTSVGCluster for helix clusters
        self.helix_cluster_dict = {} # dict { clusterid : PTSVGCluster }
        self.svg_constraint_list = [] # list of PTSVGConstraint objects
        self.svg_connector_list = []  # list of PTSVGConnector objects
                               
        # store y pos of sheet top and bottom
        # and x pos of sheet left and right in a dictionary
        # { sheet_id : (top_ypos, bottom_ypos, left_xpos, right_xpos) }
        # e.g. { 'A' : (100, 120, 30, 200) }
        # (note y before x for historical reasons, started with only ypos
        # then extended it)
        # so we can position things
        # relative to them later
        # TODO: should probably have a class to represent sheets with these
        # values as data members rather than this dict
        sheet_pos_dict = {}

        # stores above,below,left,right neighbour sheet of each sheet
        # TODO:  maybe we should use this to avoid collisions between
        #        helices and helix clusters and sheets (not just between
        #        sheets as we do currently)
        #        instead of depending on Dunnart
        #        non-overlap constraints to fix it up when this happens?
        #        (Since greedy algorithm can result in elements being placed
        #        on top of each other).
        sheet_posmap = PTPosMap()
        
        num_sheets = len(self.sheet_dict)
        num_helices = self.get_num_helices()

        # We'll make two sets: one for all the elements (sheets, helices)
        # not yet positioned, and one for those that have been positioned.
        # members of these sets are either PTNodeHelix for a helix or
        # sheet id ('A' etc.) for a sheet.
        # We also need a third set of elements to exclude from placement,
        # namely 310 and/or pi helices if we don't want to draw them.
        # (These will not be in unpositioned element (handled by iter_helices())
        # but we need a set to explicitly tell find_nearest_to_any_in_set()
        # to not return them).
        positioned_elements = Set()
        exclude_elements = Set()
        if heuristic_helix_placement:
            # only position sheets, helices placed later
            unpositioned_elements = Set(list(self.sheet_dict.keys()))
        else:
            unpositioned_elements = Set(list(self.iter_helices()) +
                                        list(self.sheet_dict.keys()))
            for node in self.iter_nodes():
                if ( isinstance(node, PTNodeHelix) and
                     ((node.get_type() == "310" and
                       not self.include_310_helices) or
                      (node.get_type() == "PI" and
                       not self.include_pi_helices) ) ):
                    exclude_elements.add(node.nodeid) #note: nodeid not node

        total_elements = len(unpositioned_elements)
        
        # get the largest element, to place first. It will be a sheet,
        # and if no sheets (alpha-only domain), the largest helix
        largest_sheet_id = self.largest_sheet()
        if largest_sheet_id != None:
            if main_sideways:
                if verbose:
                    sys.stderr.write('largest sheet ' + largest_sheet_id +
                                     ' is sideways\n')
                self.ptrelpos.set_all_sheet_strands_sideways(largest_sheet_id)
            if main_reversed:
                if verbose:
                    sys.stderr.write('largest sheet ' + largest_sheet_id +
                                     ' is reversed\n')
                self.ptrelpos.flip_all_strands_in_sheet(largest_sheet_id)
            self.build_sheet_svg(sheet_pos_dict,
                                 sse_label_scheme, seqnum_dict,
                                 largest_sheet_id,
                                 xpos, ypos, label_residue_numbers)
            unpositioned_elements.remove(largest_sheet_id)
            positioned_elements.add(largest_sheet_id)
        elif not heuristic_helix_placement:
            # no sheets, place largest helix first
            helix = self.largest_helix()
            helix.set_sideways(main_sideways)
            helix.set_reversed(main_reversed)
            self.set_helix_svginfo(helix, xpos, ypos,
                                   sse_label_scheme, seqnum_dict,
                                   label_residue_numbers)
            unpositioned_elements.remove(helix)
            positioned_elements.add(helix)

        while len(unpositioned_elements) > 0:
            assert(len(unpositioned_elements) +
                   len(positioned_elements) == total_elements)
            assert(len(unpositioned_elements & positioned_elements) == 0)
            
            # find nearest element to any already positioned element
            (positioned_element, cur_element) = \
                    self.find_nearest_to_any_in_set(positioned_elements,
                                                    exclude_elements,
                                                    heuristic_helix_placement)
            if verbose:
                sys.stderr.write('positioning ' + str(cur_element) +
                                 ' relative to ' +
                                 str(positioned_element) + '\n')

            # position the element near its closest already positioned one
            # according to position map for relative position and tableau
            # for orientation
            (relpos, pos_strand, cur_strand) = \
                     self.ptrelpos.get_relative_position(positioned_element,
                                                         cur_element)
            (xpos, ypos) = self.relpos_to_pos(positioned_element,
                                              cur_element,
                                              relpos,
                                              pos_strand, cur_strand,
                                              sheet_pos_dict)

            if isinstance(cur_element, PTNodeHelix):
                self.set_helix_svginfo(cur_element,
                                       xpos, ypos,
                                       sse_label_scheme,
                                       seqnum_dict,
                                       label_residue_numbers)
            else:
                if not isinstance(positioned_element, PTNodeHelix):
                    # positioned element is a sheet. Avoid collisions with
                    # already placed sheets by checking in posmap and
                    # positioning relative to the thing we would have collided
                    # with instead.
                    if sheet_posmap.has_key(positioned_element):
                        sheet_neighbours = sheet_posmap[positioned_element]
                        neighbour = sheet_neighbours.get_neighbour(relpos)
                        if neighbour != None:
                            if verbose:
                                sys.stderr.write("  sheet positioning: "
                                                 "can't place " +
                                                 cur_element + ' ' +
                                                 ptrelpos_to_str(relpos) + ' ' +
                                                 'sheet ' + positioned_element +
                                                 ' (occupied by ' +
                                                 neighbour +
                                                 ')\n')
                            positioned_element = neighbour
                            if verbose:
                                sys.stderr.write('  positioning ' +
                                                 str(cur_element) +
                                                 ' relative to ' +
                                                 str(positioned_element) +
                                                 ' instead\n')

                            (relpos, pos_strand, cur_strand) = \
                                     self.ptrelpos.get_relative_position(positioned_element,
                                                                         cur_element)
                            (xpos, ypos) = self.relpos_to_pos(positioned_element,
                                                              cur_element,
                                                              relpos,
                                                              pos_strand, cur_strand,
                                                              sheet_pos_dict)
                            sheet_neighbours = sheet_posmap[positioned_element]
                            neighbour = sheet_neighbours.get_neighbour(relpos)
                            if neighbour != None:
                                # note there is cut&paste code
                                # here from further up (before if isinstance...)
                                # And sometimes (e.g. 2QP2-1, the case this is
                                # intended originally to fix), we get a collision
                                # on the new positioned_element as well.
                                # So we could then use an arbtirary position
                                # (but what if none?) Don't want an lot of
                                # cut&paste trying.. maybe should loop from
                                # closest to furthest SSEs from the test element
                                # and place in relpos to one where the relpos
                                # neighbour slot is vacant.
                                # FIXME: arbitrary sheet positioning here
                                if (sheet_neighbours.west == None):
                                    relpos = RELPOS_LEFT
                                elif (sheet_neighbours.east == None):
                                    relpos = RELPOS_RIGHT
                                elif (sheet_neighbours.north == None):
                                    relpos = RELPOS_ABOVE
                                elif (sheet_neighbours.south == None):
                                    relpos = RELPOS_BELOW
                                else:
                                    sys.stderr.write('WARNING: nowhere to place ' +
                                                     cur_element +
                                                     ' relative to sheet '
                                                     + positioned_element +
                                                     '\n')
                                if verbose:
                                    sys.stderr.write('    (collision): ' + 
                                                     'positioning sheet ' +
                                                     cur_element + ' ' +
                                                     ptrelpos_to_str(relpos) +
                                                     ' sheet ' +
                                                     positioned_element +
                                                     '\n')
                                (xpos, ypos) = self.relpos_to_pos(positioned_element,
                                                                  cur_element,
                                                                  relpos,
                                                                  pos_strand, cur_strand,
                                                                  sheet_pos_dict)

                sheet_posmap.add_neighbour_obj(positioned_element, cur_element,
                                               relpos)
                self.build_sheet_svg(sheet_pos_dict,
                                     sse_label_scheme, seqnum_dict,
                                     cur_element,
                                     xpos, ypos, label_residue_numbers)

                # for strands that are both vert or both horiz and used as
                # closest elements in sheets, constrain them to be aligned
                # on their indguides, if they are above/below when vert
                # or left/right when horiz.
                # NOTE for -i option (helix placement using
                # dist matrix) pos_strand may actually be a helix,
                # so won't have an indguide). TODO: align on helix?
                if ( isinstance(pos_strand, PTNodeStrand) and
                     cur_strand.get_sideways() == pos_strand.get_sideways() and
                     ((not cur_strand.get_sideways() and
                       (relpos == RELPOS_ABOVE or relpos == RELPOS_BELOW)) or
                      (cur_strand.get_sideways() and
                       (relpos == RELPOS_LEFT or relpos == RELPOS_RIGHT))) ):
                    assert(pos_strand.indguide != None)
                    if cur_strand.get_sideways():
                        align_type = DUNNART_ALIGN_MIDDLE
                    else:
                        align_type = DUNNART_ALIGN_CENTER
                    self.svg_constraint_list.append(
                        PTSVGAlignmentConstraint(pos_strand.indguide,
                                                 cur_strand,
                                                 align_type))
            unpositioned_elements.remove(cur_element)
            positioned_elements.add(cur_element)
        # END of iteration over unpositioned_elements set
        assert(len(positioned_elements) == total_elements)

        if heuristic_helix_placement:
            # first position any helices that meet the special case of
            # being between two strands on same axis (i.e. in vertlist)
            # in sheet
            self.build_interstrand_helices_svg(sheet_pos_dict,
                                               sse_label_scheme,
                                               seqnum_dict,
                                               label_residue_numbers)

            # now position all other helices on axis of nearby strand,
            # on elsewhere according to heuristics
            self.build_helices_svg(sheet_pos_dict,
                                   sse_label_scheme, seqnum_dict,
                                   sheet_posmap, enable_sheet_gap_rule,
                                   use_helix_clustering, 
                                   helix_cluster_shading_color,
                                   color_scheme,
                                   interdomain_connectors,
                                   label_residue_numbers)

            if use_helix_clustering:
                # if using helix clustering, color nearby helix clusters 
                # the same shading color. Note that this may involve
                # making some helices a one-helix 'cluster'.
                if helix_proximity_shading_colors:
                    self.set_helix_cluster_colors(helix_proximity_shading_colors,
                                                  helix_cluster_shading_color)

                # color the helices in clusters according to the color
                # list if such is specified using 'simple' color scheme
                if color_scheme[:6] == 'simple':
                    type_color_dict = get_simple_colors(color_scheme)
                    if type_color_dict.has_key('helixcluster'):
                        helixcluster_color_list = type_color_dict['helixcluster']
                    else:
                        helixcluster_color_list = [DEFAULT_SHAPE_COLOR_HEXSTR]
                    helix_cluster_i = 0
                    for helix_cluster in self.helix_cluster_list:
                        for helix in helix_cluster.svgnodelist:
                            helix.set_color_hex(helixcluster_color_list[
                               helix_cluster_i % len(helixcluster_color_list)])
                        helix_cluster_i += 1

            # and build all the connectors
            chain_i = 0
            for nodelist in self.iter_chains():
                self.build_connectors_aligned_svg(nodelist,
                                                  chain_i,
                                                  use_connector_arrowheads,
                                                  connector_color_scheme,
                                                  interdomain_connectors)
                chain_i += 1
        else:
            # build termini
            self.build_termini_svg(interdomain_connectors)

            # build connectors
            chain_i = 0
            for nodelist in self.iter_chains():
                self.build_connectors_svg(
                    nodelist,
                    chain_i,
                    use_connector_arrowheads,
                    connector_color_scheme,
                    interdomain_connectors)
                chain_i += 1


        # build lists of residue names and ids in connectors
        for connector in self.svg_connector_list:
            connector.build_resname_sequence(self.residue_list,
                                             self.pdb_resid_dict)
            
        # build sheet clusters
        self.build_sheet_cluster_constraints_svg(sheet_shading_colors)




    def relpos_to_pos(self, reference_element, new_element, relpos,
                      ref_strand, new_strand,
                      sheet_pos_dict):
        """
        Given a reference element and an element to place, and its position
        relative to the reference element, return the absolute x and y
        co-orindates to place the new element. Uses the information from
        the nodes and the sheet_pos_dict to find the x and y
        co-ordinates of the reference element

        Parameters:
            reference_element - sheet id or PTNodeHelix of the element
                                used as reference
            new_element       - sheet id or PTNodeHelix of elmenet to place
            relpos            - position to place new_element relative to
                                reference_element. ptrelpos.RELPOS_ABOVE, etc.
            ref_strand        - strand in reference element it is relative to
                                or None if not a sheet
            new_strand        - strand in new element used for relative
                                position or None if not a sheet
                        { nodeid : (shape_xmlid, xpos, ypos, indguide_xmlid) }
            sheet_pos_dict    -  y pos of sheet top and bottom
                                and x pos of sheet left and right
                               { sheet_id : (top_y, bot_y, left_x, right_x) }
                               

        Uses data members:
            sheet_strandlists_dict  -
                dictionary of { sheet_id : list of list of nodes }
                where the list of list of nodes is
                described in build_sheet_constraints()
            ptrelpos -
                     PTRelativePosition instance to find relative
                       positions of elements 


        Return value:
            tuple (xpos, ypos) of position to place new_element
        """
        if isinstance(reference_element, PTNodeHelix):
            left_x = reference_element.xpos
            right_x = left_x
            if reference_element.get_sideways():
                right_x += reference_element.get_span() * \
                           DUNNART_HELIX_SPAN_FACTOR
            else:
                right_x += DUNNART_HELIX_WIDTH
            top_y = reference_element.ypos
            if reference_element.get_sideways():
                bot_y = top_y + DUNNART_HELIX_WIDTH
            else:
                bot_y = top_y + (reference_element.get_span() *
                                 DUNNART_HELIX_SPAN_FACTOR)
        else:
            # is a sheet id
            assert(isinstance(ref_strand, PTNodeStrand))
            left_x = sheet_pos_dict[reference_element][2]
            right_x = sheet_pos_dict[reference_element][3]
            top_y = sheet_pos_dict[reference_element][0]
            bot_y = sheet_pos_dict[reference_element][1]

        
        if relpos == RELPOS_ABOVE or relpos == RELPOS_BELOW:
            if isinstance(new_element, PTNodeHelix):
                assert(new_strand == None)
                if ref_strand != None: # position helix relative to sheet
                    xpos = ref_strand.xpos
                    if relpos == RELPOS_ABOVE:
                        if not (ref_strand.get_sideways() or
                                new_element.get_sideways()):
                            ypos = top_y - get_min_gap_size() - \
                                   new_element.get_span() * \
                                   DUNNART_HELIX_SPAN_FACTOR
                        else:
                            ypos = top_y - get_min_gap_size()
                    else:
                        ypos = bot_y + get_min_gap_size()
                else: # position helix relative to helix
                    xpos = left_x
                    if relpos == RELPOS_ABOVE:
                        if (not reference_element.get_sideways() or
                            new_element.get_sideways()):
                            ypos = top_y - get_min_gap_size() - \
                                   new_element.get_span() * \
                                   DUNNART_HELIX_SPAN_FACTOR
                        else:
                            ypos = top_y - get_min_gap_size()
                    else:
                        ypos = bot_y + get_min_gap_size()
            else:
                # new element is a sheet
                assert(new_strand != None)
                new_strand_posnum = \
                                self.ptrelpos.get_strand_posnum(new_strand)

                if (ref_strand != None and
                    isinstance(ref_strand, PTNodeStrand)):
                    # position sheet relative to sheet
                    # need to align ref_strand and new_strand on vert axis
                    # so xpos (which is pos of leftmost strand in sheet)
                    # is difference between dist of ref strand from leftmost
                    # in ref sheet and dist of new strand from leftmost in
                    # new sheet
                    ref_strand_posnum = \
                          self.ptrelpos.get_strand_posnum(ref_strand)
                    offset = ref_strand_posnum - new_strand_posnum
                    offset *= get_strand_separation()
                    xpos = left_x + offset
                    SHEET_GAP_FUDGE_FACTOR = 2.5 # see FIXME just below
                    if relpos == RELPOS_ABOVE:
#                         new_sheet_height = self.get_longest_strand_length(
#                                     self.sheet_strandlists_dict[new_element])
#                         ypos = top_y - new_sheet_height - \
#                              get_min_gap_size()
                        # FIXME really need to 'draw' (build) sheet first
                        # so we can get the height, since strand offsets
                        # mean it need bear no relation to the longest strand
                        # length, for now just something arbitrary
                        ypos = int(top_y - DUNNART_SHEET_GAP_SIZE * SHEET_GAP_FUDGE_FACTOR)
                    else:
                        ypos = int(bot_y + DUNNART_SHEET_GAP_SIZE * SHEET_GAP_FUDGE_FACTOR)

                else:
                    # we are positioning sheet relative to a helix
                    # align the correct strand on the helix vertical axis
                    offset = new_strand_posnum * get_strand_separation()
                    xpos = left_x
                    if relpos == RELPOS_ABOVE:
                        ypos = top_y - DUNNART_SHEET_GAP_SIZE
                        if not new_strand.get_sideways():
                            ypos -= new_strand.get_span() * \
                                    DUNNART_HELIX_SPAN_FACTOR
                    else:
                        ypos = bot_y + DUNNART_SHEET_GAP_SIZE
                    
        else:
            assert(relpos == RELPOS_LEFT or relpos == RELPOS_RIGHT)
            ypos = top_y
            if isinstance(new_element, PTNodeHelix):
                assert(new_strand == None)
                # positioning a helix next to a sheet or helix
                if relpos == RELPOS_LEFT:
                    if new_element.get_sideways():
                        xpos = left_x - get_min_gap_size() - \
                                   new_element.get_span() * \
                                   DUNNART_HELIX_SPAN_FACTOR
                    else:
                        xpos = left_x - get_min_gap_size()
                else:
                    xpos = right_x + get_min_gap_size()
            else:
                # new element is a sheet
                assert(new_strand != None)
                if ref_strand != None:  # positioning a sheet next to a sheet
                    if relpos == RELPOS_LEFT:
                        # first compute width of new sheet
                        new_sheet_left_strand = \
                            self.sheet_strandlists_dict[new_element][0][0]
                        new_sheet_right_strand = \
                            self.sheet_strandlists_dict[new_element][-1][0]
                        new_sheet_width = \
                             len(self.sheet_strandlists_dict[new_element]) * \
                             (DUNNART_STRAND_WIDTH + \
                              get_strand_separation())
                        # and position it starting that far to left of ref sheet
                        xpos = left_x-DUNNART_SHEET_GAP_SIZE-new_sheet_width
                                              

                    else:
                        xpos = right_x + DUNNART_SHEET_GAP_SIZE
                else:  # positioning a sheet next to a helix
                    if relpos == RELPOS_LEFT:
                        xpos = left_x - get_min_gap_size()
                    else:
                        xpos = right_x + get_min_gap_size()


        return (xpos, ypos)
        

    
    def find_nearest_to_any_in_set(self, element_set, exclude_set,
                                   sheets_only=False):
        """
        Find the nearest element (helix or sheet) to any of the elements
        in the supplied set, that is not itself in this set, or in the
        exclude_set. This is
        so we find the nearest element to an already positioned element,
        that is not itself already positioned.
        Uses the PTDistMatrix to do this.

        Parameters:
           element_set - set of PTNodeHelixs and sheet ids to find the nearest
                         element to any of them, that is not itself in this set.
           exclude_set - set of elements (PTNodeHelix obj or set id) to
                         exclude from finding (in addition to those already
                         in the element_set)
           sheets_only - (Default False) only find sheets, not helices.
           
        Uses data members (read):
           distmatrix - The PTDistMatrix that has been built already
        Return value:
           tuple (set_element, close_element)
           where set_element is an element in the supplied set and
           close_element is the
           PTNodeHelix or sheet id of the element which has minimum distance
           to set_element, and this is the minimum distance between any
           element in the element_set and any other element in the dist matrix.
           
        """
        min_dist = float("inf")
        min_dist_objid = None
        set_element = None

        # convert element (PTNodeHelix/sheet id) set to set of objids where
        # the objid is the sheet id e.g. 'A' or helix id e.g. 'HELIX_A_10'
        # FIXME: seems I got myself into a needless mess here with using
        # string ids consistently in PTDistMatrix but mixing PTNodeHelix and
        # sheet id here, necessitating this converion, should make it
        # consistent.
        element_objid_set = Set()
        for element in element_set:
            if isinstance(element, PTNodeHelix):
                objid = element.nodeid
            else:
                objid = element # sheet id
            element_objid_set.add(objid)

        for objid in element_objid_set:
            (close_objid, dist) = \
               self.distmatrix.get_min_distance_objid(
                objid, element_objid_set.union(exclude_set), sheets_only )
            if dist < min_dist:
                min_dist = dist
                set_element_objid = objid
                min_dist_objid = close_objid
#                print 'xxx',min_dist_objid,set_element_objid,dist

        # TODO: this whole objid ('A' for sheet 'HELIX_A_10' for helix etc.) is
        # pretty dodgy (see ptdistmatrix.py also), but as we ensure
        # sheet ids and helix ids are not overlapping it works
        # Should really by using classes or something.
        # Also sometimes I am using (as here) my id strings as indices for
        # dictionaries etc. and sometimes using actual object i.e. PTNode
        # etc... not very consistent.
        if len(min_dist_objid) == 1: # sheet ids are 1 char long, helix ids aren't
            close_element = min_dist_objid # just the sheet id
        else:
            close_element = self.get_node_by_id(min_dist_objid)
        if len(set_element_objid) == 1:
            set_element = set_element_objid
        else:
            set_element = self.get_node_by_id(set_element_objid)
        return (set_element, close_element)


    def aligned_strands_overlap(self, vert_list):
        """
        Return True iff the strands in the vert list (ie strands that
        are aligned on one axis) would overlap if aligned on the axis
        according to their offsets relative to neighbour strands
        (according to H bond patterns).

        Normally this doesn't happen, and in fact one of the criteria
        used to decide to put two strands on opposite sides of their common
        neighbour strand (rather than the same side, see
        ptnode.py has_strand_extent_overlap(), label_strand_bridge_sides())
        is this very criterion. It can, however happen sometimes when
        both alternatives result in overlap and geometric (dihedral angle)
        criteria are used to make the side decision. In such cases this
        function returns True and we then disable the offset alignment
        constraints so overlap can be resolved.

        Parameters:
            vert_list - list of PTNodeStrand that are aligned on an axis

        Return value:
            True if there would be overlap of strands if the strands are
            aligned on their common axis according to their H bonds to
            neighbours (offsets already set in PTNodeStrand by
            build_sheet_constraints(), acessed by get_align_pos()).
        """
        # Note we insist not only on no actual overlap but on a  minimum
        # gap of this amount (in residues, so actual position is this
        # multiplied by DUNNART_SPAN_LENGTH_FACTOR), since Dunnart will
        # mark an overlap constraint if they are even too close.
        MIN_GAP = 1
        
        if len(vert_list) < 2:
            return False
        for i in range(len(vert_list)):
            for j in range(i+1, len(vert_list)):
                strand1 = vert_list[i]
                strand2 = vert_list[j]
                if strand1.get_align_pos() < strand2.get_align_pos():
                    min_offset_strand = strand1
                    other_strand = strand2
                else:
                    min_offset_strand = strand2
                    other_strand = strand1
                if (min_offset_strand.get_align_pos() +
                    min_offset_strand.get_span() +
                    MIN_GAP > other_strand.get_align_pos()):
                    sys.stderr.write('WARNING: disabled offset constraint '
                                     'due to overlap of strands: ' +
                                     str(min_offset_strand) + ', ' +
                                     str(other_strand) +
                                     '\n')

                    return True
                

    def build_sheet_svg(self, sheet_pos_dict,
                        sse_label_scheme, seqnum_dict,
                        sheet_id, 
                        xpos, ypos,
                        label_residue_numbers):
        """
        build SVG for a sheet, with strands in left to right order
        according to the constraints previuosly calculated, including
        parallel/antiparallel relationship between strands (drawn
        as up/down arrows, same direction for parallel).
        
        Also write indguides and their alignments, and distributions.
        Each strand has a vertical indguide
        use to align vertically (required for multiple strands on
        same vertical axis) and these are also used for the
        distribution constraints (to keep uniform separation between
        strands in a sheet).

                
        Parameters:
             sheet_pos_dict - (in/out)
                              dictionary mapping from sheet id to top and
                               bottom of sheet y co-ordinate and left
                               and right x co-orindate:
                               { sheet_id : (top_ypos, bottom_ypos,
                                             left_xpos,right_xpos) }
             sse_label_scheme - 'none','sequential' or 'separate'.
                   Note that this does not affect internal number, it is
                   just for the labels on Dunnart shapes.
             seqnum_dict - dictionary { nodeid : seqnum } only used if
                            sse_label_scheme is 'sequential'
             sheet_id - id of the sheet
             xpos  - x position to start sheet
             ypos  - y position to start sheet
             label_residue_numbers - if True put start and end residue ids
                            on head and tail of helix shape

        Uses data members (readonly):
             node_list - ordered list of nodes
             ptrelpos -
                     PTRelativePosition instance to find relative
                       positions of elements
             xmlid_generator
             (write):
             svg_constraint_list

        Return value:
             None

        """
        start_xpos = xpos
        start_ypos = ypos
        
        horiz_order_list = self.sheet_strandlists_dict[sheet_id]
        prev_vert_list = None

        # NOTE: horiz/vert and x/y and height/width, top/bot
        # names may be somewhat confusing in
        # here now since sideways strands reverse the meanings of these
        # things... (i.e. height is actually width, 'horiz_indguides' are
        # vertical, etc. variable names and comments assume up/down
        # (i.e. not sideways)
        if horiz_order_list[0][0].get_sideways():  # if one strand is, all are
            sheet_is_sideways = True
        else:
            sheet_is_sideways = False

        # keep track of bottom and top of sheet so we can set sheet_pos_dict
        if sheet_is_sideways:
            top_pos = start_xpos
            bot_pos = start_xpos
        else:
            top_pos = start_ypos
            bot_pos = start_ypos

        # FIXME: shouldn't use this longest_strand_length stuff anymore,
        # need to take account of strand offsets etc. (see below too)
        (unused_longest_strand, longest_strand_length) = \
                     self.ptrelpos.get_longest_strand(horiz_order_list)
        longest_strand_length *= DUNNART_SPAN_LENGTH_FACTOR

        distribution_xmlid = None
        if len(horiz_order_list) > 1: #shouldn't get 1 strand sheets, but do
            # write distribution for separation constraints between strands
            # position it just below the bottom horizontal indguide of sheet
            if sheet_is_sideways:
                direction = DUNNART_GUIDE_TYPE_HORI
                distr_pos = xpos + longest_strand_length + 10
                        
            else:
                direction = DUNNART_GUIDE_TYPE_VERT
                distr_pos = ypos + longest_strand_length + 10

            xmlid = self.xmlid_generator.next()
            distribution = PTSVGDistribution(xmlid, direction,
                                             get_strand_separation(),
                                             distr_pos)
            self.svg_constraint_list.append(distribution)

        strand_num = 0
        for horiz_list_index in range(len(horiz_order_list)):
            vert_list = horiz_order_list[horiz_list_index]
            # write vertical indguide for strand(s) at this horizontal pos
            if sheet_is_sideways:
                direction = DUNNART_GUIDE_TYPE_HORI
                pos = ypos
            else:
                direction = DUNNART_GUIDE_TYPE_VERT
                pos = xpos
            xmlid = self.xmlid_generator.next()
            vert_indguide = PTSVGIndGuide(xmlid, pos, direction)
            self.svg_constraint_list.append(vert_indguide)

            # write a distribution constraint for this strand(s) in
            # a single vertical alignment
            # (relative to neighbour vertical alignment in horiz list)
            # note distribution constraints are to the strands'
            # corresponding vertical indguides, not strands themselves
            if strand_num > 0:
                self.svg_constraint_list.append(
                    PTSVGDistroConstraint(prev_vertlist_indguide,
                                          vert_indguide,
                                          distribution))

            else:
                # for the first strand, write a horizontal indguide to
                # be used for separation constraints which enforce the
                # vertical offset of subsequent strands relative to the
                # first one, that was computed earlier
                # (in build_sheet_constraints())
                if sheet_is_sideways:
                    direction = DUNNART_GUIDE_TYPE_VERT
                    pos = xpos
                else:
                    direction = DUNNART_GUIDE_TYPE_HORI
                    pos = ypos
                xmlid = self.xmlid_generator.next()
                first_horiz_indguide = PTSVGIndGuide(xmlid, pos, direction)
                self.svg_constraint_list.append(first_horiz_indguide)

            if sheet_is_sideways:
                sheet_start_pos = start_ypos
            else:
                sheet_start_pos = start_xpos

            # disable strand 'vertical' alignment offset constraints if
            # there would be overlap of strand shapes on the indguide
            # which would result in an unsatisfiable constraint in Dunnart
            # since we also set non-overlap constraint on.
            enable_offset_constraint = True
            if len(vert_list) > 1 and self.aligned_strands_overlap(vert_list):
                enable_offset_constraint = False
            for node in vert_list:
                pos = self.build_strand_svg(node,
                                            vert_indguide,
                                            first_horiz_indguide,
                                            xpos, ypos,
                                            sse_label_scheme,
                                            seqnum_dict,
                                            strand_num,
                                            sheet_start_pos,
                                            enable_offset_constraint,
                                            label_residue_numbers)
                if pos < top_pos:
                    top_pos = pos
                if pos + node.get_span() * DUNNART_SPAN_LENGTH_FACTOR \
                       > bot_pos:
                    bot_pos = pos + \
                              node.get_span() * DUNNART_SPAN_LENGTH_FACTOR
                strand_num += 1

            if sheet_is_sideways:
                ypos += get_strand_separation()
            else:
                xpos += get_strand_separation()
            prev_vert_list = vert_list
            prev_vertlist_indguide = vert_indguide

        if sheet_is_sideways:
            sheet_pos_dict[sheet_id] = (start_ypos, ypos,
                                        top_pos, bot_pos)
        else:
            sheet_pos_dict[sheet_id] = (top_pos, bot_pos,
                                        start_xpos, xpos)




    def build_strand_svg(self, node, vert_indguide,
                         first_horiz_indguide,
                         xpos, ypos,
                         sse_label_scheme,
                         seqnum_dict,
                         strand_num,
                         sheet_start_pos,
                         enable_offset_align_constraint,
                         label_residue_numbers):
        """
        Function used by build_sheet_svg() to build the SVG XML
        for a single strand. This consists of the actual strand shape
        along with alignment constraints.

        Parameters:
             node - node to write strand for
             vert_indguide -  vertical indguide for this strand
             first_horiz_indguide - horizontal indguide
                                   of the first strand in the sheet, which
                                   is used as the base for all separation
                                   constraints enforcing offset of strands
                                   along y axis (on vert indguide)
             xpos - current X position
             ypos - current Y position
             se_label_scheme - 'none','sequential' or 'separate'.
                    if 'sequential' number all nodes in one sequence
                   instead of sepearte sequences for strands and helices.
                   Note that this does not affect internal number, it is
                   just for the labels on Dunnart shapes.
             seqnum_dict - dictionary { nodeid : seqnum } only used if
                            sse_label_scheme is 'sequential'.
             strand_num -   Starts at 0 for first strand drawn in sheet.
                            Needed since vert alignment is relative
                            to first_horiz_indguide_xmlid which is for first
                            strand so avoid redundant indguide/align to this
                            on first strand only.
                            Also used to place the separation constraint
                            handles for each sheet at slighly different
                            positions.
             sheet_start_pos - For 'sideways' sheets, y, else x position
                               of start of sheet (first strand).
                               Used for placing the separation constraint
                               handles for relative strand offsets.
             enable_offset_align_constraint - Boolean: if True write
                             strand offset ('vertical' alignment constraint
                             (for position 'up or down' along neighbour strand)
                             else don't write the constraint.
             label_residue_numbers - Boolean: if True put start and end
                             residue ids on start and end of strand shape

        Uses data members:
             xmlid_generator

        Return value:
             the ypos of the strand just written
             (or xpos if sideways)

        """

        assert(isinstance(node, PTNodeStrand))

        # NOTE: horiz/vert and x/y and height/width
        # names may be somewhat confusing in
        # here now since sideways strands reverse the meanings of these
        # things... (i.e. height is actually width, 'horiz_indguides' are
        # vertical, etc. variable names and comments assume up/down
        # (i.e. not sideways)
        
        if sse_label_scheme == 'sequential':
            label = str(seqnum_dict[node.nodeid])
        elif sse_label_scheme == 'separate':
            label = str(node.seqnum)
            if self.num_chains() > 1:
                label = label + str(node.get_chainid()).lower()
        else:
            label = ''
            
        # mark 'edges' of barrel (the strands on end where brdige broken to
        # 'flatten' them) by putting an asterisk on the label
        # (TODO: some better way of markign this, color or something?)
        if node.get_barrel_edge():
            label += '*'
        strand_xmlid = self.xmlid_generator.next()

        if node.get_sideways():
            strand_ypos = ypos
            strand_xpos = xpos + node.get_align_pos() * DUNNART_SPAN_LENGTH_FACTOR
        else:
            strand_xpos = xpos
            strand_ypos = ypos + node.get_align_pos() * DUNNART_SPAN_LENGTH_FACTOR

        node.set_svginfo(strand_xmlid, strand_xpos, strand_ypos, label,
                         vert_indguide,
                         str(seqnum_dict[node.nodeid]))

        if label_residue_numbers:
            node.headLabel = node.resid_list[-1]
            node.tailLabel = node.resid_list[0]
        
        # align strands in vert_list on vertical alignment indguide
        if node.get_sideways():
            align_type = DUNNART_ALIGN_MIDDLE
        else:
            align_type = DUNNART_ALIGN_CENTER
        self.svg_constraint_list.append(
            PTSVGAlignmentConstraint(vert_indguide, node, align_type))


        if strand_num > 0: # not the first strand in sheet
            # write an  indguide for this strand to be used
            # for separation constraints which enforce the position
            # i.e. the align_pos
            if node.get_sideways():
                direction = DUNNART_GUIDE_TYPE_VERT
                pos = strand_xpos
            else:
                direction = DUNNART_GUIDE_TYPE_HORI
                pos = strand_ypos
            xmlid = self.xmlid_generator.next()
            this_horiz_indguide = PTSVGIndGuide(xmlid, pos, direction)
            self.svg_constraint_list.append(this_horiz_indguide)
        else:
            this_horiz_indguide = first_horiz_indguide

        # and align this strand on it, if parameter set to allow this
        if enable_offset_align_constraint:
            if node.get_sideways():
                alignment_pos = DUNNART_ALIGN_LEFT
            else:
                alignment_pos = DUNNART_ALIGN_TOP
            self.svg_constraint_list.append(
                PTSVGAlignmentConstraint(this_horiz_indguide, node,
                                         alignment_pos))


            if strand_num > 0: # not the first strand in sheet
                # write separation constraint between horiz indguide of first strand
                # (base for all offsets) and this horiz indguide.
                sephandle_pos = sheet_start_pos - strand_num * 5
                xmlid = self.xmlid_generator.next()
                sep = PTSVGDistribution(
                           xmlid, direction,
                           str(node.get_align_pos()*DUNNART_SPAN_LENGTH_FACTOR),
                           sephandle_pos)
                self.svg_constraint_list.append(sep)
                self.svg_constraint_list.append(
                    PTSVGDistroConstraint(first_horiz_indguide,
                                          this_horiz_indguide,
                                          sep))

        if node.get_sideways():
            strand_pos = strand_xpos
        else:
            strand_pos = strand_ypos
            
        return strand_pos



    def set_helix_svginfo(self, helix,  xpos, ypos,
                          sse_label_scheme,
                          seqnum_dict, label_residue_numbers):
        """
        Set the SVG info for a helix node. called by build_dunnart_svg().
        
        Parameters:
             helix   - PTSVGNodeHelix to set info in
             xpos    - x coorindate to write helix
             ypos    - y coordinate to write helix
             sse_label_scheme - if 'sequential' number all SSEs sequentially
                   instead of sepearte sequences for strands and helices.
                   Note that this does not affect internal number, it is
                   just for the labels on Dunnart shapes.
             seqnum_dict - dictionary { nodeid : seqnum } only used if
                            sse_label_scheme is 'sequential'
             label_residue_numbers - if True put start and end residue ids
                            on head and tail of helix shape

        Uses data members (readonly):
             node_list - ordered list of nodes
             include_pi_helices,include_310_helices - flag to use these or not
             (NOTE: sets is_positioned flag and svg info in helix nodes though)
             xmlid_generator

        Return value:
             None
        """
        assert isinstance(helix, PTSVGNodeHelix)

        if ( (helix.get_type() == "310" and not self.include_310_helices) or
             (helix.get_type() == "PI" and not self.include_pi_helices) ):
            return
        
        if sse_label_scheme == 'sequential':
            label = str(seqnum_dict[helix.nodeid])
        elif sse_label_scheme == 'separate':
            # helices are labelled 'A', 'B', etc. by convention
            # FIXME: should go to AA, AB, etc. if more than 26
            label = chr(ord('A')-1 + helix.seqnum)
            if self.num_chains() > 1:
                label = label + str(helix.get_chainid()).lower()
        else:
            label = ''

        xmlid = self.xmlid_generator.next()
        helix.set_svginfo(xmlid, xpos, ypos, label,
                          str(seqnum_dict[helix.nodeid]))
        if label_residue_numbers:
            helix.headLabel = helix.resid_list[-1]
            helix.tailLabel = helix.resid_list[0]
        
        helix.set_is_positioned(True)


    def get_most_nterm_visible_sse(self, nodelist):
        """
        Return the most N-terminal already placed SSE (helix or strand)
        in the supplied nodelist (ordered from N to C terminus).

        Only used when distance matrix placement is being used, depends
        on the is_positioned flag in PTNode being set.
        
        Parameters:
            nodelist - list of PTNodes ordererd from N- to C-terminus.

        Retrun value:
            PTNode that is most N-terminal which is already positioned
            or None if none found (should not happen)
        """
        i = 1 # start at second in nodelist, first is N-terminal pseudonode
        while i < len(nodelist) and not nodelist[i].get_is_positioned():
            i += 1
        if i < len(nodelist):
            return nodelist[i]
        else:
            return None

    def get_most_cterm_visible_sse(self, nodelist):
        """
        Return the most C-terminal already placed SSE (helix or strand)
        in the supplied nodelist (ordered from N to C terminus).

        Only used when distance matrix placement is being used, depends
        on the is_positioned flag in PTNode being set.
        
        Parameters:
            nodelist - list of PTNodes ordererd from N- to C-terminus.

        Retrun value:
            PTNode that is most C-terminal which is already positioned
            or None if none found (should not happen)
        """
        # start at second-last in nodelist, last is C-terminal pseudonode
        i = len(nodelist) - 2
        while i > 0 and not nodelist[i].get_is_positioned():
            i -= 1
        if i > 0:
            return nodelist[i]
        else:
            return None
        
    def build_termini_svg(self, interdomain_connectors=False):
        """
        Build SVG all the N- and C- terminus nodes (may be multiple for multi
        chains, and 'pseudo' terminus nodes for breaks in chain due to
        domain decomposition).

        Parameters:
             terminus   - PTNodeTerminus to write
              interdomain_connectors - If True, do NOT make
                       pseduo-terminus nodes
                       at domain boundaries. Instead the domain boundary
                       SSEs are left ot hve connectors to other domain
                       added later. Default False.
                       

        Uses data members (readonly):
             chain_dict - dict of chainid : node_list

        Return value:
             None
        """
        # since node list for each chain is sorted by PDB residue sequence
        # number (ascending), the N terminal node is the first and the
        # C terminal is the last.
        for nodelist in self.iter_chains():
            # position terminus symbol near corresponding
            # most N- or C- terminal element
            for term_node in [nodelist[0], nodelist[-1]]: # N-term, C-term
                # build all (include pseudo) terminus nodes if not using
                # interdomain connectors, and (even if we are using
                # interdomain conenctors), always build non-psuedo termini
                if (not interdomain_connectors or not term_node.get_pseudo()):
                    (xpos, ypos) = self.find_terminus_pos(term_node,
                                                          nodelist)
                        
                    # TODO: align on strand alignment indguide
                    
                    self.set_terminus_svginfo(term_node, xpos, ypos)


         
    def find_terminus_pos(self, term_node, nodelist):
        """
        Find the position to place the supplied (N- or C-) terminus node,
        byt finding the nearest in sequence visible SSE and returning
        a position nearby.
        
        Parameters:
             term_node   - PTNodeTerminus to write
             nodelist - list of PTNodes in this chain (containing term_node)

        Return value:
          (xpos, ypos) tuple to place terminus at.
        """
        assert(isinstance(term_node, PTNodeTerminus))
        if term_node.get_termtype() == 'N':
            is_n_term = True
        else:
            is_n_term = False
            
        if is_n_term:
            nearest_sse = self.get_most_nterm_visible_sse(nodelist)
        else: # c-term
            nearest_sse = self.get_most_cterm_visible_sse(nodelist)
            
        xpos = nearest_sse.xpos
        ypos = nearest_sse.ypos

        if ( (is_n_term and nearest_sse.get_reversed()) or 
             (not is_n_term and not nearest_sse.get_reversed()) ):
            if nearest_sse.get_sideways():
                xpos -= get_min_gap_size()
            else:
                ypos -= get_min_gap_size()
        else:
            if isinstance(nearest_sse, PTNodeStrand):
                span_length_factor = DUNNART_SPAN_LENGTH_FACTOR
            else:
                span_length_factor = DUNNART_HELIX_SPAN_FACTOR
            offset = nearest_sse.get_span() *  span_length_factor + \
                     get_min_gap_size()
            if nearest_sse.get_sideways():
                xpos += offset
            else:
                ypos += offset

        return (xpos, ypos)

        
    def set_terminus_svginfo(self, terminus, xpos, ypos):
        """
        Set the SVG infor for a N- or C- terminus node.
        
        
        Parameters:
             terminus   - PTNodeTerminus to write
             xpos    - x coorindate to write at
             ypos    - y coordinate to write at

        Uses data members (readonly):
             node_list - ordered list of nodes
             xmlid_generator

        Return value:
             None
        """
        assert isinstance(terminus, PTSVGNodeTerminus)

        xmlid = self.xmlid_generator.next()
        label = terminus.nodeid # 'C' or 'N' or 'Ca' etc.
        terminus.set_svginfo(xmlid, xpos, ypos, label)
        terminus.set_is_positioned(True)
    


    def build_sheet_cluster_constraints_svg(self,
                                            sheet_shading_colors):
        """
        Build cluster constraints for grouping
        strands into sheets. Called by build_dunanrt_svg().

        Parameters:
           sheet_shading_colors - None (use default shade for all) or
                                  'auto' (use color gradient to shade each
                                  differently) or list of colors.

        Uses data members (readonly):
             node_list - ordered list of nodes
             (write):
             sheet_cluster_list - list of SVGCluster objects, appended to here.
             xmlid_generator

        Return value:
             None
        
        """
        num_sheets = len(list(self.sheet_dict.iteritems()))
        if num_sheets < 1:
            return

        if (sheet_shading_colors):
            cluster_fill_colors = get_cluster_fill_colors(sheet_shading_colors,
                                                          num_sheets)
                    
        i = 0
        for (sheet_id, ptnode_list) in self.sheet_dict.iteritems():
            if (sheet_shading_colors):
                cluster_fill_color =  cluster_fill_colors[i]
            else:
                cluster_fill_color = DUNNART_DEFAULT_CLUSTER_FILL_COLOR
            self.sheet_cluster_list.append(PTSVGCluster(
                                               ptnode_list,
                                               self.xmlid_generator.next(),
                                               sheet_id,
                                               0,
                                               cluster_fill_color))
            i += 1


    def build_connectors_svg(self, nodelist,
                             chain_i,
                             use_connector_arrowheads=False,
                             connector_color_scheme = 'all',
                             interdomain_connectors = False):
        """
        Called by build_dunnart_svg() to 
        build SVG for  connectors for sequence in a single chain:
        a connector between each
        node (helix/strand) in sequence order.

        This new version is for use with the -i (distance matrix instead
        of heuristic helix placement) option i.e. the greedy placement
        algorithm. It places connectors on ports of helices in the same
        way as strands, i.e. using their orientation ('reversed' flag)
        as set by tableau information.
        
        Parameters:
             nodelist - list of nodes in this chain
             chain_i  - chain index (0,1,...) for selecting line color
             use_connector_arrowheads - If True make connectors directed.
             connector_color_scheme -  'all[:<color>]', 'chain[:<color_list>]',
                                       'domain[:<intra_color>,<inter_color>]',
                                       'crossing[:<color_list>]' 
             interdomain_connectors - If True, do NOT make
                       pseduo-terminus nodes
                       at domain boundaries. Instead the domain boundary
                       SSEs are left ot hve connectors to other domain
                       added later. Default False.
                       

        

        Uses data members (readonly):
             node_list - ordered list of nodes (NB sets port fields in SVGNodes)
             include_pi_helices, include_310_helices
             xmlid_generator


        Return value:
              None

        Precondition: nodelist is sorted (by start res seq ascending);
                      this is done by build_graph_from_secstruct()
                      before calling.
        """
        prevnode = None
        for nodeindex in range(len(nodelist)):
            node = nodelist[nodeindex]

            if ( isinstance(node, PTNodeHelix) and
                 ( (node.get_type() == "310" and
                    not self.include_310_helices) or
                   (node.get_type() == "PI" and
                    not self.include_pi_helices) ) ):
                continue  # skip pi/310 helix if flagged to do so

            if ( interdomain_connectors and
                 isinstance(node, PTNodeTerminus) and node.get_pseudo() ):
                continue # don't do pseudo for interdomain

            if prevnode != None:
                # add connector from prevnode to node
                if isinstance(prevnode, PTNodeTerminus):
                    srcFlags = DUNNART_DEFAULT_PORT
                else:
                    src_reversed = prevnode.get_reversed()
                    if src_reversed:
                        srcFlags = DUNNART_BOTTOM_PORT
                    else:
                        srcFlags = DUNNART_TOP_PORT
                    if prevnode.get_sideways():
                        if srcFlags == DUNNART_TOP_PORT:
                            srcFlags = DUNNART_LEFT_PORT
                        else:
                            srcFlags = DUNNART_RIGHT_PORT

                if isinstance(node, PTNodeTerminus):
                    dstFlags = DUNNART_DEFAULT_PORT
                else:
                    dst_reversed = node.get_reversed()
                    if dst_reversed:
                        dstFlags = DUNNART_TOP_PORT
                    else:
                        dstFlags = DUNNART_BOTTOM_PORT
                    if node.get_sideways():
                        if dstFlags == DUNNART_TOP_PORT:
                            dstFlags = DUNNART_LEFT_PORT
                        else:
                            dstFlags = DUNNART_RIGHT_PORT

                # line color may be overwritten later for multidomain
                linecolor = get_line_color(connector_color_scheme, chain_i)

                node.nterm_port = dstFlags
                prevnode.cterm_port = srcFlags
                    
                xmlid = self.xmlid_generator.next()
                self.svg_connector_list.append(
                    PTSVGConnector(xmlid, prevnode, node, srcFlags, dstFlags,
                                   linecolor, use_connector_arrowheads))

            prevnode = node


    def dfs_strands_seq(self, start_strand, visited, dfs_list, from_node,
                        component_id=None):
        """
        Make a depth-first search traversal of STRAND nodes
        using sequence (not bridge)
        edges starting at the specfied strand,
        returning list of (node,from_node) tuples in DFS traversal order
        where from_node is the node from which node is reached.

        Note these sequence 'edges' are simply implied by the order of
        nodes in list: a node has a sequene edge to the ones immediately
        before and after it in sequence along chain from N to C terminus.
        
        Parameters:
           start_strand - STRAND node to start at
           visited - (in/out) dictionary of {ptnode:True} visited nodes
           dfs_list - (in/out) list of (node, from_node) visited in dfs order
           from_node - node from which we are being (recursively) called
           component_id - identifier of this connected component to mark
                      each strand in it with, or None to not mark at all
                      (default).
                      

        Recursive function. call initially as
            dfslist = []
            dfs_strands_from(startnode, {}, dfslist, None)

        Return value:
            None. (output is dfs_list parameter)

        """
        visited[start_strand] = True
        if component_id != None:
            start_strand.set_seq_component_id(component_id)
        dfs_list.append((start_strand,from_node))

        # get list of strands (can only be max 2) adjacent in sequence
        # to this one and in the same sheet
        sequence_adjacent_nodes = []
        chainid = start_strand.get_chainid()
        sheetid = start_strand.get_sheet_id()
        nodelist = self.chain_dict[chainid]
        indx = nodelist.index(start_strand)
        if (indx > 1 and isinstance(nodelist[indx-1], PTNodeStrand) and
            nodelist[indx-1].get_sheet_id() == sheetid):
            sequence_adjacent_nodes.append(nodelist[indx-1])
        if (indx < len(nodelist)-2 and
            isinstance(nodelist[indx+1], PTNodeStrand) and
            nodelist[indx+1].get_sheet_id() == sheetid):
            sequence_adjacent_nodes.append(nodelist[indx+1])
        
        for node in sequence_adjacent_nodes:
            if node not in visited:
                self.dfs_strands_seq(node, visited, dfs_list, start_strand,
                                     component_id)


    def find_connected_components_seq_sheet(self):
        """
        Find the connected components (considering only STRAND nodes
        and sequence [not brdige] edges each sheet).

        This is done by a DFS traversal at every strand in the sheet
        (skipping already visited ones), giving us the partition of
        the graph of strand nodes in sheet into connected components.

        Used by set_sheet_fold_color() to color strands in sheet
        in componenets where those adjacent in sequence are colored the
        same.
        
        Parameters: None
        Return value: the number of components labelled
        Uses member data:
            sheet_dict - dict of {sheet_id : ptnode_list} represneting sheets 
            chain_dict - dict by chainid of list
                        of PTNodes in the graph (modifies PTNodes not list)

            (WRITE):
            Labels each strand node with the sheet connected component
            id it belongs to as it goes, starting from 0.
        """

        component_id = 0
        visited = {}   # dictionary of {ptnode : True} visited nodes

        for (sheet_id, ptnode_list) in self.sheet_dict.iteritems():
            for node in ptnode_list:
                if node not in visited:
                    connected_node_list = []
                    self.dfs_strands_seq(node,visited,connected_node_list,None,
                                         component_id)
                    component_id += 1
        return component_id



    def set_sheet_fold_color(self):
        """
        Set the colors of starnds in sheets for the 'fold' color scheme.
        For each sheet colors strands that are connected only by
        turns (i.e. are consecutive in sequence) one color (maybe more
        than one set of such conseuctive strands in sheet, a different
        color for each such set), and other strands (i.e. not in a sequence)
        another color(s).

        Parameters: None
        Return value: None
        Uses data members:
             nodelist, chain_dict, sheet_dict etc. 
             Modifies nodes by set_color_hex() or set_color() methods
        """

        num_components = self.find_connected_components_seq_sheet()

        if num_components == 0:
            return

        # get list of contrasting colors
        # see Glasbey et al 2007
        color_list = get_glasbey_colors_rgb()

        if len(color_list) < num_components:
            sys.stderr.write('WARNING: not enough colors in high contrast list')

        for (sheet_id, strandlist) in self.sheet_dict.iteritems():
            for strand in strandlist:
                strand.set_color(color_list[strand.get_seq_component_id() % len(color_list)])
                
            
    def set_nodelist_colorslist(self,  type_color_dict, typekey):
        """
        Utility function used by set_node_colors() to set all the nodes
        in the nodelist, that are all of the same type, with the color
        list from from type_color_dict of type typekey, which correpsonds to the
        type of nodes in nodelist (310 helix, pi helix, alpha helix).

        Nodes are colored in order of the colors in the list for the
        corresponding type.  If they run out (more elements than
        colors in list, the list is treated as circular, ie colors
        reused from the start of list again).  If no colors for that
        type, DEFAULT_SHAPE_COLOR_HEXSTR is used.

        Parameters:
            type_color_dict -  dict { type : color_list } where type is
                              'sheet','helixcluster',
                               'alpha','pi' or '310' and color_list is list
                                of color value strings where the color
                                value strings are 'rrggbb' color strings,
                                as returned by get_simple_colors()
             typekey - key to type_color_dict ie 'alpha', 'pi' etc

        Uses data members:
            Sets values in nodes in  nodelist accessed via iter_helices()
        """
        if type_color_dict.has_key(typekey):
            color_list = type_color_dict[typekey]
        else:
            color_list = [DEFAULT_SHAPE_COLOR_HEXSTR]
        node_i = 0
        for node in [n for n in self.iter_helices() if
                     n.get_type() == typekey.upper()]:
            node.set_color_hex(color_list[node_i % len(color_list)])
            node_i += 1


    def set_node_colors(self, color_scheme):
        """
        Set the color tuple in each node with a color in a gradient
        from blue at the N terminal to red at the C terminal in each
        chain.

        Parameters:
           color_scheme: string, one of the following:
               'none'     : set each shape to default color
               'simple'   : set strands in sheets from one list of colors,
                            helices in clusters (if any) to a 2nd list of colors,
                            alpha, pi and 310 helices each from their own
                            lists of colors.
               'gradient' : color from blue to red
                            along sequence from N to C terminus.
               'sheet'    : color the strands in each sheet a different
                            color (leaving helices default color)
               'fold'     : color consecutive sequence strands in sheet
                            one color, others another.
        Return value: None
        Uses data members:
             nodelist, chaindict, sheet_dict etc. (via iter_ functions)
             Modifies nodes by set_color_hex() or set_color() methods

        Raises Exceptions:
             ValueError for unknown color_scheme
         
        """
        if color_scheme[:6] == 'simple':
            type_color_dict = get_simple_colors(color_scheme)

            # first do helices

            self.set_nodelist_colorslist(type_color_dict, "alpha")
            if self.include_310_helices:
                self.set_nodelist_colorslist(type_color_dict, "310")
            if self.include_pi_helices:
                self.set_nodelist_colorslist(type_color_dict, "pi")
                
            # NB if helix is part of a cluster, will be
            # overwritten in build_helices_heuristic()
            # if helix clustering is enabled


            # now do strands in sheets

            if type_color_dict.has_key('sheet'):
                sheet_color_list = type_color_dict['sheet']
            else:
                sheet_color_list = [DEFAULT_SHAPE_COLOR_HEXSTR]
            sheet_i = 0
            for (sheet_id, ptnode_list) in self.sheet_dict.iteritems():
                for strand in ptnode_list:
                    assert(isinstance(strand, PTNodeStrand))
                    strand.set_color_hex(
                             sheet_color_list[sheet_i % len(sheet_color_list)] )
                sheet_i += 1

            # and termini
            if type_color_dict.has_key('terminus'):
                color_list = type_color_dict['terminus']
            else:
                color_list = [DEFAULT_SHAPE_COLOR_HEXSTR]
            node_i = 0
            for node in [n for n in self.iter_nodes() if
                         isinstance(n, PTNodeTerminus)]:
                node.set_color_hex(color_list[node_i % len(color_list)])
                node_i += 1

        else:
            done_color = False # used for things that don't need per-chain code
            for nodelist in self.iter_chains():
                # local_nodelist omits pi and/or 310 helices if we are not to
                # draw them
                local_nodelist = [ node for node in nodelist if
                                   not ( isinstance(node, PTNodeHelix) and
                                         ( (node.get_type() == "310" and
                                            not self.include_310_helices) or
                                           (node.get_type() == "PI" and
                                            not self.include_pi_helices) ) ) 
                                  ]
                if color_scheme == 'gradient':
                    rgb_list = list(color_gradient(len(local_nodelist)))
                    assert(len(rgb_list) == len(local_nodelist))
                    for i in range(len(local_nodelist)):
                        local_nodelist[i].set_color(rgb_list[i])
                        i += 1
                elif color_scheme == 'sheet':
                    if not done_color:
                        num_sheets = len(list(self.sheet_dict.iteritems()))
                        # set up list of strand colors, one per sheet, by using
                        # the color gradient to ensure all sheets have different
                        # color.
                        if num_sheets > 0:
                            sheet_colors = [ rgb for rgb in color_gradient(num_sheets) ]
                            i = 0
                            for (sheet_id, ptnode_list) in self.sheet_dict.iteritems():
                                for strand in ptnode_list:
                                    assert(isinstance(strand, PTNodeStrand))
                                    strand.set_color(sheet_colors[i])
                                i += 1
                        done_color = True
                    # set all helices and termini to default color
                    for node in local_nodelist:
                        if not isinstance(node, PTNodeStrand):
                            node.set_color_hex(DEFAULT_SHAPE_COLOR_HEXSTR)
                elif color_scheme == 'fold':
                    if not done_color:
                        self.set_sheet_fold_color()
                        done_color = True
                    # set all helices and termini to default color
                    for node in local_nodelist:
                        if not isinstance(node, PTNodeStrand):
                            node.set_color_hex(DEFAULT_SHAPE_COLOR_HEXSTR)
                elif color_scheme == 'none':
                    for node in local_nodelist:
                        node.set_color_hex(DEFAULT_SHAPE_COLOR_HEXSTR)
                else:
                    raise ValueError('unknown color scheme ' + color_scheme)


    ##########################################################################

    #
    # functions for heuristic helix placement etc. used when
    # not using the -i option
    #
  
    def build_interstrand_helices_svg(self, sheet_pos_dict,
                                     sse_label_scheme,
                                     seqnum_dict,
                                      label_residue_numbers):
        """
        Build SVG  helices that match the special case of being between
        two strands in the same vertilist of a sheet. We draw these
        beside the relevant sheet insead of along strands axis to
        make diagram much neater. See e.g. 2QP2-1 (using STRIDE and DDOMAIN).
        
        Called by build_dunnart_svg()

        
        Parameters:
             sheet_pos_dict - 
                              dictionary mapping from sheet id to top and
                               bottom of sheet y co-ordinate and left
                               and right x co-orindate:
                               { sheet_id : (top_ypos, bottom_ypos,
                                             left_xpos,right_xpos) }
             sse_label_scheme - if 'sequential' number all nodes in one sequence
                   instead of sepearte sequences for strands and helices.
                   Note that this does not affect internal number, it is
                   just for the labels on Dunnart shapes.
             seqnum_dict - dictionary { nodeid : seqnum } only used if
                           sse_label_scheme is 'sequential'
             label_residue_numbers - if True put start and end residue ids
                            on head and tail of helix shape

        Uses data members (readonly):
             ptrelpos -
                    PTRelativePosition instance, to get strand positions 

             node_list - ordered list of nodes
             (NOTE: sets is_positioned,is_interstrand flags in helices though)

        Return value:
             None
        """
        for nodelist in self.iter_chains():
            helix_index = 0
            while helix_index < len(nodelist):
                helix = nodelist[helix_index]
                if not isinstance(helix, PTNodeHelix):
                    helix_index += 1
                    continue  # only helices handled here

                (sheet_id, nterm_strand_index, cterm_strand_index) = \
               self.immediate_containing_sheet(helix_index, helix.get_chainid())
                if (sheet_id == None or nterm_strand_index == None or
                    cterm_strand_index == None):
                    helix_index += 1
                    continue  # not between strands in same sheet
                nterm_strand = nodelist[nterm_strand_index]
                cterm_strand = nodelist[cterm_strand_index]
                sheet_id = nterm_strand.get_sheet_id()
                assert(sheet_id == cterm_strand.get_sheet_id())
                nterm_strand_posnum = \
                               self.ptrelpos.get_strand_posnum(nterm_strand)
                cterm_strand_posnum = \
                               self.ptrelpos.get_strand_posnum(cterm_strand)
                if (nterm_strand_posnum == cterm_strand_posnum):
                    # this helix is between two strands in alignment in sheet,
                    # draw it beside sheet
                    if verbose:
                        sys.stderr.write('interstrand helices: helix '
                                         +str(helix) +
                                         ' between strands ' +
                                         str(nterm_strand)   + ',' 
                                         + str(cterm_strand) + '\n')
                    helix.set_is_interstrand(True)
                    helix.set_reversed(nterm_strand.get_reversed())
                    helix.set_sideways(nterm_strand.get_sideways())
                    if (nterm_strand_posnum >
                        len(self.sheet_strandlists_dict[sheet_id])/2):
                        relpos = RELPOS_RIGHT
                    else:
                        relpos = RELPOS_LEFT
                    if nterm_strand.get_sideways():
                        if relpos == RELPOS_RIGHT:
                            relpos = RELPOS_BELOW
                        else:
                            relpos = RELPOS_ABOVE
                    (xpos, ypos) = self.relpos_to_pos(sheet_id,
                                                      helix,
                                                      relpos,
                                                      nterm_strand,
                                                      None,
                                                      sheet_pos_dict)
                    self.set_helix_svginfo(helix,
                                           xpos, ypos,
                                           sse_label_scheme,
                                           seqnum_dict,
                                           label_residue_numbers)
                helix_index += 1



    def immediate_containing_sheet(self, node_index, chainid):
        """
        Return the sheet id of of sheet 'immediately containing' the supplied
        node. This is the sheet that the supplied node
        (usually a helix) 'interrupts', i.e. if the node is between
        two strands of the same sheet, return the id of that sheet.
        Returns None if there is no such sheet, e.g. if the node
        is between two strands which are not in the same sheet.
        Also returns the index in the node_list of the strands found
        before and after the supplied node (one or both may be None)
        that were used to determine the sheet.

        This is done simply by checking forward (towards C-terminal) and
        backward (towards N-terminal) from the node until a strand is
        found, and returning the sheet id of the strands if they are the
        from the same sheet. (Note this means that if the node to test
        is a strand in a sheet, then that sheet id will be returned, though
        this isn't the useful case of this function).

        Parameters:
            node_index - index in the node list  of the
                         PTNode to find immediate contianing sheet for.
            chainid - chainid of the node, this is the key of chain_dict
                      to get the nodelist that node_index is an index of.

        Return value:
            tuple (sheet_id, nterm_strand_index, cterm_strand_index)
            where the cterm_ and nterm_ strand indices are the indices
            in the node_list of the strands in the N- and C- terminal
            directions that were used to determine the sheet_id.
            Any (or all) of these may be None.

        Uses data members (readonly):
            chain_dict - dict by chainid of ordered list of PTNodes
            
        """
        nodelist = self.chain_dict[chainid]
        # find first strand towards N-terminal
        nterm_strand = None
        nterm_strand_index = None
        for i in range(node_index - 1, -1, -1):
            if isinstance(nodelist[i], PTNodeStrand):
                nterm_strand_index = i
                nterm_strand = nodelist[i]
                break
        # find first strand towards C-terminal
        cterm_strand = None
        cterm_strand_index = None
        for i in range(node_index + 1, len(nodelist)):
            if isinstance(nodelist[i], PTNodeStrand):
                cterm_strand_index = i
                cterm_strand = nodelist[i]
                break
        if cterm_strand != None and nterm_strand != None:
            if nterm_strand.get_sheet_id() == cterm_strand.get_sheet_id():
                return (nterm_strand.get_sheet_id(), nterm_strand_index,
                        cterm_strand_index)
        return (None, nterm_strand_index, cterm_strand_index)



    def count_helices_on_strand_axis(self, strand, strand_pos_dict):
        """
        Count the number of helices above and below (or left and right of)
        the supplied strand. Called by build_helices_svg()
        and room_for_helix_cluster().

        Parameters:
           strand - PTNodeStrand of strand to count helices in alignement
           strand_pos_dict - The strand_ypos_dict from build_helices_svg():
                    { nodeid : (ypos_above, ypos_below, sideways,
                    num_helices_aligned_above, num_helices_aligned_below) }
                    e.g. { 'STRAND_1' : (70, 330, False, 0, 1) }

        Return value:
          tuple (num_above, num_below) where num_above is number of helices
          above (or left of) strand on same axis, and num_below is the number
          below (or right of) strand on same axis.

        Uses data members: None
        """
        (ypos_above_unused, ypos_below_unused, sideways,
         num_helices_aligned_above_strand,
         num_helices_aligned_below_strand) = \
                                     strand_pos_dict[strand.nodeid]
        return (num_helices_aligned_above_strand,
                num_helices_aligned_below_strand)
        
#         xpos = strand.xpos
#         ypos = strand.ypos

#         if sideways:
#             alignpos = ypos
#             slidepos = xpos
#         else:
#             alignpos = xpos
#             slidepos = ypos
#         num_above = 0
#         num_below = 0
#         for node in self.iter_nodes():
#             nodeid = node.nodeid
#             xpos = node.xpos
#             ypos = node.ypos
#             if sideways:
#                 testpos = ypos
#                 testslidepos = xpos
#             else:
#                 testpos = xpos
#                 testslidepos = ypos
#             if (testpos == alignpos and
#                 isinstance(node, PTNodeHelix)):
#                 if testslidepos < slidepos:
#                     num_above += 1
#                 else:
#                     num_below += 1

#         return (num_above, num_below)


    def room_for_helix_cluster(self, seq_strand, seqpos,
                               strand_pos_dict):
        """
        Return True if there would be room to place the helix cluster
        before/after (seqpos) the seq_strand. Else False.
        Checks for helices aligned on strand axes on the same side of the
        sheet as the cluster would be aligned. If there are any then
        we decide that the helix cluster may not fit there.
    
        Called by build_helices_svg().

        Parameters:
           seq_strand - PTNodeStrand of strand the helix cluster would
                        be aligned (by its first helix) to
           seqpos - SEQPOS_AFTER or _BEFORE the seq_strand 
           strand_pos_dict - The strand_ypos_dict from build_helices_svg():
                    { nodeid : (ypos_above, ypos_below, sideways,
                    num_helices_aligned_above, num_helices_aligned_below) }
                    e.g. { 'STRAND_1' : (70, 330, False, 0, 1) }

        Return value:
           True if there would be room to place the helix cluster
           before/after (seqpos) the seq_strand. Else False.
           
        Uses data members:
            sheet_strandlists_dict - the dictionary (by sheet id)
                               of list of lists built by
                               build_constraints()
        """
        sheetid = seq_strand.get_sheet_id()
        if ( (seq_strand.get_reversed() and  seqpos == SEQPOS_AFTER) or
             (not seq_strand.get_reversed() and seqpos == SEQPOS_BEFORE) ):
            cluster_relpos = RELPOS_BELOW
        else:
            cluster_relpos = RELPOS_ABOVE
#        print 'ddddd',sheetid,ptrelpos_to_str(cluster_relpos)
        helices_found = 0
        for vertlist in self.sheet_strandlists_dict[sheetid]:
            for strand in vertlist: # FIXME: should not need all in vertlist
                (num_helices_above_strand,
                 num_helices_below_strand) = \
                 self.count_helices_on_strand_axis(strand,
                                                   strand_pos_dict)
#                print 'fff',strand,num_helices_above_strand,num_helices_below_strand
                if cluster_relpos == RELPOS_BELOW:
                    helices_found += num_helices_below_strand
                else:
                    helices_found += num_helices_above_strand
                if helices_found > 0:
#                    print 'eeeee',helices_found
                    return False
        return True
                    
                    
        

    def helix_is_sheet_gap(self, helix, nterm_strand, cterm_strand, seqpos,
                           seq_strand, sheet_posmap,
                           is_helix_cluster=False):
        """
        Return True if the supplied helix would be undesirably forcing
        a gap between sheets when placed on the axis using the
        build_helices_svg () method. That is, if it is to be
        positioned on an axis whic is between two sheets that have
        been placed next to each other.

        Note we do NOT return True if the helix is actually in sequence
        between strands in each of the sheets, and those strands have
        the same direction, then we DO want to place
        it between them, UNLESS the helix is part of a cluster, then
        we still don't.

        Parameters:
           helix - PTNodeHelix to test
           nterm_strand - next strand along chain in N-terminal direction
           cterm_strand - next strand along chain in C-terminal direction
           seqpos - SEQPOS_AFTER or _BEFORE as calculted in build_helices_svg()
           seq_strand - the PTNode strand whose axis the helix would be placed
                        on. Must be one of nterm_strand or cterm_strand
           sheet_posmap - the PTPosMap of sheet positionings
           is_helix_cluster - the helix is part of a helix cluster
                             default False

        Retrun value:
          True if we do not want to position the helix on the axis as we
          normally would, False otherwise.

        """
        assert(seq_strand != None)
        assert(seq_strand == nterm_strand or seq_strand == cterm_strand)
        assert(seqpos == SEQPOS_AFTER or seqpos == SEQPOS_BEFORE)
        sheet_id = seq_strand.get_sheet_id()
        if cterm_strand != None and cterm_strand != seq_strand:
            other_sheet_id = cterm_strand.get_sheet_id()
        elif nterm_strand != None and nterm_strand != seq_strand:
            other_sheet_id = nterm_strand.get_sheet_id()
        else:
            other_sheet_id = None
        sideways = seq_strand.get_sideways()
        if sideways:
            if seqpos == SEQPOS_AFTER:
                direction = RELPOS_LEFT
            else:
                direction = RELPOS_RIGHT
        else:
            if seqpos == SEQPOS_AFTER:
                direction = RELPOS_ABOVE
            else:
                direction = RELPOS_BELOW

        try:
            sheet_neighbours = sheet_posmap[sheet_id]
#            print 'yyyy',sheet_id,':',str(sheet_neighbours),sideways,ptrelpos_to_str(direction)
        except KeyError:
            return False # no neighbouring sheets

        neighbour_sheet_id = sheet_neighbours.get_neighbour(direction)

        if neighbour_sheet_id == None:
            return False # no neighbour sheet in this direction
            
        if (neighbour_sheet_id == other_sheet_id and
            nterm_strand.get_sideways() == cterm_strand.get_sideways() and
            nterm_strand.get_reversed() == cterm_strand.get_reversed() and
            not is_helix_cluster):
            return False # in seq between same direction strands in each sheet
        else:
            if verbose:
                sys.stderr.write('helix_is_sheet_gap: ' + str(helix) +
                                 ' between sheet ' + sheet_id + ' and sheet ' +
                                 neighbour_sheet_id + '\n')
            return True


    def find_helix_xypos(self, helix, seq_strand, sheet_posmap,
                        sheet_pos_dict,
                        enable_sheet_gap_rule):
        """
        Get position the helix near the sheet seq_strand is in
        according to position map for relative position and tableau
        for orientation

        Parameters:
            helix - PTNode helix to get position for
            seq_strand - strand that the helix is in sequence next to
            sheet_posmap - The PTPosMap of sheet neighbour relationships
            sheet_pos_dict - 
                         dictionary mapping from sheet id to top and
                               bottom of sheet y co-ordinate and left
                               and right x co-orindate:
                               { sheet_id : (top_ypos, bottom_ypos,
                                             left_xpos,right_xpos) }
            enable_sheet_gap_rule - If True, don't align helices on strand
                            indguides between neihbouring sheets, if the
                            total helix length is above a threshold.

        Return value:
           tuple (xpos, ypos) to place the helix at
        """
        # position the helix near the sheet seq_strand is in
        # according to position map for relative position and tableau
        # for orientation
        (relpos, pos_strand, cur_strand) = \
             self.ptrelpos.get_relative_position(seq_strand.get_sheet_id(),
                                                 helix)
        if enable_sheet_gap_rule:
            sheet_id = seq_strand.get_sheet_id()
            if sheet_posmap.has_key(sheet_id):
                # don't positino the helices where they would
                # get between neighbouring sheets, place them
                # where there is nothing beside the sheet if possible
                sheet_neighbours = sheet_posmap[sheet_id]
                if (sheet_neighbours.get_neighbour(relpos) != None):
                    if verbose:
                        sys.stderr.write('sheet gap rule: cannot place ' +
                                         str(helix) + ' ' +
                                         ptrelpos_to_str(relpos) + ' ' +
                                         'sheet ' + sheet_id + '\n')
                    # FIXME: arbitrary helix positioning here
                    if (sheet_neighbours.north == None):
                        relpos = RELPOS_ABOVE
                    elif (sheet_neighbours.south == None):
                        relpos = RELPOS_BELOW
                    elif (sheet_neighbours.east == None):
                        relpos = RELPOS_RIGHT
                    elif (sheet_neighbours.west == None):
                        relpos = RELPOS_LEFT
                    else:
                        sys.stderr.write('WARNING: nowhere to place ' +
                                         str(helix) + ' relative to sheet '
                                         + sheet_id + ' (sheet gap rule)\n')
                        # FIXME: need to do something about this...
                    if verbose:
                        sys.stderr.write('  sheet gap rule: placing ' +
                                         str(helix) + ' ' +
                                         ptrelpos_to_str(relpos) + ' of '
                                         'sheet ' + sheet_id + '\n')

        (xpos, ypos) = self.relpos_to_pos(seq_strand.get_sheet_id(),
                                          helix,
                                          relpos,
                                          pos_strand, cur_strand,
                                          sheet_pos_dict)
        return (xpos, ypos)

        
    def build_helices_svg(self, sheet_pos_dict,
                          sse_label_scheme,
                          seqnum_dict, sheet_posmap,
                          enable_sheet_gap_rule,
                          use_helix_clustering,
                          helix_cluster_shading_color,
                          color_scheme,
                          interdomain_connectors,
                          label_residue_numbers=False):
        """
        Build SVG for helices and n- and c- terminus nodes.
        Called by build_dunnart_svg()

        
        Parameters:
             sheet_pos_dict - 
                              dictionary mapping from sheet id to top and
                               bottom of sheet y co-ordinate and left
                               and right x co-orindate:
                               { sheet_id : (top_ypos, bottom_ypos,
                                             left_xpos,right_xpos) }
             sse_label_scheme - if 'sequential' number all nodes in one sequence
                   instead of sepearte sequences for strands and helices.
                   Note that this does not affect internal number, it is
                   just for the labels on Dunnart shapes.
             seqnum_dict - dictionary { nodeid : seqnum } only used if
                           sse_label_scheme is 'sequential'
             sheet_posmap - The PTPosMap of sheet neighbour relationships
             enable_sheet_gap_rule - If True, don't align helices on strand
                            indguides between neihbouring sheets, if the
                            total helix length is above a threshold.
             use_helix_clustering - If True and using heuristic
                        helix placement,
                        cluster sequential helices and place in cluster
                        with tableau and distance matrix rather than
                        aligning them all on strand axis.
            helix_cluster_shading_color - color to shade helix clusters
            color_scheme: 'none', 'simple', 'gradient', 'sheet', 'fold'.
                        Only needed because in simple scheme need to color
                        helices in clusters specially, but don't know about
                        helix clusters until build_helices_heuristic().
            interdomain_connectirs - make connectors between domains rather
                        than using pseudo-terminus nodes as normally used
                         (such as when one domain per file).

             label_residue_numbers - if True put start and end residue ids
                            on head and tail of helix shape

        Uses data members (readonly):
             node_list - ordered list of nodes
             (NOTE: sets is_positioned flag in helices though)
             include_pi_helices, include_310_helices - flags to include or not
             (write):
               helix_cluster_id_generator (initialized here)
             
        Return value:
             None
        """
        
        # NOTE: horiz/vert and x/y and height/width, top/bot
        # names may be somewhat confusing in
        # here now since sideways strands reverse the meanings of these
        # things... (i.e. height is actually width, 'horiz_indguides' are
        # vertical, etc. variable names and comments assume up/down
        # (i.e. not sideways)

        # un-reversed normally is pointing 'up' (reversed is down)
        # un-reversed when sideways is pointing 'left' (reversed is right)
        
        # dictionary mapping strand id to current y position of other nodes
        # (helices) on the vertical axis (indguide) of that strand.
        # Use to keep track of y co-ordinate to give helices when there
        # are multiple positioned above/below the strand
        # { nodeid : (ypos_above, ypos_below, sideways,
        #   num_helices_aligned_above, num_helices_aligned_below) }
        # e.g. { 'STRAND_1' : (70, 330, False, 0, 1) }
        # Note there is ypos_above and ypos_below the strand may have
        # helices before and after (above and below, or below and above,
        # respecitvely, dpeending if the strand is reversed or not) and
        # so we need to keep track of both of them seperately.
        strand_ypos_dict = {}

        # number helix clusters (if any) starting from 1
        self.helix_cluster_id_generator = GenSeqId(1)
        
        # set up the strand_ypos_dict for use in positioning helices
        # below
        for strand in self.iter_strands():
            assert(isinstance(strand, PTNodeStrand))
            nnodeid = strand.nodeid
            sheet_id = strand.get_sheet_id()
            if strand.get_sideways():
                ypos_above = sheet_pos_dict[sheet_id][2] - get_min_gap_size()
                ypos_below = sheet_pos_dict[sheet_id][3] + get_min_gap_size()
                sideways = True
            else:
                ypos_above = sheet_pos_dict[sheet_id][0] - get_min_gap_size()
                ypos_below = sheet_pos_dict[sheet_id][1] + get_min_gap_size()
                sideways = False
            strand_ypos_dict[nnodeid] = (ypos_above, ypos_below, sideways,
                                         0, 0)
                
        unpositioned_termini = [] # list of tuples (node, nodelist)
        unpositioned_helices = [] # list of tuples (node, seq_strand)
        unpositioned_helix_clusters = [] # list of tuples ([list of helices in
                                         # cluster], seq_strand, seqpos)
        
        for nodelist in self.iter_chains():
            node_index = 0
            while node_index < len(nodelist):
                node = nodelist[node_index]
                if not isinstance(node, PTNodeHelix) and \
                   not isinstance(node, PTNodeTerminus):
                    node_index += 1
                    continue  # only helices and termini handled here
                if isinstance(node, PTNodeHelix) and node.get_is_positioned():
                    node_index += 1
                    continue # this helix has already been written

                if ( isinstance(node, PTNodeHelix) and
                     ( (node.get_type() == "310" and
                        not self.include_310_helices) or
                       (node.get_type() == "PI" and
                      not self.include_pi_helices) ) ):
                    node_index += 1
                    continue # a 310 or pi helix and we're not drawing them

                if ( interdomain_connectors and
                     isinstance(node, PTNodeTerminus) and node.get_pseudo() ):
                    node_index += 1
                    continue # not drawing pseudo-termini
                     

                node_index = \
                        self.build_helices_heuristic(
                                                sheet_pos_dict,
                                                sse_label_scheme,
                                                seqnum_dict, sheet_posmap,
                                                enable_sheet_gap_rule,
                                                strand_ypos_dict,
                                                unpositioned_termini,
                                                unpositioned_helices,
                                                unpositioned_helix_clusters,
                                                nodelist,
                                                node_index,
                                                use_helix_clustering,
                                                helix_cluster_shading_color,
                                                color_scheme,
                                                interdomain_connectors,
                                                label_residue_numbers)
            # END of while loop over nodelist
        # END of iteration over chains


        # process any helices that could not be positioned earlier, now
        # by using the distance matrix placement algorithm (for unpositioned
        # helices) or helix clustering algorithm (for unpositioned helix
        # clusters)

        # first do the unpositioned helices not in a cluster
        for (helix, seq_strand) in unpositioned_helices:
            (xpos, ypos) = self.find_helix_xypos(helix, seq_strand,
                                                 sheet_posmap,
                                                 sheet_pos_dict,
                                                 enable_sheet_gap_rule)
            self.set_helix_svginfo(helix,
                                   xpos, ypos,
                                   sse_label_scheme,
                                   seqnum_dict,
                                   label_residue_numbers)

        # now do the unpositioned helix clusters
        for (helix_cluster_list, seq_strand, seqpos) in \
                unpositioned_helix_clusters:
            clusterid = self.helix_cluster_id_generator.next()
            for helix in helix_cluster_list:
                helix.set_cluster_id(clusterid)
            self.build_helix_cluster(sse_label_scheme,
                                     seqnum_dict,
                                     strand_ypos_dict,
                                     seq_strand,
                                     seqpos,
                                     helix_cluster_list,
                                     clusterid,
                                     helix_cluster_shading_color,
                                     sheet_posmap,
                                     sheet_pos_dict,
                                     label_residue_numbers)

            
        # and now process any termini nodes that we could not position earlier
        for (termnode, nodelist) in unpositioned_termini:
            if ( interdomain_connectors and termnode.get_pseudo() ):
                continue # not drawing pseudo-termini
            (node_xpos, node_ypos) =  self.find_terminus_pos(termnode,
                                                             nodelist)
            self.set_terminus_svginfo(termnode,
                                      node_xpos, node_ypos)




    def build_helices_heuristic(self,
                                sheet_pos_dict,
                                sse_label_scheme,
                                seqnum_dict, sheet_posmap,
                                enable_sheet_gap_rule,
                                strand_ypos_dict,
                                unpositioned_termini,
                                unpositioned_helices,
                                unpositioned_helix_clusters,
                                nodelist,
                                node_index,
                                use_helix_clustering,
                                helix_cluster_shading_color,
                                color_scheme,
                                interdomain_connectors,
                                label_residue_numbers = False):
        
        """
        Build SVG for helices and n- and c- terminus nodes for a single chain.
        Called by build_helices_svg()

        
        Parameters:
             sheet_pos_dict - 
                              dictionary mapping from sheet id to top and
                               bottom of sheet y co-ordinate and left
                               and right x co-orindate:
                               { sheet_id : (top_ypos, bottom_ypos,
                                             left_xpos,right_xpos) }
             sse_label_scheme - if 'sequential' number all nodes in one sequence
                   instead of sepearte sequences for strands and helices.
                   Note that this does not affect internal number, it is
                   just for the labels on Dunnart shapes.
             seqnum_dict - dictionary { nodeid : seqnum } only used if
                            sse_label_scheme is'sequential'
             sheet_posmap - The PTPosMap of sheet neighbour relationships
             enable_sheet_gap_rule - If True, don't align helices on strand
                            indguides between neihbouring sheets, if the
                            total helix length is above a threshold.
             strand_pos_dict - (IN/OUT)
                    The strand_ypos_dict from build_helices_svg():
                    { nodeid : (ypos_above, ypos_below, sideways,
                    num_helices_aligned_above, num_helices_aligned_below) }
                    e.g. { 'STRAND_1' : (70, 330, False, 0, 1) }
             unpositioned_termini - (IN/OUT) list of tuples
                                    (node, nodelist)
             unpositioned_helices - (IN/OUT) list of tuples (node, seq_strand)
             unpositioned_helix_cluster - (IN/OUT) list of ([list of helices in
                                          helix cluster], seq_strand, seqpos)
             nodelist - list of nodes in chain to process
             node_index - index of the current node to process in the nodelist
             use_helix_clustering - If True and using heuristic helix placement,
                        cluster sequential helices and place in cluster
                        with tableau and distance matrix rather than
                        aligning them all on strand axis.
             helix_cluster_shading_color - color to shade helix clusters
             color_scheme: 'none', 'simple', 'gradient', 'sheet', 'fold'.
                        Only needed because in simple scheme need to color
                        helices in clusters specially, but don't know about
                        helix clusters until build_helices_heuristic().
             interdomain_connectors - If True, do NOT make pseduo-terminus nodes
                     at domain boundaries. Instead the domain boundary
                       SSEs are left ot hve connectors to other domain
                       added later. Default False.

             label_residue_numbers - if True put start and end residue ids
                            on head and tail of helix shape


        Uses data members (readonly):
             node_list - ordered list of nodes
             (NOTE: sets is_positioned flag, cluster_id in helices though)
             include_pi_helices, include_310_helices - flags to include or not
             (read/write):
             helix_cluster_dict - dict { clusterid : PTSVGCluster)
             helix_cluster_id_generator
             xmlid_generator
             
        Return value:
          node_index where
          node_index is the index of the next node in the nodelist to be
          processed (index for while loop over nodelist in build_helices_svg())

        """
        
        # if total length of helices on an axis between two sheets is
        # greater than this then do not position the helices there,
        # instead place them later using distance matrix method
        SHEET_GAP_HELIX_LEN = 10 # 10 means ANY helix is long enough


        # when using helix 'clustering' this is the minimum number
        # helices required to make a cluster rather than just align
        # on indguide.
        MIN_HELIX_CLUSTER_SIZE = 2

        node = nodelist[node_index]
        
        # position nodes near (above or below) the sheet they
        # 'belong to' nodes after a normal (not reversed, so
        # 'pointing up') strand go below sheet if before
        # strand in seq, and above sheet if after strand in
        # seq. If strand is reversed, so is this, ie.  below
        # sheet if after in seq, above if before.
        (sheet_id, nterm_strand_index, cterm_strand_index) = \
             self.immediate_containing_sheet(node_index,
                                             node.get_chainid())

        nterm_strand = None
        cterm_strand = None
        if nterm_strand_index != None or cterm_strand_index != None:
            seqpos = None
            if cterm_strand_index == None:
                nterm_strand = nodelist[nterm_strand_index]
                seq_strand = nterm_strand
                seqpos = SEQPOS_AFTER
            elif nterm_strand_index == None:
                cterm_strand = nodelist[cterm_strand_index]
                seq_strand = cterm_strand
                seqpos = SEQPOS_BEFORE
            else:
                # we have a choice between the axis on which the
                # n-terminal or c-terminal strand is aligned.
                # Choose one with no helices already aligned on it,
                # if possible
                # but first check (even if enable_sheet_gap_rule not True)
                # if it would be a 'sheet gap' helix in one position and
                # not the other, and choose the one in which it isn't if
                # possible.
                nterm_strand = nodelist[nterm_strand_index]
                cterm_strand = nodelist[cterm_strand_index]
                (num_helices_above_nterm_strand,
                 num_helices_below_nterm_strand) = \
                 self.count_helices_on_strand_axis(nterm_strand,
                                                   strand_ypos_dict)
                (num_helices_above_cterm_strand,
                 num_helices_below_cterm_strand) = \
                 self.count_helices_on_strand_axis(cterm_strand,
                                                   strand_ypos_dict)


                seqpossg = SEQPOS_AFTER
                if nterm_strand.get_reversed():
                    if seqpossg == SEQPOS_AFTER:
                        seqpossg = SEQPOS_BEFORE
                    else:
                        seqpossg = SEQPOS_AFTER
                nterm_is_sheet_gap = self.helix_is_sheet_gap(node, nterm_strand,
                                                             cterm_strand,
                                                             seqpossg,
                                                             nterm_strand,
                                                             sheet_posmap)
                seqpossg = SEQPOS_BEFORE
                if cterm_strand.get_reversed():
                    if seqpossg == SEQPOS_AFTER:
                        seqpossg = SEQPOS_BEFORE
                    else:
                        seqpossg = SEQPOS_AFTER
                cterm_is_sheet_gap = self.helix_is_sheet_gap(node, nterm_strand,
                                                             cterm_strand,
                                                             seqpossg,
                                                             cterm_strand,
                                                             sheet_posmap)
#                print 'xxx',str(node),'nterm',nterm_is_sheet_gap,'cterm',cterm_is_sheet_gap
                if (nterm_is_sheet_gap and not cterm_is_sheet_gap):
                    seq_strand = cterm_strand
                    seqpos = SEQPOS_BEFORE
                elif (cterm_is_sheet_gap and not nterm_is_sheet_gap):
                    seq_strand = nterm_strand
                    seqpos = SEQPOS_AFTER
                elif ( (not nterm_strand.get_reversed() and
                        num_helices_above_nterm_strand == 0) or
                       (nterm_strand.get_reversed() and
                        num_helices_below_nterm_strand == 0) ):
                    seq_strand = nterm_strand
                    seqpos = SEQPOS_AFTER
                elif ( (not cterm_strand.get_reversed() and
                        num_helices_below_cterm_strand == 0) or
                       (cterm_strand.get_reversed and
                        num_helices_above_cterm_strand == 0) ):
                    seq_strand = cterm_strand
                    seqpos = SEQPOS_BEFORE
                else:
                    #use closest strand in sequence (n-terminal if tied)
                    nterm_dist = node_index - nterm_strand_index
                    cterm_dist = cterm_strand_index - node_index
                    assert(nterm_dist > 0)
                    assert(cterm_dist > 0)
                    if cterm_dist < nterm_dist:
                        seq_strand = nodelist[cterm_strand_index]
                        seqpos = SEQPOS_BEFORE
                    else:
                        seq_strand = nodelist[nterm_strand_index]
                        seqpos = SEQPOS_AFTER

            if seq_strand.get_reversed():
                if seqpos == SEQPOS_AFTER:
                    seqpos = SEQPOS_BEFORE
                else:
                    seqpos = SEQPOS_AFTER

            # for multiple helices adjacent in sequence, process them
            # all in the same way in an inner loop here.
            # We do this by making a new list of only these adjacent
            # helices (maybe with terminus), so that we can
            # reverse it if necessary.
            helix_list = []
            while not isinstance(node, PTNodeStrand): # always at least once
                if ( isinstance(node, PTNodeHelix) and
                     ( (node.get_type() == "310" and
                        not self.include_310_helices) or
                       (node.get_type() == "PI" and
                        not self.include_pi_helices) ) ):
                    node_index += 1
                    if node_index >= len(nodelist):
                        break
                    node = nodelist[node_index]
                    continue # 310 or pi helix & we're not drawing them

                if not node.get_is_positioned():
                    helix_list.append(node)
                    node_index += 1
                    if node_index >= len(nodelist):
                        break
                    node = nodelist[node_index]

            is_helix_cluster = False
            if use_helix_clustering:
                helix_cluster_list = [ helix for helix in helix_list
                                       if isinstance(helix, PTNodeHelix) ]
                if len(helix_cluster_list) > MIN_HELIX_CLUSTER_SIZE:
                    is_helix_cluster = True

            # if the list of helices would be between two sheets
            # and is long enough that we don't want it to get in
            # the way there, then put them on unpositioned list
            # to be positioned by distance matrix method later instead.
            total_helix_len = sum([helix.get_span() for helix in
                                   helix_list if
                                   isinstance(helix, PTNodeHelix)]) * \
                                         DUNNART_HELIX_SPAN_FACTOR
            if (enable_sheet_gap_rule and
                total_helix_len > SHEET_GAP_HELIX_LEN and
                self.helix_is_sheet_gap(helix_list[0], nterm_strand,
                                        cterm_strand, seqpos,
                                        seq_strand, sheet_posmap,
                                        is_helix_cluster)):
                if verbose:
                    sys.stderr.write('intersheet helices: ' +
                             str([str(helix) for helix in helix_list]) +
                             '\n')
                if is_helix_cluster:
                    unpositioned_helix_clusters.append((helix_cluster_list,
                                                        seq_strand,
                                                        None)) # no seqpos
                                                        
                    for tnode in [tn for tn in helix_list if
                                  isinstance(tn, PTNodeTerminus)]:
                        unpositioned_termini.append((tnode, nodelist))
                else:
                    for node in helix_list:
                        if isinstance(node, PTNodeHelix):
                            unpositioned_helices.append((node, seq_strand))
                            # FIXME: should change the name of the
                            #        is_interstrand flag now it is used here
                            #        also.
                            node.set_is_interstrand(True)
                        elif isinstance(node, PTNodeTerminus):
                            unpositioned_termini.append((node, nodelist))
                        else:
                            assert(False)
                # node_index already incremented in loop building
                # helix_list, so don't need to increment it to continue
                return node_index # we do not write this one now

            # reverse the list if it is below an upwards (non-reversed)
            # strand that is after (c-terminal to) it
            # or vice versa (reversed strand n-terminal to it)
            if seqpos == SEQPOS_BEFORE and not seq_strand.get_reversed() or\
               seqpos == SEQPOS_AFTER and seq_strand.get_reversed():
                helix_list.reverse()

            if ( isinstance(helix_list[-1], PTNodeHelix) and
                 (self.pdb_resid_dict[(seq_strand.get_chainid(),
                                       seq_strand.get_start_res_seq())] >
                  self.pdb_resid_dict[(helix_list[-1].get_chainid(),
                                       helix_list[-1].get_start_res_seq())]) ):
                # aligned on C-terminal strand, set last in list flag
                helix_list[-1].set_is_last_in_seq(True)
                helix_list[-1].set_reversed(seq_strand.get_reversed())


            if is_helix_cluster:
                # just put helix clusters on a list to process after all
                # aligned helices are processed, since we need the aligned
                # helices there to see if helix clusters would overlap
                # with them.
                unpositioned_helix_clusters.append((helix_cluster_list,
                                                    seq_strand,
                                                    seqpos))
                                                        
                for tnode in [tn for tn in helix_list if
                              isinstance(tn, PTNodeTerminus)]:
                    unpositioned_termini.append((tnode, nodelist))
            else:
                for node in helix_list:
                    if isinstance(node, PTNodeHelix):
                        height = node.get_span() * DUNNART_HELIX_SPAN_FACTOR
                        width = DUNNART_HELIX_WIDTH
                        node.set_reversed(seq_strand.get_reversed())
                    elif isinstance(node, PTNodeTerminus):
                        width = DUNNART_TERMINUS_SIZE
                        height = DUNNART_TERMINUS_SIZE
                    else:
                        assert(False)

                    (ypos_above, ypos_below,sideways,
                     num_helices_aligned_above,
                     num_helices_aligned_below)= \
                                 strand_ypos_dict[seq_strand.nodeid]
                    if sideways:
                        temp = height
                        height = width
                        width = temp
                        helix_extent = width
                        node.set_sideways(True)
                        node_ypos = seq_strand.ypos
                    else:
                        node_xpos = seq_strand.xpos
                        helix_extent = height

                    offset = get_min_gap_size()
                    if seqpos == SEQPOS_AFTER:
                        if sideways:
                            node_xpos = ypos_above - offset - helix_extent
                            new_above_pos = node_xpos
                        else:
                            node_ypos = ypos_above - offset - helix_extent
                            new_above_pos = node_ypos
                        strand_ypos_dict[seq_strand.nodeid] = \
                                    (new_above_pos, ypos_below, sideways,
                                     num_helices_aligned_above + 1,
                                     num_helices_aligned_below)
                    else:
                        if sideways:
                            node_xpos = ypos_below + offset
                            new_below_pos = node_xpos + helix_extent
                        else:
                            node_ypos = ypos_below + offset
                            new_below_pos = node_ypos + helix_extent
                        strand_ypos_dict[seq_strand.nodeid] = \
                                     (ypos_above, new_below_pos, sideways,
                                      num_helices_aligned_above,
                                      num_helices_aligned_below + 1)

                    if isinstance(node, PTNodeHelix):
                        self.set_helix_svginfo(node,
                                               node_xpos, node_ypos,
                                               sse_label_scheme,
                                               seqnum_dict,
                                               label_residue_numbers)
                    elif isinstance(node, PTNodeTerminus):
                        if ( interdomain_connectors and
                             isinstance(node, PTNodeTerminus) and
                             node.get_pseudo() ):
                            node_index += 1
                            continue # not drawing pseudo-termini
                        else:
                            self.set_terminus_svginfo(node,
                                                      node_xpos, node_ypos)
                    else:
                        assert(False)

                    # Also for helices (and termini)
                    # that are just before or after strand in sequence,
                    # put them on the vert indguide for that strand
                    if sideways:
                        align_type = DUNNART_ALIGN_MIDDLE
                    else:
                        align_type = DUNNART_ALIGN_CENTER
                    self.svg_constraint_list.append(
                        PTSVGAlignmentConstraint(seq_strand.indguide,
                                                 node,
                                                 align_type))

                # END of inner (for) loop over adjacent helices in helix_list
            # END of else case (not using helix clustering)
        else: # of nterm_strand_index != None or cterm_strand_index != None:
            if isinstance(node, PTNodeHelix):
                # no strand in this chain so we can't position helix
                # on a strand indguide - we will use the distance
                # matrix placement method instead:
                # find the sheet closest to this helix
                # (there must be some sheet already placed (in another
                # chain(s)) otherwise the -i option (distance matrix
                # method) would have been forced on and we wouldn't
                # be in this subroutine.
                (closest_sheetid, unused_dist) = \
                  self.distmatrix.get_min_distance_objid(node.nodeid,
                                                         Set(),
                                                         sheets_only=True)
                if verbose:
                    sys.stderr.write('positioning ' + str(node) +
                                     ' relative to ' +
                                     'sheet ' + closest_sheetid +
                                     '\n')

                # position the element near its closest already
                # positioned one according to position map for
                # relative position and tableau for orientation
                (relpos, pos_strand, cur_strand) = \
                     self.ptrelpos.get_relative_position(closest_sheetid,
                                                         node)
                (node_xpos, node_ypos) = \
                            self.relpos_to_pos(closest_sheetid,
                                               node,
                                               relpos,
                                               pos_strand, cur_strand,
                                               sheet_pos_dict)
                self.set_helix_svginfo(node, 
                                       node_xpos, node_ypos,
                                       sse_label_scheme,
                                       seqnum_dict,
                                       label_residue_numbers)
                node_index += 1 # needed to continue loop in calling routine
            elif isinstance(node, PTNodeTerminus):
                # add these to a list to process later, since we
                # want to place all helices first to ensure we
                # always have a real SSE to place terminus near to
                unpositioned_termini.append((node, nodelist))
                node_index += 1 # needed to conintue loop
                return node_index # we do not write this one now
            else:
                assert(False)

        return node_index


    def build_helix_cluster(self,
                            sse_label_scheme,
                            seqnum_dict,
                            strand_ypos_dict,
                            seq_strand,
                            seqpos,
                            helix_list,
                            helix_cluster_id,
                            helix_cluster_shading_color,
                            sheet_posmap = None,
                            sheet_pos_dict = None,
                            label_residue_numbers = False):
        """
        Build a 'cluster' of helices that are in sequence between two
        strands. Called by build_helices_heuristic().
        
        Parameters:
             sse_label_scheme - if 'sequential' number all nodes in one sequence
                   instead of sepearte sequences for strands and helices.
                   Note that this does not affect internal number, it is
                   just for the labels on Dunnart shapes.
             seqnum_dict - dictionary { nodeid : seqnum } only used if
                           sse_label_scheme is 'sequential'
             strand_ypos_dict - (IN/OUT)
                    The strand_ypos_dict from build_helices_svg():
                    { nodeid : (ypos_above, ypos_below, sideways,
                    num_helices_aligned_above, num_helices_aligned_below) }
                    e.g. { 'STRAND_1' : (70, 330, False, 0, 1) }
             seq_strand - The PTNodeStrand that first helix is list
                    is closest in sequence to, used to align that helix.
             seqpos - SEQPOS_AFTER or _BEFORE as calculted
                       in build_helices_svg()
                       or None if we cannot place cluster aligned on strand
                       (due to helix sheet gap rule), in which case
                       we place it using distance matrix placement.
             helix_list - list of helices in chain to process
             helix_cluster_id - integer id for this helix cluster
             helix_cluster_shading_color -color to shade the helix cluster

            sheet_posmap - The PTPosMap of sheet neighbour relationships
                           Only used if seqpos==None. Default None.
            sheet_pos_dict - 
                         dictionary mapping from sheet id to top and
                               bottom of sheet y co-ordinate and left
                               and right x co-orindate:
                               { sheet_id : (top_ypos, bottom_ypos,
                                             left_xpos,right_xpos) }
                           Only used if seqpos==None. Default None.
             label_residue_numbers - if True put start and end residue ids
                            on head and tail of helix shape

        Uses data members (readonly):
             node_list - ordered list of nodes
             (NOTE: sets is_positioned, is_interstrand flag in helices though)
          (write):
              helix_cluster_list - list of SVGCluster
              helix_cluster_dict - maps clutser id to SVGCluster object
              svg_constraint_list

             
        Return value:
            None
        """
        assert(seqpos in [None, SEQPOS_AFTER, SEQPOS_BEFORE])
        assert(len(helix_list) > 0)

        if not self.room_for_helix_cluster(seq_strand, seqpos,
                                           strand_ypos_dict):
            if verbose:
                sys.stderr.write('no room for helix cluster ' + str([str(h) for h in helix_list]) + '\n')
            seqpos = None # force distance matrix placement not strand alignment
        
        node = helix_list[0]
        height = node.get_span() * DUNNART_HELIX_SPAN_FACTOR
        width = DUNNART_HELIX_WIDTH
        (ypos_above, ypos_below,sideways,
         num_helices_aligned_above,
         num_helices_aligned_below)= \
                     strand_ypos_dict[seq_strand.nodeid]
        if sideways:
            temp = height
            height = width
            width = temp
            helix_extent = width
            node.set_sideways(True)
            node_ypos = seq_strand.ypos
        else:
            node_xpos = seq_strand.xpos
            helix_extent = height

        offset = get_min_gap_size()
        if seqpos == None:
            (node_xpos,node_ypos) = self.find_helix_xypos(node, seq_strand,
                                                          sheet_posmap,
                                                          sheet_pos_dict,
                                                          True) # sheet gap
        elif seqpos == SEQPOS_AFTER:
            if sideways:
                node_xpos = ypos_above - offset - helix_extent
                new_above_pos = node_xpos
            else:
                node_ypos = ypos_above - offset - helix_extent
                new_above_pos = node_ypos
            strand_ypos_dict[seq_strand.nodeid] = \
                        (new_above_pos, ypos_below, sideways,
                         num_helices_aligned_above + 1,
                         num_helices_aligned_below)
        else: # SEQPOS_BEFORE
            if sideways:
                node_xpos = ypos_below + offset
                new_below_pos = node_xpos + helix_extent
            else:
                node_ypos = ypos_below + offset
                new_below_pos = node_ypos + helix_extent
            strand_ypos_dict[seq_strand.nodeid] = \
                         (ypos_above, new_below_pos, sideways,
                          num_helices_aligned_above,
                          num_helices_aligned_below + 1)

        node.set_is_interstrand(True) # FIXME: should change name of this flag
        self.set_helix_svginfo(node, node_xpos, node_ypos,
                               sse_label_scheme,
                               seqnum_dict,label_residue_numbers)
        # Put the first helix on the list, ie the one that
        # that is just before or after strand in sequence,
        # on the vert indguide for that strand
        # Unless seqpos=None indicating helix cannot aligned on strand
        # due to helix sheet gap rule
        if seqpos != None:
            if sideways:
                align_type = DUNNART_ALIGN_MIDDLE
            else:
                align_type = DUNNART_ALIGN_CENTER
            self.svg_constraint_list.append(
                PTSVGAlignmentConstraint(seq_strand.indguide,
                                         node,
                                         align_type))

        # Now position the rest of the helices using the distance matrix
        # and tableau

        # dict { PTNodeHelix : PTNodeStrand } mapping, for a positinoed helix,
        # the strand used as as the reference to position it on its axis
        helix_refstrand_dict = {}
        
        positioned_elements = Set([node])
        unpositioned_elements = Set(helix_list[1:])
#        print 'xxx',str([str(node) for node in unpositioned_elements])
        firsthelix = True
        while (len(unpositioned_elements) > 0):
            # find nearest element to any already positioned element in cluster
            (positioned_element, cur_element, dist_unused) = \
                  self.distmatrix.find_nearest_sses_in_sets( \
                                                         positioned_elements,
                                                         unpositioned_elements)
            if verbose:
                sys.stderr.write('build_helix_cluster positioning ' +
                                 str(cur_element) +
                                 ' relative to ' +
                                 str(positioned_element) + '\n')

            # position the element near its closest already positioned one
            # according to position map for relative position and tableau
            # for orientation
            if firsthelix:
                ref_strand = seq_strand
                helix_refstrand_dict[positioned_element] = seq_strand
                firsthelix = False
            else:
                ref_strand = helix_refstrand_dict[positioned_element]
            (relpos, test_strand) = \
                     self.ptrelpos.get_helixcluster_relative_position(
                                                         positioned_element,
                                                         cur_element,
                                                         ref_strand)
            helix_refstrand_dict[cur_element] = test_strand
            assert(isinstance(cur_element, PTNodeHelix))
            cur_element.set_is_interstrand(True) # FIXME: change name of flag
            (xpos, ypos) = self.relpos_to_pos(positioned_element,
                                              cur_element,
                                              relpos,
                                              None, None,
                                              sheet_pos_dict=None)

            self.set_helix_svginfo(cur_element,
                                   xpos, ypos,
                                   sse_label_scheme,
                                   seqnum_dict,label_residue_numbers)
            positioned_elements.add(cur_element)
            unpositioned_elements.remove(cur_element)

        shading_color = get_color_list(helix_cluster_shading_color)[0] +\
                        DUNNART_CLUSTER_ALPHA_HEX
        helixcluster = PTSVGCluster(helix_list,
                                    self.xmlid_generator.next(),
                                    str(helix_cluster_id),
                                    0,  # color_num used for proximity grouping
                                    shading_color)
        self.helix_cluster_list.append(helixcluster)
        self.helix_cluster_dict[helix_cluster_id] = helixcluster
        if verbose:
            sys.stderr.write('build_helix_cluster: helices ' +
                             str([str(helix) for helix in helix_list]) + 
                             'are in cluster ' + str(helix_cluster_id) + '\n') 
            


    def build_connectors_aligned_svg(self, nodelist,
                                     chain_i,
                                     use_connector_arrowheads=False,
                                     connector_color_scheme = 'all',
                                     interdomain_connectors = False):
        """
        Called by build_dunnart_svg() to 
        build SVG for connectors for sequence in a single chain:
        a connector between each
        node (helix/strand) in sequence order.

        This is the original version for use with herustic helix
        placement (ie no -i option) along with build_helices_svg(). It puts
        connectors on ports of helices simply to try to avoid
        unnecessary crossings assuming they are aligned in order along
        strand axes, ignoring any tableau orientation they are
        supposed to have.
        
        
        Parameters:
             nodelist - list of nodes in this chain
             chain_i - chain index (0,1,...) for selecting line color
             use_connector_arrowheads - If True make connectors directed.

             connector_color_scheme -  'all[:<color>]', 'chain[:<color_list>]',
                                       'domain[:<intra_color>,<inter_color>]',
                                       'crossing[:<color_list>]' 
             
             interdomain_connectors - If True, do NOT make
                       pseduo-terminus nodes
                       at domain boundaries. Instead the domain boundary
                       SSEs are left ot hve connectors to other domain
                       added later. Default False.
                       

        

        Uses data members (readonly):
             node_list - ordered list of nodes (NB sets port fields in SVGNodes)
             include_pi_helices, include_310_helices

        Return value:
             the new current XML id

        Precondition: nodelist is sorted (by start res seq ascending);
                      this is done by build_graph_from_secstruct()
                      before calling.


        """
        prevnode = None
        prevdstFlags = None
        for nodeindex in range(len(nodelist)):
            node = nodelist[nodeindex]

            if isinstance(node, PTNodeHelix):
                if ( (node.get_type() == "310" and
                      not self.include_310_helices) or
                     (node.get_type() == "PI" and
                      not self.include_pi_helices) ):
                    continue # skip 310/pi helix if flagged to not draw

            if ( interdomain_connectors and
                 isinstance(node, PTNodeTerminus) and node.get_pseudo() ):
                continue # don't do pseudo for interdomain
            
            if prevnode != None:
                # add connector from prevnode to node
                if isinstance(node, PTNodeStrand):
                    dst_reversed = node.get_reversed()
                else:
                    dst_reversed = False
                if dst_reversed:
                    dstFlags = DUNNART_TOP_PORT
                else:
                    dstFlags = DUNNART_BOTTOM_PORT

                if prevdstFlags == None:
                    if isinstance(prevnode, PTNodeStrand):
                        src_reversed = prevnode.get_reversed()
                    else:
                        src_reversed = False
                    if src_reversed:
                        srcFlags = DUNNART_BOTTOM_PORT
                    else:
                        srcFlags = DUNNART_TOP_PORT
                else:
                    if prevdstFlags == DUNNART_TOP_PORT or \
                           prevdstFlags == DUNNART_LEFT_PORT:
                        srcFlags = DUNNART_BOTTOM_PORT
                    else:
                        srcFlags = DUNNART_TOP_PORT
                    
                if isinstance(node, PTNodeHelix):
                    if node.get_sideways():
                        prevnode_ypos = prevnode.xpos               # x coord
                        node_ypos = node.xpos                       # x coord
                    else:
                        prevnode_ypos = prevnode.ypos               # y coord
                        node_ypos = node.ypos                       # y coord
                    if node.get_is_interstrand():
                        # special case for helices between strands on same
                        # axis in sheet
                        if node.get_reversed(): # N-term strand reversed value
                            dstFlags = DUNNART_TOP_PORT
                        else:
                            dstFlags = DUNNART_BOTTOM_PORT
                    elif node.get_is_last_in_seq():
                        if node.get_reversed(): # axis strand reversed value
                            dstFlags = DUNNART_TOP_PORT
                        else:
                            dstFlags = DUNNART_BOTTOM_PORT
                    else:
                        if prevnode_ypos < node_ypos:
                            dstFlags = DUNNART_TOP_PORT
                            prevdstFlags = dstFlags
                            if not isinstance(node, PTNodeStrand):
                                srcFlags = DUNNART_BOTTOM_PORT

                if node.get_sideways():
                    if dstFlags == DUNNART_TOP_PORT:
                        dstFlags = DUNNART_LEFT_PORT
                    else:
                        dstFlags = DUNNART_RIGHT_PORT
                if prevnode.get_sideways():
                    if srcFlags == DUNNART_TOP_PORT:
                        srcFlags = DUNNART_LEFT_PORT
                    else:
                        srcFlags = DUNNART_RIGHT_PORT
                    
                if isinstance(node, PTNodeTerminus):
                    dstFlags = DUNNART_DEFAULT_PORT
                elif isinstance(prevnode, PTNodeTerminus):
                    srcFlags = DUNNART_DEFAULT_PORT

                # line color may be overwritten later for multidomain
                linecolor = get_line_color(connector_color_scheme, chain_i)

                node.nterm_port = dstFlags
                prevnode.cterm_port = srcFlags
                
                xmlid = self.xmlid_generator.next()
                self.svg_connector_list.append(
                    PTSVGConnector(xmlid, prevnode, node, srcFlags, dstFlags,
                                   linecolor, use_connector_arrowheads))
                prevdstFlags = dstFlags
            prevnode = node


    def set_helix_cluster_colors(self, helix_proximity_shading_colors,
                                 helix_cluster_shading_color):
        """
        Set shading colors for helix clusters.
        Color nearby helix clusters 
        the same shading color. Note that this may involve
        making some helices a one-helix 'cluster'.

        Called by build_dunanrt_svg().

        Parameters:
           helix_proximity_shading_colors - shade nearby helix clusters the same
                                            color: 'auto' (use color gradient
                                            to shade each differently),
                                            or list of colors.

           helix_cluster_shading_color - default color to shade helix clusters

        Uses data members (read/write):
                helix_cluster_id_generator
                helix_cluster_list
                helix_cluster_dict
                helix_cluster_dict - dict { clusterid : PTSVGCluster)
                xmlid_generator


        Return value:
             None
        
        """
        assert(helix_proximity_shading_colors != None)
        
        num_clusters = len(self.helix_cluster_list)
        if num_clusters < 1:
            return

        cluster_fill_colors = get_cluster_fill_colors(
                                              helix_proximity_shading_colors,
                                              num_clusters)

        if helix_proximity_shading_colors != 'auto':
            # put standard (-k) cluster fill color at 0 for nonassigned to gourp
            default_shading_color = \
                        get_color_list(helix_cluster_shading_color)[0] +\
                        DUNNART_CLUSTER_ALPHA_HEX
            cluster_fill_colors.insert(0, default_shading_color)

        color_num = 1 # NB start at 1 not 0; 0 used to mark non-assigned cluster
            
        self.helix_cluster_list[0].color_num = color_num
        self.helix_cluster_list[0].color = cluster_fill_colors[color_num]
        
        # for every existing helix cluster, find the nearest helix to
        # any helix in the cluster. If it is below the distance threshold,
        # then make the cluster that helix is in the same color. If it is
        # not in a cluster then make it a cluster and color it the same color
        # FIXME: this is being done very inefficiently, constructing lots
        # of sets from list of all helices and using set difference, should
        # build one set then add/subtract things from it.
        HELIX_CLUSTER_NEARBY_THRESHOLD = 7.0 # Angstroms. FIXME: adaptive?
        cluster_list = list(self.helix_cluster_list) # NB a copy not reference

        # What about transitivity? This way cluster 2
        # might be near cluster 1, so same color, then cluster
        # 3 near 2, also same color, but might not be near cluster
        # 1 at all, but ends up same color as it (via 2).
        # We resolve this by insisting all clusters are
        # pairwise 'close'.
        # TODO: FIXME: this is not so good either, as the result depends
        # on order things are processed in. E.g (7API) have cluster of helices
        # B,C,D then look at E it is close, then A, but it is not close to E
        # so not in group.
        # But for 1SNG, A  gets processed before E so get A in group with B,C,D
        # and not E.
        # In both cases, A and E are close to B,C,D cluster but not to
        # each other. What should be done? 
        # XXXX 15Feb2008 - greedy algorithm of
        # sorting nearby clusters by distance and using
        # closest first now implemented, not sure it is what we want though
        # e.g 1SNG now getting helix H close not E, 7API gets A and E.
        
        # We'll call the set of clusters that are the same color a 'group'
        group_dict = {} # dict of { color_num : list of PTSVGCluster }
        grouped_clusters = Set() # Set of PTSVGCluster that are in some group
        for helixcluster in cluster_list:
            if helixcluster in grouped_clusters:
                continue # this cluster already in a group, skip it
            cluster_helix_set = Set(helixcluster.svgnodelist)
            other_helix_set = Set(list(self.iter_helices())) \
                                  - cluster_helix_set
            close_list = self.distmatrix.find_nearby_sses_in_sets(
                cluster_helix_set, other_helix_set,
                HELIX_CLUSTER_NEARBY_THRESHOLD) # returns (dist, helix) tuples
            close_list.sort() # sort by distance (ascending)
            if verbose:
                sys.stderr.write('set_helix_cluster_colors: helices ' +
                           str([str(helix) for (dist, helix) in close_list]) +
                           ' are close to helix cluster ' +
                           helixcluster.clusterid + ' (' +
                           str([str(helix) for helix in cluster_helix_set]) +
                           ')\n')
            for (dist, nearby_helix) in close_list:
                if nearby_helix.get_cluster_id() != None:
                    nearby_cluster = \
                        self.helix_cluster_dict[nearby_helix.get_cluster_id()]
                    if verbose:
                        sys.stderr.write('  nearby helix ' +
                                         str(nearby_helix) + 
                                         ' is in cluster ' +
                                         nearby_cluster.clusterid + '\n')
                else:
                    # not in a cluster, make one just for this helix
                    helix_cluster_id = self.helix_cluster_id_generator.next()
                    nearby_cluster = PTSVGCluster([nearby_helix],
                                                  self.xmlid_generator.next(),
                                                  str(helix_cluster_id))
                    nearby_helix.set_cluster_id(helix_cluster_id)
                    self.helix_cluster_list.append(nearby_cluster)
                    self.helix_cluster_dict[helix_cluster_id] = nearby_cluster
                    if verbose:
                        sys.stderr.write('  clusters ' +
                                         helixcluster.clusterid +
                                         ' and single helix ' +
                                         str(nearby_helix) +
                                         ' (now cluster ' +
                                         str(helix_cluster_id ) +
                                         ') ' +
                                         ' are close\n')

                assert(nearby_cluster != helixcluster)
                if ( not group_dict.has_key(color_num) or
                     len(group_dict[color_num]) < 2 ):
                    all_close = True
                else:
                    all_close = True
                    test_helix_set = Set(nearby_cluster.svgnodelist)
                    for ref_cluster in group_dict[color_num]:
                        ref_helix_set = Set(ref_cluster.svgnodelist)
                        (ref_elt, test_elt, dist) = \
                           self.distmatrix.find_nearest_sses_in_sets(
                              ref_helix_set, test_helix_set)
                        if dist >= HELIX_CLUSTER_NEARBY_THRESHOLD:
                            all_close = False
                            if verbose:
                                sys.stderr.write('  cluster ' +
                                                 nearby_cluster.clusterid +
                                                 ' is not close to cluster ' +
                                                 str(ref_elt.get_cluster_id()) +
                                                 '; not adding to group\n')
                            break
                if all_close:
                    if nearby_cluster.color_num == 0:
                        thiscolornum = helixcluster.color_num
                        nearby_cluster.color_num = thiscolornum
                        nearby_cluster.color = cluster_fill_colors[thiscolornum]
                    else:
                        thiscolornum = nearby_cluster.color_num
                    if not group_dict.has_key(thiscolornum):
                        group_dict[thiscolornum]=[helixcluster, nearby_cluster]
                        grouped_clusters.add(nearby_cluster)
                        grouped_clusters.add(helixcluster)
                    else:
                        group_dict[thiscolornum].append(nearby_cluster)
                        grouped_clusters.add(nearby_cluster)
                    if verbose:
                        sys.stderr.write('  clusters ' +
                                         helixcluster.clusterid +
                                         ' and ' +
                                         nearby_cluster.clusterid +
                                         ' are close, colored with ' +
                                         'color number ' + str(thiscolornum) +
                                         '\n')
            color_num += 1

        # now remove all single-helix clusters that are not in any group
        # They are all added by this subroutine so at end of cluster_list
        if verbose:
            sys.stderr.write(str([str(c) for c in self.helix_cluster_list])
                             + '\n')
        i = len(self.helix_cluster_list) - 1
        while i > 0 and len(self.helix_cluster_list[i].svgnodelist) == 1:
            if (self.helix_cluster_list[i].color_num == 0):
                removed_cluster = self.helix_cluster_list.pop(i)
                if verbose:
                    sys.stderr.write('  helix ' +
                                 str(removed_cluster.svgnodelist[0]) +
                                 ' (cluster ' +
                                 removed_cluster.clusterid + ')' +
                                 ' not in any group so no longer a cluster\n')
            i -= 1
                
            
                    

    ##########################################################################

    #
    # functions that operate on PTSVGNodes xpos,ypos etc. having
    # already been build to build_dunnart_svg()
    #
  

    def get_bounding_box(self):
        """
        Get a bounding box for the cartoon based on all the SVG nodes
        generated already by build_dunnart_svg().

        Parameters:
           None

        Return value:
           Tuple of tuples ((x1, y1), (x2, y2)) where (x1,y1) is the top
           left corner of the bounding box and (x2,y2) is the bottom right
           corner of the bounding box.
        """
        x1 = y1 = sys.maxint
        x2 = y2 = -sys.maxint - 1
        for node in self.iter_nodes():
            if node.get_is_positioned():
                if node.xpos - node.width < x1:
                    x1 = node.xpos - node.width
                if node.xpos + node.width > x2:
                    x2 = node.xpos + node.width
                if node.ypos - node.height < y1:
                    y1 = node.ypos - node.height
                if node.ypos + node.height > y2:
                    y2 = node.ypos + node.height
        return ((x1, y1), (x2,y2))
                 

    def scale(self, scale_factor):
        """
        Scale the cartoon by uniform scaling. Can be used as primitive
        way of ensuring no overlaps.

        Parameters:
           scale_factor - factor to multiply all coordinates by

        Return value:
           None

        Updates xpos,ypos in all SVGNodes in the nodelist.
        """
        for node in self.iter_nodes():
            if node.get_is_positioned():
                node.xpos *= scale_factor
                node.ypos *= scale_factor


    def translate_relative_bbox(self, bounding_box, relpos):
        """
        Translate the cartoon to a position ABOVE, BELOW, LEFT or RIGHT of
        the supplied bounding_box (specified by tuple of
        top-left and bottom-right
        (x,y) tuples). Leaves an additional gap of DUNNART_DOMAIN_GAP_SIZE.

        Parameters:
           bounding_box - Tuple of tuples ((x1, y1), (x2, y2)) where
                           (x1,y1) is the top
                          left corner of the bounding box and (x2,y2) is
                          the bottom right corner of the bounding box.
           relpos - ptrelpos.py RELPOS_ABOVE/BELOW/LEFT/RIGHT to position
                    relative to the bounding box,
                    or None to move it to the same position as the bounding box.
        Return value:
           None.

        Updates xpos,ypos in all SVGNodes int the nodelist.
        Note must also update x/y for indguides and distribution handles

        Raises exceptions:
           ValueError for bad relpos value
        """
        (x1, y1) = bounding_box[0]
        (x2, y2) = bounding_box[1]
        this_bounding_box = self.get_bounding_box()
        (this_x1, this_y1) = this_bounding_box[0]
        (this_x2, this_y2) = this_bounding_box[1]
        xshift = 0
        yshift = 0  
        if relpos == None:
            xshift = -(this_x1 - x1)
            yshift = -(this_y1 - y1)
        elif relpos == RELPOS_ABOVE:
            yshift = -(abs(this_y2 - y1) + DUNNART_DOMAIN_GAP_SIZE)
        elif relpos == RELPOS_BELOW:
            yshift = abs(y2 - this_y1) + DUNNART_DOMAIN_GAP_SIZE
        elif relpos == RELPOS_LEFT:
            xshift = -(abs(this_x2 - x1) + DUNNART_DOMAIN_GAP_SIZE)
        elif relpos == RELPOS_RIGHT:
            xshift = abs(x2 - this_x1) + DUNNART_DOMAIN_GAP_SIZE
        else:
            raise ValueError('bad relpos value ' + str(relpos))

        # move shapes
        for node in self.iter_nodes():
            if node.get_is_positioned():
                if relpos == None:
                    node.ypos += yshift
                    node.xpos += xshift
                elif relpos == RELPOS_ABOVE or relpos == RELPOS_BELOW:
                    node.ypos += yshift
                elif relpos == RELPOS_LEFT or relpos == RELPOS_RIGHT:
                    node.xpos += xshift
                else:
                    raise ValueError('bad relpos value ' + str(relpos))
                

        # move indguides and handles
        for svgconstraint in self.svg_constraint_list:
            svgconstraint.translate(xshift, yshift)


    def write_dunnart_svg(self, fh):
        """
        Write the SVG to the supplied file using the SVG nodes built by
        build_dunnart_svg()

        Parameters:
            fh  - filehandle open for writing, to write SVG to.

        Return value:
            None
        """
        #
        # write the SVG to the file
        #

        # write out the svg that has been built for each node
        for node in self.iter_nodes():
            if node.get_is_positioned(): # e.g. 3_10, pi helices may not be
                node.write_svg(fh)

        # write out the cluster constraints (for sheets and helix clusters)
        for cluster in self.sheet_cluster_list + self.helix_cluster_list:
            cluster.write_svg(fh)

        # write out the alignment and distribution constraints
        for svgconstraint in self.svg_constraint_list:
            svgconstraint.write_svg(fh)

        # write out the connectors
        for connector in self.svg_connector_list:
            connector.write_svg(fh)


    ##########################################################################



#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------


def get_simple_colors(color_scheme):
    """
    Return the rgb color tuples for the  colors specified for the 'simple'
    color scheme, which is in the format

       'simple:sheet=<sheet_colors>.helixcluster=<helixcluster_color>.alpha=<helix_alpha_colors>.pi=<helix_pi_colors>.310=<helix_310_colors>.terminus=<terminus_colors>
        colors strands in sheets
        the sheet_color, helices in helix clusters
       (if any)  the helixcluster_color, helices the helix_color.
       the color lists are comma-delimited lists of recognized color names
       or RGB colors in hex (e.g. sheet=red,purple,#ddf02d).

    This format is assumed to have already been validated by the parameter
    parsing in main(), by regular expression (but we check for duplicate
    type names here as that is simpler than in regexp - also for valid
    colors; this functino is called as part of command line validation).
    
    Parameters:
        color_scheme: string starting 'simple:' as specfied above
    Return value:
        dict { type : color_list } where type is 'sheet','helixcluster',
        'alpha','pi' or '310' and color_list is list of color value strings
        where the color value strings are 'rrggbb' color
        hex strings (note no '#')
    Raises Excpetions:
        KeyError for unknown color names
        ValueError for dupplicate type names 
    """
    type_color_dict = {}
    if color_scheme[:7] != "simple:":
        return None # should be validated before entering this function
    type_eq_value_list = color_scheme[7:].split('.')
    for type_eq_value in type_eq_value_list:
        type_value = type_eq_value.split('=')
        typestr = type_value[0]
        if type_color_dict.has_key(typestr):
            raise ValueError("duplicate type " + typestr)
        color_list_str = type_value[1]
        color_list = get_color_list(color_list_str)
        type_color_dict[typestr] = color_list
    if verbose:
        sys.stderr.write('"simple" color scheme: ' + str(type_color_dict) + '\n')
    return type_color_dict


def get_connector_colors(connector_color_scheme):
    """
    Return the RGB color tuples for the colors specified in the
    connector_color_scheme, which is in one of the formats:

       'all[:<color>]' 
       'chain[:<color_list>]' 
       'domain[:<intra_color>,<inter_color>]'
       'crossing[:<color_list>]' 

    The color lists are comma-delimited lists of recognized color names
    or RGB colors in hex (e.g. red,purple,#ddf02d).

    It is assumed that this function is only called when a color or color_list
    acutally is present (i.e. string contains ':' followed by color(s)
    and format has already been validated by paramering parsing (by regexp
    in main()).

    Parameters:
       connector_color_scheme - as specified above
    Return value:
       list of color value strings (may have length 1 e.g. for 'all')
       where the color value strings are 'rrggbbaa' color hex strings (no '#')
       where the alpha channel value aa is always 'ff'
    Uses globals (Readonly):
        COLOR_DICT (color.py) - maps color names to RGB hex strings
    Raises Excpetions:
        KeyError for unknown color names
    """
    splitstr = connector_color_scheme.split(':')
    scheme = splitstr[0]
    color_list_str = splitstr[1]
    color_str_list = color_list_str.split(',')
    color_list = [ color + 'ff' for color in  get_color_list(color_list_str) ]
    return color_list


def get_line_color(connector_color_scheme, chain_i):
    """
    Utility function used to get connector color based on the user-specified
    color scheme

    Parmameters:
        connector_color_scheme -string as speicfied in get_connector_colors
        chain_i - index number of current chain
    Return value:
        color for this connector in this chain_i
    """

    # line color may be overwritten later for multidomain
    if connector_color_scheme[:5] == 'chain':
        if ':' in connector_color_scheme: # color list specified
            line_colors=get_connector_colors(connector_color_scheme)
            linecolor = line_colors[chain_i % len(line_colors)]
        else:   # use builtin default list
            linecolor = DUNNART_LINE_COLORS[chain_i % \
                                           len(DUNNART_LINE_COLORS)]
    elif connector_color_scheme[:3]=='all':
        if ':' in connector_color_scheme: # color list specified
            linecolor=get_connector_colors(connector_color_scheme)[0]
        else:
            linecolor = DUNNART_DEFAULT_LINE_COLOR
    elif connector_color_scheme[:6] == 'domain':
        if ':' in connector_color_scheme:
            linecolor = get_connector_colors(connector_color_scheme)[0]
        else:
            linecolor = DUNNART_DEFAULT_LINE_COLOR
    else: # "crossing" is handled in Dunnart not here, but use first color in
          # list for all connectors, crossings to be resolved in Dunnart.
        if ':' in connector_color_scheme: # color list specified
            linecolor=get_connector_colors(connector_color_scheme)[0]
        else:
            linecolor = DUNNART_DEFAULT_LINE_COLOR
    return linecolor



def write_dunnart_svg_prelude(fh, identifier,
                              pdbid,
                              num_domains, total_num_domains,
                              color_interfering_connectors=False,
                              connector_color_scheme = None,
                              use_auto_graph_layout=False):
    """
    Write the XML prelude information for the SVG.

    Parameters:
         fh  - filehandle open for writing, to write SVG to.
         identifier - string identifier (PDB id, etc) to put in comments
         pdbid - PDB id
         num_domains - the number of protein domains represneted by this file
         total_num_domains - number of domains identified in protein
         color_intefering_connectors - If True, set flags in SVG to tell
                                Dunnart to color crossing and shared path
                                connectors different colors.
                                Default False.
        connector_color_scheme -  'all[:<color>]', 'chain[:<color_list>]',
                                       'domain[:<intra_color>,<inter_color>]',
                                       'crossing[:<color_list>]' 
        use_auto_graph_layout - If True, set flags in SVG to tell Dunnart
                                   to use automatic graph layout.
                                   Default False.
    Return value:
         None.
    """
    fh.write('<?xml version="1.0" encoding="UTF-8"?>\n')
    timestamp = strftime("%d%b%Y %H:%M:%S", localtime())
    fh.write("<!-- " + identifier + " -->\n")
    fh.write("<!-- " + str(num_domains) + " of " + str(total_num_domains)  +
             " domain(s) -->\n")
    fh.write("<!-- Generated on " + timestamp
                 + " by ptgraph2 " + '\n'
                 + get_version() + '\n'
                 + " ".join(sys.argv) + " -->\n")
    fh.write('<!-- http://munk.csse.unimelb.edu.au/pro-origami -->\n')
    fh.write('<!-- for use with Dunnart (http://www.csse.monash.edu.au/~mwybrow/dunnart/) -->\n')
    fh.write('<svg ' +
         'xmlns:dunnart="http://www.csse.monash.edu.au/~mwybrow/dunnart.dtd" ' +
         'xmlns:' + PTGRAPH_NS +
             '="http://www.csse.unimelb.edu.au/~astivala/proorigami.dtd"' +
          '>\n')
    # Also write this useful information in XML to be preserved by Dunnart
    # so it stays in final version
    fh.write('<' + PTGRAPH_NS + ':identification ' + 'pdbId="' + pdbid + '" ' +
             'outputFilename="' + identifier + '" ' +
             'numberOfDomains="' + str(num_domains) + '" ' +
             'totalDomains="' + str(total_num_domains) + '" ' +
             'creationTime="' + timestamp + '" ' +
             'program="' + 'ptgraph2' + '" ' +
             'version="' + get_version() + '" ' +
             'commandLine="' + " ".join(sys.argv) + '" ' +
             '/>\n')

    connector_colors_str = ""
    if color_interfering_connectors:
        color_conn_str = "1"
        if (connector_color_scheme and connector_color_scheme[:9]=="crossing:"):
            connector_colors = get_color_list(connector_color_scheme[9:])
            connector_colors_str = ' interferingConnectorColours="' + \
                           reduce(lambda a,b : a + ',' + b, connector_colors) +\
                           '" '
    else:
        color_conn_str = "0"
    if use_auto_graph_layout:
        # FIXME - using this is probably not a good idea, should probably remove
        fh.write('  <dunnart:options automaticGraphLayout="1" nonOverlapConstraints="1" pageBoundaryConstraints="1" penaliseCrossings="1" avoidBuffer="8" colourInterferingConnectors="' + color_conn_str + '"' +
                 connector_colors_str + '/>\n')
    else:
        fh.write('  <dunnart:options automaticGraphLayout="0" nonOverlapConstraints="1" penaliseCrossings="1" avoidBuffer="10" routingBuffer="4" colourInterferingConnectors="' + color_conn_str + '"' + connector_colors_str + '/>\n')



def write_dunnart_svg_conclusion(fh):
    """
    Write the XML prelude information for the SVG.

    Parameters:
         fh  - filehandle open for writing, to write SVG to.

    Return value:
         None.
    """
    fh.write('</svg>\n')


def find_largest_domain(ptg_list):
    """
    Return the 'largest' domain in the list of PTGraph2 objects (each one
    representing one domain). 'Largest' is defined as the one with the
    largest sheet (as defined by PTGraph2.largest_sheet()). If there is
    no PTGraph2 with a sheet, use the one with the largest helix.

    There are some alternative definitions of 'largest' that also make
    sense, such as the one with the most SSEs, or with the most residues,
    but we choose this one as on the diagram (and in 3D), the largest sheet is
    what we tend to see as 'central' and orient everything around, which
    is what we want here. Note if we choose 'most SSEs' as 'largest', may
    end up with a domain with many small helices as largets, which would
    probably not be desirable, and 'most residues' may give a domain with
    large helices as largest, while we probably would prefer a sheet
    to be the 'largest' element, so that's why it is this way.

    Parameters:
       ptg_list - list of PTGraph2 objects

    Return value:
       The PTGraph2 object that is largest according to above definition.
    """
    max_sheet_size = 0
    largest_ptg = None
    for ptg in ptg_list:
        largest_sheet_id = ptg.largest_sheet()
        if largest_sheet_id != None:
            sheet_size = ptg.sheet_size(largest_sheet_id)
            if sheet_size > max_sheet_size:
                max_sheet_size = sheet_size
                largest_ptg = ptg
    if largest_ptg == None: # no largest sheet, use helices instead
        max_helix_size = 0
        for ptg in ptg_list:
            largest_helix = ptg.largest_helix()
            if largest_helix.get_span() > max_helix_size:
                max_helix_size = largest_helix.get_span()
                largest_ptg = ptg
    return largest_ptg


def domain_domain_distance(ptg1, ptg2, pdb_struct, domain_distance_dict):
    """
    Return the distance between two domains, which will be defined as
    the distance between their two closest SSEs
    (using SSE distnace defined in ptdistmatrix.py)
    
    Parameters:
        ptg1 - PTGraph2 object for one domain
        ptg2 - PTGraph2 object for the other domain
        pdb_struct - parsed PDB structure from Bio.PDB
       domain_distance_dict (In/Out) - dict { (dom1, dom2) : ret_tuple }
                     for memoizing domiain-domain distances. (dom1,dom2)
                     is tuple of two PTGraph2 objects, note both (dom1,dom2)
                     and (dom2,dom1) are always added
                     and ret_tuple is the return value tuple as defined below.

    Return value:
        tuple (dist, closest_sse1, closest_sse2, closest_res1, closest_res2)
        distance in Angstroms between the two domains, as defined above and
        closest_sse1, closest_sse2 are PTNode objects for the closest
        SSEs in ptg1 and ptg2 domains respectively and
        closest_res1 and closest_res2 are the closest residues in
        closest_sse1 and closest_sse2 respectively.
    """
    # This function is memoized by the domain_distance_dict parmeter,
    # to save recomputations of distances that are previously computed.
    if domain_distance_dict.has_key((ptg1, ptg2)):
        return domain_distance_dict[(ptg1, ptg2)]
    
    min_dist = float("inf")
    closest_sse1 = closest_sse2 = None
    closest_res1 = closest_res2 = None
    # exclude the terminus nodes
    ptg1_sses = [ node for node in ptg1.iter_nodes()
                  if not isinstance(node, PTNodeTerminus) ]
    ptg2_sses = [ node for node in ptg2.iter_nodes()
                  if not isinstance(node, PTNodeTerminus) ]
    for sse1 in ptg1_sses:
        for sse2 in ptg2_sses:
            (dist, res1, res2) = calc_sse_sse_dist(sse1, sse2, pdb_struct)
            if dist < min_dist:
                min_dist = dist
                closest_sse1 = sse1
                closest_sse2 = sse2
                closest_res1 = res1
                closest_res2 = res2
    ret_tuple12 = (min_dist,closest_sse1,closest_sse2,closest_res1,closest_res2)
    ret_tuple21 = (min_dist,closest_sse2,closest_sse1,closest_res2,closest_res1)
    domain_distance_dict[(ptg1, ptg2)] = ret_tuple12
    domain_distance_dict[(ptg2, ptg1)] = ret_tuple21
#     if verbose:
#         sys.stderr.write('dist between domain ' + ptg1.domainid + ' and ' +
#                          ptg2.domainid + ' is ' + str(min_dist) + '\n')
    return ret_tuple12


def domain_domain_orientation(ptg1, ptg2, pdb_struct):
    """
    Return the orientation (tableau code - see pttableau.py) between two
    domains. This is defined to be the orientation between the longest
    strands in the largest sheets of the two domains (or longest helix
    if no sheet).

    Parameters:
        ptg1 - PTGraph2 object for one domain
        ptg2 - PTGraph2 object for the other domain
        pdb_struct - parsed PDB structure from Bio.PDB

    Return value:
        tuple (tabcode, sse1, sse2) where tabcode
        two-character tableau code (see pttableau.py) of orientatino between
         the two domains and sse1 and sse2 are SSEs used for orientatino
         in ptg1 and pgt2 respectively
    """
    sse1 = ptg1.get_orientation_sse()
    sse2 = ptg2.get_orientation_sse()
    angle = sse1.relative_angle(sse2, pdb_struct)
    tabcode = pttableau.angle_to_tabcode(angle)
    if verbose:
        sys.stderr.write('orientation domain ' + ptg1.domainid + ',' +
                         ptg2.domainid + '; ' +
                         '  ' + str(sse1) + ',' + str(sse2) + ': ' +
                         tabcode + '\n')
    return (tabcode, sse1, sse2)

    
def find_nearest_domain(domain_set1,  domain_set2, pdb_struct,
                        domain_distance_dict):
    """
    From the PTGraph2 objects in domain_set2, find the one that is
    nearest to any of the PTGraph2 objects in domain_set1.

    Parameters:
       domain_set1 - set of PTGraph2 objects to find nearest to
       domain_set2 - set of PTGraph2 objects to find the nearest to any in
                     domain_set1
       pdb_struct - parsed PDB structure from Bio.PDB
       domain_distance_dict (In/Out) - dict { (dom1, dom2) : dd_tuple }
                     for memoizing domiain-domain distances. (dom1,dom2)
                     is tuple of two PTGraph2 objects, note both (dom1,dom2)
                     and (dom2,dom1) are always added, they are the same.
       
    Return value:
       tuple (domain1, domain2, dist_tuple)
        where domain1 is from domain_set1 and
        domain2 is from domain_set2 and the distance between the two
        is the smallest distance between any pair (a,b) with a in domain_set1
        and b inn domain_set2
        and dist_tuple is the (dist,sse1,sse2,res1,res2) tuple as defined
        by the return value of domain_domain_distance()
    """
    min_dist = float("inf")
    closest_domain1 = closest_domain2 = None
    closest_dd_tuple = None
    for domain1 in domain_set1:
        for domain2 in domain_set2:
            dd_dist_tuple = domain_domain_distance(domain1, domain2,
                                                   pdb_struct,
                                                   domain_distance_dict)
            dist = dd_dist_tuple[0]
            if (dist < min_dist):
                min_dist = dist
                closest_domain1 = domain1
                closest_domain2 = domain2
                closest_dd_tuple = dd_dist_tuple
                
    return (closest_domain1, closest_domain2, closest_dd_tuple)


def get_domain_relpos(posdom_sse, curdom_sse,
                      nearest_ref_resnum, nearest_test_resnum,
                      tabcode,
                      positioned_domain, current_domain):
    """
    Get the relative position of current_domain relative to positioned_domain
    for layout on multidomain cartoon

    Parameters:
       posdom_sse - PTNode of SSE in positioned domain for placement relative to
       curdom_sse - PTNode of SSE in current domain to place relative to posdom
       nearest_ref_resnum - residue number in reference SSE that test
                               element is closest to
       nearest_test_resnum - residue number in test SSE that is closest
                                to reference element

       tabcode - two character tableau code for orientatino of the two domains
       positioned_domain - PTGraph2 of the domain to place curdom relative to
       current_domain - PTGraph2 of the domain to place relative to posdom

    Return value:
       RELPOS_ABOVE/BELOW/etc. for placing current domain relative to
       positioned domain.
    """
    if isinstance(posdom_sse, PTNodeStrand):
        ref_element = posdom_sse.get_sheet_id()
        ref_strand = posdom_sse
    else:
        ref_element = posdom_sse
        ref_strand = None
    if isinstance(curdom_sse, PTNodeStrand):
        test_element = curdom_sse.get_sheet_id()
        test_strand = curdom_sse
    else:
        test_element = curdom_sse
        test_strand = None
    relpos = positioned_domain.ptrelpos.get_external_relpos(
        ref_element,
        test_element,
        ref_strand,
        test_strand,
        nearest_ref_resnum,
        nearest_test_resnum,
        tabcode,
        current_domain.sheet_strandlists_dict)
    return relpos
    

def write_dunnart_svg_domains(outfilehandle, outfilename,
                              ptg_list, pdb_struct,
                              sse_label_scheme,
                              use_connector_arrowheads,
                              heuristic_helix_placement,
                              sheet_shading_colors,
                              enable_sheet_gap_rule,
                              use_helix_clustering,
                              helix_cluster_shading_color,
                              connector_color_scheme,
                              color_scheme,
                              helix_proximity_shading_colors,
                              interdomain_connectors,
                              use_scaling,
                              label_residue_numbers):
    """
    Given a list of PTGraph2 objects, one per domain, write them all
    to a single SVG file so they are positioned in some sensbile way
    relative to each other on a single cartoon.
    We use a greedy algorithm similar to that used
    in positioning sheets etc. in each cartoon, i.e.

        1. Position the largest domain.
        2. Use the distance map to find closest domain to any already
           placed domain
        3. position that closest domain relative to the chosen already
           placed one, using distance/position maps and tableaux to determine
           relative position and orientation
        4. repeat from 2 until all domains placed

    Note that several fields in PTSVGNodes that are built in eacdh
    domain independently are overwritten afterwards by this function
    or things it calls for multidomains. These include the xpos and ypos
    (moving the domainds on the page), helix and strand labels
    (since these are numbered from 1 in each domain, and we may have
    to change them to be sequnetial along chains across domains),
    helix and strand colors (similarly, color gradient may now go
    across domains) and connector colors (for connectors colored by chain,
    chains now going across domains).
    
    Parameters:
       outfilehandle - open filehandle to write SVG to
       outfilename - filename (just for verbose messages)
       ptg_list - list of PTGraph2 objects, one per domain
       pdb_struct - parsed PDB structure from Bio.PDB 
       sse_label_scheme - if 'sequential' number all nodes in one sequence
                   instead of sepearte sequences for strands and helices.
                   Note that this does not affect internal number, it is
                   just for the labels on Dunnart shapes (not available to
                   be used with GraphViz either).
       use_connector_arrowheads - If True write arrowheads on connectors
                   indicating sequence direction from N- to C- terminus.
                   Only used for Dunnart.
       heuristic_helix_placement - use the original heuristic helix placement
                   instead of trying to place helices according to distance
                    matrix information.
       sheet_shading_colors - None (use default shade for all) or
                                  'auto' (use color gradient to shade each
                                  differently) or list of colors.
       enable_sheet_gap_rule - If True and using herusistic helix placement,
                       don't put 'too long' helices between sheets that are
                        neighbours.
       use_helix_clustering - If True and using heuristic helix placement,
                        cluster sequential helices and place in cluster
                        with tableau and distance matrix rather than
                        aligning them all on strand axis.
       helix_cluster_shading_color - color to shade helix clusters
       connector_color_scheme - 'all','chain','domain','crossing' (see main)
       color_scheme  - 'none', 'simple', 'gradient', 'sheet', 'fold' (see main)
       helix_proximity_shading_colors - If not None & using helix clustering,
                                          shade nearby helix clusters the same
                                          color: 'auto' (use color gradient
                                          to shade each differently),
                                          or list of colors.
       interdomain_connectors - make connectors between domains rather
                        than using pseudo-terminus nodes as normally used
                         (such as when one domain per file).
       use_scaling - if True, use scaling as primitve way to avoid overlaps
                      before dunnart processing to help avoid crashes in dunnart
      label_residue_numbers - if True put start and end residue ids
                            on head and tail of helix shape

    Return value:
       None
    """
    initial_xmlid = 1 #start XML identifiers at 1
    
    # store above,below,left,right neighbour domain of each domain, so we
    # can avoid collisions
    domain_posmap = PTPosMap()

    # dict { (dom1, dom2) : dist_tuple }
    # for memoizing domiain-domain distances. (dom1,dom2)
    # is tuple of two PTGraph2 objects, note both (dom1,dom2)
    # and (dom2,dom1) are always added, they are the same.
    # dist_tuple is defined by return value of domain_domain_distance()
    domain_distance_dict = {}

    # build svg for the largest domain as starting point
    largest_domain = find_largest_domain(ptg_list)
    current_domain = largest_domain
    build_one_dunnart_svg_domain(largest_domain, outfilename,
                                 sse_label_scheme,
                                 use_connector_arrowheads,
                                 heuristic_helix_placement,
                                 sheet_shading_colors,
                                 enable_sheet_gap_rule,
                                 use_helix_clustering,
                                 helix_cluster_shading_color,
                                 connector_color_scheme,
                                 color_scheme,
                                 helix_proximity_shading_colors,
                                 initial_xmlid,
                                 False, # main_sideways
                                 False, # main_reversed
                                 interdomain_connectors,
                                 use_scaling,
                                 label_residue_numbers)
    # ensure all XML identifiers are distinct across domains
    initial_xmlid = largest_domain.xmlid_generator.next()

    # build set of positioned and unpositioned domains. 
    positioned_domains = Set([largest_domain])
    unpositioned_domains = Set(ptg_list) - positioned_domains
    while len(unpositioned_domains) > 0:
        (positioned_domain, current_domain, dd_tuple) = \
                            find_nearest_domain(positioned_domains,
                                                unpositioned_domains,
                                                pdb_struct,
                                                domain_distance_dict)
#        (dist, posdom_sse, curdom_sse, posdom_res, curdom_res) = dd_tuple
        if verbose:
            sys.stderr.write('positioning domain ' + current_domain.domainid +
                             ' relative to domain ' + positioned_domain.domainid
                             + '\n')
#            sys.stderr.write('   ref is ' + str(posdom_sse) + ', test is ' +
#                             str(curdom_sse) + '\n')
        (tabcode, sse1, sse2) =  domain_domain_orientation(positioned_domain,
                                                           current_domain,
                                                           pdb_struct)
        # Use orientation SSEs for relative position as well
        posdom_sse = sse1
        curdom_sse = sse2
        (dist, posdom_res, curdom_res) = calc_sse_sse_dist(posdom_sse,
                                                           curdom_sse,
                                                           pdb_struct)
        (main_sideways, main_reversed) = resolve_orientation(
                                                           tabcode, sse1, sse2)
        posdom_resnum = biopdbresid_to_pdbresseq(posdom_res.get_id())
        curdom_resnum = biopdbresid_to_pdbresseq(curdom_res.get_id())
        relpos = get_domain_relpos(posdom_sse, curdom_sse,
                                   posdom_resnum, curdom_resnum,
                                   tabcode,
                                   positioned_domain, current_domain)
        if verbose:
            sys.stderr.write('positioning domain ' +
                             current_domain.domainid +
                             ' ' + ptrelpos_to_str(relpos) +
                             ' ' + positioned_domain.domainid +
                             '\n')

            
        # avoid collisoins between domains by checking in posmap and
        # positioning relative to the domain we would have collided with
        # instead
        if domain_posmap.has_key(positioned_domain):
            domain_neighbours = domain_posmap[positioned_domain]
            neighbour = domain_neighbours.get_neighbour(relpos)
            if neighbour != None:
                if verbose:
                    sys.stderr.write('  domain position: cannot place domain ' +
                                     current_domain.domainid + ' ' +
                                     ptrelpos_to_str(relpos) + ' ' +
                                     positioned_domain.domainid +
                                     ' (occupied by ' +
                                     neighbour.domainid + ')\n')
                positioned_domain = neighbour
                if verbose:
                    sys.stderr.write('   positioning domain ' +
                                     current_domain.domainid +
                                     ' relative to domain ' +
                                     positioned_domain.domainid
                                     + '\n')

                (tabcode, sse1, sse2) = \
                          domain_domain_orientation(positioned_domain,
                                                    current_domain,
                                                    pdb_struct)
                # Use orientation SSEs for relative position as well
                posdom_sse = sse1
                curdom_sse = sse2
                (dist, posdom_res, curdom_res) = calc_sse_sse_dist(posdom_sse,
                                                                   curdom_sse,
                                                                   pdb_struct)
                (main_sideways, main_reversed) = \
                                resolve_orientation(tabcode, sse1, sse2)
                posdom_resnum = biopdbresid_to_pdbresseq(posdom_res.get_id())
                curdom_resnum = biopdbresid_to_pdbresseq(curdom_res.get_id())
                relpos = get_domain_relpos(posdom_sse, curdom_sse,
                                           posdom_resnum, curdom_resnum,
                                           tabcode,
                                           positioned_domain, current_domain)
                    
                domain_neighbours = domain_posmap[positioned_domain]
                neighbour = domain_neighbours.get_neighbour(relpos)
                if neighbour != None:
                    # still a collision.
                    # FIXME: arbitrary domain positioning here
                        if (domain_neighbours.west == None):
                            relpos = RELPOS_LEFT
                        elif (domain_neighbours.east == None):
                            relpos = RELPOS_RIGHT
                        elif (domain_neighbours.north == None):
                            relpos = RELPOS_ABOVE
                        elif (domain_neighbours.south == None):
                            relpos = RELPOS_BELOW
                        else:
                            sys.stderr.write('WARNING: nowhere to'+
                                             ' place domain ' +
                                             current_domain.domainid +
                                             ' relative to domain ' +
                                             positioned_domain.domainid + 
                                             '\n')
                        if verbose:
                            sys.stderr.write('    (collision): ' + 
                                             'positioning domain ' +
                                             current_domain.domainid + ' ' +
                                             ptrelpos_to_str(relpos) +
                                             ' domain ' +
                                             positioned_domain.domainid +
                                             '\n')

        # Build SVG objects with graph and constraints for Dunnart
        build_one_dunnart_svg_domain(current_domain, outfilename,
                                 sse_label_scheme,
                                 use_connector_arrowheads,
                                 heuristic_helix_placement,
                                 sheet_shading_colors,
                                 enable_sheet_gap_rule,
                                 use_helix_clustering,
                                 helix_cluster_shading_color,
                                 connector_color_scheme,
                                 color_scheme,
                                 helix_proximity_shading_colors,
                                 initial_xmlid,
                                 main_sideways,
                                 main_reversed,
                                 interdomain_connectors,
                                 use_scaling,
                                 label_residue_numbers)
        
        # ensure all XML identifiers are distinct across domains
        initial_xmlid = current_domain.xmlid_generator.next()


        # to position each domain, we find its relative position to the
        # already positioned domain, and translate it so it is outside
        # and in the correct relative position to the bounding box of
        # the positioned one
        positioned_bbox = positioned_domain.get_bounding_box()
        # first move the current domain to same position as positioned domain
        current_domain.translate_relative_bbox(positioned_bbox, None)
        # then move it left/right/up/down relative to the positioned domain
        current_domain.translate_relative_bbox(positioned_bbox, relpos)
        domain_posmap.add_neighbour_obj(positioned_domain, current_domain,
                                        relpos)
        unpositioned_domains.remove(current_domain)
        positioned_domains.add(current_domain)
    # END while len(unpositioned_domains) > 0


    # ensure we keep unique XML ids for any SVG objects we have to add
    xmlid_generator = GenSeqId(current_domain.xmlid_generator.next())

    # convert any pseudo terminus nodes to real terminus nodes if there
    # are in fact no SSEs continuing that chain in another domain
    # (FIXME: should fix build_graph_from_secstruct() so this doesn't happen
    # in the first place but that is quite difficult as it currently stands)
    interdomain_connector_list = []
    if interdomain_connectors:
        interdomain_connector_list = fixup_terminus_nodes(
                                        list(positioned_domains),
                                         xmlid_generator,
                                         use_connector_arrowheads)

    # relabel helices and strands so not restarting in each domain
    # Note: also redoes coloring for color gradient color scheme
    keep_domain_numbering = False # TODO: make this an option
    if not keep_domain_numbering:
        relabel_nodes_multidomain(list(positioned_domains),
                                  sse_label_scheme,
                                  color_scheme)
    redo_sseseqnum_nodes_multidomain(list(positioned_domains))

    # build the connectors between domains
    if interdomain_connectors and len(positioned_domains) > 1:
        interdomain_connector_list += build_interdomain_connectors(ptg_list,
                                           xmlid_generator,
                                           use_connector_arrowheads,
                                           connector_color_scheme)

        for connector in interdomain_connector_list:
            connector.build_resname_sequence(ptg_list[0].residue_list,
                                             ptg_list[0].pdb_resid_dict)
    # For the 'chain' connector color scheme, redo the color of all connectors
    # so that each chain has a different connector color
    if connector_color_scheme[:5] == 'chain':
        recolor_connectors_chain(ptg_list, interdomain_connector_list,
                                 connector_color_scheme)
            
    # Now write out all the SVG (order doesn't matter)
    for ptg in positioned_domains:
        ptg.write_dunnart_svg(outfilehandle)
    if interdomain_connector_list:
        for conn in interdomain_connector_list:
            conn.write_svg(outfilehandle)
    

def fixup_terminus_nodes(ptg_list, xmlid_generator,
                         use_connector_arrowheads):
    """
    convert any pseudo terminus nodes to real terminus nodes if there
    are in fact no SSEs continuing that chain in another domain
    (FIXME: should fix build_graph_from_secstruct() so this doesn't happen
    in the first place but that is quite difficult as it currently stands)

    Parameters:
       ptg_list - list of PTGraph2 objects one per domain
       xmlid_generator (IN/OUT) - GenSeqId object for generating XML ids
                                   must be initizlied to start at the
                                   next unused XML id.
       use_connector_arrowheads - If True write arrowheads on connectors
                   indicating sequence direction from N- to C- terminus.
                   Only used for Dunnart.

    Return value:
       list of PTSVGConnector objects for connected added to terminus nodes
    """
    conn_list = []
    pseudoterm_list = [] # list of (node, ptg, nodelist) pseudo terminus nodes
                         # and the PTGraph2 to which it belongs, to convert
                         # to real terminus node
    for i in range(len(ptg_list)):
        ptg1 = ptg_list[i]
        for chain1 in ptg1.iter_chains():
            for termnode in [ chain1[0],  chain1[-1] ]: # N and C terminus
                if termnode.get_pseudo():
                    found_chain_other_domain = False
#                    print 'ccc checking',termnode
                    for j in range(len(ptg_list)):
                        ptg2 = ptg_list[j]
                        if ptg1 == ptg2:
                            continue
                        for chain2 in ptg2.iter_chains():
                            if (termnode.get_chainid() == chain2[0].get_chainid()
                                and ( (termnode.get_termtype() == 'C' and
                                       get_int_icode(
                                            chain2[0].get_start_res_seq())[0] >
                                        get_int_icode(
                                           termnode.get_end_res_seq())[0] - 2) or
                                      (termnode.get_termtype() == 'N' and
                                       get_int_icode(
                                              chain2[-1].get_end_res_seq())[0] <
                                       get_int_icode(
                                       termnode.get_start_res_seq())[0] + 2) ) ):
                                # NB can be +-2 as +-1 'fake' resnum
                                found_chain_other_domain = True
#                                print 'bbb found chain in domain',ptg2.domainid,'for',termnode
                                break
                    if not found_chain_other_domain:
                        pseudoterm_list.append((termnode, ptg1, chain1))

    for (pseudoterm, ptg, chain) in pseudoterm_list:
        # this chain has no continuatino in another domain, so convert
        # pseudo terminus node to a real terminus node, which involves
        # positioning it
#        print 'aaaa',pseudoterm
        (xpos, ypos) = ptg.find_terminus_pos(pseudoterm, chain)
        label = pseudoterm.get_termtype() +\
                           pseudoterm.get_chainid().lower()
        pseudoterm.set_is_positioned(True)
        pseudoterm.set_pseudo(False)
        xmlid = xmlid_generator.next()
        pseudoterm.set_svginfo(xmlid, xpos, ypos, label)

        # now add connector from last SSE in the chain to the new terminus node
        if pseudoterm.get_termtype() == 'C':
            to_node = pseudoterm
            from_node = ptg.get_most_cterm_visible_sse(chain)
#            print 'ppp',from_node
            srcFlags = from_node.get_empty_port()
            dstFlags = DUNNART_DEFAULT_PORT
        else: # 'N'
            to_node = ptg.get_most_nterm_visible_sse(chain)
            from_node = pseudoterm
#            print 'ppp2',to_node
            srcFlags = DUNNART_DEFAULT_PORT
            dstFlags = to_node.get_empty_port()
        linecolor = DUNNART_DEFAULT_LINE_COLOR
        # linecolor may be overwritten later (e.g. for line color by chain)
        xmlid = xmlid_generator.next()
        conn_list.append(PTSVGConnector(xmlid,
                                        from_node, to_node,
                                        srcFlags, dstFlags,
                                        linecolor,
                                        use_connector_arrowheads))

    return conn_list



def redo_sseseqnum_nodes_multidomain(ptg_list):
    """
    Redo the sseseqnums on the helices and strands for multidomain
    layout so that the sequential numbers are no longer restarting in
    each domain. (Same as the 'sequential' labelling scheme).

    Parameters:
       ptg_list - list of PTGraph2 objects, one for each domain

    Return value:
       None

    Raises exceptions:
       TypeError if get other than PTSVGNodeHelix or PTSVGNodeStrand
        from iter_helices() and iter_strands()
    """
    # build a list of all visible helices and strands 
    # (note: not terminus nodes) in all domains, sorted by chain id
    # and within chainid by residue sequence number (ascending)
    nodelist = []
    for ptg in ptg_list:
        nodelist += list(ptg.iter_helices()) # skips 3/5 if not selected
        nodelist += list(ptg.iter_strands())
    nodelist.sort() # depends on PTNode comparison operators 

    # just label helices and strands in same sequence from 1,
    # not restarting at chains
    seqnum = 1
    for node in nodelist:
        node.sseseqnum = str(seqnum)
        seqnum += 1


def relabel_nodes_multidomain(ptg_list, sse_label_scheme,
                              color_scheme):
    """
    Relabel the helices and strands for multidomain layout so that
    the labels are no longer restarting in each domain.

    If the color scheme is the 'gradient' color scheme, then the
    node colors are also redone so color gradient does not restart
    in each domain.

    Parameters:
       ptg_list - list of PTGraph2 objects, one for each domain
       sse_label_scheme - if 'sequential' number all nodes in one sequence
                   instead of sepearte sequences for strands and helices.
                   Note that this does not affect internal number, it is
                   just for the labels on Dunnart shapes (not available to
                   be used with GraphViz either).
       color_scheme  - 'none', 'simple', 'gradient', 'sheet', 'fold'

    Return value:
       None

    Raises exceptions:
       TypeError if get other than PTSVGNodeHelix or PTSVGNodeStrand
        from iter_helices() and iter_strands()
    """
    # build a list of all visible helices and strands 
    # (note: not terminus nodes) in all domains, sorted by chain id
    # and within chainid by residue sequence number (ascending)
    nodelist = []
    for ptg in ptg_list:
        nodelist += list(ptg.iter_helices()) # skips 3/5 if not selected
        nodelist += list(ptg.iter_strands())
    nodelist.sort() # depends on PTNode comparison operators 

    if sse_label_scheme == 'sequential':
        # just label helices and strands in same sequence from 1,
        # not restarting at chains
        seqnum = 1
        for node in nodelist:
            node.label = str(seqnum)
            if isinstance(node, PTNodeStrand) and node.get_barrel_edge():
                node.label += '*' # TODO better way to mark barrel edges
            seqnum += 1
    elif sse_label_scheme == 'separate':
        # label strands from 1, not restarting at each chain
        # label helices A,B,etc. not restarting at each chain
        # put lowercase chainid as suffix on each if multiple chains
        multiple_chains = False
        chainid_dict = {} # dict of {chainid:True} to find distinct chainids
        for ptg in ptg_list:
            if ptg.num_chains() > 1:
                multiple_chains = True
                break
            else:
                chainid = list(ptg.iter_chains())[0][0].get_chainid()
                if not chainid_dict.has_key(chainid):
                    chainid_dict[chainid] = True
        if not multiple_chains and len(chainid_dict) > 1:
            multiple_chains = True
        strand_num = 1
        helix_num = 1
        for node in nodelist:
            if isinstance(node, PTSVGNodeStrand):
                label = str(strand_num)
                if node.get_barrel_edge():
                    label += '*' # TODO better way to mark barrel edges
                strand_num += 1
            elif isinstance(node, PTSVGNodeHelix):
                # helices are labelled 'A', 'B', etc. by convention
                # FIXME: should go to AA, AB, etc. if more than 26
                label = chr(ord('A')-1 + helix_num)
                helix_num += 1
            else:
                raise TypeError('unhandled node type')
            if multiple_chains:
                label = label + str(node.get_chainid()).lower()
            node.label = label

        # if multiple chains, put chain suffix on N and C term nodes also
        if multiple_chains:
            for ptg in ptg_list:
                for chain in ptg.iter_chains():
                    assert(chain[0].get_termtype() == 'N')
                    assert(chain[-1].get_termtype() == 'C')
                    if not chain[0].get_pseudo():
                        chain[0].label = 'N' + chain[0].get_chainid().lower()
                    if not chain[-1].get_pseudo():
                        chain[-1].label = 'C' + chain[-1].get_chainid().lower()
    else: # label scheme is 'none'
        label = ''

    # Also redo the node colors, from blue to red along whole protein
    # (across domains and chains)
    # FIXME: this is not consistent with single domain mode, where
    # colr restarts for each chain - should change one or the other
    # to be consistent (probably this way, ie not restarting at each chain,
    # is best, as it is what PyMOL does?)
    if color_scheme == 'gradient':
        rgb_list = list(color_gradient(len(nodelist)))
        assert(len(rgb_list) == len(nodelist))
        for i in range(len(nodelist)):
            nodelist[i].set_color(rgb_list[i])
            i += 1
        # color each terminus node the same color as its visible neighbour
        for ptg in ptg_list:
            for chain in ptg.iter_chains():
                nterm = chain[0]
                cterm = chain[-1]
                assert(nterm.get_termtype() == 'N')
                assert(cterm.get_termtype() == 'C')
                nterm.set_color(
                    ptg.get_most_nterm_visible_sse(chain).get_color())
                cterm.set_color(
                    ptg.get_most_cterm_visible_sse(chain).get_color())

                

def build_interdomain_connectors(ptg_list,
                                 xmlid_generator,
                                 use_connector_arrowheads,
                                 connector_color_scheme):
    """
    Build connectors between domains in a multidomain cartoon.
    
    For each pair of domains, we need to go through each pair of chains
    (one from each domain), and for chains that are split between the
    domains (i.e. the same chainid is in both domains) which will have
    a pseudo-terminus that has not been drawn, then make a connector
    from the most N-terminal SSE in the one domain to the most C-terminal
    SSE in the other domain.
    
    Parameters:
        ptg_list - list of PTGraph2 objects, one per domain
        xmlid_generator (IN/OUT) - GenSeqId object for generating XML ids
                                   must be initizlied to start at the
                                   next unused XML id.
        use_connector_arrowheads - If True write arrowheads on connectors
                   indicating sequence direction from N- to C- terminus.
                   Only used for Dunnart.
        connector_color_scheme -  'all[:<color>]', 'chain[:<color_list>]',
                                  'domain[:<intra_color>,<inter_color>]',
                                  'crossing[:<color_list>]' 

    Return value:
         List of PTSVGConnector objects from interdomain connectors
    """
    # FIXME this is all rather inefficient and inelegant, should probably
    # do this better by marking things in build_graph_from_secstruct()
    # or in relabel_nodes_multidomain()
    # but ineficciency not very important due to small number of domains/chains
    
    # first build list of connected domains, which are ones that
    # have the same chainid and the pseudo-terminus nodes have the closest
    # residue sequence ids (so that if we have eg 3 domains and one chain,
    # there will only be two connections, between each pair with consecutive
    # (though note not exactly consectuive, since only SSEs have nodes)
    # residue seqwuence numbers, not between all 3 domains).
    domcon_dict = {} # dict of {chainid: (domainid1, domaind2)}
    for i in range(len(ptg_list)):
        ptg1 = ptg_list[i]
        for chain1 in ptg1.iter_chains():
            min_seqnum_diff = sys.maxint
            min_snd_domainid = None
            for j in range(i + 1, len(ptg_list)):
                ptg2 = ptg_list[j]
                for chain2 in ptg2.iter_chains():
                    if chain1[0].get_chainid() == chain2[0].get_chainid():
                        if chain1[0].get_pseudo() and chain2[-1].get_pseudo():
                            seqnum_diff = get_int_icode(
                                          chain1[0].get_start_res_seq())[0] - \
                                          get_int_icode(
                                                chain2[-1].get_end_res_seq())[0]
#                            print 'iii',chain1[0],chain2[-1]
                        elif chain1[-1].get_pseudo() and chain2[0].get_pseudo():
                            seqnum_diff = get_int_icode(
                                           chain2[0].get_start_res_seq())[0] - \
                                          get_int_icode(
                                                chain1[-1].get_end_res_seq())[0]
#                            print 'jjj',chain2[0],chain1[-1]
                        else: # cannot have both N or both C ends psuedo-termini
                            assert(False)
                        assert(seqnum_diff >= -2) # NB can be -2 as +-1 'fake' resnum
#                        print 'ttt',ptg1.domainid,ptg2.domainid,seqnum_diff
                        if seqnum_diff < min_seqnum_diff:
                            min_seqnum_diff = seqnum_diff
                            min_snd_domainid = ptg2.domainid
            if min_snd_domainid != None:
                if domcon_dict.has_key(chain1[0].get_chainid()):
                    domcon_dict[chain1[0].get_chainid()].append(
                                            (ptg1.domainid, min_snd_domainid))
                else:
                    domcon_dict[chain1[0].get_chainid()] = \
                                          [ (ptg1.domainid, min_snd_domainid) ]
#    print 'zzzzz',domcon_dict
    
    conn_list = []
    for i in range(len(ptg_list)):
        ptg1 = ptg_list[i]
        for j in range(i + 1, len(ptg_list)):
            ptg2 = ptg_list[j]
            for chain1 in ptg1.iter_chains():
                for chain2 in ptg2.iter_chains():
                    assert(chain1[0].get_termtype() == 'N')
                    assert(chain1[-1].get_termtype() == 'C')
                    assert(chain2[0].get_termtype() == 'N')
                    assert(chain2[-1].get_termtype() == 'C')
                    if (chain1[0].get_chainid() == chain2[0].get_chainid() and
                        (ptg1.domainid, ptg2.domainid) in
                        domcon_dict[chain1[0].get_chainid()]):
                        if chain1[0].get_pseudo() and chain2[-1].get_pseudo():
                            k = len(chain2)-2
                            from_node = chain2[k]
                            while k > 0 and not from_node.get_is_positioned():
                                k -= 1
                                from_node = chain2[k]
                            k = 1
                            to_node = chain2[k]
                            while (k < len(chain2) and 
                                   not to_node.get_is_positioned()):
                                k += 1
                                to_node = chain1[1]
                        elif chain1[-1].get_pseudo() and chain2[0].get_pseudo():
                            k = len(chain1) - 2
                            from_node = chain1[k]
                            while k > 0 and not from_node.get_is_positioned():
                                k -= 1
                                from_node = chain1[k]
                            k = 1
                            to_node = chain2[k]
                            while (k < len(chain2) and
                                   not to_node.get_is_positioned()):
                                k += 1
                                to_node = chain2[k]
                        else: # cannot have both N or both C ends psuedo-termini
                            assert(False)
                        srcFlags = from_node.get_empty_port()
                        dstFlags = to_node.get_empty_port()
                        if connector_color_scheme[:6] == 'domain':
                            if ':' in connector_color_scheme:
                                linecolor = get_connector_colors(
                                                     connector_color_scheme)[1]
                            else:
                                # no colors specified; make
                                # interdomain connectors the last
                                # color in list
                                linecolor = DUNNART_LINE_COLORS[-1]
                        else: # will be overwritten later for 'chain'
                            if (connector_color_scheme[:4] == 'all:' or
                                connector_color_scheme[:9] == 'crossing:'):
                                linecolor=get_connector_colors(
                                                    connector_color_scheme)[0]
                            else:
                                linecolor = DUNNART_DEFAULT_LINE_COLOR
                        xmlid = xmlid_generator.next()
                        conn_list.append(PTSVGConnector(xmlid,
                                                      from_node, to_node,
                                                      srcFlags, dstFlags,
                                                      linecolor,
                                                      use_connector_arrowheads)
                                         )
    return conn_list


                        
def recolor_connectors_chain(ptg_list, interdomain_connector_list,
                             connector_color_scheme):
    """
    Redo the connector coloring in the chase of 'chain' coloring
    so that each chain has a different color (since with multiple domains
    may have different chains in different domains, that were colored
    the same color when operating one domain at a time).

    Parameters:
        ptg_list - list of PTGraph2 objects one per domain
        interdomain_connector_list - list of PTSVGConenctor objects for
                    connectors between domains
        connector_color_scheme -  'chain[:<color_list>]'
        NOTE: modifies the color data in objects in the lists
    Return value:
        None
    """
    assert(connector_color_scheme[:5] == 'chain')

    if ':' in connector_color_scheme: # color list specified
        color_list = get_connector_colors(connector_color_scheme)
    else: # use default line colors list
        color_list = DUNNART_LINE_COLORS
    
    # build dictionary mapping each distinct chainid to index in list
    # of line colors
    chainid_colorindex_dict = {} # dict of {chainid:colorindex}
    colorindex = 0
    for ptg in ptg_list:
        for chain in ptg.iter_chains():
            chainid = chain[0].get_chainid()
            if not chainid_colorindex_dict.has_key(chainid):
                chainid_colorindex_dict[chainid] = colorindex
                colorindex = (colorindex + 1) % len(color_list)
    if verbose:
        sys.stderr.write('chain colors: ' + str(chainid_colorindex_dict) +'\n')
    # recolor each connector with its appropriate color for chain now
    # intradomain connectors
    for ptg in ptg_list:
        for conn in ptg.svg_connector_list:
            chainid = conn.src.get_chainid()
            assert(conn.dest.get_chainid() == chainid) # src,dst same chain
            conn.color = color_list[chainid_colorindex_dict[chainid]]
    # and interdomain connectors
    for conn in interdomain_connector_list:
        chainid = conn.src.get_chainid()
        assert(conn.dest.get_chainid() == chainid) # src,dst same chain
        conn.color = color_list[chainid_colorindex_dict[chainid]]


def build_one_dunnart_svg_domain(ptg, outfilename,
                             sse_label_scheme,
                             use_connector_arrowheads,
                             heuristic_helix_placement,
                             sheet_shading_colors,
                             enable_sheet_gap_rule,
                             use_helix_clustering,
                             helix_cluster_shading_color,
                             connector_color_scheme,
                             color_scheme,
                             helix_proximity_shading_colors,
                             initial_xmlid,
                             main_sideways,
                             main_reversed,
                             interdomain_connectors=False,
                             use_scaling = False,
                             label_residue_numbers=False):

    """
    Build Dunanrt SVG in the Ptgraph2 object for a single domain
    (already with constraints built, etc.)
    
    Parameters:
       ptg (read/write) - the PTGraph2 object to build the SVG in.
       outfilename - filename for verbose message text only, not opened.
       sse_label_scheme - if 'sequential' number all nodes in one sequence
                   instead of sepearte sequences for strands and helices.
                   Note that this does not affect internal number, it is
                   just for the labels on Dunnart shapes (not available to
                   be used with GraphViz either).
       use_connector_arrowheads - If True write arrowheads on connectors
                   indicating sequence direction from N- to C- terminus.
                   Only used for Dunnart.
       heuristic_helix_placement - use the original heuristic helix placement
                   instead of trying to place helices according to distance
                    matrix information.
       sheet_shading_colors - None (use default shade for all) or
                                  'auto' (use color gradient to shade each
                                  differently) or list of colors.
       enable_sheet_gap_rule - If True and using herusistic helix placement,
                       don't put 'too long' helices between sheets that are
                        neighbours.
       use_helix_clustering - If True and using heuristic helix placement,
                        cluster sequential helices and place in cluster
                        with tableau and distance matrix rather than
                        aligning them all on strand axis.
       helix_cluster_shading_color - color to shade helix clusters
       connector_color_scheme - 'none','chain','domain'
       color_scheme  - 'none', 'simple', 'gradient', 'sheet', 'fold'
       helix_proximity_shading_colors - If not None & using helix clustering,
                                          shade nearby helix clusters the same
                                          color: 'auto' (use color gradient
                                          to shade each differently),
                                          or list of colors.
       intitial_xmlid - XML identifier to start at
       main_sideways - if True, the 'main' part of the domain (largest sheet
                       or longest helix) is drawn sideways instead of vertical
       main_reversed - If True, the main part (as above) is drawn reversed
                       (down not up, or right not left when sideways)
       interdomain_connectors - If True, do NOT make pseduo-terminus nodes
                       at domain boundaries. Instead the domain boundary
                       SSEs are left ot hve connectors to other domain
                       added later. Default False.
       use_scaling - if True use uniform scaling as primitive way to
                     remove overlaps before dunnart processing.
      label_residue_numbers - if True put start and end residue ids
                            on head and tail of helix shape

    Return value:
       None.
    """
    # if this is an alpha-only domain then always use the
    # distance-matrix based helix placement algorithm as the
    # heristic ('old') algorithm requires at least one beta sheet
    domain_heuristic_helix_placement = heuristic_helix_placement
    if len(ptg.sheet_dict) == 0:
        if heuristic_helix_placement:
            sys.stderr.write('Alpha-only domain; '
                             'distance matrix helix placement used '
                             'for this domain (' + outfilename + ')\n')
            domain_heuristic_helix_placement = False

    if verbose:
        sys.stderr.write('Building SVG for ' + outfilename)
        if ptg.domainid != None:
            sys.stderr.write(' domain ' + ptg.domainid + '\n')
        else:
            sys.stderr.write('\n')
    ptg.build_dunnart_svg(sse_label_scheme,
                          use_connector_arrowheads,
                          domain_heuristic_helix_placement,
                          sheet_shading_colors,
                          enable_sheet_gap_rule,
                          use_helix_clustering,
                          helix_cluster_shading_color,
                          connector_color_scheme,
                          color_scheme,
                          helix_proximity_shading_colors,
                          initial_xmlid,
                          main_sideways,
                          main_reversed,
                          interdomain_connectors,
                          label_residue_numbers
                          )
    if use_scaling:
        ptg.scale(SCALE_FACTOR)


    
def make_graphs(pdb_filename,
                domain_program,
                secstruct_program,
                use_dot=False, use_neato=False, use_hbonds=False,
                use_dunnart_auto=False,
                sse_label_scheme = 'separate',
                use_connector_arrowheads=False,
                connector_color_scheme = 'all',
                heuristic_helix_placement=False,
                use_tableaucreator = True,
                include_310_helices = False,
                include_pi_helices = False,
                write_mfiles = False,
                use_pdb_secstruct = False,
                sheet_shading_colors = None,
                color_scheme = 'none',
                enable_sheet_gap_rule = False,
                use_helix_clustering = False,
                helix_cluster_shading_color = None,
                helix_proximity_shading_colors = None,
                multidomain_cartoon = True,
                interdomain_connectors = False,
                use_scaling = False,
                write_pmlfile = False,
                orig_pdbfilename = None,
                label_residue_numbers = False):
    """
    For the supplied filemame, read PDB format data from that file
    and create and write out a graph for the structure read.


    Paramteters:
       pdb_filename - filename of PDB file to read
       domain_program - domain parsing program ("none" or "ddomain"
                         or "cath:cath_cdf_filename") to use
       secstruct_program - secondary structure definition program
                       ('stride' or 'dssp') to use.
       use_dot    - If True use dot from GraphViz to make PostScript output
       use_neato  - If True use neato from GraphViz to make PostScript output
       use_hbonds - If True make hydrogen bond graph instead of using
                    bridge partner information to make
                    sheets from strands.
       use_dunnart_auto - If True use automatic graph layout in Dunnart
       sse_label_scheme - if 'sequential' number all nodes in one sequence
                   instead of sepearte sequences for strands and helices.
                   Note that this does not affect internal number, it is
                   just for the labels on Dunnart shapes (not available to
                   be used with GraphViz either).
       use_connector_arrowheads - If True write arrowheads on connectors
                   indicating sequence direction from N- to C- terminus.
                   Only used for Dunnart.
       connector_color_scheme - 'all','none','chain','domain','crossing'
                                     (see main)
       heuristic_helix_placement - use the original heuristic helix placement
                   instead of trying to place helices according to distance
                    matrix information.
       use_tableaucreator - if True, use TableauCreator to get angles
                            between SSEs.
       include_310_helices - if True, include 3_10 helices in the diagram
       include_pi_helices - if True, include pi helices in the diagram
       write_mfiles - if True, write MATLAB m-files to plot strands and axes.
       use_pdb_secstruct - Use HELIX and SHEET cards in PDB.
       sheet_shading_colors - None (use default shade for all) or
                                  'auto' (use color gradient to shade each
                                  differently) or list of colors.
       color_scheme  - 'none', 'simple[:sheet_color,helixcluster_color]',
                       'gradient', 'sheet', 'fold'
       enable_sheet_gap_rule - If True and using herusistic helix placement,
                       don't put 'too long' helices between sheets that are
                        neighbours.
       use_helix_clustering - If True and using heuristic helix placement,
                        cluster sequential helices and place in cluster
                        with tableau and distance matrix rather than
                        aligning them all on strand axis.
       helix_cluster_shading_color - color to shade helix clusters
       helix_proximity_shading_colors - If not None & using helix clustering,
                                          shade nearby helix clusters the same
                                          color: 'auto' (use color gradient
                                          to shade each differently),
                                          or list of colors.
       multidomain_cartoon - If True, still build PTGraph2 for each domain
                             separately, but then place all domains on the
                             same drawing rather than one per file.
       interdomain_connectors - If True and using multidomain cartoons,
                             draw connectors between domains instead of
                             the pseudo-terminus nodes use normally.
       use_scaling - If True, use uniform scaling as primitive way of
                      avoiding overlaps before dunnart processing.
       write_pmlfile - If True, write PyMOL .pml command file to load
                       structure and define SSEs according to
                       the method we used.
                       The filename of the .pml file will be PDBID.pml
                       where PDBID.pdb was the input file.
       orig_pdbfilename - input PDB filename (may be compressed; the 
                          pdb_filename is our uncompressed copy)
      label_residue_numbers - if True put start and end residue ids
                            on head and tail of helix shape

    Writes .ps or .svg files named as per description in main() below
    WARNING: overwrites these .ps or .svg files

    Note this handles PDB filenames both in the 1QLP.pdb or
    pdb1qlp.ent format, but in eithe rcase the output file is in the 1QLP.svg
    format.
    
    Return value: None
    """

    if use_pdb_secstruct:
        pdb_secstruct = ptsecstruct.read_secstruct_from_pdb_file(pdb_filename)
        if pdb_secstruct == None:
            sys.stderr.write('WARNING: error with HELIX or SHEET cards in PDB'
                             ': ' + secstruct_program +
                             ' will be used instead\n')
    else:
        pdb_secstruct = None
        
    # read secondary structure and H bond information from STRIDE or DSSP,
    # or only read hydrogen bond information from them if already have
    # secondary structure from PDB HELIX and SHEET cards.
    if secstruct_program == "stride":
        secstruct = ptsecstruct.read_secstruct_from_stride(pdb_filename,
                                                           pdb_secstruct)
    else:
        secstruct = ptsecstruct.read_secstruct_from_dssp(pdb_filename,
                                                         pdb_secstruct)

    pdb_file_basename = os.path.basename(pdb_filename)
    (name,extension) = os.path.splitext(pdb_file_basename)
    if extension.lower() == ".ent" and name[0].lower() == "d":
        # An ASTRAL/SCOP pdbstyle file such as d1apaa_.ent,
        # we will use in this case e.g. d1apaa_ as the 'pdbid'
        pdbid = name
        is_astral = True
    else: # a regular PDB file
        try:
            pdbid = secstruct.pdb_id.lstrip()
        except AttributeError:    # no PDB id, try to use filename instead
            pdbid = name
        is_astral = False

    # write PyMOL command file to load structure and show SSEs if requested
    if write_pmlfile:
        pmlfilename = pdbid + '.pml'
        sys.stdout.write('writing file ' + pmlfilename + '\n')
        pmlfile_fh = open(pmlfilename, 'w')
        secstruct.write_pymol_sse_commands(pmlfile_fh, orig_pdbfilename)
        pmlfile_fh.close()

    pdb_parser = PDBParser()
    pdb_struct = pdb_parser.get_structure(pdbid, pdb_filename) # parse PDB file

    if is_astral:
        # ASTRAL/SCOP pdbstyle file is for a single SCOP domain, so makes
        # no sense to do domain decomposition
        domain_list = [PTDomain(None, None)] # 'single domain' indicator
        if domain_program != "none":
            sys.stderr.write("WARNING: domain decomposition not run "
                             "as input file detected as an "
                             "ASTRAL/SCOP domain.\n")
    else:
        if domain_program == "none":
            domain_list = [PTDomain(None, None)] # 'single domain' indicator
        elif domain_program[0:5] == "cath:":
            cdf_filename = domain_program[5:]
            try:
                domain_list = read_domains_from_cath_cdf_file(cdf_filename, pdbid)
            except NotInCATH_Exception:
                sys.stderr.write('WARNING: PDB identifier ' + pdbid +
                                 ' not found in CDF file.')
                sys.stderr.write(' Treating as single domain.\n')
                domain_list = [PTDomain(None, None)]
        else:
            domain_list = read_domains_from_ddomain(pdb_filename,
                                                    pdb_struct[0]) # TODO: model 0
            # Sometimes DDOMAIN seems to give domain decompositions that do not
            # make sense, i.e. have domains nested one inside the other.
            # This happens for example with 2RH1. We will check for this and
            # if it happens just ignore the decomposition, making a single domain.
            domain_cd = build_domain_chaindict(domain_list)
            if not verify_domain_disjoint(domain_list, domain_cd):
                sys.stderr.write('WARNING: DDOMAIN domain decomposition is ' +
                                 'inconsistent. Treating as single domain.\n')
                domain_list = [PTDomain(None, None)]
            # NOTE: if there is only one domain, we will make it a list
            # with a single PTDomain with all data None, signifying a
            # single domain protein with no further information.  This is
            # mainly because of when there are multiple chains, in which
            # case the single domain is reported by DDOMAIN as having a
            # different chain id for start and end. If there is a single
            # domain we really don't want to do anything special, so it is
            # better to just have it as a special case where no domain
            # processing is done.
            if len(domain_list) == 1:
                domain_list = [PTDomain(None, None)]
            elif len(domain_list) == 0:
                #  This happens if DDomain crashes for example (e.g. on 1PPJ)
                sys.stderr.write("WARNING: no domain decomposition from DDOMAIN."
                                 " Treating as single domain.\n")
                domain_list = [PTDomain(None, None)]

        # output the domain decomposition in a more or less conventional format
        if len(domain_list) > 1:
            sys.stdout.write("domain decomposition: ")
            for i in range(len(domain_list)):
                sys.stdout.write(str(domain_list[i]))
                if i < len(domain_list) - 1:
                    sys.stdout.write('/')
            sys.stdout.write(' (' + domain_program + ')\n')

    # for SSEs that cross domain boundaries, move whole SSE to one of the domains
    fixup_crossdomain_sses(secstruct, domain_list)
    
    initial_xmlid = 1
    ptg_list = [] # list of PTGraph2 objects, one for each domain
    for domain in domain_list:
        ptg = PTGraph2(pdb_struct, use_hbonds,
                       include_310_helices, include_pi_helices)

        # build PTGraph2 from secondary structure and H bonds for this domain
        try:
            ptg.build_graph_from_secstruct(secstruct, domain)
        except NoSSE_Exception:
            if domain.domainid == None:
                domidext = ''
            else:
                domidext = '-' + domain.domainid
            sys.stderr.write('WARNING: No helices or strands found in domain ' +
                             pdbid + domidext +
                             ': no output written\n')
            continue

        # find sheets as connected components & label strands
        ptg.label_sheets()

        # write MATLAB m-files to plot strands in each sheet if requierd
        if write_mfiles:
            if domain.domainid == None:
                domainid = "1"
            else:
                domainid = domain.domainid
            ptg.write_sheet_mfiles(pdbid, domainid)
            ptg.write_helix_mfiles(pdbid, domainid)

        # build the PTDistMatrix distance maps for this domain
        ptg.build_dist_matrix(domain)

        # build the PTTableau tableau 
        ptg.build_tableau(pdbid, domain, use_tableaucreator)
    
        # build layout constraints 
        ptg.build_constraints()

        # label nodes with color
        ptg.set_node_colors(color_scheme)
        
        # build a graph diagram from it
        if len(domain_list) > 1:
            outfilename = pdbid + '-' + domain.domainid
        else:
            outfilename = pdbid
        if (use_dot or use_neato):
            # use GraphViz to make PostScript output file
            dg = ptgraphviz.make_graphviz_graph(ptg)
            outfilename += '.ps'
            if use_dot:
                progname = "dot"
            else:
                progname = "neato"
            sys.stdout.write('writing file ' + outfilename + '\n')
            dg.write_ps(outfilename,prog=progname) # or use ps2 not ps to make pdf later
        elif not multidomain_cartoon:
            # Build SVG objects with graph and constraints for Dunnart
            build_one_dunnart_svg_domain(ptg, outfilename,
                                     sse_label_scheme,
                                     use_connector_arrowheads,
                                     heuristic_helix_placement,
                                     sheet_shading_colors,
                                     enable_sheet_gap_rule,
                                     use_helix_clustering,
                                     helix_cluster_shading_color,
                                     connector_color_scheme,
                                     color_scheme,
                                     helix_proximity_shading_colors,
                                     initial_xmlid,
                                     False, # main_sideways
                                     False, # main_reversed
                                     False, # interdomain_connectors
                                     use_scaling,
                                     label_residue_numbers)
                
        ptg_list.append(ptg)
    # END of iteration over domain_list

    if not (use_dot or use_neato):
        # write the domains out, or build each one then write them out
        # for multidomain (need to do one at a time to use orientation relative
        # to previous one for each one)
        if connector_color_scheme[:8] == "crossing":
            color_interfering_connectors = True
        else:
            color_interfering_connectors = False
        if multidomain_cartoon:
            if len(ptg_list) > 0:
                outfilename = pdbid + '.svg'
                sys.stdout.write('writing file ' + outfilename + '\n')
                outfilehandle = open(outfilename, 'w')
                write_dunnart_svg_prelude(outfilehandle, outfilename,
                                          pdbid,
                                          len(ptg_list),
                                          len(domain_list),
                                          color_interfering_connectors,
                                          connector_color_scheme,
                                          use_dunnart_auto)
                write_dunnart_svg_domains(outfilehandle, outfilename,
                                          ptg_list, pdb_struct,
                                          sse_label_scheme,
                                          use_connector_arrowheads,
                                          heuristic_helix_placement,
                                          sheet_shading_colors,
                                          enable_sheet_gap_rule,
                                          use_helix_clustering,
                                          helix_cluster_shading_color,
                                          connector_color_scheme,
                                          color_scheme,
                                          helix_proximity_shading_colors,
                                          interdomain_connectors,
                                          use_scaling,
                                          label_residue_numbers)
                write_dunnart_svg_conclusion(outfilehandle)
                outfilehandle.close()
        else:
            for ptg in ptg_list: # one PTGraph2 object per domain
                if len(domain_list) > 1:
                    outfilename = pdbid + '-' + ptg.domainid
                else:
                    outfilename = pdbid
                outfilename += '.svg'
                sys.stdout.write('writing file ' + outfilename + '\n')
                outfilehandle = open(outfilename, 'w')
                write_dunnart_svg_prelude(outfilehandle, outfilename,
                                          pdbid,
                                          1,
                                          len(domain_list),
                                          color_interfering_connectors,
                                          connector_color_scheme,
                                          use_dunnart_auto)
                ptg.write_dunnart_svg(outfilehandle)
                write_dunnart_svg_conclusion(outfilehandle)
                outfilehandle.close()




def get_min_gap_size():
    """
    Return the value of the Dunnart minimum gap size
    """
    global DUNNART_MIN_GAP_SIZE
    return DUNNART_MIN_GAP_SIZE

def get_strand_separation():
    """
    Return the value fo the Dunnart strand separation
    """
    global DUNNART_STRAND_SEPARATION
    return DUNNART_STRAND_SEPARATION

def set_strand_gap(gapsize):
    """
    Set the value of the Dunnart strand separation and min gap size
    """
    global DUNNART_STRAND_SEPARATION
    global DUNNART_MIN_GAP_SIZE

    DUNNART_STRAND_SEPARATION = gapsize
    DUNNART_MIN_GAP_SIZE = gapsize

def usage(progname):
    """
    Print usage message and exit
    """
    sys.stderr.write("Usage: " + progname +
            " [-35acdnhrmvixgqzuwy]"
            " [ -o sse_color_scheme] [ -l connector_color_scheme ] "
            " [ -b sse_label_scheme ] "
            " [ -k <color> ]"
            " [ -g <separation> ]"
            " [ -e <color_list>|auto ] [ -f <color_list>|auto ]"
            " [-p domain_prog] [-t struct_prog] PDBfile\n")
    sys.stderr.write("  -3 include 3_10 helices in diagram\n")
    sys.stderr.write("  -5 include pi helices in diagram\n")
    sys.stderr.write("  -a use Dunnart automatic graph layout\n")
    sys.stderr.write("  -c use HELIX and SHEET cards from PDB file\n")
    sys.stderr.write("  -d use GraphViz dot instead of Dunnart SVG\n")
    sys.stderr.write("  -n use GraphViz neato instead of Dunnart SVG\n")
    sys.stderr.write("  -p domain_prog : use domain_prog to parse domains\n")
    sys.stderr.write("       supported is 'none' or 'ddomain' (default)\n")
    sys.stderr.write("       or 'cath:cdf_file_name'\n")
    sys.stderr.write("  -h graph hydrogen bonds with GraphViz\n")
    sys.stderr.write("  -b SSE labelling scheme: 'none', 'sequential', 'separate' (default)\n")
    sys.stderr.write("  -t struct_prog : use struct_prog define " \
                     "secondary structure\n")
    sys.stderr.write("       supported is 'stride' or 'dssp' (default)\n")
    sys.stderr.write("  -r compute angles internally, not with external TableauCreator\n")
    sys.stderr.write("  -m write MATLAB M-files to plot strand axes\n")
    sys.stderr.write("  -s write PyMOL .pml command file to show SSE definitions\n")
    sys.stderr.write("  -v print verbose debugging messages to stderr\n")
    sys.stderr.write("  -i use distance matrix information instead of\n"
                     "     heuristic/aesthetic algorithm for helix placement\n")
    sys.stderr.write("  -j only valid when not using -i. Don't align helices\n"
                     "     on strand axes if they would push sheets apart\n")
    sys.stderr.write("  -k <color> cluster helices, shading them all <color>\n")
    sys.stderr.write("  -e <color_list>|auto shade nearby helix clusters the same color\n")
    sys.stderr.write("  -x draw connector arrowheads\n")
    sys.stderr.write("  -f <color_list>|auto shade each sheet a different color\n")
    sys.stderr.write("  -g <separation> set the strand and minimum object separation\n")
    sys.stderr.write("  -l connector color scheme: 'all[:<color>]' (default), "
                     "'chain[:<color_list>]', "
                     "'domain[:<intra_color>,<inter_color>']"
                     ", crossing:<color_list>\n")
    sys.stderr.write("  -o SSE color scheme: 'none' (default), "
                     "'simple:sheet=<sheet_colors>.helixcluster=<helixcluster_colors>.alpha=<helix_alpha_colors>.pi=<helix_pi_colors>.310=<helix_310_colors>.terminus=<terminus_colors>', "
                     "'gradient', 'sheet', 'fold'\n")
    sys.stderr.write("  -u multidomain cartoon: place all domains in the one\n"
                     "     SVG file instead of one per file\n")
    sys.stderr.write("  -w interdomain connectors: when using multidomain\n"
                     "     cartoons, draw connectors between domains\n"
                     "     (only in conjunction with -u)\n")
    sys.stderr.write("  -q label start and end of helices and strands with\n"
                     "     first and last PDB residue id in that SSE.\n")
    sys.stderr.write("  -y use uniform scaling to try to avoid overlaps.\n"
                     "     Ugly and often does not work anyway, use only as\n"
                     "     last resort\n")
    sys.stderr.write("  -z print version information and exit\n")
    sys.exit(1)
    

#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def main():
    """
    main for ptgraph2.py

    Usage: ptgraph2 [-35acdnhrimsvxguwyz] [-p domainprog] [-t structprog]
                    [-o ssecolorscheme ] [-l connectorcolorscheme]
                    [-b sselabelsceheme ] [-f <sheet_shade_color_list>]
                    [ -k <color> ] [ -g <separation> ]
                    [-e <helixcluster_shade_color_list>]
                    PDBfile

    Output is a Dunnart file in SVG format named
    pdb-X.svg where pdb is the pdb identifier
    from the PDB file supplied and X is the domain identifier.

    -3 specifies to include 3_10 helices in the diagram. Default is only
       alpha helices.

    -5 specifies to include pi helices in the diagram. Defaul is only
       alpha helices.
       
    -a specifies to use Dunnart automatic graph layout (not compatible with
       -d or -n or -h) (Default is to use Dunnart but not auto graph layout).

    -c use the HELIX and SHEET cards from the PDB file to define secondary
       structure. DSSP or STRIDE is still needed to find hydrogen bonds,
       and will also be used if there are no HELIX or SHEET cards.

    -d specifies to use GraphViz dot instead of outputing SVG for Dunnart.
       This outputs a PostScript file.

    -n specifies to use GraphViz neato instead of outputting SVG for Dunnart.
       This outputs a PostScript file.

    -h specifies to use stride hydrogen bonds (from stride -h, unmodified
       stride can do this) instead of the default of using bridge partners
       as determined by stride (-\$ -i, requires modified stride included
       with this program) and draw an H-bond graph instead of default
       topological drawing.

    -p specifies the domain parsing program to use. Currently supported is
       "none" (no domain decomposition) or "ddomain" (default) or
       "cath:cdf_file" (CATH CDF file where cdf_file is the filename).

    -r disables the use of TableauCreator; relative angles of helices/strands
       are computed internally.

    -b specifies the SSE labelling scheme to use.
       'none' SSEs have no labels
       'sequential' SSEs are labelled 1,2,3, etc.
       'separate' (default) Strands are labelled 1,2,3,...
                  and helices are labelled A,B,C,...
       Only for Dunnart (not compatible with -d or -n or -h)

    -t specifies the secondary structure assignment program to use.
       Currently suppoed is 'dssp' and 'stride'. Default 'dssp'.

    -i specifies to use distance matrix information to place helices
       near elements they are (3d) spatially close to, rather than
       using a heurstic to place them in easy-to-read alignments on
       nearby (in sequence) strands.

    -j only valid when not using -i. Don't align helices on strand axes
       if they would exceed a length threshold 'pushing apart' sheets
       that have been placed as neighbours.

    -k <color> 
       only valid when not using -i. Draw sequential helices as a cluster
       using tableau and distance matrix rather than aligning all on strand
       axis. Shade them all <color>.

    -e <color_list> only valid when using -k.
       Color helix clusters the same shade when
       they are nearby in 3d space (although not  nearby in the diagram).
       <color_list> is as specified for -o (below). If  <color_list>
       is specified as 'auto'
       then colors are automatically selected by color gradient.
       The default is to shade them all the same  color specified by -k.

    -m writes MATLAB M-files to plot strand carbon-alpha backbone traces
       and fitted axes.

    -s writes a PyMOL .pml command file to load the structure into PyMOL
       and show cartoon with SSEs defined according to the method
       used for our cartoon (from DSSP or STRIDE or PDB) rather than
       PyMOL's internal definition.
       
    -v specifies verbose mode: debugging output is written to stderr.

    -x specifies to draw connector arrowheads indicating sequence direction
       from N- to C-terminus. Only used with Dunnart SVG.

    -f <color_list> specifies to shade each sheet a different color.
       <color_list> is as specified for -o (below). If  <color_list>
       is specified as 'auto' then
       colors are automatically selected by color gradient.
       The default is to shade them all the same (default) color.

    -g <separation> set the strand and minimum object separation value
       (defaults to 55).

    -l specifies the connector color scheme to use.
       'all[:<color>]' (default) colors all connectors same color (default black)
       'chain[:<color_list>]' colors connectors in each chain a different color
                        in order of the <color_list> or autmoatically chosen
                        if no <color_list>
       'domain[:<intra_color>,<inter_color>]' colors intra-domain connectors
                       one color (default black) and inter-domain
                       connectors another (default violet).
       'crossing[:<color_list>]' colors connectors that cross or closely
                                 follow each other
                                 different colors in order to make them easier
                                 to  distinguish. colors chosen from
                                 <color_list> in order, or automatically
                                 chosen if not specified.
        color lists are as specified for -o (below).

    -o specifies the SSE color scheme to use. 'none' is to color shapes the
       default color (default).
       'simple:sheet=<sheet_colors>.helixcluster=<helixcluster_color>.alpha=<helix_alpha_colors>.pi=<helix_pi_colors>.310=<helix_310_colors>.terminus=<terminus_colors>
        colors strands in sheets
        the sheet_color, helices in helix clusters
       (if any)  the helixcluster_color, helices the helix_color.
       the color lists are comma-delimited lists of recognized color names
       or RGB colors in hex (e.g. sheet=red,purple,#ddf02d).
       If only one color is in the list then all elements of that type are
       colored that color, otherwise the colors are used in turn.
       If they run out (more elements than colors in list, the list is
       treated as circular, ie colors reused from the start of list again).
       Not all types may be specified, those not specified will be colored
       default color. They need not be specified in order (keyword=value style
       parameters; delimited by . (no whitespace, and no ; as it is a shell
       character making it necessary to quote it).

       'gradient' colors strands and helices along a color gradient from
       N to C terminus.

       'sheet' colors the strands in each sheet a different color, leaving
        helices the default color.

       'fold' for each sheet colors strands that are connected only by
       turns (i.e. are consecutive in sequence) one color (maybe more
       than one set of such conseuctive strands in sheet, a different
       color for each such set), and other strands (i.e. not in a sequence)
       another color(s).

    -u multidomain cartoon: place all domains in the one SVG file instead of
       one domain per file.

    -w interdomain connectors: when using multidomain cartoons, draw connectors
       between domains instead of using pseudo-terminus nodes as normally
       used.

    -q residue numbers: put start and end residue numbers on start and end
       of helices and strands as labels
       
    -y use uniform scaling to try to avoid overlaps before Dunnart processing
       in order to try to reduce Dunnart crashes. Ugly, use as last resort.
       
    -z print version information and exit
    """
    global verbose
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "35acde:f:g:ijk:l:hmno:p:qrb:st:uwv?xyz")
    except getopt.GetoptError:
        usage(sys.argv[0])

    # allowed args for -p option (regexp)
    valid_domain_programs = [r"none", r"ddomain", r"cath:.*"]
    valid_domain_programs_re = [ re.compile(re_str) for re_str in
                                 valid_domain_programs ]

    valid_secstruct_programs = ["dssp", "stride"]
    valid_colorstring =  r'((#[0-9A-Fa-f]{6})|([a-zA-Z0-9]+))'
    valid_colorlist = valid_colorstring + '(,' + valid_colorstring + ')*'
    valid_type = r"((sheet)|(helixcluster)|(alpha)|(pi)|(310)|(terminus))"
    valid_typevalue = r"(" + valid_type + r"=" + valid_colorlist+ r")"
    valid_color_options = [r"none$",
                           r"simple:"+ valid_typevalue +
                                       r"(." + valid_typevalue + "){0,5}$",
                           r"gradient$",
                           r"sheet$",
                           r"fold$"]
    valid_color_options_re = [ re.compile(re_str) for re_str in
                               valid_color_options ]
    valid_connector_color_options = [r"all(:" + valid_colorstring + r")?$",
                                     r"chain(:" + valid_colorlist + r")?$",
                                     r"domain:" + valid_colorstring +
                                             r',' + valid_colorstring + '$',
                                     r"crossing(:" + valid_colorlist + r")?$"]
    valid_connector_color_options_re = [ re.compile(re_str) for re_str in
                                         valid_connector_color_options ]
    valid_sheet_shading_colors = [ r'auto$', valid_colorlist + r'$' ]
    valid_sheet_shading_colors_re = [ re.compile(re_str) for re_str in
                                      valid_sheet_shading_colors ]
    valid_helixcluster_shading_colors = [ r'auto$', valid_colorlist + r'$' ]
    valid_helixcluster_shading_colors_re = [ re.compile(re_str) for re_str in
                                             valid_helixcluster_shading_colors ]
    valid_sse_label_options = ["none", "sequential", "separate"]
    valid_k_option = valid_colorstring + '$'
    valid_k_option_re = re.compile(valid_k_option)

    
    use_dot = False
    use_neato = False
    use_hbonds = False
    verbose = False # global (python globals are only 'global' to module though)
    use_dunnart_auto = False
    sse_label_scheme = 'separate'
    use_connector_arrowheads = False
    domain_program = "ddomain"
    secstruct_program = "dssp"
    heuristic_helix_placement = True
    use_tableaucreator = True
    include_310_helices = False
    include_pi_helices = False
    write_mfiles = False
    write_pmlfile = False
    use_pdb_secstruct = False
    sheet_shading_colors = None
    color_scheme = 'none'
    enable_sheet_gap_rule = False
    use_helix_clustering = False
    helix_cluster_shading_color = None
    connector_color_scheme = 'all'
    helix_proximity_shading_colors = None
    multidomain_cartoon = False
    interdomain_connectors = False
    use_scaling = False
    label_residue_numbers = False

    for opt,arg in opts:
        if opt == "-3": # include 3_10 helices
            include_310_helices = True
        elif opt == "-5": # include pi helices
            include_pi_helices = True
        elif opt == "-a": # for Dunnart, use auto graph layout
            use_dunnart_auto = True
        elif opt == "-c": # use HELIX and SHEET cards in PDB file
            use_pdb_secstruct = True
        elif opt == "-d": # use dot
            use_dot = True
        elif opt == "-f": # shade each sheet a different color
            for valid_sheet_shading_color_re in valid_sheet_shading_colors_re:
                if valid_sheet_shading_color_re.match(arg):
                    sheet_shading_colors = arg
                    break
            if sheet_shading_colors == None:
                sys.stderr.write("valid values for -f are: " +
                                 str(valid_sheet_shading_colors) +'\n')
                usage(sys.argv[0])
            # verify that color names are known if specified
            if sheet_shading_colors != 'auto':
                try:
                    get_color_list(sheet_shading_colors)
                except KeyError:
                    sys.stderr.write("Unknown color name in -f option\n")
                    usage(sys.argv[0])
        elif opt == "-g": # set strand and minimum separation value
            set_strand_gap(int(arg))
        elif opt == "-h": # use hbonds not bridge partners
            use_hbonds = True
        elif opt == "-l": # connector color scheme
            connector_color_scheme = None
            for valid_conncolorscheme_re in valid_connector_color_options_re:
                if valid_conncolorscheme_re.match(arg):
                    connector_color_scheme = arg
                    break
            if connector_color_scheme == None:
                sys.stderr.write("valid values for -l are; " +
                                 str(valid_connector_color_options) + "\n")
                usage(sys.argv[0])
            # verify that color names are known if specified
            if ':' in connector_color_scheme:
                try:
                    get_connector_colors(connector_color_scheme)
                except KeyError:
                    sys.stderr.write("Unknown color name in -l option\n")
                    usage(sys.argv[0])
        elif opt == "-m": # write MATLAB m-files
            write_mfiles = True
        elif opt == "-n": # use neato
            use_neato = True
        elif opt == "-o": # color scheme
            color_scheme = None
            for valid_coloropt_re in valid_color_options_re:
                if valid_coloropt_re.match(arg):
                    color_scheme = arg
                    break
            if color_scheme == None:
                sys.stderr.write("valid values for -o are: " +
                                 str(valid_color_options) + "\n")
                usage(sys.argv[0])
            # verify that color names are known if specfied
            if color_scheme[:6] == "simple":
                try:
                    get_simple_colors(color_scheme)
                except KeyError:
                    sys.stderr.write("Unknown color name in -o simple option\n")
                    usage(sys.argv[0])
                except ValueError:
                    sys.stderr.write("Duplicate type name in -o simple option\n")
                    usage(sys.argv[0])
        elif opt == "-p": # domain parsing program
            domain_program = None
            for valid_domarg_re in valid_domain_programs_re:
                if valid_domarg_re.match(arg):
                    domain_program = arg
                    break
            if domain_program == None:
                sys.stderr.write("valid values for -p are: " +
                                 str(valid_domain_programs) + "\n")
                usage(sys.argv[0])
        elif opt == "-r":  # disable TableauCreator, use internal methods
            use_tableaucreator = False
        elif opt == "-b":
            if arg not in valid_sse_label_options:
                sys.stderr.write("valid options for -b are "  +
                                 str(valid_sse_label_options))
                usage(sys.argv[0])
            sse_label_scheme = arg
        elif opt == "-s": # write PyMOL .pml file
            write_pmlfile = True
        elif opt == "-t":
            if arg not in valid_secstruct_programs:
                sys.stderr.write("valid values for -t are: " +
                                 str(valid_secstruct_programs) + "\n")
                usage(sys.argv[0])
            secstruct_program = arg
        elif opt == "-i": # use distnace matrix instead of heuristic (old)
                          # helix placement algorithm
            heuristic_helix_placement = False
        elif opt == "-j": # for heuristic helix placement, not between sheets
            enable_sheet_gap_rule = True
        elif opt == "-k": # for heuristic helix placement, cluster helices
            use_helix_clustering = True
            if not valid_k_option_re.match(arg):
                sys.stderr.write("valid options for -k are: " + 
                       valid_k_option + '\n')
                usage(sys.argv[0])
            # verify that color name is known if not #RGB
            try:
                get_color_list(arg)
            except KeyError:
                sys.stderr.write("Unknown color name in -k option\n")
                usage(sys.argv[0])
            helix_cluster_shading_color = arg
        elif opt == "-e": # color nearby helix clusters with same cluster shade
            for valid_helixcluster_shading_color_re in \
                    valid_helixcluster_shading_colors_re:
                if valid_helixcluster_shading_color_re.match(arg):
                    helix_proximity_shading_colors = arg
                    break
            if helix_proximity_shading_colors == None:
                sys.stderr.write("valid values for -e are: " +
                                 str(valid_helixcluster_shading_colors) +'\n')
                usage(sys.argv[0])
            # verify that color names are known if specified
            if helix_proximity_shading_colors != 'auto':
                try:
                    get_color_list(helix_proximity_shading_colors)
                except KeyError:
                    sys.stderr.write("Unknown color name in -e option\n")
                    usage(sys.argv[0])
        elif opt == "-u": # place all domains in one cartoon SVG file
            multidomain_cartoon = True
        elif opt == "-w": # connectors between domains on multidomain cartoon
            interdomain_connectors = True
        elif opt == "-v": # verbose
            verbose = True # this module only
            ptnode_set_verbose(True) # ptnode module
            ptsecstruct.ptsecstruct_set_verbose(True) # ptsecstruct module
            ptdomain_set_verbose(True)
            ptrelpos_set_verbose(True)
            pttableau.pttableau_set_verbose(True)
        elif opt == "-q": # label residue numbers on tail and head of shapes
            label_residue_numbers = True
        elif opt == "-x":
            use_connector_arrowheads = True
        elif opt == "-y":
            use_scaling = True
        elif opt == "-z": # print version to stdout and exit
            sys.stdout.write(get_version())
            sys.stdout.write('\n')
            sys.exit(0)
        else:
            usage(sys.argv[0])

    if use_dot and use_neato:
        sys.stderr.write(
            '-d (dot) and -n (neato) options are mutually exclusive\n')
        usage(sys.argv[0])

    if use_hbonds and not (use_dot or use_neato):
        sys.stderr.write('-h (hbonds) option requires -d (dot) or -n (neato)\n')
        usage(sys.argv[0])

    if use_dunnart_auto and (use_dot or use_neato):
        sys.stderr.write('-a (Dunnart auto graph layout) cannot be used with ' +
                         'GraphViz (-d or -n)\n')
        usage(sys.argv[0])

    if sse_label_scheme == 'sequential' and (use_dot or use_neato):
        sys.stderr.write('sequential node numbering cannot be used with ' +
                         'GraphViz (-d or -n)\n')
        usage(sys.argv[0])

    if use_connector_arrowheads and (use_dot or use_neato):
        sys.stderr.write('-x (connector arrowheads) cannot be used with ' +
                         'GraphViz (-d or -n)\n')
        usage(sys.argv[0])

    if enable_sheet_gap_rule and not heuristic_helix_placement:
        sys.stderr.write("-j cannot be used with distance matrix placement (-i)\n")
        usage(sys.argv[0])

    if use_helix_clustering and not heuristic_helix_placement:
        sys.stderr.write("-k cannot be used with distance matrix placement (-i)\n")
        usage(sys.argv[0])

    if helix_proximity_shading_colors and not use_helix_clustering:
        sys.stderr.write("-e can only be used with helix clustering (-k)\n")
        usage(sys.argv[0])

    if interdomain_connectors and not multidomain_cartoon:
        sys.stderr.write("-w can only be used with multidomain cartoons (-u)\n")
        usage(sys.argv[0])

    if (connector_color_scheme[:6] == 'domain' and not interdomain_connectors):
        sys.stderr.write("WARNING: 'domain' connector color scheme only makes"
                         " sense with interdomain connectors (-w)\n")
        sys.stderr.write("         Connector color scheme reset to 'none'.\n")

        
    if len(args) != 1:
        usage(sys.argv[0])

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
        # make graph(s) from PDB file
        make_graphs(our_pdb_filename, domain_program, secstruct_program,
                    use_dot, use_neato, use_hbonds, use_dunnart_auto,
                    sse_label_scheme, use_connector_arrowheads,
                    connector_color_scheme,
                    heuristic_helix_placement, use_tableaucreator,
                    include_310_helices, include_pi_helices,
                    write_mfiles, use_pdb_secstruct, sheet_shading_colors,
                    color_scheme, enable_sheet_gap_rule,
                    use_helix_clustering, 
                    helix_cluster_shading_color,
                    helix_proximity_shading_colors,
                    multidomain_cartoon, interdomain_connectors,
                    use_scaling, write_pmlfile, pdb_filename,
                    label_residue_numbers)
    finally:
        if used_tmp_file:
            cleanup_tmpdir(TMPDIR)
    
            
if __name__ == "__main__":
    # tmpdir() annoyingly gives 'security' warning on stderr, as does
    # tmpnam(), unless we add these filterwarnings() calls. 
    warnings.filterwarnings('ignore', 'tempdir', RuntimeWarning)
    warnings.filterwarnings('ignore', 'tempnam', RuntimeWarning) 
    main()
