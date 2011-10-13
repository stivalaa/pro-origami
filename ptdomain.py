###############################################################################
#
# ptdomain.py - object to represent protein domains and functions to
#               parse into domains from external programs
#
# File:    ptdomain.py
# Author:  Alex Stivala
# Created: September 2007
#
# $Id$
#
# PTDomain is a class representing a protein domain. A domain is represented
# by a list of segments, which are contiguous subsequences of a chain.
#
# Functions are provided to parse domains using different domain parsing
# programs and return the corresponding list of PTDomain objects.
# Supported so far is:
#    . DDOMAIN (Zhou et al 2007 Protein Science 16:947-955) program
#    . CATH (CATH Domall File (CDF) 2.0) file
#
###############################################################################

import os,sys

from Bio.PDB import * # only needed for DDomain when segment spans chains
from ptutils import cleanup_tmpdir,get_int_icode
from ptsecstruct import PTSecStruct

#-----------------------------------------------------------------------------
#
# Module globals
#
#-----------------------------------------------------------------------------

verbose = False

#-----------------------------------------------------------------------------
#
# Class definitions 
#
#-----------------------------------------------------------------------------


#
# Empty classes for exceptions
#

class NotInCATH_Exception(Exception): # identifier not found in CATH CDF file
    pass

#
# Real classes
#

class PTSegment:
    """
    The PTSegment object represents a segment (contiguous subsequence of
    a chain) by the chain identifier, and start and end residue sequence
    numbers

    Note residue numbers here are integers, not PDB residue
    'numbers' which are actually strings and may have insertion
    codes, are not sequential (may have gaps, run backwards, etc.);
    unlike in PTNode, PTGraph2 etc. where we store the PDB residue numbers,
    (but use dictionary to get sequence integers when necesary).
    Using integers makes everything much
    simpler, can do calculations easily in domeval.py etc.
    For this to work strictly correctly, should use proper purely sequential
    numbering in these classes (PTSegment, PTDomain), such as that assigned
    by DSSP or STRIDE, or equivalently by using the index into list of
    residues built by Polypeptide Builder or similar from Bio.PDB.
    """
    def __init__(self, chainid, start_resnum, end_resnum):
        """
        Construct segment with supplied chain identifier and start and
        end residue sequence numbers.

        Parameters:
             chainid - PDB chain identifier (may be '-' for none)
             start_resnum - start residue sequence number
             end_resnum - end residue sequence number

        Exceptions:
           raises ValueEror if end_resnum < start_resnum

        """
        if end_resnum < start_resnum:
            raise ValueError("end residue seqnum " + str(end_resnum) + " < " +
                             "start residue seqnum " + str(start_resnum))
        self.chainid = chainid
        self.start_resnum = start_resnum
        self.end_resnum = end_resnum


    def __str__(self):
        """
        Return string representation of segment as 'chainid:start-chainid:end'
        """
        return self.chainid + ":" + str(self.start_resnum) + \
               "-" + self.chainid + ":" + str(self.end_resnum)
    

    def is_in_segment(self, res_seq_num):
        """
        Return True iff the supplied residue sequence number is in the
        interval spanned by this segment (assumed to be same chainid)

        Parameters:
           res_seq_num - PDB residue sequence number to test

        Return value:
           True if res_seq_num is in >=start_resnum and <=end_resnum else False
        """
        if res_seq_num >= self.start_resnum and \
           res_seq_num <= self.end_resnum:
            return True
        else:
            return False

    # we will define only the rich comparison operators __eq__ and __ne__
    # (not __le__, __gt__, etc.) to test for equality or non-equality
    # of segments only.

    def __eq__(self, other):
        """
        Two segments are equal if they have same chain and same start and
        end residues
        """
        if (self.chainid == other.chainid and
            self.start_resnum == other.start_resnum and
            self.end_resnum == other.end_resnum):
            return True
        else:
            return False

    def __ne__(self,other):
        """
        Two segments are '!=' (or '<>') exactly when they are not equal.
        """
        return not self.__eq__(other)



class PTDomain:
    """
    The PTDomain object represents a protein domain by a list of segments.
    Segments are contiguous subsequences of a chain, and so each is
    represented by a chain identifier, and a start and end residue sequence
    number.

    The domain consisting of domainid == None and segment_list == None
    is a special domain signifying a single-domain protein. This is used
    because we don't want to have to specify multiple segments for multiple
    chains in a single domain - a single domain protein should be treated
    just as one unit without worrying about dividing anything up.
    """
    def __init__(self, domainid, segment_list):
        """
        Create new PTDomain with supplied domain identifier and segment list.

        Parameters:
           domainid -domain identifier, string
           segment_list - list of PTSegment objects

        NOTE if both parameters are None this marks this PTDomain as the
        one used as single element in domain list to signify single-domain
        protein with no further information.
        """
        self.domainid = domainid
        self.segment_list = segment_list


    def __str__(self):
        """
        Return a representation of the domain as list of segments
        separated by ';'
        """
        if self.domainid == None and self.segment_list == None:
            return "SINGLE-DOMAIN"
        
        s = ""
        for i in range(len(self.segment_list)):
            s += str(self.segment_list[i])
            if i < len(self.segment_list) - 1:
                s += ';'
        return s

    def is_single(self):
        """
        Return True iff this is the domain with no information representing
        a single-domain protein
        """
        if self.domainid == None and self.segment_list == None:
            return True
        else:
            return False

    def is_in_domain(self, chainid, res_seq_num):
        """
        Return True iff the supplied residue specified by chainid and
        residue sequence number is in this domain.

        Parameters:
           chainid - chainid of the residue to test
           res_seq_num - PDB residue sequence number of the residue to test

        Return value:
           True if the residue is in this domain, else False.
        """
        # Note we just do a linear search, no dictionaries or anything,
        # as there is only have a maximum of maybe 5 domains, and usually
        # only 1 or 2.
        if self.is_single():
            return True # always in the special 'single domain'
        else:
            for segment in self.segment_list:
                if segment.chainid == chainid and \
                   segment.is_in_segment(res_seq_num):
                    return True
            return False

    def get_segments_for_chainid(self, chainid):
        """
        Return a list of segments of the supplied chain in this domain

        Parameters:
           chainid - id of chain to find segments of
        Return value:
           list of PTSegment objects with supplied chain id
        """
        if self.is_single():
            return []
        else:
            return [f for f in self.segment_list if f.chainid == chainid]

    def get_minmax_res_seq_in_chain(self, chainid):
        """
        Return a tuple with the lowest and highest residue sequence numbers
        in the supplied chainid in this domain.

        Parameters:
            chainid - chain id to find low and high residue sequence numbers for

        Return value:
            tuple (min_res_seq, max_res_seq) of the lowest and highest
            respectively residue sequence numbers in the supplied chainid
            in this domain.
        """
        max_res_seq = 0
        min_res_seq = sys.maxint
        for segment in self.get_segments_for_chainid(chainid):
            if segment.start_resnum < min_res_seq:
                min_res_seq = segment.start_resnum
            if segment.end_resnum > max_res_seq:
                max_res_seq = segment.end_resnum
        return (min_res_seq, max_res_seq)

    def get_chainids(self):
        """
        Return a list of all chain identifiers in this domain.

        Parameters:
           None.

        Return value:
           List of chain identifiers used by segments in this domain.
           (Each chain identifier appears only once in list).
        """
        chaindict = {}  # dict of { chainid : True } (value not used)
        if self.segment_list != None:
            for segment in self.segment_list:
                chaindict[segment.chainid] = True
        return chaindict.keys()
    

    def add_segment(self, segment):
        """
        Add a segment to this domain.

        Parameters:
           segment - PTSegement to add to the domain
        Return value: None
        Modifies member data: segment_list
        """
        if self.segment_list == None:
            self.segment_list = [segment]
        else:
            self.segment_list.append(segment)


    def remove_segment(self, segment):
        """
        Remove a segment from this domain.
        This may involve either simply removing a segment if there is one
        in the domain that corresponds exactly to the supplied segment to
        remove, otherwise the range of residues in the segment to remove
        must be deleted from some existing segment resulting in a smaller
        segment; a more complicated case can arise when the segment to
        remove spans two (or more) segments (either entirely or in part).

        Parameters:
           segment - PTSegment representing segment (continguous range of
                     residues in a chain) to remove.
           Return value: None
           Modifies member data: segment_list, and segments in the list
        """
        try:
            sindex = self.segment_list.index(segment)
            self.segment_list.pop(sindex)
        except ValueError: 
            # no segment equal to supplied one found
            # so look for a segment that entirely contains the one to remove
            found = False
            for cur_seg in self.segment_list:
                if (segment.chainid == cur_seg.chainid and
                    segment.start_resnum >= cur_seg.start_resnum and
                    segment.end_resnum <= cur_seg.end_resnum):
                    found = True
                    break
            if found:
                if (segment.start_resnum == cur_seg.start_resnum):
                    cur_seg.start_resnum = segment.end_resnum + 1
                elif (segment.end_resnum == cur_seg.end_resnum):
                    cur_seg.end_resnum = segment.start_resnum - 1
                else:
                     # need to split segment in two which we will do
                     # by shortening existing segment for first part
                     # and creating new segment for later part
                     new_seg = PTSegment(cur_seg.chainid, 
                                         segment.end_resnum + 1,
                                         cur_seg.end_resnum)
                     cur_seg.end_resnum = segment.start_resnum - 1
                     self.segment_list.append(new_seg)

            else:
                # segment is not found at all or extends ouside of a
                # segment in the list
                for cur_seg in self.segment_list:
                    if (segment.chainid == cur_seg.chainid):
                        if (segment.start_resnum >= cur_seg.start_resnum and
                            segment.start_resnum <= cur_seg.end_resnum):
                            # extends over end of cur_seg: shorten cur_seg
                            # to end at start of segment
                            cur_seg.end_resnum = segment.start_resnum - 1
                        elif (segment.end_resnum <= cur_seg.end_resnum and
                              segment.end_resnum >= cur_seg.start_resnum):
                            # extends over start of cur_seg: shorten
                            # cur_seg to start and end of segment
                            cur_seg.start_resnum = segment.end_resnum + 1
                        


#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

def read_domains_from_ddomain(pdb_filename, pdb_model, chainid=None):
    """
    Use the DDOMAIN program to parse the structure from a PDB file into
    domains and return the corresponding list of PTDomain objects.

    DDOMAIN is described in

    Zhou, Xue, Zhou 2007 'DDOMAIN: Dividing structures into domains using a
    normalized domain-domain interaction profile' Protein Science 16:947-955.

    It is available as a 64-bit linux executable and FORTRAN-77 source code
    from http://sparks.informatics.iupui.edu/Resource_files/DDOMAIN.tar.gz

    Parameters:
       pdb_filename - filename of PDB file to run DDOMAIN on
       pdb_model  - Bio.PDB model struct for this PDB entry. Note that this
                    is only needed in the case that a DDomain domain has
                    different chain identifiers for start and end and
                    is then used just to find last residue number in chain.
       chainid - (default None). If not None, only the specified chain
                 is requested.

    Return value:
       List of PTDomain objects, one for each domain.

       NOTE: if there is only one domain, we will return a list with a single
       PTDomain with all data None, signifying a single domain protein
       with no further information.
       This is mainly because of when
       there are multiple chains, in which case the single domain is reported
       by DDOMAIN as having a different chain id for start and end. If there
       is a single domain we really don't want to do anything special, so
       it is better to just have it as a special case where no domain processing
       is done.

    """
    # DDOMAIN needs the PDB file in its working directory, and it reads
    # the PDB code (e.g. 1QLP for PDB file 1QLP.pdb) from stdin
    # (optionaly with chain suffix, which we won't use)
    # Note it requires this filename format, so for format like pdb1qlp.ent
    # we need to rename the file for DDOMAIN to 1QLP.pdb

    # This is nasty, but otherwise have to modify DDOMAIN FORTRAN-77 source
    # so that's even more hassle to have to have a custom version (like we
    # did with STRIDE).
    # So we'll work in /tmp directory, make a symlink (TODO: only UNIX allows
    # this, maybe should actually copy file so works on other platforms)
    # and run DDOMAIN there.
    oldcwd = os.getcwd()
    TMPDIR = os.tempnam(None, "ptdd")
    os.mkdir(TMPDIR)
    symlink_path = None
    try:
        pdb_file_basename = os.path.basename(pdb_filename)
        (name,extension) = os.path.splitext(pdb_file_basename)
        if extension.lower() == '.pdb':   # e.g. 1QLP.pdb
            pdb_identifier = name
            pdb_file_directory = os.path.split(pdb_filename)[0]
            symlink_path = os.path.join(TMPDIR, pdb_file_basename)
            os.symlink(os.path.abspath(pdb_filename), symlink_path)
        elif extension != '.ent' or name[:3].lower() != 'pdb':
            sys.stderr.write('WARNING: unknown PDB filename format "'
                             + pdb_file_basename + '"\n')
            sys.stderr.write('  Not running DDomain\n')
            domain_list = [PTDomain(None, None)] # one-domain protein, no further info
            return domain_list
        else: # e.g. pdb1qlp.ent, make a symlink to it in format 1QLP.pdb
            pdb_identifier = name[3:7].upper()
            symlink_path = os.path.join(TMPDIR, pdb_identifier + '.pdb')
            os.symlink(os.path.abspath(pdb_filename), symlink_path)

        os.chdir(TMPDIR)
        if verbose:
            sys.stderr.write("running DDomain...")
        (ddomain_stdin, ddomain_stdout) = os.popen2("DDomain")
        if chainid != None:
            pdbchainid = pdb_identifier + chainid
        else:
            pdbchainid = pdb_identifier
        ddomain_stdin.write(pdbchainid + '\n')
        ddomain_stdin.close()
        domain_list = parse_ddomain_output(ddomain_stdout, pdb_model)
        ddomain_stdout.close()
        if verbose:
            sys.stderr.write("done\n")
    finally:
        if symlink_path:
            os.unlink(symlink_path)
        os.chdir(oldcwd)
        cleanup_tmpdir(TMPDIR)
    return domain_list
    
def parse_ddomain_output(fh, pdb_model):
    """
    Parse the output of the DDOMAIN program.

    DDOMAIN is described in

    Zhou, Xue, Zhou 2007 'DDOMAIN: Dividing structures into domains using a
    normalized domain-domain interaction profile' Protein Science 16:947-955.

    It is available as a 64-bit linux executable and FORTRAN-77 source code
    from http://sparks.informatics.iupui.edu/Resource_files/DDOMAIN.tar.gz

    Parameters:
       fh - filehandle to read DDOMAIN output from (alrady open for read)
       pdb_model  - Bio.PDB model struct for this PDB entry. Note that this
                    is only needed in the case that a DDomain domain has
                    different chain identifiers for start and end and

    Return value:
       List of PTDomain objects, one for each domain.

    """

    domain_list = []

    # Output looks like this:
    #
    # AUTHORS-trained parameters
    #   1  A     3  A   109
    #   2  B     3  B   109
    # SCOP-trained parameters
    #   1  A     3  A   109
    #   2  B     3  B   109
    # CATH-trained parameters
    #   1  A     3  A   109
    #   2  B     3  B   109
    #

    # We will use the AUTHORS-trained parameters output (see paper).

    # Note also that DDOMAIN only allows domains to be a single
    # continguous subsequence of chain anyway (i.e. not multiple
    # segments) (see paper), so we only ever have one segment in a
    # domain from this function.
    # It can, however, have a segment that includes parts of two chains
    # i.e. runs off one chain and includes (part of) another chain.
    # E.g. 1BAR:
    #
    # AUTHORS-trained parameters
    #   1  A    11  B     7
    #   2  B     8  B   138

    readout = False
    for line in fh:
        if line.strip()[0:8] == "AUTHORS-":
            readout = True
            continue
        elif line.strip()[0:5] == "SCOP-" or line.strip()[0:5] == "CATH-":
            break #finished with output
        if readout:
            splitline = line.split()
            if len(splitline) == 5:
                (domain_id,chainid_start,resnum_start,chainid_end,resnum_end) =\
                                                                    splitline
            elif len(splitline) == 3:
                # chain identifier is ' ' (space) so convert to '-' (dash)
                domain_id = splitline[0]
                chainid_start = '-'
                resnum_start = splitline[1]
                chainid_end = '-'
                resnum_end = splitline[2]
            else:
                sys.stderr.write(
                    'WARNING: error parsing line DDOMAIN output line:\n' 
                     + line)
                continue
            resnum_start = int(resnum_start)
            resnum_end = int(resnum_end)
            if resnum_start < 0:
                resnum_start = 0
                sys.stderr.write(
                   'WARNING: DDomain negative residue start number, set to 0\n')
            if resnum_end < 0:
                resnum_end = 0
                sys.stderr.write(
                   'WARNING: DDomain negative residue end number, set to 0\n')
            if chainid_start != chainid_end:
                sys.stderr.write('WARNING: DDomain different chainid in domain'+
                                 ' ' + str(domain_id) +
                                 ': splitting into two segments\n')
                # DDomain (as of Sep 2007) cannot have multiple segments in a
                # domain but it does sometimes have a different chain id for
                # start and end in a domain meaning (I think) that the
                # segment consists of the first chain from start residue up
                # to C-terminus and second chain from N-terminus up to
                # end residue (in that second chain).
                # So we make two segments, one in each of the two chains,
                # for this domain.

                start_chain = pdb_model[chainid_start]
                end_chain = pdb_model[chainid_end]
                # id of a residue in Bio.PDB is tuple (hetatm, resseqnum, icode)
                startchain_res_seqnum_list = [ res.get_id()[1] for res in
                                               start_chain.get_list()
                                               if res.get_id()[0] == ' ' ]
                max_startchain_resnum = max(startchain_res_seqnum_list)
                endchain_res_seqnum_list = [ res.get_id()[1] for res in
                                             end_chain.get_list()
                                             if res.get_id()[0] == ' ' ]
                min_endchain_resnum = min(endchain_res_seqnum_list)
                segment1 = PTSegment(chainid_start, resnum_start,
                                     max_startchain_resnum)
                segment2 = PTSegment(chainid_end, min_endchain_resnum,
                                     resnum_end)
                domain = PTDomain(domain_id, [segment1, segment2])
                domain_list.append(domain)
            else:
                if resnum_start > resnum_end:
                    # This happens e.g. on 1BMV chain 2 (only if chain 2
                    # only is requested). Don't know what it means really,
                    # but let's make sure we don't get an exception anyway
                    sys.stderr.write('WARNING: DDomain start resnum ' +
                                     str(resnum_start) + ' > end resnum ' +
                                     str(resnum_end) + ', swapping.\n')
                    tmp = resnum_start
                    resnum_start = resnum_end
                    resnum_end = tmp
                segment = PTSegment(chainid_start, resnum_start, resnum_end)
                domain = PTDomain(domain_id, [segment])
                domain_list.append(domain)
        
    return domain_list


            
def read_domains_from_cath_cdf_file(cdf_filename, pdbid, chainid=None):
    """
    Read the domain decomposition from the CATH Domall File (CDF)
    whose filename is supplied.

    These files and their description can be found at

    ftp://ftp.biochem.ucl.ac.uk/pub/cathdata/v3.1.0
    
    specifically the files README.CDF_FORMAT_2.0 for the description and
    the actual CDF file CathDomall.v3.1.0

    See http://www.cathdb.info/ for CATH in general.
    
    Parameters:
       cdf_filename - filename of the CATH CDF file to read
       pdbid - pdb identifier to read domains for
       chainid - (default None). If not None, get only this chain.

    Return value:
       List of PTDomain objects, one for each domain.

    Raises exceptions:
       NOTE: if the pdb id is not found in the file, we will return a list
       raise the NotInCATH_Exception
    
    """
    domain_list = []
    issued_warning = False
    # This is from README_CDF_FORMAT_2.0:
    #
    # KEY:
    # N  = Number of segments
    # C  = Chain character
    # I  = Insert character/code ('-' indicates no insert character)
    # S  = Start PDB number
    # E  = End PDB number
    # NR = number of residues (fragment information only)
    #
    # 1chmA  D02 F00  1  A    2 - A  156 -  1  A  157 - A  402 -
    #                 N |C    S I C    E I| N |C    S I C    E I|
    #                |<----Domain One---->|<-----Domain Two---->|
    #                   |<--Segment One-->|   |<--Segment One-->|
    #
    # This translates to:
    # 1chmA01 = Chain A; 2-156
    # 1chmA02 = Chain A; 157-402
    found = False
    cdf_fh = open(cdf_filename)
    for line in cdf_fh:
        line = line.lstrip().upper()
        if line[0] == '#':
            continue # skip comment lines
        cdf_rec = line.split()
        # CDF file appears to be sorted by PDB id ascending so we could
        # do a binary search or at least shortcut this loop but don't
        # really want to depend on it so won't bother.
        chain_name = cdf_rec[0] # 5 chars e.g. 1chmA
        if ( chain_name[:4] == pdbid.upper() and # matches our PDB identifier
             (chainid == None or chain_name[4] == chainid) ): # and chainid 
            found = True
            if cdf_rec[1][0] != 'D' or cdf_rec[2][0] != 'F':
                sys.stderr.write('WARNING: bad CDF record ignored: ' +
                                 line + '\n')
                continue
            num_domains = int(cdf_rec[1][1:])
            num_fragments = int(cdf_rec[2][1:])
            field = 3 # fields 0,1,2 where chainname, domains, fragments
            for domnum in range(num_domains):
                # CDF records actually have a different line for each
                # chain and so each chain is always a new domain in this
                # format it would appear. So we will name the domains
                # with chain identifier AND domain number so they are unique.
                # Though it seems not quite right that a new chain means
                # a new domain - domains should be able to contain multiple
                # segments of (different) chains but CDF records seem
                # to maintain chains as the higher level of hierarchy
                # (chain being part of the record identifier ie chainid
                # = pdbid + chainid, the first field).
                domain_id = (chain_name[4] +    # chain identifier
                             str(domnum + 1))   # number from one not zero
                segment_list = []
                num_segments = int(cdf_rec[field])
                field += 1
                for segnum in range(num_segments):
                    start_chainchar = cdf_rec[field]
                    field += 1
                    start_pdbnum = int(cdf_rec[field])
                    field += 1
                    insertcode1 = cdf_rec[field]
                    field += 1
                    end_chainchar = cdf_rec[field]
                    field += 1
                    end_pdbnum = int(cdf_rec[field])
                    field += 1
                    insertcode2 = cdf_rec[field]
                    field += 1

                    if start_chainchar != end_chainchar or \
                       start_chainchar != chain_name[4]:
                        # TODO: I think this never happens, but should
                        # do something with it anyway
                        sys.stderr.write('WARNING: mismatch chain characters '+
                                         'in CDF record: ' + line + '\n')
                    if start_chainchar == '0': # blank chainid in (old) PDB recs
                        if not issued_warning:
                            sys.stderr.write('WARNING: '\
                                             'CDF record ' + chain_name +\
                                             ' indicates no chain field '\
                                             'in PDB record. '\
                                             '(Will not work with'\
                                             ' remediated (2007) PDB files).\n'
                                             '  Changing chain to A.\n')
                            issued_warning = True
                        start_chainchar = 'A'
                    if start_pdbnum > end_pdbnum:
                        # This happens e.g. on 1BMV chain 2 (only if chain 2
                        # only is requested). Don't know what it means really,
                        # but let's make sure we don't get an exception anyway
                        sys.stderr.write('WARNING: CATH start resnum ' +
                                         str(start_pdbnum) + ' > end resnum ' +
                                         str(end_pdbnum) + ', swapping.\n')
                        tmp = start_pdbnum
                        start_pdbnum = end_pdbnum
                        end_pdbnum = tmp
                                         
                    segment_list.append(PTSegment(
                        cdf_chainid_to_stride_chainid(start_chainchar),
                        start_pdbnum, end_pdbnum))
                    
                domain_list.append(PTDomain(domain_id, segment_list))
            # we don't do anything with what CDF terms 'fragments'

    if not found:
        raise NotInCATH_Exception(pdbid)
    
    cdf_fh.close()
    if len(domain_list) > 0:
        return domain_list
    else:
        return [PTDomain(None, None)] # one-domain protein, no further info
    

def cdf_chainid_to_stride_chainid(cdf_chainid):
    """
    Convert a CDF (CATH Domall File) chainid to a STRIDE chainid.
    STRIDE uses '-' for a 'blank' chainid while PDB uses ' ' (space)
    and CATH (CDF) uses '0' where PDB has a blank (space) chain identfifer.
    We use the STRIDE convention ('-') in this program.
    So all this does is return the cdf_chainid unless it is '0', then
    it returns '-'.

    Parameters:
        pdb_chainid - the PDB chain identifier
    Return value:
        STRIDE chain identifier corresponding to supplied pdb_chainid
    """
    if cdf_chainid == '0':
        return '-'
    else:
        return cdf_chainid
    

def ptdomain_set_verbose(verb):
    """
    set the module global verbose flag in this module to supplied value
    Parameters: verb - True (for verbose output) or False
    Return value: None
    Uses globals: verbose (in this module)
    """
    global verbose
    verbose = verb
    

def fixup_crossdomain_sses(secstruct, domain_list):
    """
    Find any SSEs that span a domain boundary, and put each entirely
    in one domain.
    The domain is chosen as the one that contains most of the residues
    int the SSE.
    
    Parameters:
         secstruct - PTSecStruct (ptsecstruct.py) object descirbing SSEs
         domain_list - list of PTDomain objects representing all the domains
                          in this protein.
                          (in/out) NOTE: may be modified by having a segment
                           removed from a domain if SSE is only partly in 
                           the domain.
    Return value: None.
    """
    sse_list = ( [(start_chainid, start_resnum, end_chainid, end_resnum)
                  for (start_chainid, start_resnum, end_chainid, end_resnum)
                  in secstruct.strand_list] +
                 [(start_chainid, start_resnum, end_chainid, end_resnum)
                  for (start_chainid, start_resnum, end_chainid, end_resnum, helix_type)
                  in secstruct.helix_list] )
        
    for (start_chainid, start_resnum, end_chainid, end_resnum) in sse_list:
        for domain in domain_list:
            if (domain.is_in_domain(start_chainid,
                                    get_int_icode(start_resnum)[0])
                and not domain.is_in_domain(end_chainid,
                                            get_int_icode(end_resnum)[0]) ):
                # This really shouldn't happen, but does: domain
                # decomposition has determined that this SSE crosses
                # a domain boundary (really our SSE decisions don't
                # match whatever domain decomposition has done).
                # We'll have to assign the SSE to
                # a domain, and add the residues it spans into that
                # domain.

                # find domain2 as the other domain the SSE is also in
                for domain2 in domain_list:
                    if domain2 == domain:
                        continue
                    if domain2.is_in_domain(end_chainid,
                                            get_int_icode(end_resnum)[0]):
                        break
                        

                # find sse_domain as the domain with more residues of the
                # SSE in it

                domain_res_count = 0
                domain2_res_count = 0
                # FIXME: this is ignoring insertion codes etc., really
                # should convert to proper sequential residue sequence numbers
                # to do this
                start_resint = get_int_icode(start_resnum)[0]
                end_resint = get_int_icode(end_resnum)[0]
                for resint in range(start_resint, end_resint+1):
                    if domain.is_in_domain(start_chainid, resint):
                        domain_res_count += 1
                    elif domain2.is_in_domain(start_chainid, resint):
                        domain2_res_count += 1
                    else:
                        sys.stderr.write('ERROR: SSE in more than 2 domains\n')
                if domain2_res_count > domain_res_count:
                    sse_domain = domain2
                else:
                    sse_domain = domain # arbitrarily domain if equal count

                # first remove the segment from where it currently is
                seg = PTSegment(start_chainid,
                                get_int_icode(start_resnum)[0],
                                get_int_icode(end_resnum)[0])
#                print 'xxx',str(seg)
                for dom in domain_list:
#                    print 'aaa',str(dom)
                    dom.remove_segment(seg)
#                    print 'bbb',str(dom)

                    
                sys.stderr.write('WARNING: SSE ' + start_chainid + ':' +
                                 start_resnum + '-' + end_resnum +
                                 ' crosses domain boundary.\n'
                                 '  Put in domain ' + sse_domain.domainid +
                                 ' (' + str(sse_domain) + ').\n')
                sse_domain.add_segment(seg)
#                print 'zzz',str(sse_domain)

                break # no need to look at any more domains for this SSE

#     # DEBUG
#     for i in range(len(domain_list)):
#         sys.stdout.write(str(domain_list[i]))
#         if i < len(domain_list) - 1:
#             sys.stdout.write('/')
#     sys.stdout.write('\n')
#     # END DEBUG
