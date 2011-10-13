###############################################################################
#
# parsepdom.py - functions to parse the pDomains domain benchmark files
#
# File:    parsepdom.py
# Author:  Alex Stivala
# Created: December 2007
#
# $Id$
#
###############################################################################

"""
Functions to parse the pDomains protein domain decomposition benchmark
data files.

These files are available from http://pdomains.sdsc.edu

and the data and benhmarks are described in

Veretnik et al 2004 'Toward Consistent Assignment of Structural Domains in
Proteins' J. Mol. Biol. 339(3):647-678

and

Holland et al 2006 'Partitioning Protein Structures into Domains: Why is it so
Difficult?' J. Mol. Biol. 361:562-590

"""

import os,sys

from ptdomain import *

        
#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

def parse_pdomains_file(fh):
    """
    Parse a a pDomains file in the full (not raw) format, which has results
    of multiple methods for each chain. We will build a dictionary
    indexed first by chain (pdb and possibly chain identifier)
    where each entry is a dictionary indexed by method
    (PDP, DomainParser, DALI, NCBI, STERNBERG, CATH, SCOP). The values 
    for each of these methods is then the domain decomposition in the
    form a list of PTDomain (ptdomain.py) objects.

    Note that the STERNBERG method is from what is usually known as the AUTHORS
    assignment, from Islam et al 1995 'Identification and analysis of domains
    in proteins' Protein Engineering 8(6):513-525.

    This was developed with files downloaded from http://pdomains.sdsc.edu
    on 28Dec2007.

    Parameters:
       fh - filehandle open for read of a pDomains file (such as
            Benchmark_1_467)
    Return value:
       dictionary as described above
    """

    # Entries look like this:
    #
    #  Chain: 3adk
    #
    #  Method: PDP
    # Number of domains 2
    # Domain name: 3adk_a
    #  Number of fragments in  this domain: 2
    # Position of the fragment 1  start:1, end:37
    # Position of the fragment 2  start:76, end:194
    # Domain name: 3adk_b
    #  Number of fragments in  this domain: 1
    # Position of the fragment 1  start:38, end:75
    #

    # NOTE: DALI method tends to often put end residue sequence number 0
    # don't know what to do with that so we will omit DALI.
    
    chainid_dict = {} # { chainid : method_dict }
    method_dict = {} # { method string : list of PTDomain objects }
    domainlist = []
    domain = None
    segment = None
    method = None
    for line in fh:
        if line[:7] == " Chain:": # found a new entry, finish the old one
            if segment:
                domain.add_segment(segment)
                domainlist.append(domain)
                if method != "DALI": # have to omit DALI for now
                    method_dict[method] = domainlist
            if method_dict:
                chainid_dict[chainid] = dict(method_dict)
            chainid = line[8:].lstrip().rstrip('\n').rstrip().upper()
            if len(chainid) > 4:
                chainchar = chainid[4]
            else:
                chainchar = 'A' # use chain A by default (for remediated pdb)
            domainlist = []
            domain = segment = method = None
        elif line[:8] == " Method:":  # new method, finish the old one
            if segment:
                if domain:
                    domain.add_segment(segment)
                    domainlist.append(domain)
            if method != "DALI": # we'll have to omit DALI for now
                method_dict[method] = domainlist
            domainlist = []
            domain = segment = None
            method = line[9:].lstrip().rstrip('\n').rstrip()
        elif line[:12] == "Domain name:": # new domain, finish the old one
            if domain:
                if segment:
                    domain.add_segment(segment)
                domainlist.append(domain)
            domain = PTDomain(line[13:].lstrip().rstrip('\n').rstrip(), [])
            segment = None
        elif line[:24] == "Position of the fragment":
            if segment:
                domain.add_segment(segment) # new segment, finish the old one
            startline = line[line.index("start:"):]
            start_resnum = int(startline[6 : startline.index(',')])
            endline = line[line.index("end:"):]
            end_resnum = int(endline[4:].rstrip('\n'))
            try:
#                print chainid,method
                segment = PTSegment(chainchar, start_resnum, end_resnum)
            except ValueError:
                if method == "DALI":
                    pass  # DALI puts 0 as end res seq num, we'll omit it anyway
                else :
                    sys.stderr.write('WARNING: chainid ' + chainid + ' method '
                                     + method + ': end before start, '
                                     'ignoring fragment\n')
                    segment = None
                    
    if segment:
        domain.add_segment(segment)
        domainlist.append(domain)
        if method != "DALI": # have to omit DALI for now
            method_dict[method] = domainlist
    chainid_dict[chainid] = dict(method_dict)

    return chainid_dict
