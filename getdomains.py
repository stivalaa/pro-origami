###############################################################################
#
# getdomains.py - functions to get domain decomposition from multiple methods
#
# File:    getdomains.py
# Author:  Alex Stivala
# Created: December 2007
#
# $Id$
#
###############################################################################

"""
Functions to get domain decompositions via external programs or databases.

"""

from ptdomain import *
from parsepdom import *

#-----------------------------------------------------------------------------
#
# Constants
#
#-----------------------------------------------------------------------------

# regular expression to match each supported method
valid_domain_programs = [r"ddomain", r"cath:.*", r"pdomains:.*", r"none"]

#-----------------------------------------------------------------------------
#
# Module globals
#
#-----------------------------------------------------------------------------

# dictionary of domains from pDomains file built by parse_pdomains_file()
# the first time something is requested from pDomains, the file is parsed
# and dictionary stored here, subsequently it is looked up in this dictionary
pdomains_dict = None

#-----------------------------------------------------------------------------
#
# Class definitions 
#
#-----------------------------------------------------------------------------


#
# Empty classes for exceptions
#

class Unsupported_Exception(Exception): # unsupported domain method found
    pass

        
class NotInpDomains_Exception(Exception): # ident not found in pDomains file
    pass


#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------


def get_domains(domain_program, pdbid, pdb_filename, pdb_struct, chainid=None):
    """
    Get domain decomposition in the form of list of PTDomain objects
    from the nominated domain program or database.

    Parmeters:
       domain_program - 'cath:cdf_file_name' (CATH) or other supported
                          program or database
       pdbid - PDB identifier for the protein
       pdb_filename - name of PDB file
       pdb_struct - Bio.PDB parsed structure (needed for DDOMAIN)
       chainid - (default None). If not None, the chainid for single chain
                 to use.

       Return value: list of domains according to the chosen method

    Raises exceptions:
       Unsupported_Exception for unsupported domain program/db
       NotInpDomains_Exception for identifier not in pDomains database
    """
    global pdomains_dict
    
    if domain_program[0:5] == "cath:":
        cdf_filename = domain_program[5:]
        try:
            domain_list = read_domains_from_cath_cdf_file(cdf_filename, pdbid,
                                                          chainid)
        except NotInCATH_Exception:
            sys.stderr.write('WARNING: PDB identifier ' + pdbid +
                             ' not found in CDF file.')
            sys.stderr.write(' Treating as single domain.\n')
            return [PTDomain(None, None)] # one-domain protein, no further info
    elif domain_program == "ddomain":
        domain_list = read_domains_from_ddomain(pdb_filename,
                                                pdb_struct[0],
                                                chainid)
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

    elif domain_program[0:9] == "pdomains:":
        if pdomains_dict == None:
            # Build the pDomains dictionary, subsequently look up in it
            pdomains_dict = parse_pdomains_file(open(domain_program[9:]))
        # TODO: we will just always use the STERNBERG ('AUTHORS') entry for now
        if chainid != None:
            pdbid += chainid
        try:
            domain_list = pdomains_dict[pdbid]['STERNBERG']
        except KeyError:
            raise NotInpDomains_Exception(pdbid)
    elif domain_program == "none":
        # The 'none' method is the baseline method of just always assigning
        # everything to a single domain
        return [PTDomain(None, None)] # one-domain protein, no further info
    else:
        raise Unsupported_Exception("unsupported domain program/db "
                                    + domain_program)

    return domain_list


def write_domains(fh, domain_list):
    """
    output the domain decomposition in a more or less conventional format

    Parmeters:
      fh - open for write filehandle to write domain decomposition to
      domain_list - list of PTDomain objects
    Return value:
      None
    """
    fh.write(str(len(domain_list)) + ' domains:\n')
    for i in range(len(domain_list)):
        fh.write(str(domain_list[i]))
        if i < len(domain_list) - 1:
            fh.write('/')
    fh.write('\n')
    

def verify_domain_disjoint(test_domlist, chaindict):
    """
    Check that the supplied domain decomposition is valid in that
    no residue is in more than one domain.

    Parameters:
       test_domlist - list of PTDomain from the test (predicted) decomposition
       chain_dict - dict of { chainid : (min_resnum, max_resnum) }

    Return value:
       True if valid domain decomposition else False (not disjoint)
    """
    for (chain, (min_resnum, max_resnum)) in chaindict.iteritems():
        for resnum in range(min_resnum, max_resnum+1):
            num_domains_with_residue = 0
            for domain in test_domlist:
                if (domain.is_in_domain(chain, resnum)):
                    num_domains_with_residue += 1
            if (num_domains_with_residue > 1):
                sys.stderr.write('ERROR: chain ' + chain + ' residue '
                                 + str(resnum) + ' is in ' +
                                 str(num_domains_with_residue) + ' domains.\n')
                return False
    return True
    
    

def build_domain_chaindict(domlist):
    """
    Build the diction of min and max residue sequence numbers for each
    chain in the supplied domain decompositions (list of PTDomain objects).

    Parameters:
         domlist - list of PTDomain objects represetning a domain decomp.
         
    Return value:
         dict of { chainid : (min_resnum, max_resnum) }
    """
    # get the lowest and highest residue sequence number in each chain
    # of the reference domains and build a dictionary from it.
    chain_dict = {} # dict of { chainid : (min_resnum, max_resnum) }
    for domain in domlist:
        for chainid in domain.get_chainids():
            (min_resnum,max_resnum)=domain.get_minmax_res_seq_in_chain(chainid)
            if chain_dict.has_key(chainid):
                (oldmin,oldmax) = chain_dict[chainid]
                if min_resnum < chain_dict[chainid][0]:
                    chain_dict[chainid] = (min_resnum, oldmax)
                if max_resnum > chain_dict[chainid][1]:
                    chain_dict[chainid] = (oldmin, max_resnum)
            else:
                chain_dict[chainid] = (min_resnum, max_resnum)
    return chain_dict

