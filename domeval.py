###############################################################################
#
# domeval.py - functions to evaluate domain decomposition accuracy
#
# File:    domeval.py
# Author:  Alex Stivala
# Created: December 2007
#
# $Id$
#
###############################################################################

"""
Functions to evaluate domain decomposition accuracy.

The accuracy of a domain decomposition is computed as overlap of predicted
(test) and assigned (reference) residues in the domain decomposition.
If a different number of domains is assigned, the decomposition is scored
as failed an no overlap is computed.

This measure is as defined by Jones et al 1998 'Domain assignment for protein
structure using a consensus approach: Characterization and analysis',
Protein Science 7:233-242.

"""
import os,sys

from Bio.PDB import *

from ptdomain import *
from parsepdom import *
from getdomains import *


def permutations(l):
    """
    Return all the permutations of the list l.
    Obivously this should only be used for extremely small lists (like less
    than 9 members).
    Paramters:
       l - list to find permutations of
    Return value:
       list of lists, each list a permutation of l.
    """
    if len(l) > 0:
        return [ [x] + perms for x in l
                 for perms in permutations([ e for e in l if e != x]) ]
    else:
        return [[]]

def compute_overlap_score(test_domlist, ref_domlist):
    """
    Compute the overlap score between the two domain decompositions of
    the same length represented by test and ref domlist.

    For two domain assignment methods, we can't assume they order
    their domains the same way, so we have to consider every possible
    mapping between them and use the best score, as dicussed in
    Veretnik et al 2004 J. Mol. Biol. 339(3):647-678.

    
    Parameters:
       test_domlist - list of PTDomain from the test (predicted) decomposition
       ref_domlist - list of PTDomain from the reference (gold standard) decomp.

    Return value:
        The maximum overlap score over all permutations of domains
    """
    assert(len(test_domlist) == len(ref_domlist))

    # get the lowest and highest residue sequence number in each chain
    # of the reference domains and build a dictionary from it.
    chain_dict = build_domain_chaindict(ref_domlist)

    # compute the total number of residues spanned by the refernce domain list
    total_residues = 0
    for (min_resnum, max_resnum) in chain_dict.itervalues():
        total_residues += max_resnum - min_resnum + 1

    # verify that the test decompositino is valid
    if not verify_domain_disjoint(test_domlist, chain_dict):
        sys.stderr.write('ERROR: invalid domain decomposition (not disjoint)\n')
        write_domains(sys.stderr, test_domlist)
        return 0.0
    
    # now compute the score for every possible mapping of the domains
    # (so it suffices to use every permutation of one of the lists keeping
    # the other fixed). Since we don't expect ever more than 8 domains
    # (and usually 4 or fewer) this should be ok, but still expensive.
    scores = [compute_overlap_score_ordered(test_permutation, ref_domlist,
                                            chain_dict, total_residues)
              for test_permutation in permutations(test_domlist)]
    return max(scores)
    

def compute_overlap_score_ordered(test_domlist, ref_domlist,
                                  chain_dict, total_residues):
    """
    For two domain lists of the same length, ordered so that
    corresponding domains 'line up' compute the overlap score
    (discussed above) as the fraction of residues that are assigned to
    the same domain.

    Note the ordering requirement. For two domain assignment methods,
    we can't assume they order their domains the same way, so this
    function has to be called multiple times with different orderings
    to find the one with the best score.

    This is a bit more complicated as we handle multiple chains.

    Parameters:
       test_domlist - list of PTDomain from the test (predicted) decomposition
       ref_domlist - list of PTDomain from the reference (gold standard) decomp.
       chain_dict - dict of {chainid : (min_resnum, max_resnum)} built by
                    caller.
       total_resides - total number of residues in protein.

    Return value:
       overlap score in [0, 1]
    """
    assert(len(test_domlist) == len(ref_domlist))

    # calculate the overlap score by going through each residue
    # and counting 1 for overlap between the two domain decompositions.
    #
    # TODO: we could probably more efficiently compute this using the
    # cut positions and a formula like the one in Emmert-Streib and Mushegian
    # 2007 BMC Bioinformatics 8:237
    # but this is easier (if very inefficient). Howevr the Emmert-Streib
    # equation assumes domains consist of only one segment (since that's
    # how their DomainICA algorithm works) so is not general enough for
    # our purposes.
    overlap_count = 0
    for i in range(len(ref_domlist)):
        for (chain, (min_resnum, max_resnum)) in chain_dict.iteritems():
            for resnum in range(min_resnum, max_resnum+1):
                if (test_domlist[i].is_in_domain(chain, resnum)  and
                    ref_domlist[i].is_in_domain(chain, resnum)):
                    overlap_count += 1
                    
    score = float(overlap_count) / float(total_residues)
    return score

                  
def domain_eval(test_domlist, ref_domlist):
    """
    If the two domain lists are the same length, compute the overlap score
    (discussed above) as the fraction of residues that are assigned to
    the same domain.

    Otherwise, describe the test decomposition as 'undercut' (fewer
    domains than reference) or 'overcut' (more domains than reference).

    Parameters:
       test_domlist - list of PTDomain from the test (predicted) decomposition
       ref_domlist - list of PTDomain from the reference (gold standard) decomp.

    Return value:
       tuple (description, score) where description is
       'undercut', 'overcut' or 'correct'
       and score is the overlap score in [0,1] if 'correct' otherwise 0.0

    """
    if len(test_domlist) < len(ref_domlist):
        return ('undercut', 0.0)
    elif len(test_domlist) > len(ref_domlist):
        return ('overcut', 0.0)
    else:
        return ('correct', compute_overlap_score(test_domlist, ref_domlist))



def evaluate_domains(domainlist, eval_domain_program, pdbid,
                     pdb_filename, pdb_struct, chainid=None):
    """
    Evaluate the performance of the domain decmoposiotion reprresented by
    the supplied domainlist against the program or database eval_domain_program

    Parmeters:
       domain_list - list of PTDomain for our decomposition
       eval_domain_program - 'cath:cdf_file_name' (CATH) or other supported
                          program or database (see ptdomain.py)
       pdbid - PDB identifier for the protein
       pdb_filename - name of PDB file (needed for DDOMAIN)
       pdb_struct - Bio.PDB parsed structure (needed for DDOMAIN)
       chainid - (default None). If not None, only use this chain.

       Return value: tuple the (num_domains, description, score)
                     where num_domains is the number of domains
                     in the reference (ie from the eval_domain_program)
                     and decription and score are from thetuple from domain_eval
                     (see domeval.py)

    """
    ref_domain_list = get_domains(eval_domain_program,
                                  pdbid, pdb_filename, pdb_struct,
                                  chainid)
    num_domains = len(ref_domain_list)
    if verbose:
        print eval_domain_program
        write_domains(sys.stdout, ref_domain_list)
        
    (description, score) = domain_eval(domainlist, ref_domain_list)
    return (num_domains, description, score)



def run_on_pdomains_list(pdbid_list,
                         pdb_root,
                         pdomains_filename,
                         print_results,
                         get_domains_function,
                         *get_domains_args):
    """
    Run the supplied domain decomposition function get_domains_fuction
    (with args get_domains_args) over all pDomains benchmark chains in
    the specified list of PDB/chain identifiers.
    Used for training/testing/crossvalidation for tuning parameters etc.

    Parameters:
       pdbid_list - list of PDB/chain identifiers (keys in dict built by
                     parse_pdomains_file())
       pdb_root - root of the PDB divided hierarchy to find PDB files.
       pdomains_filename - the fiename of the pDomains benchmark file
       print_results - If True, write results of each chain to stdout.
       get_domains_fuction - A function that, given the following args,
                    in order:
                      pdbid - PDB identifier
                      pdb_filename - PDB file name
                      pdb_struct - Bio.PDB parsed PDB structure
                      chainid - If not None, chain idnetifier to process
                    and then those
                    in get_domains_args, returns a domain decomposition
                    in the form of a list of PTDomain objects.
       get_domains_args - variable args list for get_dmoains_function
       
    Return value:
       tuple (undercut, overcut, correct, avgscore,
       num_correct_assign, numdonmains_dict, num_processed)
       where undercut, overcut, correct are number of domains that were
       undercut (too few domains) overcut (too many domains) or had the
       correct nubmer of domains, respectively, and avgscore is the
       average score for all chains (scoring 0.0 for undercut/overcut)
       and num_correct_assign is number correctly assigned (number of
       domains correct and score over threshold)
       and numdomains_dict is a dictionary of
       { num_domains : (frequency, totalscore, avgscore,undercut,overcut,num_correct_domains,num_correct_assign)}
       mapping number of domains to scores for chains with that number of
       domains.

    Raises Exceptions:
       ValueError if bad return value from evaluate_domains()

    """
    THRESHOLD_SCORE = 0.75 # must be above this to be correct
    
    total_score = 0.0
    num_undercut = 0
    num_overcut = 0
    num_correct = 0
    num_correct_assign = 0
    num_processed = 0
    numdomains_dict = {} # described in docstring above
    for pdbchain in pdbid_list:
        if len(pdbchain) > 4:
            pdbid = pdbchain[:4]
            chainid = pdbchain[4]
        else:
            pdbid = pdbchain[:4]
            chainid = None
        pdb_dir = os.path.join(pdb_root, pdbid[1:3].lower())
        pdb_filename = os.path.join(pdb_dir, 'pdb' + pdbid.lower() + '.ent.gz')

        if not os.path.exists(pdb_filename):
            sys.stderr.write("WARNING: pdb file " + pdb_filename +
                             " not found, skipping\n")
            continue
        
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
            # parse PDB file
            pdb_parser = PDBParser()
            pdb_struct = pdb_parser.get_structure(pdbid, our_pdb_filename) 

            # run the domain decomposition method and evaluate results
            domainlist = get_domains_function(pdbid, our_pdb_filename,
                                              pdb_struct,
                                              chainid,
                                              *get_domains_args)
            if domainlist == None:
                if chainid == None:
                    chainname = ''
                else:
                    chainname = 'chain ' + chainid
                sys.stderr.write('WARNING: domain decomposition failed for '
                                 + pdbid + ' ' + chainname + '\n')
                continue
            evalresult =  evaluate_domains(domainlist,
                                           "pdomains:" + pdomains_filename,
                                           pdbid,
                                           our_pdb_filename, pdb_struct,
                                           chainid)
            (num_domains, description, score) = evalresult
            num_processed += 1
        finally:
            if used_tmp_file:
                cleanup_tmpdir(TMPDIR)

        assigndescr = 'incorrect'
        if description == 'undercut':
            num_undercut += 1
        elif description == 'overcut':
            num_overcut += 1
        elif description == 'correct':
            num_correct += 1
            if score > THRESHOLD_SCORE:
                num_correct_assign += 1
                assigndescr = 'correct'
        else:
            raise ValueError('unknown description ' + description +
                             ' from evaluate_domains\n')
        if print_results:
            sys.stdout.write(pdbchain + '\t' + str(num_domains) + '\t' +
                             description + '\t' + str(score) + ' ' +
                             assigndescr + '\n' )
            
        total_score += score
        if numdomains_dict.has_key(num_domains):
            (dfrequency, dtotalscore, davgscore,
             dundercut,dovercut,dnum_correct_domains,dnum_correct_assign) = \
                                                  numdomains_dict[num_domains]
        else:
            dfrequency = 0
            dtotalscore = 0.0
            davgscore = 0.0
            dundercut = 0
            dovercut = 0
            dnum_correct_domains = 0
            dnum_correct_assign = 0
        dfrequency += 1
        dtotalscore += score
        if description == 'undercut':
            dundercut += 1
        elif description == 'overcut':
            dovercut += 1
        elif description == 'correct':
            dnum_correct_domains += 1
            if score > THRESHOLD_SCORE:
                dnum_correct_assign += 1
        else:
            assert(False)
        numdomains_dict[num_domains] = (dfrequency, dtotalscore, davgscore,
                                        dundercut,dovercut,
                                        dnum_correct_domains,
                                        dnum_correct_assign) 

    for num_domains in numdomains_dict.iterkeys():
        (freq, total, avg, dunder,dover,dnumcd,dnumca) = numdomains_dict[num_domains]
        avg = total / float(freq)
        numdomains_dict[num_domains] = (freq,total,avg,dunder,dover,dnumcd,dnumca)
        
    avgscore = total_score / float(num_processed)
    return (num_undercut, num_overcut, num_correct, avgscore,
            num_correct_assign, numdomains_dict, num_processed)



def run_on_pdomains_file(pdb_root,
                         pdomains_filename,
                         print_results,
                         get_domains_function,
                         *get_domains_args):

    """
    Run the domain decomposition over all pDomains benchmark chains in
    the specified pDomains benchamge file.
    Used for training/testing/crossvalidation for tuning parameters etc.

    Parameters:
       pdb_root - root of the PDB divided hierarchy to find PDB files.
       pdomains_filename - the fiename of the pDomains benchmark file
       print_results - If True, print results for each chain to stdout
       get_domains_fuction - A function that, given the following args,
                    in order:
                      pdbid - PDB identifier
                      pdb_filename - PDB file name
                      pdb_struct - Bio.PDB parsed PDB structure
                      chainid - if not None, chain identnifier of chain
                    and then those
                    in get_domains_args, returns a domain decomposition
       get_domains_args - variable args for get-domains_function

    Return value:
       tuple (undercut, overcut, correct, avgscore,
       num_correct_assign, numdomain_dict)
       as described in run_on_pdomains_list()
       where undercut, overcut, correct are number of domains that were
       undercut (too few domains) overcut (too many domains) or had the
       correct nubmer of domains, respectively, and avgscore is the
       average score for all chains (scoring 0.0 for undercut/overcut)
       and num_correct_assign is number assigned correctly (correct domain
       number and score over threshold)
       and numdonains_dict maps number of domains to scores,
       as described in run_on_pdomains_list()

    """
    pdomains = parse_pdomains_file(open(pdomains_filename))
    return run_on_pdomains_list(pdomains.keys(), pdb_root, pdomains_filename,
                                print_results,
                                get_domains_function,
                                *get_domains_args)



def print_scores(num_processed,
                 num_undercut, num_overcut, num_correct, num_correct_assign,
                 avgscore, indent=0):
    """
    Neatly format the scores to stdout.
    Parameters:
       num_undercut -   number overcut
       num_overcut -    number undercut
       num_correct -    number of correctly assigned domain numbers
       num_correct_assign - number of correctly assigned domains
       avgscore -   average score in [0,1]
       indent - (default 0) number of spaces to indent
    Return value:
       None
    """
    sys.stdout.write(indent*' ')
    sys.stdout.write("number processed: %d\n" % num_processed)
    sys.stdout.write(indent*' ')
    sys.stdout.write("undercut:         %d\n" % num_undercut)
    sys.stdout.write(indent*' ')
    sys.stdout.write("overcut:          %d\n" % num_overcut)
    sys.stdout.write(indent*' ')
    sys.stdout.write("correct domains:  %d\n" % num_correct)
    sys.stdout.write(indent*' ')
    sys.stdout.write("correct assign:   %d\n" % num_correct_assign)
    sys.stdout.write(indent*' ')
    sys.stdout.write("average score:    %3.1f%%\n" % (avgscore*100.0))
    
