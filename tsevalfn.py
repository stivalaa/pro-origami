#!/usr/bin/env python
###############################################################################
#
# tsevalfn.py - Evaluate tableaux search false negative rate against SCOP
#
# File:    tsevalfn.py
# Author:  Alex Stivala
# Created: June 2008
#
# $Id$
# 
###############################################################################

"""
Evaluate false negatives using SCOP at each cutoff (rank/score),
domains that are not included above the cuttoff but are in the SOCP
related domains for the query, are false negatives.

Output is to stdout in a format that is easily usable in R with the
read.table() function, i.e. we use '#' to prefix lines not to be parsed,
and have a header line with suitable R variable names for the columns.

See usage in docstring for main()

SCOP and ASTRAL data is obtained using the Bio.SCOP library (Casbon et
al 2006 'A high level interface to SCOP and ASTRAL implemented in
Python' BMC Bioinformatics 7:10) and depends on having the data
downloaded, in SCOP_DIR (defined below).

Downloaded SCOP files from

http://scop.mrc-lmb.cam.ac.uk/scop/parse/index.html

and ASTRAL files (in scopseq-1.73) from

http://astral.berkeley.edu/scopseq-1.73.html

The files downlaoded are:

/local/charikar/SCOP/:
dir.cla.scop.txt_1.73
dir.des.scop.txt_1.73
dir.hie.scop.txt_1.73

/local/charikar/SCOP/scopseq-1.73:
astral-scopdom-seqres-all-1.73.fa
astral-scopdom-seqres-sel-gs-bib-95-1.73.id

Other files there are indices built by Bio.SCOP when first used.
"""

import warnings # so we can suppress the annoying tempnam 'security' warning
import sys,os
import getopt

 
from Bio.SCOP import *

from tsevalutils import compute_auc,parse_searchresult,eval_fn,get_betagrasp_containing_domains

from pathdefs import SCOP_DIR,SCOP_VERSION

#-----------------------------------------------------------------------------
#
# Module globals
#
#-----------------------------------------------------------------------------

verbose = False


#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

def get_scop_domains(query_sid, scop):
    """
    Get a list of SCOP domains that have the same fold as the query domain,

    Parameters:
       query_sid - SCOP domain id (eg 'd1ubia_' of the domain to query.
       scop - previously built Bio.SCOP Scop instance
    Return value:
       list of Bio.SCOP domain instances that have the same fold as the
       query_sid.
    """

    if verbose:
        sys.stderr.write('getting domains related to ' + query_sid +'...\n')
    dom =scop.getDomainBySid(query_sid)
    fold = dom.getAscendent('fold')
    related = fold.getDescendents('domain')
    if verbose:
        sys.stderr.write('found %d domains\n' % len(related))
    return related
        

def get_domains_in_same_superfamily(sid, scop):
    """
    Return a list of SCOP domain instances that are in the same 
    superfamily as the supplied SCOP domain id.

    Parameters:
       sid - (string) SCOP domain identifier as sid e.g. d1ubia_
       scop - previuosly created Bio.SCOP Scop object

    Return value:
       list of SCOP domains that are in the same superfamily as the
       supplied sid (including the domain for the sid itself).
    """
    return scop.getDomainBySid(sid).getAscendent(
                                       'superfamily').getDescendents('domain')




def filter_domains_astral95(scop_domain_list, scop, astral):
   """
   Given list of Bio.SCOP domain objects, return list of those domains
   that are in the ASTRAL 95% sequence identity nonredundant subset.

   Parameters:
      scop_domain_list - list of Bio.SCOP domain objects
      scop - previously built Bio.SCOP Scop instance
      astral - previously build Bio.SCOP Astral instance
     
   """
   related95 = [ r for r in scop_domain_list if astral.isDomainInId(r, 95) ]
   return related95


def get_domains_not_in_searchresult(searchresult, scop, astral,
                                    use_nonredundant):
    """
    Return a list of domain identifiers that are in SCOP but are not
    present in the search result.

    Parameters:
      searchresult - list of (score, domainid)
      scop - previously built Bio.SCOP Scop instance
      astral - previously build Bio.SCOP Astral instance
      use_nonredundant - If True, filter out domains not in ASTRAL 95% nr subset

    Return value:
      list of domainid where each domainid is in SCOP but not in searchresult
    """
    scoproot = scop.getRoot()
    all_domains = scoproot.getDescendents('domain')
    if verbose:
        sys.stderr.write('got %d domains total\n' % len(all_domains))
    if use_nonredundant:
        all_domains = filter_domains_astral95(all_domains, scop, astral)
        if verbose:
            sys.stderr.write('filtered to %d domains in ASTRAL 95\n' %
                             len(all_domains))
    searchresult_dict = dict([(sid, True) 
                              for sid in 
                              [sid for (score,sid) in searchresult] ])
    domainid_list = [scopdom.sid for scopdom in all_domains 
                     if scopdom.sid not in searchresult_dict]
    return domainid_list


def get_goldstd_domains(query_sid, use_prosmos_dataset, use_superfamily,
                        use_nonredundant, 
                        scop, astral):
    """
    Geth the "gold standard" list of domains. This is the domains that are
    in the same fold or superfamily (according to use_superfamily option)
    as the supplied domain specified by query_sid (or the beta grasp containg
    domains if use_prosmos_dataset is True).

    Parameters:
      query_sid - SCOP domain identifier (sid) of query
      use_prosmos_dataset - Bool gold standard is the ProSMos beta-grasp data set
      use_superfamily - Bool use domains in same superfamily rather than same fold
      use_nonredundant - Bool if True use ASTRAL 95% nr subset
      scop - previously built Bio.SCOP Scop instance
      astral - previously build Bio.SCOP Astral instance

    Return value:
      list of sids for domains in same fold/superfamily as query_sid
    """
    if use_prosmos_dataset:
        goldstd_domains = get_betagrasp_containing_domains(scop)
    elif use_superfamily:
        if verbose:
            sys.stderr.write('getting domains in same superfamily as ' + query_sid
                             + '...\n')
        goldstd_domains = get_domains_in_same_superfamily(query_sid, scop)
        if verbose:
            sys.stderr.write('found ' + str(len(goldstd_domains)) + ' domains\n')
    else:
        goldstd_domains = get_scop_domains(query_sid,scop)

    if use_nonredundant:
        goldstd_domains = filter_domains_astral95(goldstd_domains, scop, astral)

    if verbose:
        sys.stderr.write('got ' + str(len(goldstd_domains)) + ' domains in ASTRAL.\n')

    return goldstd_domains


def write_slrtab(searchresult, goldstd_domains, bottom_score,
                 use_nonredundant, 
                 scop, astral):
    """
    Write the slrtab output to stdout. This is two columms with score in
    first column and label (0 or 1) in second column, for use with R
    ROCR package.
    
    Parameters:
       searchresult -list of (score, domainid) tuples parsed from pogram to eval
       goldstd_domains - list of gold standard domains to eval against
       bottom_score - if not None, a float for score to give any domain
                      not assigned a score in the searchresult
       use_nonredundant - Bool if True use ASTRAL 95% nr subset
       scop - previously built Bio.SCOP Scop instance
       astral - previously build Bio.SCOP Astral instance
    """
    # convert goldstd list of domainids to dictionary keyed by domainid
    # for fast lookup as we iterate through search results.
    # The dictionary is { domainid : True } (we don't have a value,
    # just need to quickly test for presence of domainid in gold std
    # postive list).
    goldstd_pos_dict = dict([(sid, True) for
                       sid in [scopdom.sid for scopdom in goldstd_domains]])
    for (score, domainid) in searchresult:
        if goldstd_pos_dict.has_key(domainid):
            label = 1
        else:
            label = 0
        sys.stdout.write('%20.8f  %d\n' % (score,label))

    if bottom_score != None:
        # write out an entry with the 'lowest' score for all domains
        # that are in the gold standard data but not given a score
        # in the parsed output from the program being evaluated.
        lowscore_domains = get_domains_not_in_searchresult(searchresult,
                                                           scop, astral,
                                                           use_nonredundant)
        for domainid in lowscore_domains:
            if goldstd_pos_dict.has_key(domainid):
                label = 1
            else:
                label = 0
            sys.stdout.write('%20.8f  %d\n' % (bottom_score, label))
        if verbose:
            sys.stderr.write('set score to %f for %d domains\n' %
                             (bottom_score, len(lowscore_domains)))




def tsevalfn_set_verbose(verb):
    """
    set the module global verbose flag in this module to supplied value
    Parameters: verb - True (for verbose output) or False
    Return value: None
    Uses globals: verbose (in this module)
    """
    global verbose
    verbose = verb

#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname + " [-xvrbfluno] [-z score] "
                     " <query_sid> <tabsearch.outputfile>\n")
    sys.stderr.write('  -v verbose messages to stderr\n')
    sys.stderr.write('  -x output is from external database, ignore key erors\n')
    sys.stderr.write('  -b use the ProSMoS data set for validation of beta-grasp query\n')
    sys.stderr.write('  -n negate scores\n')
    sys.stderr.write('  -o take log10 of scores\n')
    sys.stderr.write('  -f use full SCOP not ASTRAL 95% nonredundant subset\n')
    sys.stderr.write('  -l write table of values and labels intead of TPR,FPR\n')
    sys.stderr.write('  -u evaluate at superfamily rather than fold level\n')
    sys.stderr.write('  -z score : assign identifiers not present in the output a score of score\n')
    sys.stderr.write('  -e version: SCOP version (default %s)\n' %  SCOP_VERSION)
    sys.exit(1)

    
def main():
    """
    main for tsevalfn.py

    Usage: tsevalfn.py [-vbxnoflu] [-e version] [-z score]  <query_sid> <tabsearch.outputfile>

    
    -v turns on debug output to stderr

    -n negate all scores (so that most -ve score becomes highest score)
       (if used with -o then log10 is taken and then negated)

    -o take log10 of all scores

    -b use the PRoSMoS data set (Shi et al 2007) data set for validation
       of the beta-grasp substructure query

    -x results are from external program eg TableauSearch webserver so
       ignore SCOP sids that are in our database but not in search results
       (since they may be from older version of database that does not
       have these domains).  This is also needed when evaluating output
       from a program that does not output a score for every query
       (many only output scores for queries with score above some threshold
       e.g. VAST or top n scores), otherwise invalid rates and AUC are
       calculated. (This latter is not applicable to -l option).

    -b use the PRoSMoS data set (Shi et al 2007) data set for validation
       of the beta-grasp substructure query

    -f use the full SCOP data set not the ASTRAL 95% sequence identity subset

    -l instead of computing TPR/FPR table and AUC, just
       write table with one column of values (scores) and other of 
       corresponding true labels (0/1 for same/different fold as query)
       for use with the R ROCR packge.

    -u use domains from the same superfamily as the query rather than from
       the same fold as the query as the gold standard.

    -z score : any identifiers that are not present in the output but
               are in the gold standard data are given the specified score.
               This is for programs that do not assign a score to all
               domains, but only those above some threshold or the top n,
               or just cannot assign a score to some domains for some reason.
               This would generally be specified as some score lower 
               than all other scores.

    -e version: use SCOP version specified e.g. -e1.73

    <query_sid> is the SCOP id (e.g. 'd1ubia_') of the query domain
    
    <tabsearch.outputfile> is the output from the tabsearchqpml.file,
    which is a text file where each line is identifier then whitespace
    then score

    d1xksa_ -35.99999999
    d3sila_ -35.99999999
    ....
    d2mhua_ -0.499999999

    May be specified as - for stdin.

    The table of positive and false negative rates is printed to stdout.
    """
    negateflag = False
    ignore_search_keyerror = False
    use_prosmos_dataset = False
    use_nonredundant = True
    compute_rates = True
    use_superfamily = False
    bottom_score = None
    logflag = False
    scop_version = SCOP_VERSION

    try:
        opts,args = getopt.getopt(sys.argv[1:], "e:lbfxrnos:uz:v?")
    except:
        usage(os.path.basename(sys.argv[0]))
    for opt,arg in opts:
        if opt == "-s" : #self match score
            sys.stderr.write("-s option is obsolete; ignored\n")
        elif opt == "-r": # reversed: higher scores better
            sys.stderr.write("-r option is obsolete, fix caller: exiting\n")
            sys.exit(1)
        elif opt == "-n": # negate scores
            negateflag = True
        elif opt == "-o": # take log10 of scores
            logflag = True
        elif opt == "-v": # verbose
            tsevalfn_set_verbose(True)
        elif opt == "-x": # results from search of external db
            ignore_search_keyerror = True
        elif opt == "-b": # use the ProSMoS dataset to validate beta-grasp query
            use_prosmos_dataset = True
        elif opt == "-f": # use full SCOP not nonredundant subset
            use_nonredundant = False
        elif opt == "-l": # write scores and labels for ROCR, don't compute
            compute_rates = False
        elif opt == "-u": # evaluate at superfamily not fold level
            use_superfamily = True
        elif opt == "-z": # score to give to domains that have no score
            bottom_score = float(arg)
        elif opt == "-e": # specify SCOP version
            scop_version = float(arg)
        else:
            usage(os.path.basename(sys.argv[0]))

    if len(args) != 2:
        usage(os.path.basename(sys.argv[0]))

    query_sid = args[0]
    tabsearch_file = args[1]

    # read SCOP and ASTRAL data
    if verbose:
        sys.stderr.write('Reading SCOP data...\n')
    scop = Scop(dir_path=SCOP_DIR,version=scop_version)
    astral = Astral(dir_path=SCOP_DIR,version=scop_version,scop=scop)

    goldstd_domains = get_goldstd_domains(query_sid, use_prosmos_dataset,
                                          use_superfamily, 
                                          use_nonredundant, scop, astral)

    if verbose:
        sys.stderr.write('parsing search results...\n')
    if tabsearch_file == '-':
        tabsearch_fh = sys.stdin
    else:
        tabsearch_fh = open(tabsearch_file)
    (searchresult,commentlist) = parse_searchresult(tabsearch_fh, negateflag,
                                                    logflag)
    if tabsearch_file != '-':
        tabsearch_fh.close()

    sys.stdout.write('#' + ' '.join(sys.argv) + '\n') #identifying info about us
    sys.stdout.write('# results from:\n')
    for line in commentlist:
        sys.stdout.write('#  ')
        sys.stdout.write(line) # identifying information about search run
    sys.stdout.write('\n')
    if compute_rates:
        if bottom_score != None:
            lowscore_domains = get_domains_not_in_searchresult(searchresult,
                                                               scop, astral,
                                                               use_nonredundant)
            searchresult += [(bottom_score, sid) for sid in lowscore_domains]
            if verbose:
                sys.stderr.write('set score to %f for %d domains\n' %
                                 (bottom_score, len(lowscore_domains)))
            searchresult.sort() # sort by ascending score

        eval_fn([scopdom.sid for scopdom in goldstd_domains], 
                searchresult,
                ignore_search_keyerror)
    else:
        sys.stdout.write('score  label\n')
        sys.stdout.write('#-----------\n')
        write_slrtab(searchresult, goldstd_domains, bottom_score, 
                     use_nonredundant,
                     scop, astral)
            
if __name__ == "__main__":
    warnings.filterwarnings('ignore', 'tempnam', RuntimeWarning) 
    main()
