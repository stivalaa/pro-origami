###############################################################################
#
# tsevalutils.py - Functions for tableau search evaluation
#
# File:    tsevalutils.py
# Author:  Alex Stivala
# Created: June 2008
#
# $Id$
# 
###############################################################################

"""
Functions for protein db search program evaluation scripts tsevalfn.py
and tsevalrmsd.py and mkroc50tab, and others.

Includes functions to parse output, compue AUC and ROC50 values,
get lists of true hits from SCOP, filter down to nonredundant ASTRAL
subsets, etc.
"""

import sys
from math import log10

from Bio.SCOP import *


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


def compute_auc(fprlist, tprlist):
    """
    Compute area under ROC curve (AUC) given a list of false postive
    rates and corresponding list of true postive rates.
    Sums the areas of trapezoids formed by each adjacent entry in
    the lists.

    Parameters:
        fprlist - list of false positive rates (FPR), ordered from 0 to 1
        tprlist - corresponding list of true positive rates (TPR)
    Return value:
        AUC (area under [ROC] curve)
    """
    # in R: auc <- sum(diff(fprlist)*(tprlist[-1]+tprlist[-length(tprlist)]))/2

    n = len(fprlist)
    assert(len(tprlist) == n)
    widths = [b - a for (a,b) in zip(fprlist[:n-1], fprlist[1:])]
    avgheights = [(a+b)/2.0 for (a,b) in zip(tprlist[1:], tprlist[:n-1])]
    auc = sum([x*y for (x,y) in zip(widths, avgheights)])
    return auc




def parse_searchresult(fh, negateflag=False, logflag=False, sortflag=True):
    """
    Parse a text file where each line is identifier then whitespace
    then score,  e.g.

    d1xksa_ -35.99999999
    d3sila_ -35.99999999
    ....
    d2mhua_ -0.499999999

    Then sort it by score from lowest to highest
    ie this means the best hit as the bottom of the file, and worst at top.
    If the negateflag is True then the scores are negated before sorting
    (so in the example above where -36 is the 'best' score, then it becomes
    36 and so is last in the list after sorting).

    Parameters:
        fh - open (read) filehandle to parse from
        negateflag - if True, negate scores before sorting
        logflag - if True take log10 of scores. If used with negateflag
                  then log10 is taken first then negated
        sortflag - If True sort results by score ascending. Default True.

    Return value:
         tuple (scorelist, commentlines) where
         scorelist is list of (score, domainid), sorted by score (if 
           sortflag=TRUE)
         commentlines is list of strings, each one a comment line from file
    """
    reslist = []
    commentlist = []
    for line in fh:
        if line[0] == '#':
            commentlist.append(line)
            continue
        splitline = line.split()
        if len(splitline) < 2:
            sys.stderr.write('bad line: ' + line + '\n')
            continue
        domainid = splitline[0]
        score_str = splitline[1]
        if score_str.lower() == 'nan' or score_str == '********':
            # if we get any NaN values then then sort() gets completely
            # screwed up and the evaluation is all wrong, so skip them.
            sys.stderr.write('skipping NaN: ' + line + '\n')
            continue
        try:
            score = float(score_str)
        except ValueError:
            sys.stderr.write('skipping invalid score ' + line + '\n')
            continue
        if logflag:
            score = log10(score)
        if negateflag:
            score = -score
        reslist.append((score, domainid))
    if sortflag:
        reslist.sort()  # sort ascending
    return (reslist, commentlist)

    

def eval_fn(goldstd_positive_domains, searchresult, 
            ignore_search_keyerror=False):
    """
    Evaluate the false negative rate for different score cutoffs
    and print table to stdout.
    
    Parameters:
       goldstd_positive_domains - 
                       list of identifiers that the query
                       should match i.e. the 'gold standard' list of 
                       hits.
       searchresult -  list of (score, domainid), sorted by ascending score
        ignore_search_keyerror - If True, ignore key errors in search result.
                      Usually this should not happen, but when using external
                      results i.e. those not from same database (ASTRAL subset)
                      used here, it can. Eg. using TableauSearch webserver
                      results which have older SCOP as database,
                      so there are  SCOP sids in our database that are
                      not in the serach results at all. Then this option
                      just ignores them rather than raising exception KeyError.

    Return value:
       None - table written to stdout.
       
    """
        
    sys.stdout.write('score  fp_count     tpr     fpr\n')
    sys.stdout.write('#------------------------------\n')

    fprlist = []
    tprlist = []

    resultlist = searchresult

    # remove gold standard SCOP domains that are not in the search
    # result, which can happen from using external results, also
    # happens for some domains where for some reason the Bio.SCOP
    # isDomainInId(95) does not contain some domains that actually
    # are in the downlaoded ASTRAL SCOP 95% set (e.g. d2aeph1)
    # This can also happen when the search got an error for a domain
    # so it is not in the search results.
    # Use the -x option (ignore_search_keyerror) to handle these cases
    # (not always on by default since these things "shouldn't" happen).
    if ignore_search_keyerror:
        search_dict = dict([(pdbid, (score, rank)) for (rank, (score, pdbid))
                            in enumerate(searchresult)])
        our_goldstd_positive_domains = [scopdom for scopdom in goldstd_positive_domains if search_dict.has_key(scopdom) ]
    else:
        our_goldstd_positive_domains = goldstd_positive_domains

#    sys.stderr.write('original len, reduced len ' +str(len(goldstd_positive_domains)) +  ' , ' + str(len(our_goldstd_positive_domains)) + '\n')

    
    # convert goldstd list of domainids to dictionary
    # keyed by domainid for fast lookup as we iterate through search results.
    # The dictionary is { domainid : True } (we don't have a value,
    # just need to quickly test for presendce of domainid in gold stad
    # postive list0

    goldstd_pos_dict = dict([(scopdom, True) for
                             scopdom in our_goldstd_positive_domains])
    
    # start at classifying all as true (TPR=FPR=1)
    tp_count = len(our_goldstd_positive_domains)
    fp_count = len(resultlist) - len(our_goldstd_positive_domains)
    for cutoff_rank in xrange(len(resultlist)):
        cutoff_score = resultlist[cutoff_rank][0]
        if cutoff_rank > 0:
            # we are now classifying the previous one and all below as negative
            prev_scopsid =  resultlist[cutoff_rank - 1][1]
            if goldstd_pos_dict.has_key(prev_scopsid):
                tp_count -= 1
            else:
                fp_count -= 1

        tpr = float(tp_count) / float(len(our_goldstd_positive_domains)) #sensitivity = true pos rate
        fpr = float(fp_count) / float(len(resultlist) - len(our_goldstd_positive_domains)) #FP rate
        specificity = 1.0 - fpr

        fprlist.append(fpr)
        tprlist.append(tpr)
            
        sys.stdout.write('%5.1f %8d    %5.3f   %5.3f\n' %
                         (cutoff_score, fp_count, tpr, fpr))

    fprlist.reverse()
    tprlist.reverse()
    auc = compute_auc(fprlist, tprlist)
    sys.stdout.write('\n')
    sys.stdout.write('# AUC = %5.3f\n' % auc)


def iter_searchresult(fh, multiquery = False, skip_self_query = False,
                      negateflag = False, logflag = False):
    """
    This is a generator function version of parse_searchresult(),
    that yields one searh result at a time from the file, instead
    of parsing the whole thing and returning a list. Note because of this
    it cannot sort the result, the file must already be sorted.

    Also, if multiquery=True then this parses a format with multiple
    queryies in the one file. In this format, a special 'comment' line
    delimits results from each queryu. This line has the format:

    # QUERYID = querid

    or
    
    # QUERY ID = querid
    
    e.g.
    
    # QUERYID = d1ubia_

    to delimit results for each separate query. The blastout2col.sh etc.
    scripts create this format. Note that the ordering of results 
    (scores) is only within each individual query.
    
    The format of the file is identifier then whitespace then score e.g.:

    d2mhua_ 0.499999999
    ....
    d1xksa_ 35.99999999
    d3sila_ 35.99999999


    Parameters:
        fh - open (read) filehandle to parse from
        multiquery - Boolean. If True multiquery format as described above.
        skip_self_query - Boolean. If True, omit match of query against itself
                          (only for multiquery=True).
        negateflag - if True, negate scores before sorting
        logflag - if True take log10 of scores. If used with negateflag
                  then log10 is taken first then negated

    Return value:
         Yields a tuple (score, domainid) for each search result line in the
         file. Ordered according to the file order, i.e. must be from
         worst to best hit.
         If multiquery  = True then the tuples yielded are 
         (queryid, score, domainid). All same queryid must be consecutive
         and ordering is only within each block of consecutive identical
         queryid.
    """
    queryid = None
    for line in fh:
        if line[0] == '#':
            if multiquery:
                splitline = line.split()
                if ( (len(splitline) > 1 and splitline[1] == 'QUERYID') or
                     (len(splitline) > 2 and splitline[1] == 'QUERY' and
                      splitline[2] == 'ID') ):
                    queryid = line.split('=')[1].lstrip().rstrip()
            continue
        splitline = line.split()
        if len(splitline) != 2:
            sys.stderr.write('bad line: ' + line + '\n')
            continue
        domainid = splitline[0]
        score_str = splitline[1]
        if score_str.lower() == 'nan' or score_str == '********':
            # if we get any NaN values then then sort() gets completely
            # screwed up and the evaluation is all wrong, so skip them.
            sys.stderr.write('skipping NaN: ' + line + '\n')
            continue
        score = float(score_str)

        if logflag:
            if (score == 0.0):
                score = -1e308  # dodgy: very small number
            else:
                score = log10(score)
        if negateflag:
            score = -score

        if multiquery:
            if skip_self_query and queryid.lower() == domainid.lower():
                continue
            else:
                yield((queryid, score, domainid))
        else:
            yield((score, domainid))


def get_searchresult_info(fh):
    """
    Used in conjunction with iter_searchresult() to get comments and 
    number of search result lines from the search results file; cannot
    combine this functionality in iter_searchresult() as that is a generator
    function unlike parse_searchresult() - this function and iter_searchresult()
    together form the generator function replacement for parse_searchresult()
    for very large files so we don't have to read the whole file at once
    into core.

    Parameters:
        fh - open (read) filehandle to parse from

    Return value:
       tuple (num_results, commentlines) 
         num_results is number of results lines in file, i.e. the number
         of result tuples that will be returned by iter_searchresult()
         commentlines is list of strings, each one a comment line from file

    """
    num_results = 0
    commentlist = []
    for line in fh:
        if line[0] == '#':
            commentlist.append(line)
            continue
        splitline = line.split()
        if len(splitline) != 2:
            sys.stderr.write('bad line: ' + line + '\n')
            continue
        domainid = splitline[0]
        score_str = splitline[1]
        if score_str.lower() == 'nan' or score_str == '********':
            # if we get any NaN values then then sort() gets completely
            # screwed up and the evaluation is all wrong, so skip them.
            sys.stderr.write('skipping NaN: ' + line + '\n')
            continue
        num_results += 1
    return (num_results, commentlist)


def iter_slrtab(fh):
    """
    This is a generator function
    that yields one (score, label) tuple at a time from the slrtab
    file, instead of parsing the whole thing and returning a list.
    Note because of this
    it cannot sort the result, the file must already be sorted, from
    'lowest' (i.e. least likely to be a positive instance, i.e. 'worst hit')
    to 'highest' (least likely to be a postivie instance).

    The format of the slrtab file is score then class label (0 or 1) e.g.

    0.499999999  0
    ...
    35.99999999  1
    35.99999999  1


    Parameters:
        fh - open (read) filehandle to parse from

    Return value:
         Yields a tuple (score, label) for each search result line in the
         file. Ordered according to the file order, i.e. must be from
         best to worst hit.
    """
    for line in fh:
        if line[0] == '#':
            continue
        splitline = line.split()
        if len(splitline) == 0:
            continue # skip blank lines
        if splitline[0] == "score":
            continue
        if len(splitline) != 2:
            sys.stderr.write('bad line: ' + line + '\n')
            continue
        score_str = splitline[0]
        label_str = splitline[1]
        score = float(score_str)
        label = int(label_str)
        yield((score, label))


def get_slrtab_info(fh):
    """
    Used in conjunction with iter_slrtab() to get comments and 
    number of (Score, label) lines from the slrtab file; cannot
    combine this functionality in iter_slrtab() as that is a generator.

    Parameters:
        fh - open (read) filehandle to parse from

    Return value:
       tuple (num_results, num_positive, commentlines, total_tp_possible,
              num_queries) 
         num_results is number of results lines in file, i.e. the number
         of (score, label) tuples that will be returned by iter_slrtab()
         num_positive is number of results with class label 1
         commentlines is list of strings, each one a comment line from file
         total_tp_possible is the total number of true positives possible in the database
         num_queries is the number of queries in the slrtab

           total_tp_possible and num_queries are parsed from the 
           TOTAL_TP_POSSIBLE and NUM_QUERIES comment lines generated
           by mkslrtabmultiquery.py or mkslrtab.py

    """
    num_results = 0
    num_positive = 0
    total_tp_possible = None
    num_queries = None
    commentlist = []
    for line in fh:
        if line[0] == '#':
            commentlist.append(line)
            splitline = line.split('=')
            if len(splitline) == 2:
                if splitline[0] == '# TOTAL_TP_POSSIBLE ':
                    total_tp_possible = int(splitline[1])
                elif splitline[0] == '# NUM_QUERIES ':
                    num_queries = int(splitline[1])
            continue
        splitline = line.split()
        if len(splitline) != 2:
            sys.stderr.write('bad line: ' + line + '\n')
            continue
        num_results += 1
        score_str = splitline[0]
        label_str = splitline[1]
        score = float(score_str)
        label = int(label_str)
        if (label == 1):
            num_positive += 1
        elif (label != 0):
            sys.stderr.write('bad label "%d"\n' % label)

    return (num_results, num_positive, commentlist,
            total_tp_possible, num_queries)



def get_domain_by_sid(query_sid, scop):
    """
    Get a Bio.SCOP domain id for supplied sid. This could be done
    directly by scop.getDomainBySid(query_sid) except the various 
    exceptions and error conditions that can arise and we need to deal
    with (see comments in function).

    Parameters:
        query_sid - SCOP domain id (eg 'd1ubia_' of the domain to query.    
        scop - previously built Bio.SCOP Scop instance

    Return value:
        Bio.SCOP domain instance fo the query_sid
    """
    dom =scop.getDomainBySid(query_sid)
    if dom == None:
        # some domains seem to start with g rather than d
        # e.g. g1avo.4 or g1dk9.1 
        # they are in the astral-scopdom-seqres-all-1.73.fa file
        # with the g at the start, but for some reason get None when
        # looking up that sid, need to replace g with d instead.
        if query_sid[0] == 'g':
            query_sid = 'd' + query_sid[1:]
        # when using sids from older versions of ASTRAL SCOP
        # (e.g. 1.65 as used in Hulsen et al. (2006) sometimes
        # we get identifiers with _ as chainid e.g. d1ab4__ which
        # have been changed to have a as chainid instead as per recent PDB
        # remediation. So we'll try chainid a instead of _
        if query_sid[5] == '_':
            query_sid = query_sid[:5] + 'a' + query_sid[6]

        dom = scop.getDomainBySid(query_sid)

    return dom


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
        sys.stderr.write('getting domains in same fold as ' + query_sid +'...\n')
    dom = get_domain_by_sid(query_sid, scop)
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
    dom = get_domain_by_sid(sid, scop)
    return dom.getAscendent('superfamily').getDescendents('domain')


def get_domains_in_same_family(sid, scop):
    """
    Return a list of SCOP domain instances that are in the same 
    family as the supplied SCOP domain id.

    Parameters:
       sid - (string) SCOP domain identifier as sid e.g. d1ubia_
       scop - previuosly created Bio.SCOP Scop object

    Return value:
       list of SCOP domains that are in the same superfamily as the
       supplied sid (including the domain for the sid itself).
    """
    dom = get_domain_by_sid(sid, scop)
    return dom.getAscendent('family').getDescendents('domain')


def filter_domains_astral_nrpercent(scop_domain_list, scop, astral, nrpercent):
   """
   Given list of Bio.SCOP domain objects, return list of those domains
   that are in the ASTRAL nr% sequence identity nonredundant subset.

   Parameters:
      scop_domain_list - list of Bio.SCOP domain objects
      scop - previously built Bio.SCOP Scop instance
      astral - previously build Bio.SCOP Astral instance
      nrpercent - (integer) percent sequence identity subset to use (e.g. 95)
     
   """
   related = [r for r in scop_domain_list if astral.isDomainInId(r, nrpercent)]
   return related


def get_db_size(scop, astral, nrpercent):
    """
    Return number of domains in supplied SCOP. If nrpercent is not None
    then return nubmer of domains in specified ASTRAL sequence identity
    nonredundant subset.

    Parameters:
      scop - previously built Bio.SCOP Scop instance
      astral - previously build Bio.SCOP Astral instance
      nrpercent - (integer) percent sequence identity subset to use (e.g. 95)
                  or None for whole SCOP.

    Return value:
      Number of domains in SCOP or ASTRAL nonredundant subset                  
    """
    all_domains = scop.getRoot().getDescendents('domain')
    if nrpercent != None:
        all_domains = filter_domains_astral_nrpercent(all_domains,
                                                      scop, astral,
                                                      nrpercent)
    return len(all_domains)
                      
                      


def get_seq_db_size(scop, astral, nrpercent):
    """
    Return number of sequences in supplied SCOP. If nrpercent is not None
    then return nubmer of domains in specified ASTRAL sequence identity
    nonredundant subset.
    This is not the same as get_db_size(), which gives all the domains
    in SCOP, here we try to find the number of domains that have
    sequences in the FASTA (or actually, due to the way Bio.SCOP works,
    in the .id file, which is sometimes not the same for reasons
    I don't understand).

    Parameters:
      scop - previously built Bio.SCOP Scop instance
      astral - previously build Bio.SCOP Astral instance
      nrpercent - (integer) percent sequence identity subset to use (e.g. 95)
                  or None for whole SCOP.

    Return value:
      Number of domains in SCOP or ASTRAL nonredundant subset                  
    """
    if nrpercent == None:
        db_size = len(astral.fasta_dict) # dodgy direct access, but works
    else:
        db_size = len(astral.domainsClusteredById(nrpercent))
    return db_size
                      
                     
 
def get_goldstd_domains(query_sid, level,
                        use_nonredundant, 
                        scop, astral, nrpercent):
    """
    Geth the "gold standard" list of domains. This is the domains that are
    in the same fold or superfamily or family (according to level option)
    as the supplied domain specified by query_sid.

    Parameters:
      query_sid - SCOP domain identifier (sid) of query
      level - 'fold' or 'superfamily' or 'family' - the level of SCOP
              hierarchy to use to define gold standard.
      use_nonredundant - Bool if True use ASTRAL 95% nr subset
      scop - previously built Bio.SCOP Scop instance
      astral - previously build Bio.SCOP Astral instance
      nrpercent - (integer) percent sequence identity subset to use (e.g. 95)
                  Only used if use_nonredundant=True

    Return value:
      list of sids for domains in same fold/superfamily as query_sid
    """
    if level == 'superfamily':
        if verbose:
            sys.stderr.write('getting domains in same superfamily as ' + query_sid
                             + '...\n')
        goldstd_domains = get_domains_in_same_superfamily(query_sid, scop)
        if verbose:
            sys.stderr.write('found ' + str(len(goldstd_domains)) + ' domains\n')
    elif level == 'family':
        if verbose:
            sys.stderr.write('getting domains in same family as ' + query_sid
                             + '...\n')
        goldstd_domains = get_domains_in_same_family(query_sid, scop)
        if verbose:
            sys.stderr.write('found ' + str(len(goldstd_domains)) + ' domains\n')
    elif level == 'fold':
        goldstd_domains = get_scop_domains(query_sid,scop)
    else:
        raise ValueError('unknown level: ' + level)

    if use_nonredundant:
        goldstd_domains = filter_domains_astral_nrpercent(goldstd_domains, scop, astral, nrpercent)
        if verbose:
            sys.stderr.write('got ' + str(len(goldstd_domains)) + ' domains in ASTRAL ' + str(nrpercent) + '% sequence nr subset.\n' )

    return goldstd_domains


def get_betagrasp_containing_domains(scop):
    """
    Get a list of SCOP domains that contain the beta-grasp motif.
    Note that these are not just the domains in the SCOP beta-grasp
    (ubiquitin-like) fold, but the three categories described in

    Shi et al (2007) 'Searching for three-dimensional secondary structure
    patterns in proteins with ProSMoS' Bioinformatics 23(11):1331-1338

    Refer to Table 1 in the above for the list of supefamilies built
    by this subroutine.

    Parameters:
        scop - previously built Bio.SCOP Scop instance

    Return value:
        list of Bio.SCOP domain instances that contain the beta-grasp
        fold, either as (1) the core (2) a gregarious fold or (3) strutural
        drift (see Shi et al 2007 and refernces therein).

    Implemented with SCOP 1.73 (so note some sids are different due to
    remediation of chainid e.g. d1ubq__ (SCOP 1.69 as used in Shi et al)
    is now d1ubqa_ in SCOP 1.73).
    """
    #
    # Category 1: beta-grasp core
    #

    # ubiquitin-like
    
#     core = get_domains_in_same_superfamily('d1ubqa_', scop)
#     core += get_domains_in_same_superfamily('d1c9fa_', scop)
#     core += get_domains_in_same_superfamily('d1fm0d_', scop)
#     core += get_domains_in_same_superfamily('d1frda_', scop)
#     core += get_domains_in_same_superfamily('d2saka_', scop)
#     core += get_domains_in_same_superfamily('d1an8a2', scop)
#     core += get_domains_in_same_superfamily('d1pgxa_', scop)
#     core += get_domains_in_same_superfamily('d1tifa_', scop)
#     core += get_domains_in_same_superfamily('d1f1hl1', scop)
#     core += get_domains_in_same_superfamily('d1qf6a2', scop)
#     core += get_domains_in_same_superfamily('d1mfwa_', scop)
#     core += get_domains_in_same_superfamily('d1t0qc_', scop)
#     core += get_domains_in_same_superfamily('d2fug13', scop) # Nqo1 like
#     core += get_domains_in_same_superfamily('d2gria1', scop) # NSP3A like

    core = get_scop_domains('d1ubia_', scop) # ubiquitin-like fold
        
    # Nudix
    core += get_domains_in_same_superfamily('d1muta_', scop)

    # Anthrax protective antigen
    core += get_domains_in_same_superfamily('d1acca_', scop)

    # Oxidoreductase domain
    core += get_domains_in_same_superfamily('d1soxa3', scop)

    # AF oxidoreductase domain
    core += get_domains_in_same_superfamily('d1aorb2', scop)

    #
    # Category 2: gregarous folds
    #
    gregarious = []
    # Ribosomal protein L25-like
    gregarious += get_domains_in_same_superfamily('d1d6ka_', scop)

    # BtrG-like
    gregarious += get_domains_in_same_superfamily('d1vkba_', scop)
    
    # RNA-polymerase
    gregarious += get_domains_in_same_superfamily('d1i3qb_', scop)

    # Hypothetical protein HP Ym1108w
    gregarious += get_domains_in_same_superfamily('d1n6za_',scop)

    # QueA-like
    gregarious += get_domains_in_same_superfamily('d1vkyb_',scop)

    #
    # Category 3: structural drift
    #
    drift = [] 

    # NB as per notes to Table 1 in Shi et al (2007), note that in Category 3
    # (structural drift) proteins some SSEs of the beta-grasp motif are not
    # part of the domain core, so not all superfamily members may hvae the
    # beta-grasp motif. Hence we do NOT get all domains in the superfamily
    # of the respresentative in the table, only the representative itself.


    # 4Fe-4S ferredoxins
    drift.append(scop.getDomainBySid('d1h0hb_'))
    drift.append(scop.getDomainBySid('d1vlen2')) # prosmos found, manually checked
    drift.append(scop.getDomainBySid('d1kqfb1')) # prosmos found, manually checked
    

    # TIM barrel - Enolase C-domain-like
    drift.append(scop.getDomainBySid('d1e9ia1'))
    drift.append(scop.getDomainBySid('d1yela1')) # prosmos, manually checked
    drift.append(scop.getDomainBySid('d1iyxa1')) # prosmos, manually checked
    drift.append(scop.getDomainBySid('d2fyma1')) # prosmos, manually checked


    # TIM barrel - (Trans)glycosidases
    drift.append(scop.getDomainBySid('d1fhla_'))
    drift.append(scop.getDomainBySid('d1hjsa_')) # prosmos, manually checked
    drift.append(scop.getDomainBySid('d1ur4a_')) # prosmos, manually checked
    drift.append(scop.getDomainBySid('d1hjqa_')) # prosmos, manually checked


    # alpha/alpha-Toroid - six-hairpin glycosyltransferases
    drift.append(scop.getDomainBySid('d1ut9a1'))

    # Metal ATPase domain - Metal cation-transporting ATPase
    drift.append(scop.getDomainBySid('d1su4a3'))
    drift.append(scop.getDomainBySid('d1q3ia_')) # prosmos, manually checked
    drift.append(scop.getDomainBySid('d1wpga3')) # prosmos, manually checked (dubious on this one though)

    
    # Cystatin-like - NTF2-like
    drift.append(scop.getDomainBySid('d1e3va_'))
    drift.append(scop.getDomainBySid('d1q40a_')) # prosmos, manually checked
    drift.append(scop.getDomainBySid('d1q40b_')) # prosmos, manually checked
    drift.append(scop.getDomainBySid('d1jkgb_')) # prosmos, manually checked
    drift.append(scop.getDomainBySid('d1nwwa_')) # prosmos, manually checked
    drift.append(scop.getDomainBySid('d1tuha_')) # prosmos, manually checked
    drift.append(scop.getDomainBySid('d1s5aa_')) # prosmos, manually checked
    drift.append(scop.getDomainBySid('d2a15a1')) # prosmos, manually checked


    # Knottins - Growth factor receptor domain
    drift.append(scop.getDomainBySid('d1igra3'))

    return core + gregarious + drift


def is_true_positive(domainid, goldstd_pos_dict):
    """
    Return true if domainid is in positive class according to supplied
    dict (i.e. in the dict) - need a function not just simple dict
    lookup due to dodgy stuff with domains not starting with 'd' when
    using Bio.SCOP (seem comments in code for details).
    
    Parameters:
        domainid -domain identifier (sid) to test
        goldstd_pos_dict - dict { domainid : True } for domains in pos class
    Return value:
        True if domainid is in positive class according to supplied dict
    """

    if not goldstd_pos_dict.has_key(domainid):
        # some domains seem to start with g rather than d
        # e.g. g1avo.4 or g1dk9.1 
        # they are in the astral-scopdom-seqres-all-1.73.fa file
        # with the g at the start, but for some reason get None when
        # looking up that sid, need to replace g with d instead.
        if domainid[0] == 'g':
            domainid = 'd' + domainid[1:]
        # when using sids from older versions of ASTRAL SCOP
        # (e.g. 1.65 as used in Hulsen et al. (2006) sometimes
        # we get identifiers with _ as chainid e.g. d1ab4__ which
        # have been changed to have a as chainid instead as per recent PDB
        # remediation. So we'll try chainid a instead of _
        if domainid[5] == '_':
            if len(domainid) < 7:
                sys.stderr.write('WARNING: Bad domain id %s\n' % domainid)
                return False # bizarrely, we get 'd1orf_' from swsse2 sometimes
            domainid = domainid[:5] + 'a' + domainid[6]
        return goldstd_pos_dict.has_key(domainid)
    else:
        return True
    


def tsevalutils_set_verbose(verb):
    """
    set the module global verbose flag in this module to supplied value
    Parameters: verb - True (for verbose output) or False
    Return value: None
    Uses globals: verbose (in this module)
    """
    global verbose
    verbose = verb
