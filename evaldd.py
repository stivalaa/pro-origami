#!/usr/bin/env python
###############################################################################
#
# evaldd.py - Evaluate a domain decomposition method against a reference method
#
# File:    evaldd.py
# Author:  Alex Stivala
# Created: December 2007
#
# $Id$
# 
###############################################################################

"""
This script evaluates one domain decomposition method using another
one (usually a database rather than a 'method' e.g. CATH) as the
reference or 'gold standard'.

The actual processing/parsing of the domain method/db is in ptdomain.py
(which is also used by ptgraph2.py and domainfc.py (so far)
and the domain overlap measure is implemented in domeval.py

See usage in docstring for main()
"""
import warnings # so we can suppress the annoying tempnam 'security' warning
import sys,os
import getopt
import re
from time import strftime,localtime
import resource # for getrusage()

from Bio.PDB import *  # only needed for DDomain when segment spans chains

from ptdomain import *
from ptutils import cleanup_tmpdir
import domeval
from getdomains import *

#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

    
def get_domains_wrapper(pdbid, pdb_filename, pdb_struct, chainid,
                        domain_program):
    """
    Get domain decomposition in the form of list of PTDomain objects
    from the nominated domain program or database.

    This is just a wrapper for the getdomains.py get_domains() function
    that makes domain_program last arg not first so it matches the
    prototype for passing to run_on_pdomains_file (domeval.py).
    TODO: should just change the get_domains() function to have the
    right parameter order instead.
    
    Parmeters:
       pdbid - PDB identifier for the protein
       pdb_filename - name of PDB file
       pdb_struct - Bio.PDB parsed structure (needed for DDOMAIN)
       chainid - (default None). If not None, the chainid for single chain
                 to use.
       domain_program - 'cath:cdf_file_name' (CATH) or other supported
                          program or database

    Return value: list of domains according to the chosen method

    """
    return get_domains(domain_program, pdbid, pdb_filename, pdb_struct, chainid)


def usage(progname):
    """
    Print usage message and exit
    """
    
    sys.stderr.write("Usage: " +progname + " [-vs] <test_domain_method> <ref_domain_method> PDBfile OR\n")
    sys.stderr.write("       " +progname + " [-vs] -b <pdomains_file> <test_domain_method> <PDBroot>\n")
    sys.stderr.write("  -v print verbose debugging messages to stderr\n")
    sys.stderr.write("  -s show both test and ref domain decompositions\n")
    sys.stderr.write(" valid domain programs are: " + str(valid_domain_programs)
                     + "\n")
    sys.stderr.write("  -b is an alternate mode that runs over all chains in \n"
                     "     the specified pDomains file with PDBroot as the \n"
                     "     root of the divided PDB hierarchy\n")

    sys.exit(1)
    

#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def main():
    """
    main for evaldd.py

    Usage: evaldd [-vs] <test_domain_method> <ref_domain_method> PDFile OR
           evaldd [-v] -b <pdomains_file> <test_domain_method> <PDBroot>

    test and ref domain methods are the domain programs or files supported
    by ptdomain.py (ddomain, cath:cdf_filename, etc.)

    -s prints the domain decompositions to stdout along with the scores
    -v turns on debug output to stderr
    -b Instead of specifying PDBfile and optional chainid, if the -b mode
       is used the process is run over all chains in the specified pDomains
       benchmark file using the PDBroot as the root of the PDB divided
       hierarchy to get PDB files.

    The score (as per domeval.py) is printed to stdout.
    """
    global verbose

    print_domains = False
    
    valid_domain_programs_re = [ re.compile(re_str) for re_str in
                                 valid_domain_programs ]

    benchmark_mode = False
        
    try:
        opts,args = getopt.getopt(sys.argv[1:], "b:sv?")
    except:
        usage(os.path.basename(sys.argv[0]))
    for opt,arg in opts:
        if opt == "-b":  # use all chains in pDomains benchmark file
            benchmark_mode = True
            pdomains_filename = arg
        elif opt == "-v": # verbose
            verbose = True # this module only
            ptdomain_set_verbose(True) # ptdomain module
        elif opt == "-s": # show domain decompositions
            print_domains = True
        else:
            usage(os.path.basename(sys.argv[0]))


    if benchmark_mode:
        if len(args) != 2:
            usage(os.path.basename(sys.argv[0]))
        if print_domains:
            sys.stderr.write("-s option not valid with -b\n")
            usage(os.path.basename(sys.argv[0]))
        test_domain_method = args[0]
        pdbroot = args[1]
        found_valid_method = False
        for valid_domarg_re in valid_domain_programs_re:
            if valid_domarg_re.match(test_domain_method):
                found_valid_method = True
                break
        if not found_valid_method:
            sys.stderr.write("valid values for domain parsing methods are: " +
                             str(valid_domain_programs) + "\n")
            usage(os.path.basename(sys.argv[0]))
        timestamp = strftime("%d%b%Y %H:%M:%S", localtime())
        sys.stdout.write("# Generated on " + timestamp
                         + " by evaldd.py $Revision$\n"
                         + "# " + " ".join(sys.argv) + "\n")
        (num_undercut, num_overcut, num_correct, avgscore,
         num_correct_assign, numdomains_dict, num_processed) = \
                       domeval.run_on_pdomains_file(pdbroot, pdomains_filename,
                                            True,
                                            get_domains_wrapper,
                                            test_domain_method)
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
             domeval.print_scores(freq, undercut, overcut, domain_num_correct,
                                  domain_assign_correct, avg, indent=2)
             sys.stdout.write('\n')
        sys.stdout.write("\nTotals:\n")
        domeval.print_scores(num_processed,num_undercut, num_overcut,
                             num_correct,
                             num_correct_assign, avgscore, indent=0)
        exit_status = 0
    else:
        if len(args) != 3:
            usage(os.path.basename(sys.argv[0]))
        test_domain_method = args[0]
        ref_domain_method = args[1]
        found_valid_method = False
        for domprog in [test_domain_method, ref_domain_method]:
            for valid_domarg_re in valid_domain_programs_re:
                if valid_domarg_re.match(domprog):
                    found_valid_method = True
                    break
        if not found_valid_method:
            sys.stderr.write("valid values for domain parsing methods are: " +
                             str(valid_domain_programs) + "\n")
            usage(os.path.basename(sys.argv[0]))

        pdb_filename = args[2]
        # check for compressed files. We only support gzip (.gz)
        # Note we are not using the zlib or GzipFile python modules
        # since we are calling to external programs which require the
        # file uncompressed themsevles anyway so we'll just run gzip
        # to uncompress the file to a temporary directory.
        pdb_file_basename = os.path.basename(pdb_filename)
        (name,extension) = os.path.splitext(pdb_file_basename)
        pdbid = name.upper()
        if len(pdbid) >= 6 and pdbid[:3] == "PDB":
            pdbid = pdbid[3:7]
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
            exit_status = 0
            # parse PDB file - only needed for DDOMAIN when segment spans chains
            pdb_parser = PDBParser()
            pdb_struct = pdb_parser.get_structure(pdbid, our_pdb_filename)

            test_domains = get_domains(test_domain_method, pdbid,
                                       our_pdb_filename, pdb_struct)
            ref_domains = get_domains(ref_domain_method, pdbid,
                                      our_pdb_filename, pdb_struct)
            if print_domains:
                print test_domain_method
                write_domains(sys.stdout, test_domains)
                print ref_domain_method
                write_domains(sys.stdout, ref_domains)
            print domeval.domain_eval(test_domains, ref_domains)
        except NotInCATH_Exception,ex_pdbid:
            sys.stderr.write(str(ex_pdbid) + " not found in CATH CDF file\n")
            exit_status = 1
        finally:
            if used_tmp_file:
                cleanup_tmpdir(TMPDIR)

    sys.exit(exit_status)

            
if __name__ == "__main__":
    warnings.filterwarnings('ignore', 'tempnam', RuntimeWarning) 
    main()
