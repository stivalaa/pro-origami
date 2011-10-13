#!/usr/bin/env python
###############################################################################
#
# traindomfc.py - Run domainfc.py with varying threshold to find optimum value
#
# File:    traindomfc.py
# Author:  Alex Stivala
# Created: December 2007
#
# $Id$
#
###############################################################################

"""
Run domainfc.py with the distance threshold at varying values
in order to find the best value on the training dataset.

A table is written to stdout showing the average score overall and for
each number of domains for all the thresholds tested.
This table has a header and is in a format so it can be read in R
with e.g.
   trainscores <- read.table('traindomfc_m_467.out',header=TRUE)
"""

import warnings # so we can suppress the annoying tempnam 'security' warning
import sys,os

from parsepdom import *
from domeval import *
from domainfc import get_domainfc_domains


#-----------------------------------------------------------------------------
#
# Constants
#
#-----------------------------------------------------------------------------

PDBROOT = "/local/charikar/pdb/pdb"

#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------

def main():
    """
    main for traindomfc.py
    """
    
    if len(sys.argv) != 2:
        sys.stderr.write('Usage: ' + os.path.basename(sys.argv[0])
                         + ' pdomains_benchmark_file\n')
        sys.exit(1)

    pdbroot = PDBROOT
    pdomains_filename = sys.argv[1]
    MAX_NUM_DOMAINS = 8
    sys.stdout.write(" t    all   ")
    for num_domains in range(1, MAX_NUM_DOMAINS+1):
        sys.stdout.write("dom%d  " % num_domains)
    sys.stdout.write("\n")

    for dist_threshold in range(3, 30):
        (num_undercut, num_overcut, num_correct, avgscore,
         num_correct_assign, numdomains_dict, num_processed) = \
                       run_on_pdomains_file(pdbroot, pdomains_filename,
                                            False,
                                            get_domainfc_domains,
                                            dist_threshold,
                                            True,   # use MCL
                                            'pdb',  # secondary structure
                                            False, False, # use dot, use neato
                                            True,   # include 310 helices
                                            True,   # include pi helices
                                            False,  # use edgeweights
                                            True,   # add loop nodes
                                            True)  # use modularity cut
        sys.stdout.write("%5.2f %4.3f " % (dist_threshold, avgscore))
        for num_domains in range(1, MAX_NUM_DOMAINS+1):
            try:
                (freq,total,avg,undercut,overcut,
                 domain_num_correct, domain_assign_correct) = \
                 numdomains_dict[num_domains]
                sys.stdout.write("%4.3f " % avg)
            except:
                sys.stdout.write("NA    ")
        sys.stdout.write("\n")

            
if __name__ == "__main__":
    warnings.filterwarnings('ignore', 'tempnam', RuntimeWarning) 
    main()
