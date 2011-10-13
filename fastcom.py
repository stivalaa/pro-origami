#!/usr/bin/env python
###############################################################################
#
# fastcom.py - Wrapper to run FastCommunity and get resulting communities
#
# File:    fastcom.py
# Author:  Alex Stivala
# Created: December 2007
#
# $Id$
#
# Usage: fastcom.py [pairsfile]
#
#    pairsfile is the input .pairs file describing the graph in terms
#    of (one per line) tab-delimited pairs of node numbers (from 0)
#    indicating undirected edges, as described in the FastCommunity
#    documentation (see website below).
#    If not specified, then input is from stdin
#
#    Output is to stdout.
#  
# This is a wrapper around the FastCommunityFH program from Aaron Clauset
# that hides the necessity to run it twice, the first time to get the step
# number with max modularity and second to actually get the groups.
#
# FastCommunity is available from
#
# http://cs.unm.edu/~aaron/research/fastmodularity.htm
#
# where the usage and file formats are also described.
# The relevant publication to cite is
#
# A. Clauset, M.E.J. Newman and C. Moore, "Finding community structure
# in very large networks."  Phys. Rev. E 70, 066111 (2004).
#
###############################################################################


import warnings # so we can suppress the annoying tempnam 'security' warning
import sys,os,glob
import getopt



#-----------------------------------------------------------------------------
#
# Constants
#
#-----------------------------------------------------------------------------

FASTCOMMUNITY = "FastCommunityMH" # binary to run (in PATH)

#-----------------------------------------------------------------------------
#
# Functions
#
#-----------------------------------------------------------------------------

def cleanup_tmpdir(tmpdir):
    """
    Remove a temporary directory and its contents
    Parameters:
       tmpdir - temporary directory to remove
    Return value: None
    """
    try:
        for filename in glob.glob(os.path.join(tmpdir, "*")):
            os.remove(filename)
        os.rmdir(tmpdir)
    except OSError, inst:
        sys.stderr.write('WARNING: could not remove temp files'
                         ' in ' + tmpdir + '\n' + str(inst) + '\n')

def usage(progname):
    """
    Print usage message and exit
    """
    sys.stderr.write("Usage: " + progname + " <pairsfile>\n")
    sys.exit(1)

#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------
    
def main():
    """
    main for fastcom.py
    """
    if len(sys.argv) > 2:
        usage(os.path.basename(sys.argv[0]))
    elif len(sys.argv) == 2:
        pairsfilename = sys.argv[1]
        use_stdin = False
    else:
        use_stdin = True
    
    # We need to make a temp directory to work in, 
    # as FastCommunityMH insists on putting all
    # work files in directory of input file (-o option although documented
    # is not in code).
    # If input is a file then symlink to it (note, must end in .pairs)
    # else write our stdin to a file (ending with .pairs)
    TMPDIR = os.tempnam(None, "fc")
    os.mkdir(TMPDIR)
    os.chdir(TMPDIR)
    try:
        oldcwd = os.getcwd()
        if use_stdin:
            pairsfile_basename = 'mystdin.pairs'
            (pairsfile_prefix, suffix) = os.path.splitext(pairsfile_basename)
            fh = open(pairsfile_basename, 'w')
            for line in sys.stdin:
                fh.write(line)
            fh.close()
        else:
            pairsfile_basename = os.path.basename(pairsfilename)
            (pairsfile_prefix, suffix) = os.path.splitext(pairsfile_basename)
            symlink_path = os.path.join(TMPDIR, pairsfile_basename)
            os.symlink(os.path.abspath(pairsfilename), symlink_path)

        firstlabel = "first"   # label for first run
        secondlabel = "second" # and for second run

        # Run FastCommunity the first time to .info file and to get max Q step
        os.system(FASTCOMMUNITY + " -f " + pairsfile_basename +
                  " -l " + firstlabel + " >/dev/null")
        info_filename = pairsfile_prefix + "-fc_" + firstlabel + ".info"
        for line in open(info_filename):
            if line[:11] == "STEP------:":
                step = int(line.split()[1])
                break

        # Now we have the max Q step, use it on the -c option on second run
        os.system(FASTCOMMUNITY + " -f " + pairsfile_basename +
                  " -l " + secondlabel + " -c " + str(step) + " >/dev/null")
        groups_filename = pairsfile_prefix + "-fc_" + secondlabel + ".groups"
        # and just output the groups file to stdout as the result
        for line in open(groups_filename):
            sys.stdout.write(line)
        
    finally:
        if not use_stdin:
            os.unlink(symlink_path)
        os.chdir(oldcwd)
        cleanup_tmpdir(TMPDIR)


if __name__ == "__main__":
    # tmpdir() annoyingly gives 'security' warning on stderr, as does
    # tmpnam(), unless we add these filterwarnings() calls.
    warnings.filterwarnings('ignore', 'tempdir', RuntimeWarning)
    warnings.filterwarnings('ignore', 'tempnam', RuntimeWarning) 
    main()
