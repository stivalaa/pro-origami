#!/usr/bin/env python
###############################################################################
#
# buildtableauxdb.py - build a database of protein tableaux
#
# File:    buildtableauxdb.py
# Author:  Alex Stivala
# Created: May 2008
#
# $Id$
#
###############################################################################

"""
Build a database of protein tableaux from either the divided PDB hierarchy
(compressed files e.g. pdb1qlp.ent.gz or uncompressed)
or the ASTRAL pdb-style hierarchy
(uncompressed files e.g. d1qlpa_.ent).

Also used to build a database of SSE axis midpoint distance matrices.

See docstring for main() for usage.

The database format is a hash table mapping PDB identifiers to a list
of tableaux in packed format (only storing one triangle, only tableau
codes, indexed by SSE sequential number so no PTNode objects etc).
that is saved using Python pickle module so can be loaded straight into
another Python program.
Ie the intention is that the whole db sits in memory (as it does when
built here).

Refer to tableaubuild.py and pttableau.py for details about tableaux.

It is written in Python and depends on some Python libraries:

. BioPython (including Bio.PDB)
  http://www.biopython.org

  Reference for Bio.PDB is:
  Hamelryck and Manderick 2003 "PDB parser and structure class implemented
  in Python" Bioinformatics 19:2308-2310

  which in turn depends on Numeric
  http://sourceforge.net/projects/numpy


Developed on Linux 2.6.9 (x86_64) with Python 2.5.1
and BioPython 1.43 with Numeric 24.2

"""

import warnings # so we can suppress the annoying tempnam 'security' warning
import sys,os
import re
import getopt
import pickle

import numpy.oldnumeric as Numeric

import ptsecstruct
from ptnode import ptnode_set_verbose
from ptdomain import *
from ptutils import cleanup_tmpdir
import getdomains
from tableaubuild import get_tableaux
from pttableau import PTTableauPacked
from ptversion import get_version


#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------




def build_db(input_root,
             secstruct_program = 'dssp',
             domain_program = 'none',
             include_310_helices = True,
             include_pi_helices = True,
             min_sse_len = None,
             use_numeric = False,
             use_hk  = False,
             build_dist_matrices = False):
    """
    Build the tableaux db in memory.

    Parameters:
       input_root - root of PDB or ASTRAL pdbstyle divided hierarchy
       secstruct_program - secondary structure definition program
                       ('stride' or 'dssp' or 'pdb') to use.
       domain_progam - domain decompositino method ('ddomain','cath', etc.)
       include_310_helices - if True, include 3_10 helices in the graph
       include_pi_helices - if True, include pi helices in the graph
       min_sse_len - if not None, minimum SSE length
       use_numeric - if True build database of Numeric array (Omega matrix)
                      rather than PTTableauPacked
       use_hk - If True build database with HH and KK codes for strands in
                same sheet.
       build_dist_matrices - If True build database of SSE axis midpoint
                distance matrices rather than tableaux.

    Return value:
       dict of { pdbid : [PTTableauPacked list] } for use_numeric=False OR
               { pdbid : [Numeric.array list ]}   for use_numeric=True OR
               { pdbid : [Numeric.array list ]}   for build_dist_matrices=True

                pdbid is e.g. 1QLP when built from PDB or
                         e.g. 1qlp1a when built from ASTRAL pdbstyle
    
    """
    tableau_db = {} # dict of { pdbid : [PTTableauPacked list] }, OR
                    # { pdbid : [Numeric.array list ]}  for use_numeric=True    
                    #                               or build_dist_matrices=True
                    # pdbid is e.g. 1QLP when built from PDB or
                    #          e.g. 1qlp1a when built from ASTRAL pdbstyle

    key_count = 0
    tableaux_count = 0
    keyerror_count = 0
    file_count = 0
    for root,dirs,files in os.walk(input_root):
        for pdb_filename in [filename for filename in files
                             if (os.path.splitext(filename)[1] == '.ent'
                                 or os.path.splitext(filename)[1] == '.pdb'
                                 or os.path.splitext(filename)[1] == '.gz')]:
            sys.stderr.write('processing ' + pdb_filename + '\n')
            file_count += 1
            (pdbid, tableaux_list, sse_string_list) = \
                    get_tableaux(os.path.join(root, pdb_filename),
                                    secstruct_program, domain_program,
                                    include_310_helices, include_pi_helices,
                                    None,  # sse_id_list
                                    min_sse_len,
                                    use_numeric,
                                    use_hk,
                                    build_dist_matrices)
            if not (use_numeric or build_dist_matrices):
                tableaux_list = [PTTableauPacked(tableau) for tableau
                                 in tableaux_list]
            if tableau_db.has_key(pdbid):
                sys.stderr.write("ERROR: duplicate key " + pdbid + "\n")
                keyerror_count += 1
            else:
                tableau_db[pdbid] = tableaux_list
                key_count += 1
                tableaux_count += len(tableaux_list)
    sys.stdout.write("processed %d files\n" % file_count)
    sys.stdout.write("resulting in %d db entries\n" % key_count)
    if build_dist_matrices:
        sys.stdout.write("             %d SSE distance matrices\n" % tableaux_count)
    else:
        sys.stdout.write("             %d tableaux\n" % tableaux_count)
    sys.stdout.write("with %d duplicate key errors\n" % keyerror_count)

    return tableau_db


#-----------------------------------------------------------------------------
#
# Main
#
#-----------------------------------------------------------------------------


def usage(progname):
    """
    Print usage message and exit
    """
    sys.stderr.write("Usage: " + progname +
            " [-35knv] [-d] [-m min_sse_len] [-t struct_prog] "
            "[-p domain_prog]  <pdbroot> <dbname>\n")
    sys.stderr.write("  -3 include 3_10 helices\n")
    sys.stderr.write("  -5 include pi helices\n")
    sys.stderr.write("  -d build SSE distance matrices not tableaux\n")
    sys.stderr.write("  -k use HH and KK codes for anti/parallel strands in same sheet\n")
    sys.stderr.write("  -m  min_sse_len : minimum SSE length\n")
    sys.stderr.write("  -n use numeric values (Omega matrix) rather than tableua\n")
    sys.stderr.write("  -p domain decomposition method/db\n"
                     "     valid values are none (default), "
                     "ddomain, cath:cdffile, pdomains:pdomainsfile\n")
    sys.stderr.write("  -t struct_prog : use struct_prog define " \
                     "secondary structure\n")
    sys.stderr.write("       supported is 'pdb' or 'stride' or 'dssp' (default)\n")
    sys.stderr.write("  -v print verbose debugging messages to stderr\n")
    sys.stderr.write("\nWARNING: dbname is overwritten if it exists\n")
    sys.exit(1)


def main():
    """
    main for buildtableauxdb.py

    Usage: pytableaucreate [-35knv] [-d] [-m min_sse_len ]
             [-t structprog] [-p domainprog]
             <pdbroot> <dbname>


    WARNING: dbname is overwritten if it exists
    
    -3 specifies to include 3_10 helices in the diagram. Default is only
       alpha helices.

    -5 specifies to include pi helices in the diagram. Defaul is only
       alpha helices.

    -d build SSE distance matrices not tableaux

    -k use the HH and KK codes for respectively antiparallel and parallel
       strands in the same sheet, rather than the O, P etc. codes.

    -m min_sse_len :  specifies the minimum SSE length to include in tableaux.

    -n use numeric values (Omega matrix) rather than tableau.

    -p specify the domain decomposition method.
       Valid values are 'none' (default), 'ddomain', 'cath:cdf_filename'.

    -t specifies the secondary structure assignment program to use.
       Currently suppoed is 'pdb' and 'dssp' and 'stride' and 'pmml'.
       Default 'pdb'.

    -v specifies verbose mode: debugging output is written to stderr.
    """
    global verbose
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "35dknm:p:t:v?")
    except getopt.GetoptError:
        usage(os.path.basename(sys.argv[0]))

    valid_secstruct_programs = ["dssp", "stride", "pdb", "pmml"]
    valid_domain_programs = getdomains.valid_domain_programs + [r"none"]
    valid_domain_programs_re = [ re.compile(re_str) for re_str in
                                 valid_domain_programs ]

    verbose = False # global (python globals are only 'global' to module though)
    secstruct_program = "dssp"
    include_310_helices = False
    include_pi_helices = False
    use_hk = False
    domain_program = "none"
    min_sse_len = None
    use_numeric = False
    build_dist_matrices = False

    for opt,arg in opts:
        if opt == "-3":   # include 3_10 helices
            include_310_helices = True
        elif opt == "-5": # include pi helices
            include_pi_helices = True
        elif opt == "-d": # build distance matrices not tableaux
            build_dist_matrices = True
        elif opt == "-k": # use HH and KK codes
            use_hk = True
        elif opt == "-m": # min sse length
            min_sse_len = int(arg)
        elif opt == "-n": # use numeric values (Omega matrix)
            use_numeric = True
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
        elif opt == "-t":
            if arg not in valid_secstruct_programs:
                sys.stderr.write("valid values for -t are: " +
                                 str(valid_secstruct_programs) + "\n")
                usage(sys.argv[0])
            secstruct_program = arg
        elif opt == "-v": # verbose
            verbose = True # this module only
            ptnode_set_verbose(True) # ptnode module
            ptsecstruct.ptsecstruct_set_verbose(True) # ptsecstruct module
            ptdomain_set_verbose(True) # ptdomain module
        else:
            usage(sys.argv[0])
        
    if len(args) != 2:
        usage(os.path.basename(sys.argv[0]))

    if build_dist_matrices:
        if use_numeric:
            use_numeric = False
            sys.stderr.write("WARNING: -n (numeric) ignored for -d (distance matrix)\n")
        if use_hk:
            sys.stderr.write("-k (use HH and KK) invalid for -d (distance matrix)\n");
            usage(sys.argv[0])

    if use_numeric and use_hk:
        sys.stderr.write("-n (numeric) and -k (use HH and KK codes) are "
                         "mutually exlusive\n")
        usage(sys.argv[0])

    sys.stdout.write(sys.argv[0] + ': version is: ' + get_version() + '\n')
    sys.stdout.write(sys.argv[0] + ': options are: ' +  str(sys.argv[1:]) + '\n')
                     
    input_root = args[0]
    output_filename = args[1]

    fh = open(output_filename, 'w')
    tableau_db = build_db(input_root, secstruct_program, domain_program,
                          include_310_helices, include_pi_helices,
                          min_sse_len, use_numeric, use_hk, build_dist_matrices)
    if build_dist_matrices:
        sys.stdout.write('writing SSE distance matrix db to ' + output_filename + '...\n')
    else:
        sys.stdout.write('writing tableaux db to ' + output_filename + '...\n')
    pickle.dump(tableau_db, fh)
    fh.close()
    sys.stdout.write('done.\n')

if __name__ == "__main__":
    warnings.filterwarnings('ignore', 'tempnam', RuntimeWarning) 
    main()


