#!/usr/bin/env python
###############################################################################
#
# pytableaucreate - Python implementation of protein Tableau creator
#
# File:    pytableaucreate.py
# Author:  Alex Stivala
# Created: February 2008
#
# $Id$
#
#
# Create a protein tableau and write it to stdout.
# The implemntation is actually in pttableau.py which is used by ptgraph2.py
# (Pro-Origami), this is basically just a wrapper for testing / standalone
# tableau creation (see pttableau.py).
#
# Also used to create SSE midpoint distance matrix.
#
# Tableaux are described by Kamat and Lesk 2007
# 'Contact Patterns Between Helices and Strands of Sheet Define Protein
#  Folding Patterns' Proteins 66:869-876
# and Lesk 2003 'From Electrons to Proteins and Back Again'
# Int. J. Quant. Chem. 95:678-682
# and Lesk 1995 'Systematic representation of folding patterns'
# J. Mol. Graph. 13:159-164.
#
# The implementation is based on Arun Konagurthu's TableauCreator program, see
# Konagurthu, Stuckey and Lesk 2008 'Structural search and retrieval using
# a tableau representation of protein folding patterns' Bioinformatics
# (advance access, to be published Jan 5 2008).
#
# Example usage:
#
#  pytableaucreate.py 1QLP.pdb
#
# Filenames may be either in the format above or the pdbq1lp.pdb format.
# Compressed pdb files are supported (gzip) (e.g. pdb1qlp.ent.gz).
#
# It is written in Python and depends on some Python libraries:
#
# . BioPython (including Bio.PDB)
#   http://www.biopython.org
#
#   Reference for Bio.PDB is:
#   Hamelryck and Manderick 2003 "PDB parser and structure class implemented
#   in Python" Bioinformatics 19:2308-2310
#
#   which in turn depends on Numeric
#   http://sourceforge.net/projects/numpy
#
#
# Developed on Linux 2.6.9 (x86_64) with Python 2.5.1
# and BioPython 1.43 with Numeric 24.2
#
###############################################################################


import warnings # so we can suppress the annoying tempnam 'security' warning
import sys,os
import getopt
import re
import pickle
import random
import copy
from math import degrees
import numpy.oldnumeric as Numeric
from Bio.PDB import *

import ptsecstruct
from ptnode import ptnode_set_verbose
from ptdomain import *
from ptutils import cleanup_tmpdir,isNaN
import getdomains
from tableaubuild import TableauBuild,make_tableaux
from pttableau import PTTableauPacked

#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

def write_tableau(n, tableau, permutation,  use_numeric,
                  fortran_format, build_distance_matrix):
    """
    Write tableau or distance matrix to stdout.

    n - order of tableau or distance matrix (n by n)
    tableau - PTTableau object for tableau or Numeric matrix for
              Omega matrix or Numeric matrix for distance matrix
    permutation - permuted list of ingeters in interval [0, n-1] to
                  permute the rows+cols of the tableau/matrix by
                  (so [0,1,2,...n-1] for no permutation.
    use_numeric - boolean. If true, tableau is a Numeric Omega matrix
                  not a tableau.
    fortran_format - boolean. If True, put in lower triangle format
                     for FORTRAN programs tsrchd etc.
    build_distance_matrix - boolean. If True is a distance matrix not a
                      tableau or Omega matrix.

    """
    if build_distance_matrix:
        distmatrix = tableau
        if fortran_format:
            for k in range(n):
                for l in range(k+1):
                    kprime = permutation[k]
                    lprime = permutation[l]
                    if isNaN(distmatrix[kprime,lprime]):
                        dist = 0.0
                    else:
                        dist = distmatrix[kprime,lprime]
                    if dist > 99.9:
                        sys.stderr.write('WARNING: distance %f at (%d,%d) truncated to 99.9 for fortran format\n' % (dist,kprime,lprime))
                        dist = 99.9
                    sys.stdout.write("%6.3f " % dist)
                sys.stdout.write("\n")
        else:
            for k in range(n):
                for l in range(n):
                    kprime = permutation[k]
                    lprime = permutation[l]
                    sys.stdout.write("% 6.2f " % distmatrix[kprime,lprime])
                sys.stdout.write("\n")
    elif use_numeric:
        Omega = tableau
        if fortran_format:
            for k in range(n):
                for l in range(k+1):
                    kprime = permutation[k]
                    lprime = permutation[l]
                    if isNaN(Omega[kprime,lprime]):
                        angle = 0.0
                    else:
                        angle = Omega[kprime,lprime]
                    sys.stdout.write("%6.3f " % angle)
                sys.stdout.write("\n")
        else:
            for k in range(n):
                for l in range(n):
                    kprime = permutation[k]
                    lprime = permutation[l]
                    sys.stdout.write("% 4.3f " % Omega[kprime,lprime])
                sys.stdout.write("\n")
    else:
        if fortran_format:
            for k in range(n):
                for l in range(k+1):
                    kprime = permutation[k]
                    lprime = permutation[l]
                    sys.stdout.write(tableau[(kprime,lprime)] + " ")
                sys.stdout.write("\n")
        else:
            # can't just sys.stdout.write(str(tableau)) if shuffled...
            for k in range(n):
                for l in range(n):
                    kprime = permutation[k]
                    lprime = permutation[l]
                    sys.stdout.write(tableau[(kprime,lprime)] + " ")
                sys.stdout.write('\n')


def write_tableau_old_format(n, Omega, ssestr):
    """
    Write tableau to stdout in the original
    (Arun) TableauCreator format, with angles in degrees,
    full matrix, number of SSEs on first line and SSE sequence
    (DSSP codes E,H) on second line.

    n - order of tableau  matrix (n by n)
    Omega - Numeric matrix for Omega matrix
    sse_str - SSE string correspdonding to the Omega matrix
    """
    sys.stdout.write(str(len(Omega)) + '\n')
    sys.stdout.write(ssestr + '\n')
    for k in range(n):
        for l in range(n):
            angle = degrees(Omega[k,l])
            if isNaN(angle) or k == l:
                angle = -999.0
            sys.stdout.write("% 7.1f " % angle)
        sys.stdout.write("\n")

def write_distmatrix_old_format(n, dmat, ssestr):
    """
    Write  distance matrix to stdout in format for TableauComparer
    (Arun)
    full matrix, number of SSEs on first line and SSE sequence
    (DSSP codes E,H) on second line.

    n - order of  distance matrix (n by n)
    dmat - Numeric matrix for distance matrix
    sse_str - SSE string correspdonding to the distance matrix
    """
    sys.stdout.write(str(len(dmat)) + '\n')
    sys.stdout.write(ssestr + '\n')
    for k in range(n):
        for l in range(n):
            d = dmat[k,l]
            if isNaN(d) or k == l:
                d = -999.0
            sys.stdout.write("% 7.1f " % d)
        sys.stdout.write("\n")

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
            " [-35knuvfe] [-d|-b] [-t struct_prog] "
            "[-p domain_prog] [-a domainid] [-s sse_num_list] [-c chainid] "
            "[-m min_sse_len] [-o <savefile>] [-i identifier] "
            "<PDBfile>\n")
    sys.stderr.write("  -3 include 3_10 helices\n")
    sys.stderr.write("  -5 include pi helices\n")
    sys.stderr.write("  -k use HH and KK codes for anti/parallel strands in same sheet\n")
    sys.stderr.write("  -n output numeric matrix rather than tableau\n")
    sys.stderr.write("  -e output numeric tableau angles in degrees, in original TableauCreator .angles file format\n")
    sys.stderr.write("  -f output in FORTRAN style format for TSRCHN\n")
    sys.stderr.write("  -d build SSE axis midpoint distance matrix not tableau\n")
    sys.stderr.write("  -b build both tableau and distance matrix\n")
    sys.stderr.write("  -p domain decomposition method/db\n"
                     "     valid values are none (default), "
                     "ddomain, cath:cdffile, pdomains:pdomainsfile\n")
    sys.stderr.write("  -a domainid : only output for specified domain\n")
    sys.stderr.write("  -t struct_prog : use struct_prog define " \
                     "secondary structure\n")
    sys.stderr.write("       supported is 'pdb' (default) or 'stride' or 'dssp' or 'pmml'\n")
    sys.stderr.write("  -s sse_num_list : specifies comma-separated list of "
                     "SSE sequential numbers to include in the tableau\n")
    sys.stderr.write("  -m min_sse_len : minimum number of residues in SSE to "
                     "be included in tableau\n")
    sys.stderr.write("  -c chainid : specify chain identifier; only build"
                     "tableau for that chain\n")
    sys.stderr.write("  -i identifier : when using -f, specify identifier "
                     " to use rather than deriving from filename\n")
    sys.stderr.write("  -o savefile : save tableau in packed format for use "
                     "in other programs such as tabsearchqpml.py\n"
                     "  WARNING: savefile is overwritten if it exists.\n")
    sys.stderr.write("  -u randomly permute the tableau/distance matrix\n")
    sys.stderr.write("  -v print verbose debugging messages to stderr\n")
    sys.exit(1)


def main():
    """
    main for pytableaucreate.py

    Usage: pytableaucreate [-35nefuv] [-d|-b] [-t structprog] [-p domainprog]
                [-a domainid]
                [-s sse_num_list] [-c chainid] [-m min_sse_len]
                [-o savefile] <PDBfile>


    -3 specifies to include 3_10 helices in the diagram. Default is only
       alpha helices.

    -5 specifies to include pi helices in the diagram. Defaul is only
       alpha helices.

    -k use the HH and KK codes for respectively antiparallel and parallel
       strands in the same sheet, rather than the O, P etc. codes.
       
    -n output a numeric omega matrix instead of tableau.

    -e output numeric tableau angles in degrees, in the original
       TableauCreator .angles file format, with number of entries on
       first line, SSE sequence description on second line (E/H), then
       (full) matrix with angles in degrees (rather than radians).
       For distance matrix, same format with distances between SSEs
       in Angstroms.

    -f output the matrix in 'FORTRAN style' lower triangle with
       header line suitable for input to TMATN.

    -d build SSE axis midpoint distance matrix rather than tableau.

    -b build both the tableau and distance matrix and output together,
       for use with tsrchd etc. for example. If -u is used to permute
       the matrices, they are permuted the same way so they are still
       consistent.

    -p specify the domain decomposition method.
       Valid values are 'none' (default), 'ddomain', 'cath:cdf_filename'.

    -a domainid : only output specified domain

    -t specifies the secondary structure assignment program to use.
       Currently suppoed is 'pdb' and 'dfh,ssp' and 'stride' or 'pmml'.
       Default 'pdb'.

    -s sse_num_list specifies a comman-separated
       list of SSE sequential ids to build the
       tableau for. SSE sequential id's start at 1 and go from N to C
       terminus. E.g. -s1,5,8 includes only the 1st, 5th and 8ths SSEs.
       Numbers do not restart at chains (but do restart in each domain).
       These nubmers are those assigned by 'ptgraph2 -b sequential' option.

       TODO: this currently does not make sense when multiple domains
       are being procssed, this option applies to each domain.
    
   -c chainid : specify chain identifier; only build tableau for that chain

   -m min_sse_len : minimum nubmer of residues in SSE for it to be included

   -i identifier : when using fortran format (-f), specify the identifier
      to use in the output rather than deriving it from the filename

    -o savefile : save tableau in packed format for use in other
       programs, such as tabsearchqpml.py
       WARNING: savefile is overwritten if it exists

       TODO: this currently does not make sense when multiple domains
       are being procssed, this option only saves first domain.
       
    -u randomly pemute the rows+cols (symmetric) of the tableau/distance matrix.
       writes the permutation vector in form 
       permutation = i,j,..,m
       e.g. 
       permutation = 3,1,2,4
       as first line of output before identifier information and tableau

    -v specifies verbose mode: debugging output is written to stderr.
    """
    global verbose
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "35bdfknep:a:t:s:c:m:i:o:uv?")
    except getopt.GetoptError:
        usage(os.path.basename(sys.argv[0]))

    valid_secstruct_programs = ["dssp", "stride", "pdb", "pmml"]
    valid_domain_programs = getdomains.valid_domain_programs + [r"none"]
    valid_domain_programs_re = [ re.compile(re_str) for re_str in
                                 valid_domain_programs ]

    verbose = False # global (python globals are only 'global' to module though)
    secstruct_program = "pdb"
    include_310_helices = False
    include_pi_helices = False
    domain_program = "none"
    sse_id_list = None
    use_numeric = False
    use_hk = False
    savefilename = None
    min_sse_len = None
    fortran_format = False
    build_distance_matrix = False
    chainid = None
    fident = None
    do_shuffle = False
    build_both = False # both tableau and dist matrix
    use_old_format = False # size + SSE chain + degrees omega matrix
    domainid = None

    for opt,arg in opts:
        if opt == "-3":   # include 3_10 helices
            include_310_helices = True
        elif opt == "-5": # include pi helices
            include_pi_helices = True
        elif opt == "-d":  # build SSE midpoint distance matrix not tableau
            build_distance_matrix = True
        elif opt == "-b": # build both tableau and distance matrix
            build_both = True
        elif opt == "-k": # use HH and KK codes
            use_hk = True
        elif opt == "-n": # output numeric matrix not tableau
            use_numeric = True
        elif opt == "-e": # use TableauCreator .angles file format
            use_old_format = True
        elif opt == "-f":  # FORTRAN style format for TMATN
            fortran_format = True
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
        elif opt == "-a":  # only output tableau for specified domain id
            domainid = arg
        elif opt == "-t":
            if arg not in valid_secstruct_programs:
                sys.stderr.write("valid values for -t are: " +
                                 str(valid_secstruct_programs) + "\n")
                usage(sys.argv[0])
            secstruct_program = arg
        elif opt == "-s":
            sse_id_list_str = arg.split(',')
            sse_id_list = []
            sse_id_uniq_dict = {} # { id : True } just for checking all unique
            for sse_id_str in sse_id_list_str:
                if sse_id_str.isdigit():
                    if sse_id_uniq_dict.has_key(int(sse_id_str)):
                        sys.stderr.write("duplicate SSE sequential number "  +
                                         sse_id_str + "\n")
                        usage(sys.argv[0])
                    sse_id_uniq_dict[int(sse_id_str)] = True
                    sse_id_list.append(int(sse_id_str))
                else:
                    sys.stderr.write("not a valid SSE sequential number '" +
                                     sse_id_str + "'\n")
                    usage(sys.argv[0])
            sse_id_list.sort() # ensure SSEs are in order
        elif opt == "-c": # chain identifier
            if len(arg) != 1:
                sys.stderr.write("invalid chain identifier for -c option\n")
                usage(sys.argv[0])
            chainid = arg.upper()
        elif opt == "-m": # min sse len
            min_sse_len = int(arg)
        elif opt == "-i": # identifier to use for fortran format
            fident = arg
        elif opt == "-o": # save tableau in packed format
            savefilename = arg
        elif opt == "-u": # randomly permute the tableau/matrix
            do_shuffle = True
        elif opt == "-v": # verbose
            verbose = True # this module only
            ptnode_set_verbose(True) # ptnode module
            ptsecstruct.ptsecstruct_set_verbose(True) # ptsecstruct module
            ptdomain_set_verbose(True) # ptdomain module
        else:
            usage(sys.argv[0])

    if use_numeric and use_hk:
        sys.stderr.write("-n (numeric) and -k (use HH and KK codes) are "
                         "mutually exlusive\n")
        usage(sys.argv[0])

    if build_distance_matrix and build_both:
        sys.stderr.write("WARNING: both -d (build dist matrix) and -b "
                         "(build both) specified, ignoring -d\n")
        build_distance_matrix = False

    if savefilename and do_shuffle:
        sys.stderr.write('WARNING: saved tableau will not be shuffled\n')

    if build_distance_matrix:
        if use_numeric:
            use_numeric = False
            sys.stderr.write("WARNING: -n (numeric) ignored for -d (distance matrix)\n")
        if use_hk:
            sys.stderr.write("-k (use HH and KK) invalid for -d (distance matrix)\n");
            usage(sys.argv[0])

    if (secstruct_program == "pmml" and 
        (min_sse_len == None  or min_sse_len < 3)):
        sys.stderr.write("WARNING: PMML can give SSEs of length 1 or 2 causing axis fitting to fail, setting minimum length to 3 as if -m3 were specfified\n")
        min_sse_len = 3

    if fident:
        if not fortran_format:
            sys.stderr.write("-i is only valid with -f\n")
            usage(sys.argv[0])
        elif len(fident) > 8:
            sys.stderr.write("identifier must be 8 chars or less\n")
            usage(sys.argv[0])

    if use_old_format and (build_both or
                           use_hk or use_numeric or fortran_format or
                           do_shuffle or savefilename):
        sys.stderr.write("-e (use old .angles format) is not compatible "
                         "with -b -k or -n or -f or -u or -o\n")
        usage(os.path.basename(sys.argv[0]))
              
    if len(args) != 1:
        usage(os.path.basename(sys.argv[0]))

    pdb_filename = args[0]

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
        if fortran_format and fident:
            pdbid = fident
        else:
            pdbid = name.upper()
            if len(pdbid) >= 6 and pdbid[:3] == "PDB":
                pdbid = pdbid[3:7]
            if chainid:
                pdbid += '_' + chainid

        # parse PDB file
        pdb_parser = PDBParser()
        pdb_struct = pdb_parser.get_structure(pdbid, our_pdb_filename)
        # create the Tableaux and output them
        (tableaux_list, ssestr_list) = make_tableaux(our_pdb_filename,
                                      pdb_struct,
                                      secstruct_program,
                                      domain_program,
                                      include_310_helices,
                                      include_pi_helices,
                                      (use_numeric or use_old_format),
                                      sse_id_list,
                                      use_hk,
                                      min_sse_len,
                                      build_distance_matrix,
                                      chainid,
                                      domainid)
        if build_both:
            (distmatrix_list, ssestr_list) = make_tableaux(our_pdb_filename,
                                            pdb_struct,
                                            secstruct_program,
                                            domain_program,
                                            include_310_helices,
                                            include_pi_helices,
                                            use_numeric,
                                            sse_id_list,
                                            use_hk,
                                            min_sse_len,
                                            True, # build_distance_matrix
                                            chainid,
                                            domainid)
        i = 1
        for tableau in tableaux_list:
            n = len(tableau)
            permutation = range(n) # used to permute rows/cols: null permutation
            if do_shuffle:
                random.shuffle(permutation) # actually permute for shuffle mode
                if verbose:
                    sys.stderr.write('permutation is: ' + str(permutation)+'\n')
                sys.stdout.write('permutation = ' + ','.join([str(x+1) for x in permutation]) + '\n')
            if i > 1:
                sys.stdout.write('\ndomain ' + str(i) + ':\n')

            if fortran_format:
                sys.stdout.write("%7s %4d\n" % (pdbid.upper(), n))

            if use_old_format:
                if build_distance_matrix:
                    write_distmatrix_old_format(n, tableau, ssestr_list[i-1])
                else:
                    write_tableau_old_format(n, tableau, ssestr_list[i-1])
            else:
                write_tableau(n, tableau, permutation, use_numeric,
                              fortran_format, build_distance_matrix)

            if build_both:
                write_tableau(n, distmatrix_list[i-1],
                              permutation, use_numeric,
                              fortran_format, True)
                
            i += 1
    finally:
        if used_tmp_file:
            cleanup_tmpdir(TMPDIR)


    if savefilename:
        if verbose:
            sys.stderr.write('writing tableau to ' + savefilename +'\n')
        fh = open(savefilename, "w")
        if len(tableaux_list) > 1:
            sys.stderr.write('WARNING: only saving first tableau in list\n')
        if build_distance_matrix:
            pickle.dump(distmatrix, fh)
        elif use_numeric:
            # Numeric/numpy seems to have no 'packed' format for symmetric
            # matrices, so we just have to dump the whole thing.
            pickle.dump(Omega, fh)
        else:
            pickle.dump(PTTableauPacked(tableaux_list[0]), fh)
        fh.close()
        
if __name__ == "__main__":
    warnings.filterwarnings('ignore', 'tempnam', RuntimeWarning) 
    main()
