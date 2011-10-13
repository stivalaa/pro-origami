###############################################################################
#
# ptutils.py - Miscellaneous utility functions
#
#
# File:    ptutils.py
# Author:  Alex Stivala
# Created: October 2007
#
# $Id$
#
###############################################################################

import os,sys
import glob

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
    


def get_int_icode(res_seq):
    """
    Return tuple (int, icode) with integer residue sequence number and
    single char icode from PDB residue sequence string such as '60A' or '61'
    etc.

    Parameters:
       res_seq - PDB resisue sequence number string, with or without icode
    Return value:
       tuple (int, icode) where
         int is integer part of res_seq, and icode is char insertino code or None
    """
    if not res_seq[-1].isdigit():
        int1 = int(res_seq[0:len(res_seq)-1])
        icode = res_seq[-1]
    else:
        int1 = int(res_seq)
        icode = None
    return (int1, icode)


def biopdbresid_to_pdbresseq(biopdb_residueid,ignore_insertion_codes=False):
    """
    Give a Bio.PDB Residue id tupe (hetatm, resseqnum, icode), return
    the PDB residue sequence number string consisting of the sequence
    number and the insertion code, if not blank.

    Parameters:
       biopdb_residueid - tuple (hetatm, resseqnum, icode) from Residue.get_id()
            ignore_insertion_codes - If True, a hack to make it work with
                      PMML (only) which does not report insertion codes
                      unlike DSSP and STRIDE
    Return value:
       string residue PDB sequence number e.g. '60' or '60A'.
    """
    # Residue.get_id() gives tuple (hetatm, resseqnum, icode)
    res_seq = str(biopdb_residueid[1])
    if not ignore_insertion_codes:
        if biopdb_residueid[2] != ' ':
            res_seq +=  biopdb_residueid[2]
    return res_seq


def pdb_res_seq_cmp(res_seq1, res_seq2):
    """
    Comparison function for PDB residue sequence numbers, which
    are strings that may have insertion code on end e.g. '60'
    or '60A'. Compare in integer part as integers, if same then
    use insertion code e.g. '60' < '61' regardless of what insertion
    codes we could put on the end, but also '60' < '60A' < '60B' < '61'.

    DANGER: this comparison is not always correct for all PDB files:
            is assumes residues are ordered as above (60, 60A, 60B, 61, etc.)
            but some PDB files do NOT work this way, e.g. 1HVC, where the
            first part of the chain numbered 1 up to 99 all has insertion code
            B, then there are some more residues, then numbering starts
            again at 1 but with insertion code A on all residues!
            
    Parameters:
       res_seq1 - PDB residue sequence string as above
       res_seq2 - PDB resiue seqwuence string as above
    Return value:
       -1 if res_seq1 < res_seq2
        0 if res_seq1 = res_seq2
        1 if res_seq1 > res_seq2,
      according to the ordering defined above.
    Uses no data members.
    """
    (int1, icode1) = get_int_icode(res_seq1)
    (int2, icode2) = get_int_icode(res_seq2)

    if int1 < int2:
        return -1
    elif int1 > int2:
        return 1
    else:  # compare icodes, note None < x for any x
        if icode1 < icode2: 
            return -1
        elif icode1 > icode2:
            return 1
        else:
            return 0

def char_if_not_blank(char):
    """
    IF supplied character is not space, return it, else return empty
    string.
    Paramaters:
       char - char to test
    Return value:
       char if char is not space, else ''.
    """
    if char == ' ':
        return ''
    else:
        return char



def isNaN(x):
    """
    Test if supplied float is an IEEE not-a-number (NaN).
    For some reason Python does not hav a function to do this,
    and nor does Numeric (although numpy and scipy have support for it).
    
    Parameters:
        x - float to test for NaN

    Return value:
        True if x is NaN, else False.
    """
    # NaN is the only float value that is not equal to itself (IEEE
    # standard)
    if x != x:
        return True
    else:
        return False
 

