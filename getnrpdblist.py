#!/usr/bin/env python
###############################################################################
#
# getnrpdblist.py - get a list of nonredundant PDB identifiers
#
# File:    getnrpdblist.py
# Author:  Alex Stivala
# Created: May 2008
#
# See usage in docstring for main()
#
# $Id$
#
###############################################################################

"""
This program is used to get a list of nonredundant PDB identifiers
for eg building a cartoon database with ptgraph2.py, in order not to
have to run on entire PDB.

It uses the cluster file output of CD-HIT, which is used to find
clusters of PDB entries above a certain sequence identity (very quickly).

The references for CD-HIT, which is available from
http://bioinformatics.ljcrf.edu/cd-hi are:

1. Clustering of highly homologous sequences to reduce the size of
   large protein databases. Weizhong Li, Lukasz Jaroszewski & Adam
   Godzik. Bioinformatics (2001) 17:282-283

2. Tolerating some redundancy significantly speeds up clustering of
   large protein databases. Weizhong Li, Lukasz Jaroszewski & Adam
   Godzik. Bioinformatics (2002) 18: 77-82

3. Cd-hit: a fast program for clustering and comparing large sets of
   protein or nucleotide sequences. Weizhong Li & Adam
   Godzik. . Bioinformatics (2006) 22:1658-1659


We parse the output of running CD-HIT over the pdb_seqres.txt file
obtained from ftp://ftp.wwpdb.org/pub/pdb/derived_data
with a command like:

cd-hit -i pdb_seqres.txt -o pdb_seqres.nr95 -c 0.95 -n 5 -s 0.9 -d 80

resulting in the pdb_seqres_nr95.clstr file with the clusters, which
is what this program parses. This file looks like:

>Cluster 0
0	2999aa, >1s1i_3 mol:na length:2999  5.8S/25S ribosomal RNA... *
>Cluster 1
0	2922aa, >1ffk_0 mol:na length:2922  23S RRNA... *
1	2922aa, >1jj2_0 mol:na length:2922  23S RRNA... at 100%

>Cluster 2919
0	378aa, >1k9o_I mol:protein length:378  ALASERPIN... *

where each '>' at start of line is a new cluster, and each entry
in a cluster has fields: sequence length, FASTA identifier/description,
then '*' for representative sequence else 'at x%' showing perc. identity
with representative sequence. (See CD-HIT documentation for more details).

Note we are depending on the format of entries in pdb_seqres.txt (and
hence the .clstr file) here, namely that the first thing on the FASTA
id line is the PDB identifier with underscore and chainid appended.
This also allows us to use the mol: field in the FASTA line to only
select mol:protein since we don't want RNA etc.

CD-HIT version 2007-0131 was used.

Given this clustering, we choose one sequence from each cluster as
the representative (note we do not just use the CD-HIT representative,
we want to choose highest resolution structure, sequence clustering knows
nothing about that):

1. Sort by sequence length and choose the set of highest sequence length
   entries (this will always inlcude the represnetative). This is since
   shorter ones are generally either fragments or are missing residues.

2. Among this reduced set, look in the PDB files and (try to, since this
   information is not in a standardized format) find the method (X-ray
   crystallography or NMR). Choose X-ray crystallography over if possible.

3. Among the remaining entries, choose the one (or several) with highest
   resolution.

4. If still more than one candidate, choose that with highest percent
   identify with CD-HIT representative.

5. If still more than one, choose arbitrary one [TODO: maybe use most
   recent date, but problems with getting that info from header sometimes]

TODO: maybe should try something more sophisticated like the SPACI in
ASTRAL (Brenner et al (2000) Nucleic Acids Res. 28(1):254-256)

Bio.PDB from BioPython (http://www.biopython.org) is required for parsing
PDB files.
"""
import warnings # so we can suppress the annoying tempnam 'security' warning
import os,sys
from Bio import PDB
from ptutils import cleanup_tmpdir


#-----------------------------------------------------------------------------
#
# Constants
#
#-----------------------------------------------------------------------------

# location of divided PDB files (hierarcy e.g. ql/pdb1qlp.ent.gz)
PDBDIV_ROOT = "/local/charikar/pdb/pdb"


#-----------------------------------------------------------------------------
#
# Class definitions
#
#-----------------------------------------------------------------------------

class CDHIT_ParseError(Exception):
    """
    Exception class for parse errors.
    """
    pass

class MethodResolution:
    """
    struct for entries in dictionary keyed by PDB identifier
    """
    is_xray    = False # boolean: EXPDTA is X-RAY DIFFRACTION
    resolution = float("inf")  # float resolution in Ansgstroms if is_xray

class CDHIT_sequence:
    """
    The CDHIT_sequence class represents one entry (sequence) in a CD-HIT cluster
    file.
    """
    def __init__(self, sequencenum, seqlen, descr, is_repr, percent_identity):
        """
        Construct a CDHIT_sequence from information in a single line in .clstr file

        Parameters:
            sequencenum - int sequential id in cluster for this sequence
            seqlen - int sequence length for this sequence
            descr - description string (FASTA descr line) for this sequence
            is_repr - bool True if sequence is the representative for the
                      cluster this sequence is in, else False.
            percent_identity - int percentage sequence identity with
                      reprsentative for this cluster.
        """
        self.sequencenum = sequencenum # int sequential id in cluster for this
                                 # sequence
        self.seqlen = seqlen     # int sequence length for this sequence
        self.descr = descr       # description string (FASTA descr line) for
                                 # this sequence
        self.is_repr = is_repr   # bool True if sequence is the
                                 # representative for the cluster this
                                 # sequence is in, else False.
        self.percent_identity = percent_identity # int percentage
                                                 # sequence identity
                                                 # with representative
                                                 # for this cluster.

        # following are filled in later with PDB information, access directly
        self.is_xray = False     # bool True for exp method X-ray diffraction
        self.resolution = float("inf")   # float resolution in Angstroms if is_xray

    def __str__(self):
        """
        Return string rep of sequence 
        """
        if self.is_repr:
            isrep = "*"
        else:
            isrep = "."
        if self.is_xray:
            method = "X"
        else:
            method = "?"
        if self.resolution < float("inf"):
            resolution = "%3.1f" % self.resolution
        else:
            resolution = "???"
        return "%6d %6d %s %3d %s %s %s" % (self.sequencenum, self.seqlen,
                                      isrep, self.percent_identity,
                                      method, resolution,
                                      self.descr)

    def is_protein(self):
        """
        Return True if this is a protein sequence, else False.
        Note we are depending on the format of entries in pdb_seqres.txt (and
        hence the .clstr file) here, namely the
        the mol: field in the FASTA line to only
        select mol:protein since we don't want RNA etc.

        Parameters:
            None.
        Return value:
            True if sequence is marked as protein in pdb_seqres FASTA descr line
            else False.
        """
        if self.descr.split()[1] == "mol:protein":
            return True
        else:
            return False


class CDHIT_cluster:
    """
    The CDHIT_cluster class represents one cluster from CD-HIT in the
    CD-HIT .clstr file. It is a list of CDHIT_sequence objects, each one
    being a sequence in the cluster.
    """
    def __init__(self, clusternum):
        """
        Construct an empty CDHIT_cluster object from information on the
        >Cluster line, ie nothing byt cluster id yet. Add entries
        for each sequence in the cluster with add_sequence.
        """
        self.clusternum = clusternum
        self.seqlist = []  # list of CDHIT_sequence for this cluster

    def add_sequence(self, cdhit_sequence):
        """
        Add a CDHIT_sequence to this cluster.

        Parameters:
           cdhit_sequence - CDHIT_sequence object to add to list
        Return value:
            None
        Modifies member data:
           seqlist
        """
        self.seqlist.append(cdhit_sequence)

    def __str__(self):
        """
        Return string rep of cluster, as cluster id with each sequence
        one per line.
        """
        s = "cluster " + str(self.clusternum) + " (" + str(len(self.seqlist)) +\
            " sequences)\n"
        for seq in self.seqlist:
            s += str(seq) + "\n"
        return s

    def __len__(self):
        """
        Return the number of sequences in the cluster
        """
        return len(self.seqlist)

    def __getitem__(self, index):
        """
        Return the index'th sequence in this cluster
        """
        return self.seqlist[index]


    def discard_short_seqs(self):
        """
        Discard sequences that are shorter than the longest sequence.

        Parameters:
            None
        Return value:
            None:
        Modifies member data:
            seqlist
        """
        maxseqlen = 0
        for seq in self.seqlist:
            if seq.seqlen > maxseqlen:
                maxseqlen = seq.seqlen
        self.seqlist = [seq for seq in self.seqlist if seq.seqlen == maxseqlen]

    def discard_non_xray(self):
        """
        Discard sequences that are are not EXPDTA X-RAY DIFFRACTION,
        but only if there are any EXPDTA X-RAY DIFFRACTION
        sequences, i.e. don't discard non x-ray if they are the only
        available.

        Parameters:
            None
        Return value:
            None:
        Modifies member data:
            seqlist

        """
        found_xray = False
        for seq in self.seqlist:
            if seq.is_xray:
                found_xray = True
                break
        if found_xray:
            self.seqlist = [seq for seq in self.seqlist if seq.is_xray]


    def discard_lower_resolution(self):
        """
        Discard seqwuences that are lower than the highest resolution.

        Parameters:
            None
        Return value:
            None:
        Modifies member data:
            seqlist
        """
        minres = float("inf")
        for seq in self.seqlist:
            if seq.resolution < minres:
                minres = seq.resolution
        if minres < float("inf"):
            self.seqlist = [seq for seq in self.seqlist if seq.resolution <= minres]

    def discard_lower_similarity(self):
        """
        Discard seqwuences that are lower than the highest similarity
        to (CD-HIT) representative sequence.

        Parameters:
            None
        Return value:
            None:
        Modifies member data:
            seqlist

        """
        maxsim = 0
        for seq in self.seqlist:
            if seq.percent_identity > maxsim:
                maxsim = seq.percent_identity
        self.seqlist = [seq for seq in self.seqlist if seq.percent_identity == maxsim]
            
#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

    
def parse_cdhit_sequence_line(line):
    """
    Given a single line from the CD-HIT .clstr file, such as:

    0	378aa, >1sek_A mol:protein length:378  SERPIN K... *

    or

    1	378aa, >1vgm_B mol:protein length:378  378aa long hypothetical citrate synthase... at 100%

    parse it and return a CDHIT_seq object.
    Note there is a tab after the entry number

    Parameters:
        line - sequence entry line from .clstr file to parse
    Return value:
        CDHIT_sequence object parsed from line
    Raises exceptions:
        CDHIT_ParseError if line cannot be parsed
    """
    splittab = line.split('\t')
    sequencenum = int(splittab[0])
    seqlen_str = splittab[1].split(',')[0]
    if seqlen_str[len(seqlen_str)-2:len(seqlen_str)] != "aa":
        raise CDHIT_ParseError("could not parse seqlen on line: " + line)
    seqlen = int(seqlen_str[:len(seqlen_str)-2])
    split_descr = splittab[1].split(", >")[1].split("... ")
    descr = split_descr[0]
    percident_str = split_descr[1].strip()
    if percident_str == "*":
        is_repr = True
        percident = 100
    else:
        is_repr = False
        percstr = percident_str.split()[1]
        if not percstr[:len(percstr)-1].isdigit():
            raise CDHIT_ParseError("could not parse percent identity on line: " + line)
        percident = int(percstr[:len(percstr)-1])
    return CDHIT_sequence(sequencenum, seqlen, descr, is_repr, percident)


def yield_cluster_from_file(fh):
    """
    Generator funcion which Parses the CD-HIT cluster file, yielding
    a CDHIT_cluster objects for each cluster in the file.

    Parameters:
        fh - open (read) filehandle of CD-HIT .clstr file

    Return value:
        list of CDHIT_cluster objects 
    """
    cluster = None
    for line in fh:
        if line[0] == ">": # New cluster
            if cluster:
                yield cluster
            clusternum = int(line.split()[1])
            cluster = CDHIT_cluster(clusternum)
        elif line.split()[0].isdigit():
            sequence = parse_cdhit_sequence_line(line)
            cluster.add_sequence(sequence)
        else:
            sys.stderr.write("unrecognized line: " + line)
            cluster = None
    if cluster:
        yield cluster


def get_nr_pdb_list(TMPDIR):
    """
    The main program logic to get the nonredundant list of pdb identifiers,
    selecting the highest resolution as representative.
    See module docstring at top of file for description

    Parameters:
       TMPDIR - name of temp directory to use
    Return value:
       None.

       Output is to stdout:
       list of list of pdb ids, each entry in list (line) is a list of pdb ids
       reprsenting a cluster; first in the inner (cluster) list is
       the chosen represenstative.
    """
    pdb_dict = {} # dict of {pdbid : MethodResolution} to cache info from PDB
    
    for cluster in yield_cluster_from_file(sys.stdin):
        if not cluster[0].is_protein(): # since clustered, if one not, all not
            continue # discard non-protein sequences
        orig_seqlist = list(cluster.seqlist) # keep copy before deleting some
        cluster.discard_short_seqs()
        if len(cluster) > 1:
            # now we need to look in PDB files to find highest res X-ray struct
            for seq in cluster.seqlist:
                pdbid = seq.descr[:4].lower()
                if pdb_dict.has_key(pdbid):
                    methres = pdb_dict[pdbid]
                    seq.is_xray = methres.is_xray
                    seq.resolution = methres.resolution
                else:
                    name = "pdb" + pdbid
                    pdbfile = os.path.join(PDBDIV_ROOT,
                                           os.path.join(pdbid[1:3], name + ".ent.gz"))
                    tmp_pdbfilename = os.path.join(TMPDIR, name)
                    os.system("gzip " + pdbfile + " -d -c > " + tmp_pdbfilename)
                    pdbheader = PDB.parse_pdb_header(tmp_pdbfilename)
                    if 'x-ray' in pdbheader['structure_method'].lower():
                        seq.is_xray = True
                        seq.resolution = float(pdbheader['resolution'])
                    methres = MethodResolution()
                    methres.is_xray = seq.is_xray
                    methres.resolution = seq.resolution
                    pdb_dict[pdbid] = methres
                    os.unlink(tmp_pdbfilename)
            cluster.discard_non_xray()
        if len(cluster) > 1:
            cluster.discard_lower_resolution()
        if len(cluster) > 1:
            cluster.discard_lower_similarity()
        if len(cluster) > 1:
            cluster.seqlist = [cluster.seqlist[0]] # arbitrary: use first seq

        repr_id =  cluster.seqlist[0].descr[:6].lower()
        sys.stdout.write(repr_id + ": ")
        for seq in orig_seqlist:
            other_id = seq.descr[:6].lower()
            if other_id != repr_id:
                sys.stdout.write(other_id + " ")
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
    sys.stderr.write("Usage: " + progname + "\n"
                     "  Input is CD-HIT cluster file on stdin\n"
                     "  Output is list of PDB identifiers on stdout\n"\
                     "    Each has format: \n"
                     "      <repr-id>: <list of ids in cluster>\n")
    sys.exit(1)


def main():
    """
    main forgetnrpdblist.py

    Usage: getnrpdblist.py 
              Input is CD-HIT cluster file on stdin
              Output is list of PDB identifiers on stdout
              Each has format:
              <repr-id>: <list of ids in cluster>

           where <repr-id> is the PDB (+ chain e.g. 2ssp_B) id in pdb_seqres
           format that has been chosen as the representative for the
           cluster and following is list of pdb_seqres identifiers in that
           cluster. Only proteins are output, not DNA, RNA, etc.
           (see module docstring at top of file).
    """
    if len(sys.argv) > 1:
        usage(sys.argv[0])

    TMPDIR = os.tempnam(None, "pdbgz")
    os.mkdir(TMPDIR)
    try:
        get_nr_pdb_list(TMPDIR)
    finally:
        cleanup_tmpdir(TMPDIR)
    
        
if __name__ == "__main__":
    warnings.filterwarnings('ignore', 'tempnam', RuntimeWarning) 
    main()
