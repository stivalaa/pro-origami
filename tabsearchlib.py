###############################################################################
#
# tabsearchlib - functions used in tableau searching
#
# File:    tabsearchlib.py
# Author:  Alex Stivala
# Created: June 2008
#
# $Id$
#
#
###############################################################################
"""
Library of functions used in different tableau searching methods, such
as score functions for tableaux entries.
"""

import pickle
import sys
from math import pi,degrees

from pttableau import PTTableauPacked

#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

def score_discrete(tabcode_a, tabcode_b):
    """
    Return the tableau matching score between the two tableaux entries
    tabcode_a amd tabcode_b, as per Kamat et al (2008) (s5.1).
    This score is 2 if the tableau entries are equal, 1 if they are equal
    in only one position, else -2.
    """
    if tabcode_a[0] == tabcode_b[0]:
        if tabcode_a[1] == tabcode_b[1]:
            return 2
        else:
            return 1
    elif tabcode_a[1] == tabcode_b[1]:
        return 1
    else:
        return -2

def score_numeric(omega_a, omega_b):
    """
    Return the tableau matching score between two Omega matrix entries
    omega_a and omega_b, as per Kamat et al (2008)

    Parameters:
        omega_a - angle in (-pi, pi]
        omega_b - angle in (-pi, pi]
    Return value:
        score betweem omega_a and omega_b
    """
    delta = min( abs(omega_a - omega_b), 2*pi - abs(omega_a - omega_b) )
    score = pi/4 - delta;
    return score

def score_numeric_deg(omega_a, omega_b):
    """
    Return the tableau matching score between two Omega matrix entries
    omega_a and omega_b, as per Kamat et al (2008),
    but this time exactly as in the paper (so results for sequence
    alignment style dp are the same) - although parameters are still
    in radians convert to degrees in range (0,360] and calculate score
    from that, so values are larger than in score_numeric().

    Parameters:
        omega_a - angle in (-pi, pi]
        omega_b - angle in (-pi, pi]
    Return value:
        score betweem omega_a and omega_b
    """
    omega_a_deg = degrees(omega_a) + 180.0
    omega_b_deg = degrees(omega_b) + 180.0
    delta = min( abs(omega_a_deg - omega_b_deg),
                 360.0 - abs(omega_a_deg - omega_b_deg) )
    score = 45.0 - delta;
    return score



def load_db(db_filename, use_numeric, verbose=False):
    """
    Load the tableaux database into memory.

    Parameters:
       db_filename - filename of a tableaux database created with
                     buldtableauxdb.py
       use_numeric - if True, entries are Numeric.array omega matrices
                     not PTTableauPacked objects
       verbose - If True, write information to stderr

    Return value:
           dict of { pdbid : [PTTableauPacked list] } or
                   { pdbid : [Numeric.array list] } if use_numeric
                pdbid is e.g. 1QLP when built from PDB or
                         e.g. 1qlp1a when built from ASTRAL pdbstyle
    
    """
    if verbose:
        sys.stderr.write('loading tableaux database from ' + db_filename +
                         '...\n')
    db = pickle.load(open(db_filename))
    if not isinstance(db, dict):
        return None
    if not isinstance(db.iterkeys().next(), str):
        return None
    if not isinstance(db.itervalues().next(), list):
        return None
    if not use_numeric:
        if not isinstance(db.itervalues().next()[0], PTTableauPacked):
            return None
    if verbose:
        sys.stderr.write('loaded ' + str(len(db)) + ' entries.\n')
    return db

