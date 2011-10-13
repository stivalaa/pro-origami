###############################################################################
#
# geometry.py - geometry (3D) functions.
#
# File:    geometry
# Author:  Alex Stivala
# Created: November 2007
#
# $Id$
#
# Utility functions for geometry in 3D space.
#
###############################################################################

from numpy import * 
from Bio.PDB import Vector

def LineLineIntersect(p1, p2, p3, p4):
    """
    Calculate the line segment PaPb that is the shortest route between
    two lines P1P2 and P3P4. This line is perpendicular to both P1P2 and P3P4.
    Calculate also the values of mua and mub
    where
      PA = P1 + mua (P2 - P1)
      Pb = P3 + mub (P4 - p3)

    This function is just a simple Python implementation of the
    algorithm described and implememted in C by Paul Bourke:

      http://local.wasp.uwa.edu.au/~pbourke/geometry/lineline3d/

    Parameters:
        p1 - Vector of first point for first line
        p2 - Vector of second point for first line
        p3 - Vector of first point of second line
        p4 - Vector of second point of second line
        

    Return value:
        tuple (pa, pb, mua, mub) where
          pa is  Vector of Pa, first point on shortest line
          pb is Vector of Pb, second point on shortest line
          mua is the value s.t. Pa = P1 + mua (P2 - P1)
          mub is the value s.t. Pb = P3 + mub (P4 - P4)
        or None if there is no solution
    """
    EPS = 1.0e-08
    
    p13 = p1 - p3
    p43 = p4 - p3

    if alltrue(less(abs(p43.get_array()), EPS)):
        return None
    p21 = p2 - p1
    if alltrue(less(abs(p21.get_array()), EPS)):
        return None

    d1343 = p13[0] * p43[0] + p13[1] * p43[1] + p13[2] * p43[2]
    d4321 = p43[0] * p21[0] + p43[1] * p21[1] + p43[2] * p21[2]
    d1321 = p13[0] * p21[0] + p13[1] * p21[1] + p13[2] * p21[2]
    d4343 = p43[0] * p43[0] + p43[1] * p43[1] + p43[2] * p43[2]
    d2121 = p21[0] * p21[0] + p21[1] * p21[1] + p21[2] * p21[2]

    denom = d2121 * d4343 - d4321 * d4321
    if abs(denom) < EPS:
        return None

    numer = d1343 * d4321 - d1321 * d4343
    mua = numer / denom
    mub = (d1343 + d4321 * mua) / d4343

    # Note using Numeric.array '*' operator here for element-wise multiplication
    # as Bio.PDB.Vector '*' operator is vector dot product.
    # (Note also must have Vector + int * array and NOT
    # int * array + Vector due to Python type coercion rules).
    pa = p1 + mua * p21.get_array()
    pb = p3 + mub * p43.get_array()

    return (pa, pb, mua, mub)


def ProjectPointOntoLine(A, B, P):
    """
    Project a point (in 3D), P, onto a line (in 3D) defined by the two
    points A and B.
    This point Q is the point on line AB such that the line PQ is orthogonal
    to AB (PQ is the shortest line between the point P and the line AB)

    Parameters:
       A - Vector for point A on the line AB
       B - Vector for point B on the line AB
       P - Vector for point P to project onto the line AB

    Return value:
       Vector representing the point Q on line AB such that PQ is the shortest
       line from P to line AB.
    """
    # this basically involves finding the point Q such that
    # (P - Q) * (B - A) = 0
    # where * is dot product.
    # This is done by solving for u after substituting the equation for line AB:
    # Q = A + u(B - A)
    # into the first equation.
    # See  http://local.wasp.uwa.edu.au/~pbourke/geometry/pointline/

    # Note * here is just scalar multiplication (float)
    u =((P[0]-A[0])*(B[0]-A[0])+(P[1]-A[1])*(B[1]-A[1])+(P[2]-A[2])*(B[2]-A[2]))\
        / (B - A).normsq()
    Q = A + u*(B - A).get_array() # this * is scalar*vector (Numeric.array *)
    return Q
