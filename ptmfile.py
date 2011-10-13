###############################################################################
#
# ptmfile.py - Functions to write MATLAB .m file for drawing 3D strand axes
#
# File:    ptmfile.py
# Author:  Alex Stivala
# Created: November 2007
#
# $Id$
#
# Utility functions for writing MATLAB .m files to plot 3D strand backbone
# carbon-alpha traces and the axes fit to them.
#
###############################################################################

import sys
from time import strftime,localtime

from numpy.oldnumeric import array
from Bio.PDB import Vector

def mfile_write_prelude(fh):
    """
    Write the start of the m-file, which is just intructions to MATLAB
    to hold plot etc.

    Parameters:
        fh - open filehandle of m-file to write to
    Return value: None
    """
    timestamp = strftime("%d%b%Y %H:%M:%S", localtime())
    fh.write("% Generated on " + timestamp
             + " by ptmfile.py $Revision$: "
             + " ".join(sys.argv) + "\n")

    fh.write('hold on\n')
    
    
def mfile_write_strand(fh, c_alpha_veclist, midpoints, centroid, dircos,
                       nterm_point, cterm_point, label):
    """
    Write the data for a single strand to the open m-file. This information
    is an array of the carbon-alpha co-ordinates of residues in the strand,
    for which a 3D line is plotted, a vector of the centroid, for which
    an asterisk is plotted, and a direction cosines vector for the axis,
    along which a line is drawn (through the centroid).
    A red asterisk is plotted at one end of the axis line segment to indicate
    the direction along the line of the direction cosines vector.
    A blue plus is plotted at the projection of the most nterminal midpoitn
    on the axis, and a red circle at the projection of the most cterminal
    midpoint on the axis.

    Parameters:
       fh - open filehandle of m-file to write to
       c_alpha_veclist - list of Bio.PDB.Vector for each C-alpha atom co-ord
       midpoints - array of Bio.PDB.Vectors of each pair of consecutive
                   C_alpha atoms, or None.
       centroid - Bio.PDB.Vector of centroid
       dircos - Bio.PDB.Vector of direction cosine of computed strand axis

    Return value: None
    """
    coords_str = "A = "
    coords_str += str(array([list(vec) for vec in c_alpha_veclist]))
    coords_str += ";\n"

    if midpoints != None:
        midpoints_str = "M = "
        midpoints_str += str(array([list(vec) for vec in midpoints]))
        midpoints_str += ";\n"

    centroid_str = "c = " + str(list(centroid)) + ";\n"
    dircos_str = "v = " + str(list(dircos)) + ";\n"
    
    fh.write(coords_str)
    fh.write(centroid_str)
    fh.write(dircos_str)

    fh.write("plot3(A(:,1),A(:,2),A(:,3),'r');\n") # red c-alpha backbone line

    # cyan line along axis through centroid in direction of direction cosines vector
    fh.write("l = [c(1)-30*v(1), c(2)-30*v(2), c(3)-30*v(3); ")
    fh.write("     c(1)+30*v(1), c(2)+30*v(2), c(3)+30*v(3)];\n")
    fh.write("plot3(l(:,1), l(:,2), l(:,3),'c');\n")

    fh.write("plot3(c(1),c(2),c(3),'b*');\n") # blue asterisk for centroid

    # plot red asterisk at end of axis line segment to indicate direction
    fh.write("d = c + 30 * v;\n")
    fh.write("plot3(d(1), d(2), d(3), 'r*');\n")

    # plot blue plus at projection of n-terminal midpoint on axis, and
    # red circle at projection of c-terminal midpoint on axis
    nterm_str = "np = " + str(list(nterm_point)) + ";\n"
    cterm_str = "cp = " + str(list(cterm_point)) + ";\n"
    fh.write(cterm_str)
    fh.write(nterm_str)
    fh.write("plot3(np(1),np(2),np(3),'b+');\n")
    fh.write("plot3(cp(1),cp(2),cp(3),'ro');\n")
    
    # write label at centroid
    fh.write("text(c(1), c(2), c(3), '" + " " + label + "')\n")


def mfile_write_helix(fh, c_alpha_veclist, midpoints, centroid,
                      dircos, nterm_point, cterm_point, label):
    """
    Write the data for a single helix to the open m-file. This information
    is an array of the carbon-alpha co-ordinates of residues in the helix,
    for which a 3D line is plotted, a vector of midpoints of consectuive
    C_alpha triples, for wihch red crosses are plotted,
    a vector of the centroid, for which
    an asterisk is plotted, and a direction cosines vector for the axis,
    along which a line is drawn (through the centroid).
    A red asterisk is plotted at one end of the axis line segment to indicate
    the direction along the line of the direction cosines vector.
    A blue plus is plotted at the projection of the most nterminal midpoitn
    on the axis, and a red circle at the projection of the most cterminal
    midpoint on the axis.
    
    Parameters:
       fh - open filehandle of m-file to write to
       c_alpha_veclist - list of Bio.PDB.Vector for each C-alpha atom co-ord
       midpoints - list of Bio.PDB.Vector of midpoints of consectuive
                   C_alpha triples used to fit axis to.
       centroid - Bio.PDB.Vector of centroid
       dircos - Bio.PDB.Vector of direction cosine of computed helix axis
       nterm_point - Bio.PDB.Vector of project of most N-terminal midpoint onto
                     axis
       cterm_point - Bio.PDB.Vector of projection of most C-terminal midpoint
                     onto axis.
       label - string to write as label for the helix

    Return value: None
    """
    coords_str = "A = "
    coords_str += str(array([list(vec) for vec in c_alpha_veclist]))
    coords_str += ";\n"

    midpoints_str = "M = "
    midpoints_str += str(array([list(vec) for vec in midpoints]))
    midpoints_str += ";\n"
    
    centroid_str = "c = " + str(list(centroid)) + ";\n"
    dircos_str = "v = " + str(list(dircos)) + ";\n"
    
    fh.write(coords_str)
    fh.write(midpoints_str)
    fh.write(centroid_str)
    fh.write(dircos_str)

    fh.write("plot3(A(:,1),A(:,2),A(:,3),'m');\n") # magenta c-alpha backbone line
    # plot red cross for each midpiont of C_alpha triple
    fh.write("for i = 1:size(M,1)\n")
    fh.write("   plot3(M(i,1), M(i,2), M(i,3), 'rx');\n")
    fh.write("end;\n")
    
    fh.write("plot3(c(1),c(2),c(3),'k*');\n") # black asterisk for centroid

    # cyan line along axis through centroid in direction of direction cosines vector
    fh.write("l = [c(1)-30*v(1), c(2)-30*v(2), c(3)-30*v(3); ")
    fh.write("     c(1)+30*v(1), c(2)+30*v(2), c(3)+30*v(3)];\n")
    fh.write("plot3(l(:,1), l(:,2), l(:,3),'c');\n")

    # plot red asterisk at end of axis line segment to indicate direction
    fh.write("d = c + 30 * v;\n")
    fh.write("plot3(d(1), d(2), d(3), 'r*');\n")

    # plot blue plus at projection of n-terminal midpoint on axis, and
    # red circle at projection of c-terminal midpoint on axis
    nterm_str = "np = " + str(list(nterm_point)) + ";\n"
    cterm_str = "cp = " + str(list(cterm_point)) + ";\n"
    fh.write(cterm_str)
    fh.write(nterm_str)
    fh.write("plot3(np(1),np(2),np(3),'b+');\n")
    fh.write("plot3(cp(1),cp(2),cp(3),'ro');\n")
    
    # write label at centroid
    fh.write("text(c(1), c(2), c(3), '" + label + "')\n")
        
    
def mfile_write_conclusion(fh):
    """
    Write the concluding commands to the open m-file after all strands
    have been writtein with mfile_write_strand().
    This is instructions to turn on grid, etc.

    Parameters:
       fh - open filehandle of m-file to write to

    Return value: None
    """
    fh.write("grid on\n")
    
