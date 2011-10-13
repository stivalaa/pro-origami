###############################################################################
#
# ptposmap.py - class to keep track of neighbours in a 2d grid 
#
# File:    ptrelpos.py
# Author:  Alex Stivala
# Created: October 2007
#
# The PTPosMap class keeps track of items in a 2d grid, ie each item
# can have only four immediate neighbours (above, below, left, right
# or (equivalently) north, south, east, west).
#
# $Id$
#
#
###############################################################################

from ptrelpos import RELPOS_ABOVE, RELPOS_BELOW, RELPOS_LEFT, RELPOS_RIGHT

#-----------------------------------------------------------------------------
#
# Module globals
#
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#
# Class definitions 
#
#-----------------------------------------------------------------------------
    
class PTNeighbours:
    """
    The PTNeighbours class is just a struct with 4 members, one for each
    of the directions we use. Access directly or via get_neighbour(),
    set to the object type
    we are keeping track of, or None.
    """
    def __init__(self):
        """
        Initialize PTNeighbours struct with all neighbours None.
        """
        self.north = None
        self.south = None
        self.east = None
        self.west = None

    def __str__(self):
        """
        Return string represetnation of the struct
        """
        return str(self.north) + "\t" + str(self.south) + "\t" + \
               str(self.east) + "\t" + str(self.west)
        

    def get_neighbour(self, relpos):
        """
        Given a relpos (RELPOS_ABOVE etc. (ptrelpos.py), return the
        neighbour in that direction.

        Parameters:
           relpos - RELPOS_ABOVE/BELOW/LEFT/RIGHT to get neighbour in
                    that direction
        Return value:
           the value of the correpsonding north/south/east/west neighbour
           (may be None).

        Raises exceptions:
           TypeError if direction is not one of RELPOS_ABOVE, etc.

        """
        # FIXME: we could have avoided all this stuffing about with
        # RELPOS_ABOVE meaning north etc. by just indexing a list by
        # RELPOS_ values instead of having north,south,etc. data members
        if relpos == RELPOS_ABOVE:
            return self.north
        elif relpos == RELPOS_BELOW:
            return self.south
        elif relpos == RELPOS_LEFT:
            return self.east
        elif relpos == RELPOS_RIGHT:
            return self.west
        else:
            raise TypeError('bad direction value ' +str(relpos) + '\n')

    
class PTPosMap:
    """
    The PTPosMap class keeps track of items in a 2d grid, ie each item
    can have only four immediate neighbours (above, below, left, right
    or (equivalently) north, south, east, west).

    
    We implement it as a mapping container, i.e. using __getitem__ and
    __setitem__ so that elements can bet get/set with dictionary/array
    type syntax e.g. posmap[objid].

    """
    def __init__(self):
        """
        Initialize empty PTPosMap
        """
        self.pos_dict = {}  # dictionary of { object : PTNeighbours }

    def __str__(self):
        """
        Return string representation of the PTPosMap
        """
        s = ""
        for obj,neighbours in self.pos_dict.iteritems():
            s += str(obj) + ": " + str(neighbours) + "\n"
        return s
    

    #
    # methods defined to implement container type
    #

    def __len__(self):
        """
        Return number of entries set in the posmap

        Parameters: None
        Return value: Number of entries set in posmap
        """
        return len(self.pos_dict)


    def __getitem__(self, obj):
        """
        Return the entry in the posmap for the supplied object.

        Parameters:
           obj - object (any type) to look up in posmap
        Return value:
           PTNeighbours object with the neighbours of the obj
        
        """
        return self.pos_dict[obj]


    def __setitem__(self, obj, neighbours):
        """
        Set the entry in the posmap for the supplied obj to the supplied
        neighbours struct.
        
        Parameters:
           obj - object (any type) to set entry for
           neighbours - PTNeighbours struct to set entry to

        Return value: None

        Raises exceptions:
           TypeError  if neighbours is not a PTNeighbours instance
        """
        if not isinstance(neighbours, PTNeighbours):
            raise TypeError("bad neighbours parameter\n")

        self.posdict[obj] = neighbours


    def has_key(self, obj):
        """
        Return True iff the obj has an entry in the posmap
        Paratmers: obj - object (any type) to check for in posmap
        Return value: True if obj is a key in the posmap else False
        """
        return self.pos_dict.has_key(obj)

    # have not implemented: __delitem__, __iter__, __contains__
    
    def add_neighbour_obj(self, obj, neighbour_obj, direction):
        """
        Update the entry for obj with a new neighbour in the specified
        direction. Overwrites existing neighbour in that direction if any
        (does not check for this).  Entry will be created if not
        already there. Also adds the symetric entries i.e. if
        A is added as LEFT of B, then B is also updated as RIGHT of A

        Parameters:
           obj - object (any type) to update entry for (need not exist)
           neighbour_obj - object (any type) that is the neighbour
           direction - RELPOS_ABOVE, etc direction of neighbour

        Return value: None

        Raises exceptions:
           TypeError if direction is not one of RELPOS_ABOVE, etc.
        """
        if self.pos_dict.has_key(obj):
            neighbours = self.pos_dict[obj]
        else:
            neighbours = PTNeighbours()

        if self.pos_dict.has_key(neighbour_obj):
            sym_neighbours = self.pos_dict[neighbour_obj]
        else:
            sym_neighbours = PTNeighbours()
            
        if direction == RELPOS_ABOVE:
            neighbours.north = neighbour_obj
            sym_neighbours.south = obj
        elif direction == RELPOS_BELOW:
            neighbours.south = neighbour_obj
            sym_neighbours.north = obj
        elif direction == RELPOS_LEFT:
            neighbours.east = neighbour_obj
            sym_neighbours.west = obj
        elif direction == RELPOS_RIGHT:
            neighbours.west = neighbour_obj
            sym_neighbours.east = obj
        else:
            raise TypeError('bad direction value ' +str(direction) + '\n')

        self.pos_dict[obj] = neighbours
        self.pos_dict[neighbour_obj] = sym_neighbours
        
        
