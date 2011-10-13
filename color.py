###############################################################################
#
# color.py - color space (HSV, RGB) functions.
#
# File:    color.py
# Author:  Alex Stivala
# Created: December 2007
#
# $Id$
#
# Utility functions for color space operations.
#
###############################################################################

import re

#-----------------------------------------------------------------------------
#
# Constants
#
#-----------------------------------------------------------------------------

# Dictionary mapping color names to hex color strings as per HTML 4.01
# specifiction (http://www.w3.org/TR/REC-html40/types.html#h-6.5)
# all names here in lowercase, should lowercase names for lookup as color
# names should be case insensitive.
HTML_COLOR_DICT = dict([('black',   '000000'),
                        ('silver',  'c0c0c0'),
                        ('gray',    '808080'),
                        ('white',   'ffffff'),
                        ('maroon',  '800000'),
                        ('red',     'ff0000'),
                        ('purple',  '800080'),
                        ('fuchsia', 'ff00ff'),
                        ('green',   '008000'),
                        ('lime',    '00ff00'),
                        ('olive',   '808000'),
                        ('yellow',  'ffff00'),
                        ('navy',    '000080'),
                        ('blue',    '0000ff'),
                        ('teal',    '008080'),
                        ('aqua',    '00ffff'),

                        ('beige',   'f0f0d2')  # Dunnart default shape color
                        ])

# TODO: there are only 16 colors in the HTML spec, should change to using
# SVG colors (http://www.w3.org/TR/SVG/types.html#ColorKeywords).
COLOR_DICT = HTML_COLOR_DICT


DUNNART_CLUSTER_ALPHA = 0.33 # alpha channel (opactiy) value for clusters
DUNNART_CLUSTER_ALPHA_HEX = '55' # alpha channel (opacity) value in hex

# list of cluster shading colors, in format
# rrggbbaa (hex: red, green, blue, alpha) for Dunnart.
# Note we still use this if there are few enough clusters since
# although the color gradient looks fine on Dunnart, when using
# Inkscape (eg converting to PNG) or Firefox, most of the shading ends up being
# rendered as the same default color for some reason, eg green (00ff0055)
# does not work at all, ends up being rendered same as blue (!?) 
# (This was a Dunnart bug, which I've now fixed (1/7/08), but might as
# well leave this default selection of colors).
DUNNART_DEFAULT_CLUSTER_FILL_COLOR = "60cdf355"
DUNNART_CLUSTER_FILL_COLORS = [ DUNNART_DEFAULT_CLUSTER_FILL_COLOR,
                                "ff000055",   # red
                                "551a8b55",   # purple4
                                "0000ff55"    # blue
                              ]


# Data frm Table II of Glasbey et al 2007 "Colour Displays for Categorical
# Images" Color Research and Application 32(4):304-309
# These are the first 32 colours found by sequential search to solve 
# their formulation of maximizing the minimimum distance between colours
# (in the CIELAB colour space).
# These are the RGB tuple values of those colors (each in [0,255]).
# Use to obtain a list of up to 32 colors that are "maximally distinguishable"
# ie to human perception are least likely to be confused.
GLASBEY_COLORS     = [(255, 255, 255),
                      (  0,   0, 255),
                      (255,   0,   0),
                      (  0, 255,   0),
                      (  0,   0,  51),
                      (255,   0, 182),
                      (  0,  83,   0),
                      (255, 211,   0),
                      (  0, 159, 255),
                      (154,  77,  66),
                      (  0, 255, 190),
                      (120,  63, 193),
                      ( 31, 150, 152),
                      (255, 172, 253),
                      (177, 204, 113),
                      (241,   8,  92),
                      (254, 143,  66),
                      (221,   0, 255),
                      ( 32,  26,   1),
                      (114,   0,  85),
                      (118, 108, 149),
                      (  2, 173,  36),
                      (200, 255,   0),
                      (136, 108,   0),
                      (255, 183, 159),
                      (133, 133, 103),
                      (161,   3,   0),
                      ( 20, 249, 255),
                      (  0,  71, 158),
                      (220,  94, 147),
                      (147, 212, 255),
                      (  0,  76, 255)
                     ]

#-----------------------------------------------------------------------------
#
# Function definitions
#
#-----------------------------------------------------------------------------

def get_glasbey_colors_rgb():
   """
   Get data from Table II of Glasbey et al 2007 "Colour Displays for Categorical
   Images" Color Research and Application 32(4):304-309
   These are the first 32 colours found by sequential search to solve 
   their formulation of maximizing the minimimum distance between colours
   (in the CIELAB colour space).
   These are the RGB tuple values of those colors (each in [0,1]),
   omitting the first which is white.
   Use to obtain a list of up to 31 colors that are "maximally distinguishable"
   ie to human perception are least likely to be confused.

   Praameters: None
   Return value: list of (r,g,b) tuples; each r,g,b in [0,1]
   """
   return [(r/255.0, g/255.0, b/255.0) for (r,g,b) in GLASBEY_COLORS[1:]]


def hsv_to_rgb(h, s, v):
    """
    Convert a color specified by a point in HSV color space to the
    corresponding point in RGB color space.

    Parameters:
       h - Hue value in [0,360] (ie degrees), may be None (when s is 0)
       s - Saturation value in [0,1]
       v - Value value in [0,1]

    Retrun value:
       (r, g, b) tuple corresponding to the h,s,v supplied.
       Each of r,g,b is in [0,1]

    Raises exceptions:
       ValueError if s is 0 but h is not None

    This is an implementation of the procedure given in
    Foley and van Dam (1982) "Fundamentals of Interative Computer Graphics"
    (Addison-Wesley) p. 616.
    See also the wikipedia "HSL Color Space" entry (which cites the above book).
    """
    assert(0.0 <= h <= 360.0)
    assert(0.0 <= s <= 1.0)
    assert(0.0 <= v <= 1.0)
    
    if s == 0:
        if h == None:
            r = v
            g = v
            b = v
        else:
            raise ValueError('hsv_to_rgb: if S is 0 then H is undefined')
    else:
        if h == 360:
            h = 0
        h = h / 60.0
        i = int(h)
        f = h - i
        p = v * (1.0 - s)
        q = v * (1.0 - s*f)
        t = v * (1.0 - s*(1.0 - f))
        if i == 0:
            (r, g, b) = (v, t, p)
        elif i == 1:
            (r, g, b) = (q, v, p)
        elif i == 2:
            (r, g, b) = (p, v, t)
        elif i == 3:
            (r, g, b) = (p, q, v)
        elif i == 4:
            (r, g, b) = (t, p, v)
        elif i == 5:
            (r, g, b) = (v, p, q)
        else:
            assert(False)
    return (r, g, b)


def color_gradient(n):
    """
    Return a color gradient from blue to red consisting of n (parameter)
    different colors (RGB).

    This is done by varying Hue in HSV color space while holding
    Saturation and Value constant, and converting the each to RGB.

    This is a generator function.

    Parameters:
       n - The number of different colors to return.

    Return value:
       YIELDs an (r,g,b) tuple (n times).
    """
    saturation = 1.0  
    value = 1.0       # saturation and value are kept constant
    final_hue = 0.0   # red is at 0 degrees.
    initial_hue = 240.0  # initial hue is at 240 degrees, which is blue.
                      # (240 degrees is 4 60degree steps around the HSV hexcone)
    stepsize = initial_hue / n
    hue = initial_hue
    i = 0
    while i < n:
        (r, g, b) = hsv_to_rgb(hue, saturation, value)
        yield (r, g, b)
        hue -= stepsize
        i += 1
        


def rgb_tuple_to_hex_str(rgb, a=None):
    """
    Return a hex string representation of the supplied RGB tuple
    (as used in SVG etc.) with alpha, rescaling all values from [0,1]
    into [0,255] (ie hex x00 to xff).

    Parameters:
       rgb - (r,g,b) tuple (each of r,g,b is in [0,1])
       a - (default None) alpha channel (transparancy) value or None

    Return value:
       string 'rrggbb' where rr, gg, bb are hex representations of r,g,b resp.
       or 'rrggbbaa' if a!=None where aa is hex repr. of alpha value a.
    """
    xrgb = tuple([int(round(x*255)) for x in rgb])
    if a == None:
        return '%02x%02x%02x' % xrgb
    else:
        a = int(round(a*255))
        return '%02x%02x%02x%02x' % (xrgb[0],xrgb[1],xrgb[2],a)
    
    
    
def hex_str_to_rgb_tuple(hexstr):
    """
    Return RGBA or RGBA tuple (all values in [0,1]) correspdonding to supplied
    hex color string.

    Parameters:
        hexstr - 'rrggbb' or 'rrggbbaa' hex color string
    Return value:
        (r,g,b) or (r,g,b,a) tuple (all values in [0,1])
        correspdonding to hex color string.
    Raises exceptions:
        ValueError if invalid hexstring format.
    """
    rgbstring_re = re.compile(r'[0-9A-Fa-f]{6}')
    rgbastring_re = re.compile(r'[0-9A-Fa-f]{8}')

    if rgbstring_re.match(hexstr):     # RRGGBB
        r = int(hexstr[0:2], 16) / 255.0
        g = int(hexstr[2:4], 16) / 255.0
        b = int(hexstr[4:6], 16) / 255.0
        return (r,g,b)
    elif rgbastring_re.match(hexstr):  # RRGGBBAA
        r = int(hexstr[0:2], 16) / 255.0
        g = int(hexstr[2:4], 16) / 255.0
        b = int(hexstr[4:6], 16) / 255.0
        a = int(hexstr[6:8], 16) / 255.0
        return (r,g,b,a)
    else:
        raise(ValueError("Invalid hex color string '" + hexstr + "'"))



def get_color_list(color_list_str):
    """
    Return the RGB color tuples for the colors specified 
    as comma-delimited lists of recognized color names
    or RGB colors in hex (e.g. red,purple,#ddf02d).

    Parameters:
       color_list_str - list of colors as above
    Return value:
       list of color value strings
       where the color value strings are 'rrggbb' color hex strings (no '#')
    Uses globals (Readonly):
        COLOR_DICT (color.py) - maps color names to RGB hex strings
    Raises Excpetions:
        KeyError for unknown color names
    """
    color_str_list = color_list_str.split(',')
    color_list = []
    for color_str in color_str_list:
        if color_str[0] == '#':
            color_hex = color_str[1:]
        else:
            color_hex = COLOR_DICT[color_str.lower()]
        color_list.append(color_hex)
    return color_list



def get_cluster_fill_colors(color_list_str, num_clusters):
    """
    Get list of cluster fill (shading) colros from supplied list of
    colors, or if 'auto', generate them from default colors or color
    gradient.
    If list of colors is supplied, but num_clusters > length of that list
    then treat list as circular (i.e. when run out go back and resuse
    colors from start).
    

    Parameters:
         num_clusters - number of clusters we need to shade
         color_list_str -  'auto' (use color gradient to shade each
                            differently) or list of colors.

    Return value:
       list of color value strings
       where the color value strings are 'rrggbbaa' color hex strings (no '#')
       with alpha channel (opacity) aa set to the default value '55'.
       length of list is num_clusters.
       
    Uses globals (Readonly):
        DUNNART_CLUSTER_FILL_COLORS, DUNNART_CLUSTER_ALPHA
    
    """
    if color_list_str == 'auto':
        # set up list of cluster fill colors, one per sheet,
        # by using color gradient to ensure all sheets
        # different color Note we still use
        # DUNNART_CLUSTER_FILL_COLORS if there are few enough
        # sheets since although the color gradient looks fine
        # on Dunnart, when using Inkscape (eg converting to
        # PNG), most of the shading ends up being rendered as
        # the same default color for some reason. (FIXME)
        if num_clusters <= len(DUNNART_CLUSTER_FILL_COLORS):
            cluster_fill_colors = list(DUNNART_CLUSTER_FILL_COLORS)
        else:
            cluster_fill_colors = [ rgb_tuple_to_hex_str(rgb,
                                                  DUNNART_CLUSTER_ALPHA)
                                    for rgb in color_gradient(num_clusters)
                                  ]
    else:
        color_list = [ color + DUNNART_CLUSTER_ALPHA_HEX for color in
                       get_color_list(color_list_str) ]
        i = 0
        cluster_fill_colors = []
        while (i < num_clusters):
            cluster_fill_colors.append(color_list[i % len(color_list)])
            i += 1
    return cluster_fill_colors


