import numpy as np
import matplotlib.pyplot as pl
import itertools

import LT.box as B

# create a limiter data file for orbit

# grouping functions
def group0(lst, n):
      for i in range(0, len(lst), n):
            val = lst[i:i+n]
            yield tuple(val)

def group3(lst, n, defaultvalue = 0):
        """group([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 3) => iterator
    
        Group an iterable into an n-tuples iterable. Incomplete tuples
        are padded with defaultvalue
    
        >>> list(group(range(10), 3))
        [(0, 1, 2), (3, 4, 5), (6, 7, 8), (9, 0, 0)]
        """
        # array of slices starting at range(n) and stepping over n elements
        # first element: [0, 3, 6, 9]
        # second       : [1, 4, 7]
        # third        : [2, 5, 8]
        slices_array = [itertools.islice(lst, i, None, n) for i in range(n)]
        # zip the arrays together fill missing elements with fill value
        # *slices_array passes the arrays as an arbitrary number of argument
	if defaultvalue == None:
		return group0(lst, n)
	else:
		return itertools.izip_longest( *slices_array, fillvalue = defaultvalue)


def print_group5(x, out):
    for g in group3(x, 5):
        out.write(std_line_f5(g[0], g[1], g[2], g[3], g[4]))

def print_group6(x, out, skip_last_nl = False):
      lines = []
      # for g in (group3(x, 6, defaultvalue = None)):
      for g in (group3(x, 6, defaultvalue = None)):
	      # lines.append(std_line_f6(g[0], g[1], g[2], g[3], g[4], g[5]))
	      lines.append(std_line_f(*g))
      for l in lines[:-1]:
            o.write(l)
      last_line = lines[-1]
      if skip_last_nl:
            o.write(last_line[:-1])
      else:
            o.write(last_line)

# Formatting according to fortran


def std_line_f5(x0,x1,x2,x3,x4):
	return "{0:16.9e},{1:16.9e},{2:16.9e},{3:16.9e},{4:16.9e},\n".format(x0,x1,x2,x3,x4)

def std_line_f6(x0,x1,x2,x3,x4,x5):
	return "{0:16.9e},{1:16.9e},{2:16.9e},{3:16.9e},{4:16.9e},{4:16.9e},\n".format(x0,x1,x2,x3,x4,x5)

def std_line_i2(x0,x1):
	return "{0:5d} {1:5d}".format(x0,x1)

# variable number of floats
def std_line_f(*x):
    fmt0 = "{"
    fmt1 = "{0:d}:16.9"
    fmt2 = "}, "
    string = ""
    for i,xx in enumerate(x):
        fmt = fmt0 + fmt1.format(i) + fmt2
        string += fmt
    string += "\n"
    return string.format(*x)




fi = B.get_file("mast_limiter_inner.data")
fo = B.get_file("mast_limiter_outer.data")

ri = B.get_data(fi,'r'); zi = B.get_data(fi,'z')
ro = B.get_data(fo,'r'); zo = B.get_data(fo,'z')

# now write the limiter data in a poper format
# there is only one toroidal region

o = open("MASTLIMIT00.DAT",'w')

o.write(" 1                                       ! only one toroidal region\n")
o.write("MAST\n")
o.write("  0.0                                    ! start of region\n")
o.write(std_line_i2( len(ri), len(ro) )+"                    ! number of points on small-R & large R \n")
print_group6(ri[::-1], o)  # reverset the ordering of the innter data
o.write("\n")
print_group6(zi[::-1], o)
o.write("\n")
print_group6(ro, o)
o.write("\n")
print_group6(zo, o, skip_last_nl = True)
end_txt =\
"""
;				Below is info on how to display limiter in plot
 2				; # limiter profiles to plot
 1     0			;1st profile is toroidal segment 1, linetype 0
 4     3			;2nd profile: toroidal seg 4 (HHFW ant) dashed

"""
o.write(end_txt)
o.close()
