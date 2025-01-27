"""Main command-line interface for the eccen_adjust library.

This script requires a file that contains V1 labels. That file
can be an MGZ or label file, but it needs to have the value 1 for
V1 and not anywhere else.
"""

import sys
import argparse

import neuropythy as ny

from . import adjust_eccen_in_v1

# This needs to do a few things:
# (1) parse the arguments for a FreeSurfer directory and
#     eccentricity and v1 label files (varea); also the
#     hemlabel (which hemisphere) option.
#     (optionally) let the user set shape and min_eccen
#     and max_eccen with options.
# (2) load all of these into memory
# (3) call adjust_eccen_in_v1 with the above parameters

# First use argparse:
print("IT WORKS")
sys.exit(0)

# Load in data:
sub = ny.freesurfer_subject(input_directory)
eccen = ny.load(input_eccen)

# Once we have loaded the subject (sub) and hemisphere name (hemlabel) and the
# eccentricity map (eccen) and the labels or "varea" (label):
mask = (label == 1)
hem = sub.hemis[hemlabel]
#mesh = hem.surface('midgray') # could be an option, but nbd
v1_eccen = eccen[mask]
v1_sarea = hem.prop('midgray_surface_area')[mask]
adjusted_eccen = adjust_eccen_in_v1(v1_eccen, v1_sarea)
# ^^^ optionally also include shape, min_eccen etc.

# Write out adjusted_eccen file.
# (This will need a command-line option to know the filename)
ny.save(output_filename, adjusted_eccen)


