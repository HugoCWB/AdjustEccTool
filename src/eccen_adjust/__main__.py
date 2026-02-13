"""Main command-line interface for the eccen_adjust library.

This script requires a file that contains V1 eccentricity and area labels. 
That file can be an MGZ or label file, but it needs to have the value 1 for
V1 and not anywhere else. Additionally, this file needs to be placed in the 
participant's FreeSurfer surf (if .mgz/mgh) or label (if .label) folder.

The script also retrieves information about the participant from the FreeSurfer 
subjects directory - Please make sure to set SUBJECTS_DIR directory accordingly 
beforehand.
"""

import sys
import argparse

import neuropythy as ny
import numpy as np

from . import adjust_eccen_in_v1
from . import plot_originalvsadjusted, plot_distributionECCvalues, plot_comparisonCorticalECC,plot_CMF

# This needs to do a few things:
# (1) parse the arguments for a subject ID in a FreeSurfer directory and
#     eccentricity (eccen) and v1 label files (varea); also the
#     hemlabel (which hemisphere) option.
#     (optionally) let the user set shape and min_eccen
#     and max_eccen with options.
# (2) load all of these into memory
# (3) call adjust_eccen_in_v1 with the above parameters
# (4) print out useful information
# (5) write out adjusted eccentricity output
# (6) write out vertex labels into a txt file for use in Matlab
# (7) Plot the figures for:
#           - Benson pre- versus post-eccentricity values
#           - Distribution of Benson eccentricity values for pre- and post-eccentricity adjustment 
#           - Eccentricity maps on cortical surface for pre- and post-eccentricity adjustment 
#           - Cortical magnification factor curve for pre-, post-eccentricity adjustment, H&H model

print("IT WORKS")

# Initialize parser
parser = argparse.ArgumentParser(description="Adjusting Eccentricity in Benson eccentricity map")
# Add mandatory arguments 
parser.add_argument('-i1', '--input1', type=str, required=True, help='FreeSurfer subject ID') 
parser.add_argument('-i2', '--input2', type=str, required=True, help='Hemisphere label (e.g., lh)') 
parser.add_argument('-i3', '--input3', type=str, required=True, help='Eccentricity file and location within the FreeSurfer directory (e.g., surf/lh.benson14_eccen.mgz (or .mgh))') 
parser.add_argument('-i4', '--input4', type=str, required=True, help='Area file and location within the FreeSurfer directory (e.g., surf/lh.benson14_v_area.mgz (or .mgh))') 
parser.add_argument('-i5', '--input5', type=str, required=True, help='Output Directory')
parser.add_argument('-i6', '--input6', type=str, required=True, default='.png', help='Format for saving figures (default = .png)')
# Add optional arguments 
parser.add_argument('-i7', '--input7', type=float, required=False, default=0.75, help='Shape (default = 0.75)')
parser.add_argument('-i8', '--input8', type=float, required=False, default=0, help='Minimum eccentricity (default = 0)')
parser.add_argument('-i9', '--input9', type=float, required=False,default=90, help='Maximum eccentricity (default = 90)')
parser.add_argument('-i10', '--input10', type=float, required=False, default=0, help='Existing pRF maps (default = 0)')

# Retrieve arguments
args = parser.parse_args()
input_subject_id = args.input1
hemlabel = args.input2
input_eccen = args.input3
input_varea = args.input4
output_dir = args.input5
figformat = args.input6
shape = args.input7
min_ecc = args.input8
max_ecc = args.input9
pRFexist = args.input10

# Load in data:
sub = ny.freesurfer_subject(input_subject_id) # e.g., <path-to-freesurfer-subjects-folder>
eccen = ny.load(sub.path_join(input_eccen))# e.g., input_eccen = surf/lh.benson14_eccen.mgz (or .mgh)
label = ny.load(sub.path_join(input_varea)) # e.g., input_varea = surf/lh.benson14_varea.mgz (or .mgh)

# Filter vertices of interest once we have loaded the subject (sub), hemisphere name (hemlabel)
# the eccentricity map (eccen) and the labels or "varea" (label):
hem = sub.hemis[hemlabel.lower()]
mask = (label == 1)
#mesh = hem.surface('midgray') # could be an option, but nbd
v1_eccen = eccen[mask]
v1_sarea = hem.prop('midgray_surface_area')[mask]

# Perform the eccentricity adjustment
print("\nSubject ID:", sub.name)
adjusted_eccen,scale = adjust_eccen_in_v1(v1_eccen, v1_sarea, shape, min_ecc, max_ecc)
# ^^^ optionally also include shape, min_eccen, max_eccen

# Print out number of vertices in V1 mask, V1 surface area, initial and adjusted eccentricity range
print(hemlabel.upper(), " – V1 Surface Area (in mm\u00b2): ", v1_sarea.sum(), " (scale = ", scale, "; number of vertices = ", v1_sarea.size, ")", sep='')
print(hemlabel.upper(), ' - Benson14 Eccentricity (in degrees): min = ', v1_eccen.min(), ', max = ', v1_eccen.max(), sep='')
print(hemlabel.upper(), ' - Adjusted Eccentricity (in degrees): min = ', adjusted_eccen.min(), ', max = ', adjusted_eccen.max(), sep='') 

# Write out adjusted_eccen output as MGH file
# (This will need a command-line option to know the output directory)
output_fname = output_dir + '/' + hemlabel.lower() + '.V1-adjusted-eccen.mgh'
ny.save(output_fname, adjusted_eccen)

# Write out vertex labels of selected vertices on cortical surface into a txt file
# This is used to link back to vertices indices in SamSrf Srf mat file
# (This will need a command-line option to know the output directory)
output_fname = output_dir + '/' + hemlabel.lower() + '-V1-adjusted-eccen-VtxLabels.txt'
v1_vtxidx = hem.labels[mask] # Vertex Index of selected vertices on cortical surface
np.savetxt(output_fname, v1_vtxidx, fmt="%d")

# Plot the Benson original versus Benson adjusted eccentricity values 
plot_fname = output_dir + '/' + hemlabel.lower() + '-V1-originalvsajusted' + figformat
plot_originalvsadjusted(v1_eccen, adjusted_eccen, fname=plot_fname)

# Plot the Benson original versus Benson adjusted eccentricity distributions 
plot_fname = output_dir + '/' + hemlabel.lower() + '-V1-distributionECCvalues' + figformat
plot_distributionECCvalues(v1_eccen, adjusted_eccen, v1_sarea, scale,shape=shape,min_eccen=min_ecc,max_eccen=max_ecc,fname=plot_fname)

# Plot the eccentricity maps on cortical surface for Benson original, Benson adjusted, and pRF (if existpRF=1) 
plot_fname = output_dir + '/' + hemlabel.lower() + '-V1-comparisonCorticalECC' + figformat
plot_comparisonCorticalECC(hem,v1_vtxidx,v1_eccen, adjusted_eccen,existpRF=pRFexist,fname=plot_fname)

# Plot the cortical magnification factor curve for Benson original, Benson adjusted, H&H model, and pRF (if existpRF=1)
plot_fname = output_dir + '/' + hemlabel.lower() + '-V1-CMF' + figformat
plot_CMF(v1_sarea,v1_eccen,adjusted_eccen,scale,shape=shape,existpRF=pRFexist,fname=plot_fname)


