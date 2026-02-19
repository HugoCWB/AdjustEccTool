# AdjustEccTool - Benson Map Eccentricity Adjustment Tool (Docker)

This repository provides a Docker-based tool for correcting the eccentricity distribution derived from the **Benson atlas** (via *neuropythy*). The method implements the procedure described in:

**Chow‑Wing‑Bom et al. (2025). _Mapping Visual Contrast Sensitivity and Vision Loss Across the Visual Field with Model‑Based fMRI_. eLife.**

The tool loads the neuropythy outputs from a participant’s FreeSurfer folder, extracts V1 eccentricity and area maps, applies the eccentricity‑adjustment algorithm, and saves corrected maps, vertex labels, and diagnostic plots.

--- 

## TODO / Future Extensions

In `__init__.py` and `__main__.py`: The current plotting functions assume Benson atlas inputs. Consider adding support for pRF-derived eccentricity maps, allowing users to overlay or compare pRF estimates with Benson-adjusted values in the diagnostic plots.

---

## Requirements

### 1. Benson atlas outputs
The script requires two files per hemisphere:

- **Eccentricity map:**  
  `?h.benson14_eccen.mgh` (or .mgz)
- **V1 area label map:**  
  `?h.benson14_varea.mgh` (or .mgz)

These must:
- be located in the participant’s FreeSurfer directory  
  (`surf/` for `.mgz/.mgh`, or `label/` for `.label`)
- contain **value 1 for V1 vertices only** in the V1 area map (i.e., `?h.benson14_varea.mgh` (or .mgz))

### 2. FreeSurfer SUBJECTS_DIR
The script reads the participant’s FreeSurfer reconstruction.  
Before running the container, ensure to set SUBJECTS_DIR directory accordingly beforehand (see below).

---

## Before You Start
Ensure that **Docker Desktop** is running.

---

## Building the Docker Image (~9-10 minutes)

Neuropythy’s base image is **AMD64**, so the image must be built using that platform. To do so, use the following command:

```
docker build --platform=linux/amd64 -t <DOCKER PROJECT IMAGE NAME> "<PATH TO DOCKER SCRIPT FOLDER>"
```

For instance, if you want to name your docker image smriproject and you downloaded this toolbox in /Users/.../Desktop/StructuralMRI: 
```
docker build --platform=linux/amd64 -t smriproject "/Users/.../Desktop/StructuralMRI"
```

Alternatively, you can go to the downloaded directory (e.g., cd /Users/.../Desktop/StructuralMRI) and execute the following: 
```
docker build --platform=linux/amd64 -t smriproject .
```
(**Note: the final . is important**)

To verify your image installation:
```
docker run --platform=linux/amd64 smriproject --help
```

---

## Running the Tool on a Participant

```
docker run --platform=linux/amd64 \
  -v <FREESURFER SUBJECTS DIRECTORY PATH>:/subjects \
  -v <PARTICIPANT RESULT FOLDER DIRECTORY PATH>:/out \
  -e SUBJECTS_DIR=/subjects \
  smriproject \
  -i1 <PARTICIPANT FREESURFER ID> \
  -i2 ?h \
  -i3 surf/?h.benson14_eccen.mgh \
  -i4 surf/?h.benson14_varea.mgh \
  -i5 /out \
  -i6 '.png'
```

Hemisphere notation: Replace ?h with the appropriate hemisphere label (Left hemisphere: use lh; Right hemisphere: use rh)
For example, for the left hemisphere, the relevant arguments become:
`-i2 lh`, `-i3 surf/lh.benson14_eccen.mgh`, and `-i4 surf/lh.benson14_varea.mgh`.

All outputs (adjusted maps, vertex labels, and plots) will be written to the mounted results folder.

Example for left hemisphere from participant C001:
```
docker run --platform=linux/amd64 \
  -v /Users/hugo/freesurfer/subjects:/subjects \
  -v /Users/hugo/results/C001:/out \
  -e SUBJECTS_DIR=/subjects \
  smriproject \
  -i1 C001 \
  -i2 lh \
  -i3 surf/lh.benson14_eccen.mgh \
  -i4 surf/lh.benson14_varea.mgh \
  -i5 /out \
  -i6 '.png'
```

---

## Cleaning Up Docker Images
You can remove the image via Docker Desktop or:

```
docker images
docker rmi smriproject
```
