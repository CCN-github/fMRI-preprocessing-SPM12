# fMRI-preprocessing-SPM12
This is a SPM12 batch script that runs a standard fMRI preprocessing pipeline on a BIDS formatted data-set.

## Included preprocessing steps 
- defacing
- field map esimtation
- realignment and unwarping
- slice-time correction
- normalization
- smoothing

## Requirements
- a BIDS formatted fMRI data-set See http://bids.neuroimaging.io/ for more information.
- a T1 anatomical image
- a field map
- at least 1 functional run (saved as a 4D nifti file)

## Bugs / Improvements
If you find bugs in this script or have suggestions for improvement, please report both here https://github.com/CCN-github-beta/fMRI-preprocessing-SPM12/issues

## Contact
david.wisniewski@ugent.be

