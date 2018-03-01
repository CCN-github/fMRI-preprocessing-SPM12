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
The scripts are written with the imaging sequences used at the GIfMI (http://gifmi.ugent.be/) in mind. 

## Bugs / Improvements
If you find bugs in this script or have suggestions for improvement, please report both here https://github.com/CCN-github/fMRI-preprocessing-SPM12/issues

## Import raw data in the BIDS format
If you are interested in importing your data according to BIDS standards, go check out the scripts here https://github.com/NeuroStat/CustomFormatBIDS

## Contact
david.wisniewski@ugent.be

