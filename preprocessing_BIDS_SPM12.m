function preprocessing_BIDS (sub)

%{
This functions runs a standard preprocessing pipeline in SPM12. Make sure SPM12 is in your matlab path. 
 
The scripts requires a BIDS formatted dataset. See http://bids.neuroimaging.io/
It also assumes you acquired a T1 analtomical image, at least 1 functional run saved as a 4D nifti file, and a field map (which I highly recommend).
 
In order to run this script, call the preprocessing_BIDS function, giving it one subject numer.
The script will then run the preprocessing for that subejct.   

23/01/2018 v.0.9.0: - stable version

David Wisniewski (david.wisniewski@ugent.be)
 %}


%% FLAGS 
%{
The script is partly modular. You can switch specific preprocessing step on (=1) or off (=0) here. 
The order of the steps is fixed at the moment though. It is assumed that you run your preprocessing
in the order: deface -> fieldmap -> realign and unwarp -> slice time correction -> coregistration -> normalization -> smoothing
%}

do_deface = 1;          % deface the anatomical image to anonymize the data
do_fieldmap = 1;        % calculate fieldmap (vdm file)
do_realign_unwarp = 1;  % realignment (get rid of movement artifacts) + unwarping (account for disortions in the magnetic field using the fieldmap)
do_slice_time = 1;      % slice time correction
do_coregister = 1;      % coregistration of T1 to EPI images
do_norm=1;              % normalize images to MNI space
do_smooth=1;            % smooth images    

%% PARAMETERS 
%{
specify the paramters of your dataset here
%}
% specify the number of functional runs acquired
nruns=5;                
% slice time corrections
param.st.nslices = 33; % number of slices per volume
param.st.tr = 2; %TR in sec
param.st.ta = 1.93939393939394; %time of acquisition for one slice, i.e. TR./nslices
param.st.so = [33:-1:1]; % slice acquisition order
param.st.refslice = 17; % reference slice, for fMRI usually the middle slice
% smoothing
param.smooth.fwhm = [8 8 8]; % define the smoothing kernel in mm

%% DIRECTORIES
% convert subject number into a string for easier handling
substr=num2str(sub,'%.2d');
% specify the directory in which your BIDS compatible raw data are localized
basedir = '\\folder\BIDS';
% specify the directory to which all derivatives of the raw data are to be saved
% side note: BIDS splits your data into the raw unprocessed images, and derivatives. The latter are all images derived from the raw data, e.g. preprocessed images.
derivativesdir = '\\folder\derivatives';
subdir=['\sub-' substr];


%% RUN THE ANALYSIS

%% DEFACE T1
% in order to guarantee anaonymity, we deface the anatomical images
if do_deface
    % run the defacing step
    matlabbatch{1}.spm.util.deface.images = cellstr(spm_select('FPList',[basedir subdir '\anat\'],'^sub.*.nii$'));
    spm_jobman('run', matlabbatch);
    clear matlabbatch;
    % file management
    % delete original file with the face still intact
    delete (spm_select('FPList',[basedir subdir '\anat\'],'^sub.*.nii$'));
    % rename anonimized file
    fname=spm_select('FPList',[basedir subdir '\anat\'],'^anon.*.nii$');
    newname= regexprep(fname,'anon_','');
    movefile (fname,newname);
    % copy anonymized file to \derivates folder
    destination = regexprep(newname,'BIDS','derivatives');
    copyfile (newname,destination) ;
end


%% FIELDMAP
% estimate a field map 
if do_fieldmap
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.phase =  cellstr(spm_select('FPList',[basedir subdir '\fmap\'],'^sub.*phase2.nii$'));
        % of the two magnitude images, use the one with the shorter TE,which will have greater signal
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.data.presubphasemag.magnitude = cellstr(spm_select('FPList',[basedir subdir '\fmap\'],'^sub.*magnitude2.nii$'));
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.et = [5.19 7.56];
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.maskbrain = 0;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.blipdir = -1;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.tert = 32.64;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.epifm = 0;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.ajm = 0;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.method = 'Mark3D';
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.fwhm = 10;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.pad = 0;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.uflags.ws = 1;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.template = {'F:\Data\Programs Scripts\spm12\toolbox\FieldMap\T1.nii'};
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.fwhm = 5;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.nerode = 2;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.ndilate = 4;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.thresh = 0.5;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.defaults.defaultsval.mflags.reg = 0.02;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.session.epi = cellstr(spm_select('FPList',[basedir subdir '\func'],'^sub.*run-1_bold.nii$'));
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchvdm = 1;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.sessname = 'session';
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.writeunwarped = 0;
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.anat = '';
        matlabbatch{1}.spm.tools.fieldmap.calculatevdm.subj.matchanat = 0;
        spm_jobman('run', matlabbatch);
        clear matlabbatch;
        
        %file management: move newly generated files into derivatives
        %folder
        fname=spm_select('FPList',[basedir subdir '\fmap\'],'^fpm.*.nii$');
        destination = regexprep(fname,'BIDS','derivatives');
        movefile (fname,destination)  ;
        fname=spm_select('FPList',[basedir subdir '\fmap\'],'^scsub.*.nii$');
        destination = regexprep(fname,'BIDS','derivatives');
        movefile (fname,destination)  ;
        fname=spm_select('FPList',[basedir subdir '\fmap\'],'^vdm.*.nii$');
        destination = regexprep(fname,'BIDS','derivatives');
        movefile (fname,destination)  ;
        fname=spm_select('FPList',[basedir subdir '\func\'],'^usub.*.nii$');
        destination = regexprep(fname,'BIDS','derivatives');
        movefile (fname,destination)  ;
        fname=spm_select('FPList',[basedir subdir '\func\'],'^wfmag.*.nii$');
        destination = regexprep(fname,'BIDS','derivatives');
        movefile (fname,destination);
    
end

%% REALIGN & UNWARP
% correct for movement artifacts and magnetic field inhomogeneities. The latter requires the estimated field map from the previous step.
if do_realign_unwarp
    clear matlabbatch;
    
    vdmfile = cellstr(spm_select('FPList',[derivativesdir subdir '\fmap\'],'^vdm.*\phase2.nii$'));
    for run=1:nruns
        filter = ['^sub.*run-' num2str(run) '_bold.nii$'] ;
        matlabbatch{1}.spm.spatial.realignunwarp.data(run).scans = cellstr(spm_select('FPList',[basedir subdir '\func\'],filter));
        matlabbatch{1}.spm.spatial.realignunwarp.data(run).pmscan = vdmfile;
    end
       
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.quality = 0.9;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.sep = 3;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.fwhm = 5;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.rtm = 0;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.einterp = 4;
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.ewrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realignunwarp.eoptions.weight = '';
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.basfcn = [12 12];
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.regorder = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.lambda = 100000;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.jm = 0;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.fot = [4 5];
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.sot = [];
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.uwfwhm = 4;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.rem = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.noi = 5;
    matlabbatch{1}.spm.spatial.realignunwarp.uweoptions.expround = 'Average';
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.uwwhich = [2 1];
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.rinterp = 4;
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.mask = 1;
    matlabbatch{1}.spm.spatial.realignunwarp.uwroptions.prefix = 'u';

    spm_jobman('run', matlabbatch);
    clear matlabbatch;
    
    % file management
	% move files to /derivatives folder
    fname=cellstr(spm_select('FPList',[basedir subdir '\func\'],'^usub.*.nii$'));
    destination = regexprep(fname,'BIDS','derivatives');
    for i = 1:length(fname)
        movefile (fname{i},destination{i});
    end
    fname=spm_select('FPList',[basedir subdir '\func\'],'^mean.*.nii$');
    destination = regexprep(fname,'BIDS','derivatives');
    movefile (fname,destination);
    fname=cellstr(spm_select('FPList',[basedir subdir '\func\'],'^*.mat$'));
    destination = regexprep(fname,'BIDS','derivatives');
    for i = 1:length(fname)
        movefile (fname{i},destination{i});
    end
    fname=cellstr(spm_select('FPList',[basedir subdir '\func\'],'^*.txt$'));
    destination = regexprep(fname,'BIDS','derivatives');
    for i = 1:length(fname)
        movefile (fname{i},destination{i});
    end
end

%% SLICE TIME CORRECTION
% corrects for slice acquisition order
if do_slice_time
    clear matlabbatch img_files;
    img_files = [{cellstr(spm_select('FPList', [derivativesdir subdir '\func\'], '^usub.*.nii$'))}];        
    matlabbatch{1}.spm.temporal.st.scans = img_files;
    matlabbatch{1}.spm.temporal.st.nslices = param.st.nslices;    % check this again
    matlabbatch{1}.spm.temporal.st.tr = param.st.tr;
    matlabbatch{1}.spm.temporal.st.ta = param.st.ta;
    matlabbatch{1}.spm.temporal.st.so = param.st.so; % check this again whther this is actually true. also see 
    % https://static.healthcare.siemens.com/siemens_hwem-hwem_ssxa_websites-context-root/wcm/idc/groups/public/@global/@imaging/@mri/documents/download/mdaz/nzmy/~edisp/mri_60_graessner-01646277.pdf
    matlabbatch{1}.spm.temporal.st.refslice = param.st.refslice;
    matlabbatch{1}.spm.temporal.st.prefix = 'a';
    spm_jobman('run', matlabbatch);
    clear matlabbatch;
end


%% COREGISTRATION
% in order to normalize the functional images, first coregister the T1 to the mean EPI image.
if do_coregister
    clear matlabbatch;
    % use the mean EPI as the stationary reference image
    matlabbatch{1}.spm.spatial.coreg.estwrite.ref = cellstr(spm_select('FPList', [derivativesdir subdir '\func\'], '^mean.*.nii$')); 
    % use the MPRAGE as the source image that will be warped to match the mean EPI
    matlabbatch{1}.spm.spatial.coreg.estwrite.source = cellstr(spm_select('FPList', [basedir subdir '\anat\'], '^sub.*.nii$'));
    matlabbatch{1}.spm.spatial.coreg.estwrite.other = {''};
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.cost_fun = 'nmi';
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.sep = [4 2];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.tol = [0.02 0.02 0.02 0.001 0.001 0.001 0.01 0.01 0.01 0.001 0.001 0.001];
    matlabbatch{1}.spm.spatial.coreg.estwrite.eoptions.fwhm = [7 7];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.interp = 4;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.wrap = [0 0 0];
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.mask = 0;
    matlabbatch{1}.spm.spatial.coreg.estwrite.roptions.prefix = 'r';
    spm_jobman('run', matlabbatch);
    clear matlabbatch;
    
    %file management: move newly generated files into derivatives
    %folder
    fname=spm_select('FPList',[basedir subdir '\anat\'],'^rsub.*.nii$');
    destination = regexprep(fname,'BIDS','derivatives');
    movefile (fname,destination)  ;
        
end


%% NORMALIZATION
%normalize the images to MNI space
if do_norm
    clear matlabbatch;
    % pick the images to be normalized
	img_files = [cellstr(spm_select('FPList', [derivativesdir subdir '\func\'], '^ausub.*\.nii$'))];
    % use the coregistered T1 as a reference image      
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.vol = cellstr(spm_select('FPList', [derivativesdir subdir '\anat\'], '^rsub.*.nii$'));
    % define files to normalize
    matlabbatch{1}.spm.spatial.normalise.estwrite.subj.resample = img_files;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasreg = 0.0001;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.biasfwhm = 60;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.tpm = {'F:\Data\Programs Scripts\spm12\tpm\TPM.nii'};
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.affreg = 'mni';
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.reg = [0 0.001 0.5 0.05 0.2];
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.fwhm = 0;
    matlabbatch{1}.spm.spatial.normalise.estwrite.eoptions.samp = 3;
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.bb = [-78 -112 -70
                                                                 78 76 85];
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.vox = [2 2 2];
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.interp = 4;
    matlabbatch{1}.spm.spatial.normalise.estwrite.woptions.prefix = 'w';
    spm_jobman('run', matlabbatch);
    
    
end

%% SMOOTHING
% smooth the images
if do_smooth
    clear matlabbatch;
    img_files = [cellstr(spm_select('FPList', [derivativesdir subdir '\func\'], '^wausub.*.nii$'))];    
    matlabbatch{1}.spm.spatial.smooth.data = img_files;
    matlabbatch{1}.spm.spatial.smooth.fwhm = param.smooth.fwhm;
    matlabbatch{1}.spm.spatial.smooth.dtype = 0;
    matlabbatch{1}.spm.spatial.smooth.im = 0;
    matlabbatch{1}.spm.spatial.smooth.prefix = 's';
    spm_jobman('run', matlabbatch);
    clear matlabbatch;
end

