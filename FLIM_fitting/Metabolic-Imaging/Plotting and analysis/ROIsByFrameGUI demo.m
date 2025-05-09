daypath = 'Z:\Lab\Tim\Boston IVF\Discarded_study\2019-07-05 40X zscans, 20X overnight\';
cd(daypath)

acqpath = [daypath 's1_a1_zscan'];

%% Standard functions performed as usual:

MultiD_sdt_sort_Lab(acqpath);
MultiD_tiff_convert(acqpath);
Ilastik_TrainImGen;
IlastikProbMap_shell(acqpath,'Mod_N_nrm');
IlastikProbMap_shell(acqpath,'Mod_F_nrm');



%% to run ICM + troph together, first run the analysis as usual
Make_Masks_from_IlastikProb(acqpath);
FLIM_decay_from_probmap(acqpath); 
FLIM_batch_fitting_probmaps(acqpath);

% Will default to 'JointMasks' and save outputs with that label.

% Proceed to rerun with separate ICM and troph ROIs...

%% Run ROIsByFrameGUI
ROIsByFrameGUI
% ManROIs saves to acq path. Drawn ICM ROIs are auto-saved as you go.
% When you are done drawing all the ICMs, click 'Save Troph' to save and
% additional set of ROIs that are just the inverse of the ICMs.

% Run 'Make_Masks_from_IlastikProb' with a new input, 'ManROILab'. First
% run on ICM:
Make_Masks_from_IlastikProb(acqpath,-1,-1,-1,-1,-1,-1,-1,-1,-1,'ICM');

% MaskLab becomes 'MasksICM', and since the mask is 'joint', label is
% 'JointMasksICM'
FLIM_decay_from_probmap(acqpath,-1,'JointICMMasks'); 
% FLIMGetFixedShiftsBGs(acqpath) % Calculate average shifts and BG in case you want to fix them
FLIM_batch_fitting_probmaps(acqpath,-1,'JointICMMasks');

% Repeat the analysis for the trophectoderm regions. Results will be saved
% with separate labels.
% NOTE: we have to reset the multiD_indices.mat file by running
% 'MultiD_sdt_sort_Lab.m' again. This is because the previous Make_Masks
% call would have found a lot of empty frames and flagged those frames for
% exclusion.
MultiD_sdt_sort_Lab(acqpath);
MultiD_tiff_convert(acqpath);
Make_Masks_from_IlastikProb(acqpath,-1,-1,-1,-1,-1,-1,-1,-1,-1,'Troph');
FLIM_decay_from_probmap(acqpath,-1,'JointTrophMasks'); 
% FLIMGetFixedShiftsBGs(acqpath) % Calculate average shifts and BG in case you want to fix them
FLIM_batch_fitting_probmaps(acqpath,-1,'JointTrophMasks');

%% Proceed with MDpars and comparative plots
FLIM_master_list_and_MDpars(daypath)
