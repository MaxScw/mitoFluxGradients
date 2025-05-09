
% On Metascope
FLIM_IllProfCal(file);
HWPPowerCalMaiTai(file);
IntensitymatrixTiffs(acqpath,WrVal,Ch2Bool,Bit16Bool); % Generic intensity converstion

CopyName2MultiDIndices % For reanalyzing older data sets

%% sdt sort and tiff convert
MultiD_sdt_sort_Lab(acqpath);
MultiD_tiff_convert(acqpath, SHGBool, WellCropBool, bit16Bool);

% Frame exclude
MultiD_exclude_frames(acqpath,exframes,posnum);

%% Segmentation and probability map options
Ilastik_TrainImGen;
IlastikProbMap_shell(acqpath,ModName,poss,Frs,thresh);
IlastikProbMapCons(acqpath,SaveBool) % Optional: if you want to try to constrain the prob masses, e.g. for perturbations

Make_Polygon_Crop(acqpath) % Optional: if there is some stationary crap you want to crop out
ROIsByFrameGUI % Optional: if we want to go frame-by-frame and manually set regions of interest

Make_Masks_from_IlastikProb(acqpath,poss,area_cuts,thresh,Gblur,SingMaskBool,MaskLab,MergeBool,RemDeadCells);
% NOTE: also calculates background intensities. To plot, use 'BG_ints_plot(acqpath,pos)'
% Alt simple segmentation, e.g. for including all pixels from an FAD sample... just use kmeans and have 1 group.
Make_Mask_Generic_WoW(acqpath,poss,area_cuts,Gblur,nClusts,KeepClusts,thresh,NumWells,frames)

%% Decay extraction and fitting
FLIM_decay_from_probmap(acqpath,poss,MaskLab,MinFrames,ProbThreshMeth,MasksToDo); 
FLIMGetFixedPars([daypath 'sample_subfolder']) % Calculate average shifts and BG in case you want to fix them
FLIM_batch_fitting_probmaps(acqpath,poss,MaskLab,Bin,fxshftBool,fxBgBool,IRFname,nexpo);
% Fitting calls 
FLIMAcqParamPlot(acqpath)


%% High level plotting and analysis
FLIM_master_list_and_MDpars(daypath,Zrange,trange,ExcludeFolders,Label) % Create master list and MDpars for whole day
% FLIM_batch_averages % previous, obsolete
% Use '-\MATLAB PROGRAMS\Metabolic Imaging\Plotting and analysis' routines to make plots from Masterlist and MDpars
FLIM_mdpars_4Plots_Bswarm_Ttests_3dPlots
FLIM_mdpars_TimePlots
FLIM_mdpars_ZPlots

% To do pixel-by-pixel fits to generate FLIM param images:
FLIM_param_imgs % Unfinished... but would be useful.


load('C:\Google Drive\MATLAB PROGRAMS\Metabolic Imaging\Plotting and analysis\colorblind_colormap.mat')
%% Other 
% Division time GUI
divisions_GUI

FLIM_CatT_mdpars(acqs); % For concatenating separate acquisitions that are really the same experiment
LaserOPplot(acqs); % For verifying that laser output didn't drift during acquisition. (logged in AcqSdtCorr.txt)

% Additional sometimes-useful functions
IlastikProbMapCons(acqpath,SaveBool); % Useful for perturbation experiments. 

% Picking color schemes or gradients for plots
SelectingColors4Colorblind_wGradients

%% Useful snippets



% Find samples and time points with particular values 
Frs = squeeze(MDpars(4,1,:,:,1,1,1));
[rw,T] = find(Frs>.5); % <- Condition
[ml(rw,1) num2cell(T)]

% Label cell for doing plots, etc, in loops
Labs = {'Nirr','Nbound','Ntau1','Ntau2','Firr','Fbound','Ftau1','Ftau2'};

% Find plotyy handles after the fact.
AX=findall(0,'type','axes'); AX(1).YLim = [.17 .3]; AX(2).YLim = [.19 .27];










% %%%%%%% Other older stuff %%%%%%%
% 
% % Segmentation 
% Make_Mask_Generic_WoW(acqpath,poss,area_cuts,Gblur,nClusts,KeepClusts,thresh,NumWells,frames);
% Make_Mask_Generic_Kmeans(acqpath) % All embs in one big mask
% 
% % Prev from masks
% FLIM_decay_from_mask([acqpath 's1_a1'],-1,-1,10);
% % FLIMGetFixedShifts([acqpath 's1_a1']); % <- GET SHIFTS! <- NOPE! FIT SHIFTS
% FLIM_batch_fitting([acqpath 's1_a1'],-1,0,-1,3);
% 
% % Weka (Obsolete)
% MatIJdatpathconv(acqpath)
% Make_Masks_from_IJWeka(acqpath,poss,area_cuts,thresh,Gblur,SingMaskBool,MaskLab,RemDeadCells); % Using Weka seg in ImageJ
% 
% 
% FLIM_decay_from_mask(acqpath);
% FLIM_batch_fitting(acqpath);
% % C elegans:
% MultiD_sdt_sort(acqpath);
% MultiD_sdt_sort_from_meta(acqpath); 
% MultiD_sdt_sort_dualTL(acqpath,tl);
% 
% MultiD_tiff_convert(acqpath);
% MultiD_exclude_frames(acqpath,exframes,posnum);
% 
% CalBeadProc_Cembs(acqpath); % If cal beads were used. OBSOLETE
% UserEggCenters(acqpath)
% Make_Masks_TouchingEggs_Drift(acqpath,MasksFrame)
% % or Make_Masks_TouchingEggs (sans drift)
% % or Make_Mask_Generic(acqpath) for zoom in egg
% 
% FLIM_decay_from_mask(acqpath);
% % FLIM_decay_from_thresh(acqpath,threshN,threshF); % Obsolete
% 
% FLIM_batch_fitting(acqpath);
% 
% % Collective analysis
% FLIM_batch_averages(acqpath,Zrange,trange,ExcludeFolders,Label)
% FLIM_batch_TimeArrs(acqpath,Zrange,trange,ExcludeFolders,Label,tstmpBool)
% 
% 
% % Cumulus cells
% MultiD_sdt_sort(acqpath);
% MultiD_sdt_sort_from_meta(acqpath);
% 
% MultiD_tiff_stk_convert(acqpath,AniBool);
% 
% Make_Mask_Generic(acqpath)
% % Can play with 'area_cuts,level_fact', but defaults seem to work
% 
% FLIM_decay_from_mask(acqpath);
% 
% FLIM_batch_fitting(acqpath);
% 
% FLIM_batch_averages(acqpath,Zrange,trange)