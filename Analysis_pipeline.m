%Analysis pipeline.

% this script simply steps through the analyses in order, with some
% descriptions of output. See the individual scripts for more details.

% add path with analysis/plot scripts to the search directory:


%laptop:
setmydirs_detectv3;
% addpath(genpath([pwd filesep 'support_functions']));


%% Preprocessing:
% Import and wrangle Unity output for Matlab Analysis:
j0_VR_import_plot;

%Find the peaks in head position data to epoch per step/stride
j1_findpeaksHeadPos;
%%
% detect and interpolate blinks from eye-movement data
j1A_cleanEyeMovementpertrial;

% detect onsets of saccades per trial.
j1B_timestamp_Saccades_pertrial;

%% step/stride analysis
% classify events (target onsets, RTs, etc), with respect to percentage
% within a step...
j2_appendgaitcycleData;
% ... add whether step was left/right foot
j2B_appendgaitcycleData_LRft;
% ... and percentage with respect to linked steps (stride)
j2C_append_doubGCData;

% perform the same analysis (step /stride allocation) for the saccade table
j2D_appendSaccades_toGait;
j2E_appendSaccades_doubGC;

% with the step start/fin times clearly defined, epoch all time-series
% of head position for later plots:
j3_epochGait_timeseries;

%% Participant and group-level effects
% calculate participant level effects (DVs per step/stride %),
% included binned versions 
% aggregate into a group data structure.
j4_cleanPFX_concatGFX;
% same for saccades (separate table so need a separate script).
j4A_crunchPFX_concatGFX_saccades;

% Supplementary Figure: DVs relative to saccade onset.
j4B_crunchPFX_concatGFX_relativetoSaccades; 

%% Permutation tests
% create null effects by shufflinf the event position uniformly per
% permutation.
j5_createNullGFX;
% same for saccade data.
j5_createNullGFX_sacc;

% repeat main analysis (Fourier Fits), on each permutation of our null
% data (group-level effects). Also test observed data over a wide range of
% frequencies (in cycles per stride).
j8_testfourier_Obs_Null_GFX;

%perform same, testing each individual participant against their shuffled
%likelihood:
j10_testfourier_Obs_Null_PFX;
%% >>>>>>>                  Plotting:
%basic data inspection (not used in MS).
% j6__Plot_Data_withinGait_basic;

% this performs individual level psychometric fits, and group level results
% for MS figure:
j7__Plot_Data_withinGaits_PFits_v3;  % (Manuscript Figure 2)

% plot main results (observed data, Fourier fits, shuffled results, etc).
j9__Plot_Data_FourierFits_GFX;       % (Manuscript Figure 3)

%plot main results (participant level effects), includes prevalence
%calculations and plots.
j11__Plot_Data_FourierFits_PFX;       % (Manuscript Figures 4 and 5)
