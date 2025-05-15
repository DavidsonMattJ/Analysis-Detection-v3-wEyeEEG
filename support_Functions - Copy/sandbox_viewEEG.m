%first attempt at viewing walking EEG data.
clear all ; close all;
%reorient to datadir:
fncdir = pwd; 
cd ../data_Raw/DSI
dsidatadir=pwd;
%%
%%
%list available files:
filesavail = dir([pwd filesep '*raw.csv']);
disp({filesavail(:).name}');
%%
% for ifile = 1%:length(filesavail)
eeglab;
cd(dsidatadir)
EEG=pop_WearableSensing_ljb(filesavail(1).name);
eeglab redraw

% end %% end ifile 