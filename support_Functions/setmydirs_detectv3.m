%setmydirs_detectv3

% PC / laptop
datadir= 'C:\Users\mdav0285\Documents\GitHub\Analysis-Detection-v3-wEyeEEG\data_Raw\Unity';
%Mac:
% datadir= '/Users/matthewdavidson/Documents/GitHub/Analysis-Detection-v3-wEyeEEG/data_Raw/Unity'

% work computer:
% datadir='C:\Users\mdav0285\Documents\GitHub\Analysis-Detection-v3-wEyeEEG\data_Raw\Unity'
cd(datadir)
cd ../../data_Processed
procdatadir = pwd;
cd ../
addpath([pwd filesep 'support_Functions']);
cd 'Figures'; figdir= pwd;

% cd(procdatadir)
% pfols= dir([pwd  filesep '*summary_data.mat']);
% nsubs= length(pfols);
% %show ppant list:
% tr= table((1:length(pfols))',{pfols(:).name}' );
% disp(tr)
