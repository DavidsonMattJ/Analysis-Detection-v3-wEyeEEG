% rejTrials_detectv3%

% certain bad trials (identified in figure folder).
% look at 'TrialHeadYPeaks', to see if the tracking dropped out on any
% trials. Or if participants were not walking smoothly.

% abbreviations for rejected trials:
% s = poor signal quality (drop-outs/discontinuities)
% g= poor gait extraction (head tracking unclear).

% step in to reject particular ppant+trial combos.
badtrials=[];
switch subjID
    case 'p01_2022-08-12-01-51'
        badtrials= [];
 
    case 'p02_2022-09-22-08-57'
        badtrials=[ 65];
    case 'p03_2022-09-21-11-22'
        badtrials=[103,138];
    case 'p04_2022-09-20-11-26'; %nice example trials of a clean gait
         badtrials=[]; % no bad
    case 'p05_2022-08-09-01-41'
        badtrials=62;
    
    case 'p06_2022-09-21-09-14'; % new
badtrials=[]; % none
    case 'p07_2022-07-01-02-24'
        badtrials=59; % g.

    case 'p08_2022-09-22-02-50'
        badtrials=[];
    case 'p09_2022-09-21-12-35'
        badtrials=[]; % no bad
    case 'p10_2022-09-23-03-41'
        badtrials = 21;
    case 'p11_2022-09-22-10-11'
        badtrials=[21,180]; % signal, gait.
    case 'p12_2022-08-02-02-09'
        badtrials = [];
    case 'p13_2022-09-22-12-35'
        badtrials=[21,48,49,69,171]; %s/g,g,g,g,s

    case 'p14_2022-08-23-01-07'
        badtrials=[];
    case 'p15_2022-08-08-09-55'
        badtrials=[];
    case 'p16_2022-08-23-11-18'
        badtrials=[]; %none  
    case 'p17_2022-08-02-03-45' %
        badtrials = []; % no bad trials

    case 'p18_2022-08-23-10-09'
        badtrials=155;% s
    case 'p19_2022-08-12-11-25'
        badtrials=[]; % no bad trials.
    case 'p20_2022-07-25-02-11'
        badtrials=[21,151,152]; % g,g,g,
    
    case 'p21_2022-09-23-10-19' % 
        badtrials=[34,118,119]; % s G G

    case 'p22_2022-08-23-03-21'
        badtrials =[]; %s
    case 'p23_2022-07-01-05-06'
        badtrials=[];
    case 'p24_2022-08-09-11-45' % concatenated data set (208 trials total).
        badtrials=[];% no bad trials

    case 'p25_2022-08-12-08-59'
        badtrials= [95,97,136]; % s s s
    case 'p26_2022-07-22-09-54'
        badtrials =[24,54,55,60]; % s, s,s,s

    case 'p27_2022-08-09-02-44'
        badtrials=189;
    case 'p28_2022-08-19-02-04'
        badtrials=23;        % g

    case 'p29_2022-09-23-09-00'
        badtrials=[];
    case 'p30_2022-09-21-04-00'
        badtrials=[21,91,176];
    case 'p31_2022-08-09-04-08' % 
        badtrials= 165;
    case 'p32_2022-08-09-09-12'
        badtrials=[];
    case 'p33_2022-08-12-02-45'
        badtrials=[]; % no bad
    case 'p34_2022-08-19-10-06'
        badtrials=[25,30]; % s s
    case 'p35_2022-09-22-03-46'
        badtrials=[21,133,178];
    case 'p36_2022-08-26-02-20' % % 
        badtrials=[];
end

    
%%
if ismember(itrial,badtrials)
    disp(['Skipping bad trial ' num2str(itrial) ' for ' subjID]);
    skip=1;
end
%%