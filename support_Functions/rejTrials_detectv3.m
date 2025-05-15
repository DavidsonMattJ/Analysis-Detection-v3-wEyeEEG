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
    case 'ABB_2022-08-12-01-51'
        badtrials= [];
     case 'AC_2022-08-23-02-16' %rejected
        badtrials=[21,61,62];
 
    case 'AH_2022-09-22-08-57'
        badtrials=[ 65];
    case 'AS_2022-09-21-11-22'
        badtrials=[103,138];
    case 'BL_2022-09-20-11-26'; %nice example trials of a clean gait
         badtrials=[]; % no bad
    case 'BT_2022-08-09-01-41'
        badtrials=62;
    
    case 'CC_R_2022-09-21-09-14'; % new
badtrials=[]; % none
    case 'CP_2022-07-01-02-24'
        badtrials=59; % g.

    case 'DA_2022-09-22-02-50'
        badtrials=[];
    case 'EJ_2022-09-21-12-35'
        badtrials=[]; % no bad
    case 'EK_2022-09-23-03-41'
        badtrials = 21;
    case 'EPM_2022-09-22-10-11'
        badtrials=[21,180]; % signal, gait.
    case 'GC_2022-08-02-02-09'
        badtrials = [];
    case 'IRK_2022-09-22-12-35'
        badtrials=[21,48,49,69,171]; %s/g,g,g,g,s

    case 'JC_2022-07-25-04-50' % ?? REJECT 143 trials only
        badtrials = 23; % s.
    case 'JEPM_2022-09-20-02-34'
        badtrials=21;
    case 'JH_2022-08-23-01-07'
        badtrials=[];
    case 'JT_2022-08-08-09-55'
        badtrials=[];
    case 'JW_2022-09-23-02-24'
        badtrials=[];% none
    case 'JZ_2022-08-26-01-09'
        badtrials=157; % s

    case 'KK_2022-08-23-11-18'
        badtrials=[]; %none  
    case 'KP_2022-08-02-03-45' % Some drift across the room.
        badtrials = []; % no bad trials

    case 'LJTN_2022-08-23-10-09'
        badtrials=155;% s
    case 'LJ_2022-08-12-11-25'
        badtrials=[]; % no bad trials.
    case 'LL1_2022-07-25-02-11'
        badtrials=[21,151,152]; % g,g,g,
    case 'LL2_2022-08-25-11-03'
        badtrials=[149,180];% s s

    case 'LYC_2022-09-23-10-19' % Large slant, check overall performance.
        badtrials=[34,118,119]; % s G G

    case 'MCC_2022-08-23-03-21'
        badtrials =[]; %s
    case 'MD1_2022-07-01-05-06'
        badtrials=[];
    case 'MD2_2022-08-09-11-45' % concatenated data set (208 trials total).
        badtrials=[];% no bad trials

    case 'ML_2022-08-12-08-59'
        badtrials= [95,97,136]; % s s s
    case 'NL_2022-07-22-09-54'
        badtrials =[24,54,55,60]; % s, s,s,s

    case 'NT_2022-08-09-02-44'
        badtrials=189;
    case 'QL_2022-08-19-02-04'
        badtrials=23;        % g

    case 'RJS_2022-09-23-09-00'
        badtrials=[];
    case 'SC_2022-09-21-04-00'
        badtrials=[21,91,176];
    case 'SF_2022-08-09-04-08' % Some drift present
        badtrials= 165;
    case 'SO_2022-08-09-09-12'
        badtrials=[];
    case 'TMO_2022-08-12-02-45'
        badtrials=[]; % no bad
    case 'WB_2022-08-19-10-06'
        badtrials=[25,30]; % s s
    case 'YL_2022-09-22-03-46'
        badtrials=[21,133,178];
    case 'YW_2022-08-19-01-04'
        badtrials=[];
    case 'ZT_2022-08-26-02-20' % % some drift
        badtrials=[];
end

    
%%
if ismember(itrial,badtrials)
    disp(['Skipping bad trial ' num2str(itrial) ' for ' subjID]);
    skip=1;
end
%%