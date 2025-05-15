%adjustGait_subjlist_detectver3
%this script simply holds the hard-coded adjustments per individual

usetype=0; % default 
% usetype=1% details below. catches reduced step lengths for shorter

% participants.
if strcmp(subjID, 'ABB')
    usetype=1;
elseif strcmp(subjID,'AC') && itrial==28
    usetype=1;
elseif strcmp(subjID,'AH') && itrial== 61
    usetype=1;
elseif strcmp(subjID, 'AS')
    usetype=1; 
elseif strcmp(subjID, 'BT')
    usetype=1;
elseif strcmp(subjID, 'CP') 

elseif strcmp(subjID, 'DA')
    usetype=1;
elseif strcmp(subjID, 'GC') && ismember(itrial, [41,58,189])
    usetype=1;
elseif strcmp(subjID, 'KP') && ismember(itrial, [22,28])
    usetype=1;
    pkdist= 20;
elseif strcmp(subjID, 'JC') && ismember(itrial, [96,111])
    pkdist= 15;
    pkheight = .01;
elseif strcmp(subjID, 'JH') && itrial==103
    pkdist= 15;
    pkheight = .01;
elseif strcmp(subjID, 'JT') && ismember(itrial,[104,141])
    pkdist= 15;
    pkheight = .01;

elseif strcmp(subjID, 'JZ')
    pkdist= 15;
    pkheight = .01;
elseif strcmp(subjID, 'LJTN') && itrial==101
    pkdist= 15;
    pkheight = .01;
elseif strcmp(subjID,'LJ')
    usetype=1;
elseif strcmp(subjID, 'LL2')
    usetype=1;
elseif strcmp(subjID, 'LYC')
    usetype=1;

elseif strcmp(subjID, 'MCC') && ismember(itrial ,[67,75,107])
   usetype=1;
elseif strcmp(subjID, 'MD2')
    usetype=1;
elseif strcmp(subjID, 'ML')
    usetype=1;
elseif strcmp(subjID, 'NL')
    usetype=1;

elseif strcmp(subjID, 'NT')
   usetype=1;

elseif strcmp(subjID, 'SF')
    usetype=1;
elseif strcmp(subjID, 'SO')
    usetype=1;
elseif strcmp(subjID, 'YL')
    usetype=1;
elseif strcmp(subjID, 'YW')
    usetype=1;
elseif strcmp(subjID, 'ZT')
    usetype=1;
else % back to original:
    pkdist = ceil(pkduration*Fs);
    pkheight = .02;
end

%% change search params:
if usetype==1
    pkdist=15;
    pkheight=.005;
end