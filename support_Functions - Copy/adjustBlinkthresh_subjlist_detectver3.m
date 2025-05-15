function blinkthresh =adjustBlinkthresh_subjlist_detectver3(subjID)

usetype=0; % default 
% usetype=1% details below. catches reduced step lengths for shorter

% participants.

if strcmp(subjID, 'ABB')
elseif strcmp(subjID,'AH') 
elseif strcmp(subjID, 'AS')
elseif strcmp(subjID, 'BL')

elseif strcmp(subjID, 'BT')
elseif strcmp(subjID, 'CC')
    usetype=2;

elseif strcmp(subjID, 'CP') 
elseif strcmp(subjID, 'DA')
elseif strcmp(subjID, 'GC')
elseif strcmp(subjID, 'KP') 
elseif strcmp(subjID, 'JC') 
elseif strcmp(subjID, 'JH') 
elseif strcmp(subjID, 'JT') 
    usetype=1; % lower threshold

elseif strcmp(subjID, 'JZ')
elseif strcmp(subjID, 'LJTN') 

elseif strcmp(subjID,'LJ')
elseif strcmp(subjID, 'LL2')
elseif strcmp(subjID, 'LYC')
elseif strcmp(subjID, 'MCC') 
elseif strcmp(subjID, 'MD2')
elseif strcmp(subjID, 'ML')
elseif strcmp(subjID, 'NL')

elseif strcmp(subjID, 'NT')

elseif strcmp(subjID, 'SF')
elseif strcmp(subjID, 'SO')
elseif strcmp(subjID, 'YL')
elseif strcmp(subjID, 'YW')
elseif strcmp(subjID, 'ZT')
else % back to original:    
end



%% change search params:
switch usetype
    case 0
    blinkthresh=0.4; %default

    case 1
    blinkthresh= 0.3;
    case 2
    blinkthresh= 0.5;
        
end
end
