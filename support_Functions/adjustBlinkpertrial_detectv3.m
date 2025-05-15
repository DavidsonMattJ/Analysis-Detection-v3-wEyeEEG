% adjustBlinkpertrial_detectv3

% this script just outsources some hardcoded trial adjustments, necessary
% to clean blink data.


if max(blinksAt)> max(blinksEnd)

    blinksAt(end)= []; % remove last

elseif strcmp(subjID,'AS') && ismember(itrial, [12,58,111,121]);

    if itrial==121 % special case (low thresh).
        blinksAt= [13; blinksAt];

    end

elseif strcmp(subjID, 'CC') && ismember(itrial, [1,68,159,164,177]);

    if ismember(itrial,[1,68])
        blinksAt(1)= []; % remove first

    elseif itrial ==159
        % weird half blink.
        blinksAt= blinksAt([1,2,3,5]);

        blinksEnd= blinksEnd([1,2,4,6]);
    elseif itrial==164
        blinksAt= [blinksAt, 51];

    end

elseif strcmp(subjID, 'CP') && ismember(itrial,[45,104,122,137,165]);
    if ismember(itrial,[137,165])
        blinksAt(1)=[];
    else
        blinksEnd= maxl;
    end

elseif strcmp(subjID,'DA') &&  ismember(itrial, [172])

    if itrial==172
        blinksAt= [19; blinksAt];
    end
elseif strcmp(subjID, 'EK') && ismember(itrial, [36,170,172,177]);

    if itrial==177
        blinksAt(2)=[];

    else
        blinksAt(1) = [];% 36, 170
    end
elseif strcmp(subjID, 'JEPM') && ismember(itrial, [1,74])
    if itrial==1
        blinksEnd= [87; blinksEnd];
    elseif itrial==74
        blinksEnd= [55; blinksEnd];
    end
elseif strcmp(subjID, 'JH') && ismember(itrial, [46])
    blinksEnd= maxl;

elseif strcmp(subjID, 'JT') && ismember(itrial, [1,14,21,28,48,82,86,91,93,107,110,143,148, 165,195])
    if itrial==1
        blinksEnd= 30; % trial 1
    elseif itrial==14
        blinksEnd(1)=[];
    elseif ismember(itrial, [21,28,82,86,91,93,107,110,195])
        blinksAt(1)=[];
    elseif itrial==48
        blinksAt(end)=[];
    elseif itrial==143
        blinksAt=[1,277, 543, 657];
        blinksEnd=[21,292, 552,711];
    elseif itrial==148
        blinksAt=[ 84, 110, 146, 211,450, 613 ];
        blinksEnd=[ 92, 123, 182, 223,462, 622 ];
    elseif itrial==165
        blinksAt=[29,64,532, 575, 745];
        blinksEnd=[42,90,554, 590, 763];
    end
elseif strcmp(subjID, 'KK') && ismember(itrial,[1])
    if itrial==1
        blinksEnd= [43,726];
    end
elseif strcmp(subjID, 'LJTN') && ismember(itrial, [165])
    blinksAt=28; blinksEnd= 47;
elseif strcmp(subjID, 'LL1') && ismember(itrial, [180])
    blinksEnd = [9,66];
elseif strcmp(subjID, 'LYC') && ismember(itrial, [10,158,193])
    blinksAt(1) = [];
elseif strcmp(subjID, 'MD1') && ismember(itrial, [116,194])
    blinksAt(1) = [];
elseif strcmp(subjID, 'MD2') && ismember(itrial, [92,98,146,157,191])
    if ismember(itrial,92)
        blinksAt(1) = [];
    elseif itrial==98
        blinksAt = [134, 757];
    elseif itrial==146 || itrial==157
        blinksEnd(1)=[];
    elseif itrial==191
        blinksAt= [77, 668]
    end
elseif strcmp(subjID, 'NL') && ismember(itrial, [49,51])
    if itrial==49
        blinksAt(1)=[];
    elseif itrial==51
        blinksEnd(1)=[];
    end
elseif strcmp(subjID, 'QL') && ismember(itrial, [20])

    if itrial==20
        blinksAt(1)=[];
    end
elseif strcmp(subjID, 'RJS') && ismember(itrial, [33,44,49,80,85,120,125])

    if itrial==49
        blinksAt = [179, 245, 314, 460, 496];
        blinksEnd = [196,252, 328, 476, 512];
    elseif itrial == 80
        blinksEnd=[7,53,83];
    elseif itrial==85
        blinksAt = [66, 148, 543,676,784];
    elseif itrial==120
        blinksEnd = [72,576,660,773];
    else
        blinksAt(1)=[];
    end
elseif strcmp(subjID, 'WB') && ismember(itrial, [122])
    blinksAt(1)=[];
elseif strcmp(subjID, 'YL') && ismember(itrial, [114])
    blinksAt(1)=[];
else %
    %%

    figure(10); clf;
    plot(trial_EyeData.Y); hold on;
    %                     plot(Eyetrack);
    shg
    error(['check code: ' subjID ' ippant ' num2str(ippant) ' trial ' num2str(itrial)]);
end % query subj.
