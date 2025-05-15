% j7_Plot_Data_withinGait_PFits_v3

% 

%plots the psychometric fits per condition (walk, stand), as well as fits 
% as a function of gait cycle (splitting into thirds?).


setmydirs_detectv3;
cd(procdatadir)
pfols= dir([pwd  filesep '*summary_data.mat']);
nsubs= length(pfols);
%
%show ppant list:
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)

%%
job.concat_GFX=0;

%participant level effects:
job.plot_PsychmFits_pfx=0; % compare types (walking, standing, L/R ft).

% using the psignifit function library:
job.calc_Psychmfits_altversion =0; % trying out the psignifit function
job.plot_PsychmFits_pfx_alt=0;% follows the figure design from above, for comparison.

% using Palamaedes function library:
job.plot_PsychmFits_gaitquantile_pfx=0; % compare within Gait Cycle
% job.plot_PsychmFits_gaitquantile_pfx_alt=1; % compare within Gait Cycle

%group level effects:
% job.plot_Acc_labelled = 1; % scatter of PFX accuracy, per condition, labelled by subj ID.
job.plot_PsychmFits_GFX=0;
job.plot_PsychmFits_gait_GFX=0; % compare within Gait Cycle (LR, RL, combined).


job.plot_ConditionGFX_MSver1= 0; % GFX above.
job.plot_PsychmFits_gait_MSver1=0; %fits over the gait.. combined single steps.


job.plot_ConditionGFX_MSver2_altFits=1; % MSver with alternate psych fit routine. (precomputed in calc..)
job.plot_PsychmFits_gait_MSver2_AltFits=0; % MS version of fits over the gait.. combined single steps.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%data wrangling first: Concat, and create null distribution.

pidx1=ceil(linspace(1,100,11)); % length n-1
pidx2= ceil(linspace(1,200,21));% 
%%
nquants= 5; % how many quantiles to subdivide gait cycle?
gaittypes = {'single gait' , 'double gait'};
%%
if job.concat_GFX
    %preallocate storage:
    dataINrespPos=[];

    GFX_headY=[];
    GFX_TargPosData=[];
    subjIDs=cell(1,nsubs);

    for ippant =1:nsubs
        cd(procdatadir)    %%load data from import job.
        load(pfols(ippant).name, ...
            'trial_summaryTable', 'subjID', 'gait_ts_resamp', 'gait_ts_gData');

        subjIDs{ippant} = subjID;


        % first retrieve index for separate gaits (L/Right feet).
        Ltrials= strcmp(trial_summaryTable.trgO_gFoot, 'LR');
        Rtrials= strcmp(trial_summaryTable.trgO_gFoot, 'RL');

        % mean head pos:
        GFX_headY(ippant).gc = mean(gait_ts_resamp,1, 'omitnan');

        %same for ts data
        Ltrials_ts= strcmp(gait_ts_gData.gaitFeet, 'LR');
        Rtrials_ts= strcmp(gait_ts_gData.gaitFeet, 'RL');

        %can also create a doubled version: LRL
        GFX_headY(ippant).doubgc_LRL= [mean(gait_ts_resamp(Ltrials_ts,:),1,'omitnan'), mean(gait_ts_resamp(Rtrials_ts,:),1,'omitnan')];

        %% for standing still data, prepare for PF fits.
        T=trial_summaryTable;

        %note that per ppant, we may also need to reject trials from
        %         rejTrials_questv3;

        exprows = find(T.isPrac==0);
        exprows=exprows(2:end); % remove first target out of staircase.
        strows= find(T.isStationary==1);
        wkrows= find(T.isStationary==0);


        %intersect of experiment, and relevant condition:
        standingrows = intersect(strows, exprows);
        walkingrows = intersect(wkrows, exprows);

        % for entire condition (wlk. stand. without split by gait details),
        % extract the data (for accuracy and RT calcs).
        useRows=[];
        useRows{1} = standingrows;
        useRows{2} = walkingrows;
        %% REGION Condition level analysis:
        for iStWlk=1:2
            extractrows = useRows{iStWlk};

            % collect data for standing portion:

            trgContrIDX = T.targContrastPosIdx(extractrows);
            trgCorr= T.targCor(extractrows);


            OutOfNum = ones(1, size(trgContrIDX,1));
            [~, NumPer, TotalPer] = PAL_MLDS_GroupTrialsbyX(trgContrIDX, trgCorr,...
                OutOfNum);

            % also collect RTs per contr.
            trgRT = T.clickRT(extractrows);
            trgContr = T.targContrast(extractrows);

            [meanRTs, meanContr]=deal(zeros(1,7));

            for itargC=1:7
                tmpr=find(trgContrIDX==(itargC-1));

                % avoid negative (these were missed targs).

                meanRTs(itargC) = nanmean(trgRT(tmpr));
                meanContr(itargC)= nanmean(trgContr(tmpr));

            end

            %store data:
            if iStWlk==1
                NumPerContr_standing = NumPer;
                TotalPerContr_standing = TotalPer;
               
                meanRT_standing = meanRTs;
                meanContr_standing= meanContr;
            else
                NumPerContr_wlking = NumPer;
                TotalPerContr_wlking=TotalPer;

                meanRT_wlking = meanRTs;
                meanContr_wlking= meanContr;
            end
        end
        %end region condition level:
        %% REGION:  gait split analysis

        %% also precompute for different gait sizes (1-2), and L/R leading feet:
        nGait=1;
        pidx=pidx1; % indices for binning (single or double gc).
        useL=Ltrials;
        useR=Rtrials;
        alltrials = 1:length(Ltrials);
        trialstoIndex ={useL, useR, alltrials};

        % split by gait/ foot (L/R)
        for iLR=1:3
            uset = trialstoIndex{iLR};

            %% %%%%%%%%%%%%%%%%
            % Step through different data types :
            %%%%%%%%%%%%%%%%%%

            %% store the histcounts, for target contrast level per gaitsize,
            % and step (L/R).
            %Target onset data:
            targContrIDX = T.targContrastPosIdx(uset);
            targPos= T.trgO_gPcnt(uset); % gait pcnt for target onset.
            targCorr = T.targCor(uset);

            targContr = T.targContrast(uset);
            targRT = T.clickRT(uset);
            %Each entry of StimLevelsall corresponds to single trial
            OutOfNum = ones(1,size(targContr,1));

            % pre compute result per subj
            [ValPerContr_iLR, NumPerContr_iLR, TotalPerContr_iLR] = PAL_MLDS_GroupTrialsbyX(targContrIDX, targCorr,...
                OutOfNum);
            RTPerContr_iLR=zeros(1,7);
            %overwrite value per (currently an index), with mean contr per:
            for itargC=1:7
                tmpr=find(targContrIDX==(itargC-1));

                % avoid negative (these were missed targs).
                ValPerContr_iLR(itargC)= mean(targContr(tmpr), 'omitnan');
                RTPerContr_iLR(itargC) = mean(targRT(tmpr), 'omitnan');
            end

            % end region (iGC)

            %% Now also compute, but restrict to gait quantiles.

            [ValPerContr_gaitQntl, NumPerContr_gaitQntl, TotalPerContr_gaitQntl,RTPerContr_gaitQntl] = deal([]);

            qntlBounds = round(quantile(0:pidx(end),nquants-1)); %
            pcntBounds = [1, qntlBounds, pidx(end)];
            %n bounds,
            for iq=1:length(pcntBounds)-1

                tmpA = find(targPos>pcntBounds(iq));
                tmpB = find(targPos<=pcntBounds(iq+1));
                useC = intersect(tmpA,tmpB);

                trgContrIDX_in = targContrIDX(useC);
                trgContr_in = targContr(useC);
                trgCorr_in = targCorr(useC);
                trgRT_in = targRT(useC);
                OutOfNum= ones(1,size(trgContrIDX_in,1));

                [ValPerContr_gaitQntl(iq,:),...
                    NumPerContr_gaitQntl(iq,:), ...
                    TotalPerContr_gaitQntl(iq,:)] = PAL_MLDS_GroupTrialsbyX(trgContrIDX_in, trgCorr_in,...
                    OutOfNum);

                for itargC=1:7
                    tmpr=find(trgContrIDX_in==(itargC-1));

                    % avoid negative (these were missed targs).
                    ValPerContr_gaitQntl(iq,itargC)= mean(trgContr_in(tmpr), 'omitnan');
                    RTPerContr_gaitQntl(iq,itargC) = mean(trgRT_in(tmpr), 'omitnan');
                end

            end

            %% store


            %using Targ pos as index:            
            %stationary:
            GFX_TargPosData(ippant,iLR).gc_ContrVals_allstatnry= meanContr_standing;
            GFX_TargPosData(ippant,iLR).gc_NumPerContr_allstatnry= NumPerContr_standing;
            GFX_TargPosData(ippant,iLR).gc_TotalPerContr_allstatnry= TotalPerContr_standing;
            GFX_TargPosData(ippant,iLR).gc_RTperContr_allstatnry=meanRT_standing;
            
            %walking (all)
            %the counts per gait "%"
            GFX_TargPosData(ippant,iLR).gc_ContrVals_allwlking= meanContr_wlking;
            GFX_TargPosData(ippant,iLR).gc_NumPerContr_allwlking = NumPerContr_wlking;
            GFX_TargPosData(ippant,iLR).gc_TotalPerContr_allwlking = TotalPerContr_wlking;
            GFX_TargPosData(ippant,iLR).gc_RTperContr_allwlking =meanRT_wlking;
            
            %walking (split by LR)
            GFX_TargPosData(ippant,iLR).gc_ContrVals_LRwlking= ValPerContr_iLR;
            GFX_TargPosData(ippant,iLR).gc_NumPerContr_LRwlking = NumPerContr_iLR;
            GFX_TargPosData(ippant,iLR).gc_TotalPerContr_LRwlking = TotalPerContr_iLR;
            GFX_TargPosData(ippant,iLR).gc_RTperContr_LRwlking =RTPerContr_iLR;

            %walking (split by LR, and  % gait cycle)
             GFX_TargPosData(ippant,iLR).gc_ContrVals_qntlwlking= ValPerContr_gaitQntl;
            GFX_TargPosData(ippant,iLR).gc_NumPerContr_qntlwlking = NumPerContr_gaitQntl;
            GFX_TargPosData(ippant,iLR).gc_TotalPerContr_qntlwlking = TotalPerContr_gaitQntl;
            GFX_TargPosData(ippant,iLR).gc_RTperContr_qntlwlking =RTPerContr_gaitQntl;

            GFX_TargPosData(ippant,iLR).gc_qntl_bounds= pcntBounds;


        end % iLR


disp(['fin concat for ppant ' subjID]);
    end % ppant


    cd([procdatadir filesep 'GFX']);
    save('GFX_Data_inGaits_PFfits', ...
        'GFX_headY', 'GFX_TargPosData',...
        'subjIDs');%, '-append');
else
    
%%
pidx1=ceil(linspace(1,100,11)); % length n-1
pidx2= ceil(linspace(1,200,21));% 
%
nquants= 5; % how many quantiles to subdivide gait cycle?
gaittypes = {'single gait' , 'double gait'};
    cd([procdatadir filesep 'GFX']);
    
    load('GFX_Data_inGaits_PFfits');
    load('GFX_Data_inGaits_SigniFITs.mat')

end
%% now for plotting jobs
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% plots at Participant level:

%%
%PFX
if job.plot_PsychmFits_pfx% target onset relative to gait.
 %% for each participant, plot the PM fits for standing and walking
 % compares L and R ft also.
 cfg=[];
 cfg.subjIDs = subjIDs;
 cfg.type = 'Target';
 cfg.figdir= figdir; % for orienting to figures folder
 cfg.HeadData= GFX_headY;
 cfg.pidx1= pidx1;
 cfg.pidx2= pidx2;
 cfg.plotlevel='PFX';
 % cycles through ppants, plots with correct labels.
 plot_PFfits_comparison(GFX_TargPosData, cfg);
 
end

if job.plot_PsychmFits_gaitquantile_pfx
    %% 
    cfg=[];
 cfg.subjIDs = subjIDs;
 cfg.type = 'Target';
 cfg.datadir= datadir; % for orienting to figures folder
 cfg.HeadData= GFX_headY;
 cfg.pidx1= pidx1;
 cfg.pidx2= pidx2;
 cfg.plotlevel='PFX';
 cfg.figdir= figdir;
 % cycles through ppants, plots with correct labels.
 plot_PFfits_gaitcycle(GFX_TargPosData, cfg);
 
end

% % % alternate version using the Psignifit functions:

%PFX
if job.calc_Psychmfits_altversion% target onset relative to gait.
 %% for each participant, plot the PM fits for standing and walking
 % compares L and R ft also.
 cfg=[];
 cfg.subjIDs = subjIDs;
 cfg.type = 'Target';
 cfg.figdir= figdir; % for orienting to figures folder
 cfg.HeadData= GFX_headY;
 cfg.pidx1= pidx1;
 cfg.pidx2= pidx2;
 cfg.plotlevel='PFX';
 % cycles through ppants, saves at GFX level to be called for plots.
 calc_PFfits_psignifit(GFX_TargPosData, cfg);
 
end
if job.plot_PsychmFits_pfx_alt% target onset relative to gait.
 %% for each participant, plot the PM fits for standing and walking
 % compares L and R ft also.
 cfg=[];
 cfg.subjIDs = subjIDs;
 cfg.type = 'Target';
 cfg.figdir= figdir; % for orienting to figures folder
 cfg.HeadData= GFX_headY;
 cfg.pidx1= pidx1;
 cfg.pidx2= pidx2;
 cfg.plotlevel='PFX';
 % cycles through ppants, plots with correct labels.
 plot_PFfits_comparison_altver(GFX_signifit, cfg);
 
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GFX : group effects plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if job.plot_PsychmFits_GFX% target onset relative to gait.
 %% for each participant, plot the PM fits for standing and walking
 % compares L and R ft also. (slow)
 cfg=[];
 cfg.subjIDs = subjIDs;
 cfg.type = 'Target';
 cfg.figdir= figdir; % for orienting to figures folder
 cfg.HeadData= GFX_headY;
 cfg.pidx1= pidx1;
 cfg.pidx2= pidx2;
 cfg.plotlevel='GFX';
 % cycles through ppants, plots with correct labels.
 plot_PFfits_comparison(GFX_TargPosData, cfg);
 
end

%%

if job.plot_PsychmFits_gait_GFX
    %% 
    cfg=[];
 cfg.subjIDs = subjIDs;
 cfg.type = 'Target';
 cfg.datadir= datadir; % for orienting to figures folder
 cfg.HeadData= GFX_headY;
 cfg.pidx1= pidx1;
 cfg.pidx2= pidx2;
 cfg.plotlevel='GFX';
 cfg.figdir=figdir;
 % cycles through ppants, plots with correct labels.
 plot_PFfits_gaitcycle(GFX_TargPosData, cfg);
 
end
%%
if job.plot_ConditionGFX_MSver1==1; % modelled off the above.
%% for each participant, plot the PM fits for standing and walking
 % compares L and R ft also. (slow)
 cfg=[];
 cfg.subjIDs = subjIDs;
 cfg.type = 'Target';
 cfg.figdir= figdir; % for orienting to figures folder
 cfg.HeadData= GFX_headY;
 cfg.pidx1= pidx1;
 cfg.pidx2= pidx2;
 cfg.plotlevel='GFX';
 % cycles through ppants, plots with correct labels.
 plot_PFfits_comparison_MSVer1(GFX_TargPosData, cfg);
 
end



if job.plot_ConditionGFX_MSver2_altFits==1; % modelled off the above.
%% for each participant, plot the PM fits for standing and walking
 % compares L and R ft also. (slow)
 cfg=[];
 cfg.subjIDs = subjIDs;
 cfg.type = 'Target';
 cfg.figdir= figdir; % for orienting to figures folder
 cfg.HeadData= GFX_headY;
 cfg.pidx1= pidx1;
 cfg.pidx2= pidx2;
 cfg.plotlevel='GFX';
 % cycles through ppants, plots with correct labels.
 plot_PFfits_comparison_MSVer1_altFits(GFX_signifit, cfg);
 




end
if job.plot_PsychmFits_gait_MSver1==1; % MS version. combined single.
    %%
    cfg=[];
 cfg.subjIDs = subjIDs;
 cfg.type = 'Target';
 cfg.datadir= datadir; % for orienting to figures folder
 cfg.HeadData= GFX_headY;
 cfg.pidx1= pidx1;
 cfg.pidx2= pidx2;
 cfg.plotlevel='GFX';
 cfg.figdir=figdir;
 % cycles through ppants, plots with correct labels.
 plot_PFfits_gaitcycle_GFX_MSver1(GFX_TargPosData, cfg);
end

if job.plot_PsychmFits_gait_MSver2_AltFits==1; % MS version. combined single.
    %%
    cfg=[];
 cfg.subjIDs = subjIDs;
 cfg.type = 'Target';
 cfg.datadir= datadir; % for orienting to figures folder
 cfg.HeadData= GFX_headY;
 cfg.pidx1= pidx1;
 cfg.pidx2= pidx2;
 cfg.plotlevel='GFX';
 cfg.figdir=figdir;

 cfg.debug=0; % produce secondary plots of individual fits.
 cfg.FitatGroup= 0; % apply fitting procedure to group level proportion correct
 cfg.shareX =1; % combine on the same points of the  X axis.
 cfg.cleanOnly=1; % perform participant outlier rejection.
 cfg.applyOffset=0;
 cfg.FitsfromSharedX=1;
 cfg.fontsize= 15;
 % cycles through ppants, plots with correct labels.
 plot_PFfits_gaitcycle_GFX_MSver2_AltFits(GFX_signifit_byGait, cfg);
end
