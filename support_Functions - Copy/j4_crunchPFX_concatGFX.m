% j4_crunchPFX_concatGFX


% Here we load all data from previous analyses, preparing to save before
% plotting in next round of jobs.

%%%%%% QUEST DETECT version 3 (wEEG and eye) %%%%%%



setmydirs_detectv3;
% show ppant numbers:
cd(procdatadir)
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%   


%% concat data:
%preallocate storage:

pidx1=ceil(linspace(1,100,21)); % length n-1
pidx2=ceil(linspace(1,100,41));% 
useindexes = {pidx1, pidx2};
gaittypes = {'single gait' , 'double gait'};
%

GFX_headY=[];
GFX_EyeD=[];
GFX_TargPosData=[];
GFX_RespPosData=[];
subjIDs=cell([1, length(pfols)]);

avWalkParams=[] % added post hoc to calculate average stride length (duration sec).

for ippant =1:nsubs
    cd(procdatadir)    
    load(pfols(ippant).name, 'subjID', 'trial_summaryTable',...
        'gait_ts_gData', 'gait_ts_resamp','gait_ts_eyedirY_resamp', ...
        'gait_ts_eyedirZ_resamp', 'gait_ts_eyeposY_resamp', 'gait_ts_eyeposZ_resamp',...
        'doubgait_ts_resamp');
  
    
    subjIDs{ippant} = subjID;
    
    disp(['concatenating subject ' subjID]);
    
    % first retrieve index for separate gaits (L/Right feet).
    Ltrials= strcmp(trial_summaryTable.trgO_gFoot, 'LR');
    Rtrials= strcmp(trial_summaryTable.trgO_gFoot, 'RL');

    % mean head pos:
    GFX_headY(ippant).gc = mean(gait_ts_resamp,1, 'omitnan');
     %can also create a doubled version: 
    GFX_headY(ippant).doubgc= mean(doubgait_ts_resamp,1, 'omitnan');

    GFX_EyeD(ippant).gc_dirY = mean(gait_ts_eyedirY_resamp,1, 'omitnan');
    GFX_EyeD(ippant).gc_dirZ = mean(gait_ts_eyedirZ_resamp,1, 'omitnan');
    GFX_EyeD(ippant).gc_posY = mean(gait_ts_eyeposY_resamp,1, 'omitnan');
    GFX_EyeD(ippant).gc_posZ = mean(gait_ts_eyeposZ_resamp,1, 'omitnan');


    Ltrials_ts= strcmp(gait_ts_gData.gaitFeet, 'LR');
    Rtrials_ts= strcmp(gait_ts_gData.gaitFeet, 'RL');
    avWalkParams(ippant)= mean(gait_ts_gData.gaitDuration(Ltrials_ts)) + mean(gait_ts_gData.gaitDuration(Rtrials_ts)); 
   
    
    %% also calculate binned versions
    % These take the mean over an index range, per gait cycle position
    
    for nGait=1:2  % use to have option for 1 or 2.
  
        pidx=useindexes{nGait}; % indices for binning (single or double gc).
        ppantData= trial_summaryTable;
    
        if nGait==1
        useL=Ltrials;   
        useR=Rtrials;    
        useAll = 1:length(Ltrials);

        % split by foot (LR, RL, all combined)
        else


            %%% % % % ! ! ! ! ! ! ! UNFINISHED
        % note that here we have data in separate columns (so find non-nan
        % rows.)
        useL =~isnan(trial_summaryTable.trgO_gPcnt_LRL);
        useR =~isnan(trial_summaryTable.trgO_gPcnt_RLR);
        useall= 1:length(Ltrials);


        end 
        trialstoIndex ={useL, useR, useAll};

        for iLR=1:3
            uset= trialstoIndex{iLR};
            
            %% %%%%%%%%%%%%%%%%
            % Step through different data types :
            
            % note that now, we will sample randomly from all
            % positions, based on the number in real data.
            
            % 3 main analysis (for behaviour).
            %- RT relative to target gait position
            %- RT relative to response gait position
            %- Acc relative to target gait position.
            % nb: acc relative to response, is moot, since incorrect responses
            %(False Alarms) were so few.
            %%%%%%%%%%%%%%%%%%
            % for each of these three types, define the search values
            % to shuffle over:
            
            if nGait==1
                searchPosIDX = {ppantData.trgO_gPcnt(uset), ... % trg onset as %
                    ppantData.respO_gPcnt(uset), ... % resp O as pcnt
                    ppantData.trgO_gPcnt(uset)}; % targ onset as %
            else %(Much messier, but due to overlap in events occurring across 2 gaits.)
                if iLR==1
                    searchPosIDX = {ppantData.trgO_gPcnt_LRL(uset), ... % trg onset as %
                        ppantData.respO_gPcnt_LRL(uset), ... % resp O as pcnt
                        ppantData.trgO_gPcnt_LRL(uset)}; % targ onset as %
                elseif iLR==2
                    searchPosIDX = {ppantData.trgO_gPcnt_RLR(uset), ... % trg onset as %
                        ppantData.respO_gPcnt_RLR(uset), ... % resp O as pcnt
                        ppantData.trgO_gPcnt_RLR(uset)}; % targ onset as %


                elseif iLR==3
                    % both. Problem is they can have identical row info.

                    searchPosIDX = {[ppantData.trgO_gPcnt_LRL(useL); ppantData.trgO_gPcnt_RLR(useR)], ... % trg onset as %
                        [ppantData.respO_gPcnt_LRL(useL); ppantData.respO_gPcnt_RLR(useR)], ... % resp O as pcnt
                        [ppantData.trgO_gPcnt_LRL(useL);ppantData.trgO_gPcnt_RLR(useR)]}; % targ onset as %


                end
            end


            % normally, we can match the gait percent (event) to the DV:
            searchPosDVs = {ppantData.clickRT(uset), ... %RT, RT, Acc.
                ppantData.clickRT(uset), ...
                ppantData.targCor(uset)};


            % %except in this particular case:
            if iLR==3 && nGait==2
           
             %and define the DV of interest:
            searchPosDVs = {[ppantData.clickRT(useL);ppantData.clickRT(useR)], ... %RT, RT, Acc.
                [ppantData.clickRT(useL);ppantData.clickRT(useR)],  ...
                [ppantData.targCor(useL); ppantData.targCor(useR)]};

            end




            for itype=1:3
                %preallocate:
                outgoing = nan(1,length(pidx)); % 
                allpos = searchPosIDX{itype};
                allDVs = searchPosDVs{itype};
                
                
                % for actual gait pos, select reldata
                for ip=1:pidx(end)
                    useb = find(allpos==ip);
                    
                    
                    %% take mean for these  trials (if RT):
                    if itype<=2
                        outgoing(ip) = mean(allDVs(useb),1, 'omitnan');
                        
                    elseif itype==3
                        % else compute accuracy.
                        
                        % remove nans from length calcs.
                        tmpResp = allDVs(useb);
                        useResp = tmpResp(~isnan(tmpResp));
                        outgoing(ip)= sum(useResp, 'omitnan') / length(useResp);
                    end
                    
                end
                
                switch itype
                    case 1
                        targOns_RTs = outgoing;
                    case 2
                        respOns_RTs= outgoing;
                    case 3
                        targOns_Acc= outgoing;
                end
            end % itype
            
            
          
            %% store the histcounts, for target contrast level per pos.
            tPos= ppantData.trgO_gPcnt(uset); % gait pcnt for target onset.
            targContrIDX = ppantData.targContrastPosIdx(uset);
            targContr = ppantData.targContrast(uset);
            
            % % again, the special case:
            if nGait==2
                if iLR==1
                    tPos= ppantData.trgO_gPcnt_LRL(uset);
                elseif iLR==2
                    tPos= ppantData.trgO_gPcnt_RLR(uset);
                elseif iLR==3
                    tPos= [ppantData.trgO_gPcnt_LRL(useL);ppantData.trgO_gPcnt_RLR(useR)];
                      targContrIDX = [ppantData.targContrastPosIdx(useL);ppantData.targContrastPosIdx(useResp)];
                      targContr = [ppantData.targContrast(useL);ppantData.targContrast(useR)];
                end

            end
            
            [tgcMat]= deal(nan(7,pidx(end)));
            tcontrasts = nan(1,7);
            for  icontr= 1:7
                tindx = find(targContrIDX==(icontr-1));
                tgcMat(icontr,:) = histcounts(tPos(tindx), pidx(end));
                %store the mean of all contrasts per indx:
                tcontrasts(icontr) = mean(targContr(tindx));
            end
            
            %% store histcounts, per target onset pos, and resp onset pos.
            clkPos= ppantData.respO_gPcnt(uset); % gait pcnt for response onset.

            % % correct for special case: 
            if nGait==2
                if iLR==1
                    clkPos= ppantData.respO_gPcnt_LRL(uset);
                elseif iLR==2
                    clkPos= ppantData.respO_gPcnt_RLR(uset);
                elseif iLR==3
                    clkPos= [ppantData.respO_gPcnt_LRL(useL);ppantData.respO_gPcnt_RLR(useR)];
                end

            end
             % % % !!

            usetypes = {tPos, clkPos};
            for itype=1:2
                datain = usetypes{itype};
                % remove zero pos points before counting:
                datain(datain==0)=[];
                datain(isnan(datain))=[];
                % raw hit counts, per gaitsample:
                tmp_counts = histcounts(datain, pidx(end));
                
                %critical step -convert 0 to nans. this avoids averaging
                %problems.
                tmp= tmp_counts; % trgO_counts, / clkOcounts saved below.
                tmp(tmp==0)=NaN;
                tmptoAvg=tmp;
                
                % rename
                if itype==1
                    trgO_counts = tmp_counts;
                    trgOtoAvg = tmptoAvg; % with nans removed.
                elseif itype==2
                    respO_counts = tmp_counts;
                    respOtoAvg = tmptoAvg;                
                end
            end
            %% %%%%%%%%%%%%%%%%%%%
            % Step through all data types, performing binning based on
            % length of pidx.
            %%%%%%%%%%%%%%%%%%%
            %trg onset counts
            %resp onset counts
            %acc per trg onset
            %acc per resp onset
            %rt per trg onset
            %rt per resp onset
            usedata= {trgOtoAvg, respOtoAvg, ... % using ver with NaNs removed
                targOns_Acc,targOns_RTs, respOns_RTs};
            
            for itype= 1:5
                
                %% for all types, take mean over bin indices.
                out_bin=nan(1,length(pidx)-1);
                datatobin = usedata{itype};
                for ibin=1:length(pidx)-1
                    idx = pidx(ibin):pidx(ibin+1);
                    
                    out_bin(ibin) = mean(datatobin(idx), 'omitnan');
                    
                end
                
                switch itype
                    case 1
                        trgCounts_bin = out_bin;
                    case 2
                        respCounts_bin = out_bin;
                    case 3
                        targOns_Acc_bin = out_bin;                   
                        
                    case 4
                        targOns_RT_bin= out_bin;
                    case 5
                        respOns_RT_bin = out_bin;
                end
                
            end
            %% %%%%%%%%%%%%%%%%%%%%%%%%% now save:
            gaitnames = {'gc_', 'doubgc_'}; % change fieldname based on nGaits:
            
            %add targ data to one matrix, resp data to another:
            %trg contrast info:
            GFX_TargPosData(ippant,iLR).([gaitnames{nGait} 'trgcontr'])= tcontrasts;
            GFX_TargPosData(ippant,iLR).([gaitnames{nGait} 'trgcontrIDXmatrx']) = tgcMat;
            
            % rest of the relevant (behavioural outcomes) data:
            GFX_TargPosData(ippant,iLR).([gaitnames{nGait} 'counts']) = trgO_counts;
            GFX_TargPosData(ippant,iLR).([gaitnames{nGait} 'Acc']) = targOns_Acc;
            GFX_TargPosData(ippant,iLR).([gaitnames{nGait} 'rts']) = targOns_RTs;
            
            
            GFX_TargPosData(ippant,iLR).([gaitnames{nGait} 'binned_counts']) = trgCounts_bin;
            GFX_TargPosData(ippant,iLR).([gaitnames{nGait} 'binned_Acc']) = targOns_Acc_bin;
            GFX_TargPosData(ippant,iLR).([gaitnames{nGait} 'binned_rts']) = targOns_RT_bin;
            
            %also using response position as index:
            GFX_RespPosData(ippant,iLR).([gaitnames{nGait} 'counts']) = respO_counts;
%             GFX_RespPosData(ippant,iLR).([gaitnames{nGait} 'FAs']) = FAPos_counts;
%             GFX_RespPosData(ippant,iLR).([gaitnames{nGait} 'Acc']) = respOns_Acc;
            GFX_RespPosData(ippant,iLR).([gaitnames{nGait} 'rts']) = respOns_RTs;
            
            GFX_RespPosData(ippant,iLR).([gaitnames{nGait} 'binned_counts']) = respCounts_bin;
%             GFX_RespPosData(ippant,iLR).([gaitnames{nGait} 'binned_Acc']) = respOns_Acc_bin;
            GFX_RespPosData(ippant,iLR).([gaitnames{nGait} 'binned_rts']) = respOns_RT_bin;
            
            
        end % iLR
    end % nGaits.
    
% % after all said and done, concatenate the two single gaits to make a
% % stride (for plots).
% % store as a new double gc (i.e. a stride).
%     
% useD = {GFX_TargPosData, GFX_RespPosData};
% %store field names, but only first time through, otherwise we will append
% %too many times.
% if ippant==1 
% allfieldsTarg = fieldnames(GFX_TargPosData);
% allfieldsResp= fieldnames(GFX_RespPosData);
% end
% 
% for idatatype = 1:2
%     
%     if idatatype==1
%         allfields= allfieldsTarg;
%     else
%         allfields= allfieldsResp;
%     end
%   
% 
%     for ifield = 1:length(allfields)
%         % per field, store with new field name prefix (doubgc_) and
%         % horizonrally concatenate the data
%         if ifield==4
%             pause(.1)
%         end
%         useD{idatatype}(ippant,4).(['doub' allfields{ifield}]) = ...
%             [useD{idatatype}(ippant,1).([allfields{ifield}]), ... %LR
%             useD{idatatype}(ippant,2).([allfields{ifield}])]; %RL
%     end
% 
% 
%    
% end
%  %rename to save  with new info   
% GFX_TargPosData = useD{1};
% GFX_RespPosData= useD{2};
%     
end % ppant
%% save GFX
cd([procdatadir filesep 'GFX']);
disp('Saving GFX');
save('GFX_Data_inGaits', ...
    'GFX_headY', 'GFX_TargPosData','GFX_RespPosData',...
    'subjIDs', 'pidx1', 'pidx2', 'gaittypes','GFX_EyeD');%, '-append');

