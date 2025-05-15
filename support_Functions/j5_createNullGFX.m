% j5_createNullGFX
%% similar to the concat job (j4), except we disregard gait position information, 
 % and instead sample uniformly from all positions to create a null distribution.
%  [nbins x nPerm] per subject and datatype.



setmydirs_detectv3;
%% show ppant numbers:
cd(procdatadir);
pfols = dir([pwd filesep '*summary_data.mat']);
nsubs= length(pfols);
tr= table((1:length(pfols))',{pfols(:).name}' );
disp(tr)
%

%% remove based on gazeDVA
remGaze=1;


%try randomly shuffling all stride-cycle onsets across the experiment.

%% concat data:
%preallocate storage:

gaittypes = {'single gait' , 'double gait'};
%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   %preallocate storage:
   
     %use same spacing as observed data, for the bins:
     %add extra fields for the shuffled data.
     cd([procdatadir filesep 'GFX'])
    load('GFX_Data_inGaits', 'GFX_TargPosData', 'GFX_RespPosData','pidx1', 'pidx2');
    
%     GFX_TargPos_nullData=[]; % store the results per subj, for stats.
%     GFX_RespPos_nullData=[];
    
    nPerm = 1000; % how many times to perform the resampling?
    
    for ippant =1:length(pfols)
        
        cd(procdatadir)    %%load data from import job.
        load(pfols(ippant).name, ...
            'trial_summaryTable','subjID');
        
        if remGaze==1
            %remtrials
            remtrials= find(trial_summaryTable.targOutsideDVA==1);
            trial_summaryTable(remtrials,:) = [];
            
            
        end
        
        Ppant_gaitData= trial_summaryTable;
    
        Ltrials= strcmp(trial_summaryTable.trgO_gFoot, 'LR');
    
        Rtrials= strcmp(trial_summaryTable.trgO_gFoot, 'RL');
       
        disp(['constructing null distributions for ' subjID]);
 
        %%

        useindexes= {pidx1, pidx2};
        for nGait=1:2
        pidx=useindexes{nGait}; % indices for binning (single or double gc).
        ppantData= Ppant_gaitData;
        useL=Ltrials;
        useR=Rtrials;
        useAll = 1:length(Ltrials);

              
             
        trialstoIndex_gait ={useL, useR, useAll};
        
        % split by gait/ foot (L/R)
        for iLR=3%1:3
           
               uset=find(trialstoIndex_gait{iLR});
                
                
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

                searchPosIDX = {ppantData.trgO_gPcnt(uset), ... % trg onset as %
                    ppantData.respO_gPcnt(uset), ... % resp O as pcnt
                    ppantData.trgO_gPcnt(uset)}; % targ onset as %

                %and define the DV of interest:
                searchPosDVs = {ppantData.clickRT(uset), ... %RT, RT, Acc.
                    ppantData.clickRT(uset), ...
                    ppantData.targCor(uset)};
                
                for itype=1:3
                    %preallocate:
                    outgoingnull = nan(nPerm, pidx(end));
                    null_CVs = nan(3,nPerm); % 5, median, and 95% intervals.
                    
                    allpos = searchPosIDX{itype};
                    allDVs = searchPosDVs{itype};
                    allpos_real = allpos(~isnan(allpos));
                    allDVs_real = allDVs(~isnan(allpos));
                
                % for actual gait pos, select another at random:
                for ip=1:pidx(end)     
                     
%                     % we may want to match the N datapoints:
%                     % how many observed data points this index?
%                      actualIDX=find(allpos_real==ip);
%                      actualN = length(actualIDX);
%                      actualObs = nanmean(allRTs_real(actualIDX));

                        tmpstor=nan(1, nPerm);                     
                    for iperm = 1:nPerm
                        % instead of the actual index (ip), replace with an indices at
                        % random from all possible:
                        useIDX = randi(pidx(end),1);
                        
%                         thisrand = randi(length(allpos_real),[1,actualN]);
                        thisrand= find(allpos_real==useIDX);
                        
                        %% take mean for these random trials (if RT):                       
                        if itype<=2
                        tmpstor(iperm) = nanmean(allDVs_real(thisrand));
                        
                        elseif itype==3 % else compute accuracy
                            
                        tmpstor(iperm)= nansum(allDVs_real(thisrand) / length(allDVs_real(thisrand))); 
                        end
                    end
%                    
                    %store null
                    outgoingnull(:, ip) = tmpstor;
                end
                    %could also store 95% CritValues at this stage:                    
%                     null_CVs(1:3,ip) =quantile(tmpstor, [.05, .5, .95]);

                    switch itype
                        case 1
                            targOns_RTs_null = outgoingnull;
                        case 2
                            respOns_RTs_null = outgoingnull;
                        case 3
                            targOns_Acc_null = outgoingnull;
                    end
                end % itype
                
                %% %%%%%%%%%%%%%%%% 
                % Step through slightly different analysis for shuffling
                % counts per timepoint.
                
                
                % 2 main analysis (for counts). 
                
                %-targOnset 
                %-clickOnset
                
                %work out what the count per actual IDX was, then shuffle
                %these counts. 
                
                %%%%%%%%%%%%%%%%%% 
                
                % now adding saccade counts.
               
                %% store the histcounts, for target contrast level per pos.
                usePos= {ppantData.trgO_gPcnt(uset),...% gait pcnt for target onset.                
                ppantData.respO_gPcnt(uset)}; % gaitpcnt for resp onset.
            
            for itype=1:2
                useposdata = usePos{itype};
                
                % raw hit counts, per gaitsample:
                real_counts= histcounts(useposdata, pidx(end));
                
                %critical step -convert 0 to nans. this avoids averaging
                %problems (as above, in data anaysis)
                
                real_counts(real_counts==0)=NaN;
                %now create a null, by shuffling these counts per perm.
                
                % preallocate:
                outgoingnull = nan(nPerm, pidx(end));
                for iperm = 1:nPerm
                    
                    
                    a = rand(1,length(real_counts));
                    [~, i] = sort(a);
                    outgoingnull(iperm,:) =real_counts(i);
                end
                
                switch itype
                    case 1
                        targOns_Counts_null = outgoingnull;
                    case 2
                        respOns_Counts_null = outgoingnull;
                end
            end % both types.
            
                %% %%%%%%%%%%%%%%%%%%%
                % Step through all shuffled data, storing binned versions, per
                % permutation previously performed.
                %%%%%%%%%%%%%%%%%%%
               
                binmedata = {targOns_RTs_null, respOns_RTs_null, targOns_Acc_null,...
                    targOns_Counts_null, respOns_Counts_null};
                
                for itype = 1:length(binmedata)
                   
                   usedata = binmedata{itype};                   
                   outgoingbinned=nan(nPerm, length(pidx)-1);
                   
                   for ibin=1:length(pidx)-1
                       idx = pidx(ibin):pidx(ibin+1);
                       
                       outgoingbinned(:,ibin) = nanmean(usedata(:,idx),2);
                       
                   end
                    
                   %rename this time. to avoid headaches
                   switch itype
                       case 1
                           trgRT_bin_shuff= outgoingbinned;
                       case 2
                           respRT_bin_shuff= outgoingbinned;
                       case 3
                           trgAcc_bin_shuff= outgoingbinned;
                       case 4
                           trgCounts_bin_shuff = outgoingbinned;
                       case 5
                           respCounts_bin_shuff = outgoingbinned;
                   end
                end % per type
                
               %% store for all ppants, preparing to save:
                gaitnames = {'gc_', 'doubgc_'}; % change fieldname based on nGaits:
                
                %add binned, shuffled vers:
                GFX_TargPosData(ippant,iLR).([gaitnames{nGait} 'binned_counts_shuff']) = trgCounts_bin_shuff;
                GFX_TargPosData(ippant,iLR).([gaitnames{nGait} 'binned_Acc_shuff']) = trgAcc_bin_shuff;
                GFX_TargPosData(ippant,iLR).([gaitnames{nGait} 'binned_rts_shuff']) = trgRT_bin_shuff;
                
                %also using click pos as index:                
                GFX_RespPosData(ippant,iLR).([gaitnames{nGait} 'binned_counts_shuff']) = respCounts_bin_shuff;           
                GFX_RespPosData(ippant,iLR).([gaitnames{nGait} 'binned_rts_shuff']) = respRT_bin_shuff;
                
                
            end % iLR
        end
 

    
    end % ppant
    
    %%
    cd([procdatadir filesep 'GFX']);
    
    disp('resaving GFX data with shuffle attached...');
    save('GFX_Data_inGaits', ...
      'GFX_TargPosData','GFX_RespPosData',...
        '-append');