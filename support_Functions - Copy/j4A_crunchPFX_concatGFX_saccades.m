% j4A_crunchPFX_concatGFX_saccades


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

GFX_SaccadePosData=[];
subjIDs=cell([1, length(pfols)]);


for ippant =1:nsubs
    cd(procdatadir)    
    load(pfols(ippant).name, 'subjID', 'saccadeSummaryTable');
    
    subjIDs{ippant} = subjID;
    
    disp(['concatenating subject ' subjID]);
    
    
  %% NOW  we will process the saccade onset data (where saccades occurred during the gait). 
    % also calculate binned versions
    % 
    ppantData= saccadeSummaryTable;
    
    
    for nGait=1:2  % use to have option for 1 or 2.
  
        pidx=useindexes{nGait}; % indices for binning (single or double gc).
    
        if nGait==1 % single step.
      
            % first retrieve index for separate gaits (L/Right feet).
            useL= strcmp(ppantData.saccOns_gFoot, 'LR');
            useR= strcmp(ppantData.saccOns_gFoot, 'RL');
            useAll = 1:length(useL);

        % split by foot (LR, RL, all combined)
        else

% As above, just use all trials per speed.
        useall= 1:length(useL);


        end 
        trialstoIndex ={find(useL), find(useR), useAll};
       


        for iLR=3% 1:3 % Left Right, both
            
            % use the index of both:
            uset=  trialstoIndex{iLR};
           
            
          %% store the saccade count per gPcnt (binned later).
           %% %%%%%%%%%%%%%%%%
            % Step through different data types :

            % 1 main analysis (for behaviour).
            %-saccade onset
           
            %%%%%%%%%%%%%%%%%%
           
            if nGait==1
            searchPosIDX = ppantData.saccOnset_gPcnt(uset); % targ onset as %

            else 

                searchPosIDX = [ppantData.saccOnset_gPcnt_step1inStride(uset); ppantData.saccOnset_gPcnt_step2inStride(uset)];

            end

            %and define the DV of interest:
            %note that in each case, its the saccade RT. we will subselect
            %rows based on saccade direction. 
        

            
                %preallocate:
                outgoing = nan(1,pidx(end));
                allpos = searchPosIDX;

               

            %% store the histcounts, for saccadess (all), toward, and away
            % 
            usetypes= {searchPosIDX};
            for itype=1%:length(usetypes)
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
            switch itype
                case 1
                    sacc_counts_all = tmp_counts;
                    sacc_toAvg_all = tmptoAvg; % with nans removed.
%                 case 2
%                     sacc_counts_toward = tmp_counts;
%                     sacc_toAvg_toward = tmptoAvg; % with nans removed.
% 
%                 case 3
%                     sacc_counts_away = tmp_counts;
%                     sacc_toAvg_away = tmptoAvg; % with nans removed.
            end

            end
            %% %%%%%%%%%%%%%%%%%%%
            % Step through all data types, performing binning based on
            % length of pidx.
            %%%%%%%%%%%%%%%%%%%
            % sacc all (counts)
            % sacc toward
            % sacc away

            % sacc all (rts)
            % sacc toward
            % sacc away

            usedata= {sacc_toAvg_all};%,...
%                 sacc_toAvg_toward,...
%                sacc_toAvg_away ,...
%                sacc_rts_all,...
%                sacc_rts_toward,...
%                sacc_rts_away};%
            
            for itype= 1%:6
                
                %% for all types, take mean over bin indices.
                out_bin=nan(1,length(pidx)-1);
                datatobin = usedata{itype};
                for ibin=1:length(pidx)-1
                    idx = pidx(ibin):pidx(ibin+1);
                    
                    out_bin(ibin) = mean(datatobin(idx), 'omitnan');
                    
                end
                
                switch itype
                    %counts
                    case 1
                        sacc_counts_all_bin = out_bin;
%                     case 2
%                         sacc_counts_toward_bin = out_bin;
%                     case 3
%                         sacc_counts_away_bin = out_bin;                   
%                   %RTs
%                     case 4
%                         sacc_rts_all_bin = out_bin;
%                     case 5
%                         sacc_rts_toward_bin = out_bin;
%                     case 6
%                         sacc_rts_away_bin = out_bin;                   
% %                   
                end
                
            end
            %% %%%%%%%%%%%%%%%%%%%%%%%%% now save:
            gaitnames = {'gc_', 'doubgc_'}; % change fieldname based on nGaits:
            
            %add sacc data
            % counts (relative to gait).
           GFX_SaccadePosData(ippant,iLR).([gaitnames{nGait} 'sacc_all_counts']) =sacc_counts_all;
           GFX_SaccadePosData(ippant,iLR).([gaitnames{nGait} 'sacc_all_binned_counts']) =sacc_counts_all_bin;

%           
          
            
        end % iLR
        
      



    end % nGaits.
    disp(['fin ippant ' num2str(ippant)])
end % ppant
%% save GFX
cd([procdatadir filesep 'GFX']);
disp('Saving GFX');
save('GFX_Data_inGaits', ...
   'GFX_SaccadePosData',...
   'pidx1', 'pidx2', 'gaittypes','-append');%, '-append');
